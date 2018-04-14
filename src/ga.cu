#include "kp.h"
#include "common.h"
#include <memory.h>
#include <cuda.h>
#include <curand_kernel.h>

#define get_bit(arr, i) (arr[i / 32] & (1 << (i % 32)))
#define set_bit(arr, i) (arr[i / 32] |= (1 << (i % 32)))
#define unset_bit(arr, i) (arr[i / 32] &= ~(1 << (i % 32)))

// Generate roughly 50/50 distribution of enabled and disabled gene-bits.
__device__ void
gen_population_nums(uint32_t *pop, uint32_t len, uint32_t genes, curandState_t *state)
{
  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < genes; j++) {
      pop[i * genes + j] = curand(state) % 2;
    }
  }
}

// Generate genes with a specified number of bits enabled.
__device__ void
gen_population_bits(uint32_t *pop, uint32_t len, uint32_t genes, uint32_t bits, curandState_t *state, int seed_adjust)
{
  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < bits; j++) {
      pop[i * genes + (curand(state) % genes)] = 1;
    }
  }
}

void
print_cell(uint32_t *cell, uint32_t len)
{
    for (uint32_t i = 0; i < len; i++) {
        printf("%u", cell[i]);
    }
    printf("\n");
}

void
print_population(uint32_t **pop, uint32_t len, uint32_t genes)
{
  for (uint32_t i = 0; i < len; i++) {
      print_cell(pop[i], genes);
  }
}

void
print_cell_nodes(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
{
    uint32_t value = 0, weight = 0;
    for (uint32_t i = 0; i < len; i++) {
        if (cell[i] == 1) {
            value += knapsack[i].value;
            weight += knapsack[i].weight;
            printf("{ value: %u\tweight: %u }\n", knapsack[i].value, knapsack[i].weight);
        }
    }
    printf("value: %u\t weight: %u\n", value, weight);
}

knapsack_node_t
test_fitness(uint32_t *cell, uint32_t len,
             knapsack_node_t *knapsack, uint64_t max_weight)
{
  knapsack_node_t score = {0, 0};

  for (uint32_t i = 0; i < len; i++) {
    if (cell[i] == 1) {
      score.weight += knapsack[i].weight;
      score.value += knapsack[i].value;

      // Note that breaking once we hit the max weight only counts the values up
      // to that point, and so might induce bias.
      if (score.weight > max_weight) {
        score.value = score.value / 2;
        break;
      }
    }
  }

  if (score.weight == 0) {
      return score;
  }

  // Prioritize the value, but give a bump for a better value-to-weight ratio.
  //return value_total + (uint64_t) (value_total / weight_total);
  return score;
}

__device__ knapsack_node_t
test_fitness_kernel(uint32_t *cell, uint32_t len,
                    knapsack_node_t *knapsack, uint64_t max_weight)
{
  knapsack_node_t score = {0, 0};

  for (uint32_t i = 0; i < len; i++) {
    if (cell[i] == 1) {
      score.weight += knapsack[i].weight;
      score.value += knapsack[i].value;

      // Note that breaking once we hit the max weight only counts the values up
      // to that point, and so might induce bias.
      if (score.weight > max_weight) {
        score.value = score.value / 2;
        break;
      }
    }
  }

  if (score.weight == 0) {
      return score;
  }

  // Prioritize the value, but give a bump for a better value-to-weight ratio.
  //return value_total + (uint64_t) (value_total / weight_total);
  return score;
}

__device__ uint32_t
get_weight_kernel(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
{
    uint32_t weight = 0;
    for (uint32_t i = 0; i < len; i++) {
        if (cell[i] == 1) {
            weight += knapsack[i].weight;
        }
    }

    return weight;
}

__device__ uint32_t
get_value_kernel(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
{
    uint32_t value = 0;
    for (uint32_t i = 0; i < len; i++) {
        if (cell[i] == 1) {
            value += knapsack[i].value;
        }
    }

    return value;
}

uint32_t
get_weight(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
{
    uint32_t weight = 0;
    for (uint32_t i = 0; i < len; i++) {
        if (cell[i] == 1) {
            weight += knapsack[i].weight;
        }
    }

    return weight;
}

uint32_t
get_value(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
{
    uint32_t value = 0;
    for (uint32_t i = 0; i < len; i++) {
        if (cell[i] == 1) {
            value += knapsack[i].value;
        }
    }

    return value;
}

__device__ void
breed(uint32_t *child1, uint32_t *child2, uint32_t *parent1, uint32_t *parent2, uint32_t half)
{
  for (uint32_t i = 0; i < half; i++) {
    child1[i] = parent1[i];
    child1[i + half] = parent2[half - i];
    child2[i] = parent1[half - i];
    child2[i + half] = parent2[i];
  }
}

__device__ void
mutate(uint32_t *cell, uint32_t len, curandState_t *state)
{
    cell[curand(state) % len] = !cell[curand(state) % len];
}

__global__ void
cycle(knapsack_node_t *knapsack, parameters_t settings,
           uint32_t *dMigration, uint32_t *dResults, uint32_t islands)
{
  int id = threadIdx.x;
  uint32_t genes = settings.item_count;
  uint32_t* pop0 = (uint32_t *) malloc(settings.pop_size * genes * sizeof(uint32_t));
  memset(pop0, 0, settings.pop_size * genes * sizeof(uint32_t));
  uint32_t* pop1 = (uint32_t *) malloc(settings.pop_size * genes * sizeof(uint32_t));
  memset(pop1, 0, settings.pop_size * genes * sizeof(uint32_t));
  uint32_t* best_cell = (uint32_t *) malloc(genes * sizeof(uint32_t));
  uint32_t* migr_indices = (uint32_t *) malloc(settings.migration_size * sizeof(uint32_t));

  uint32_t* pop = pop0;
  uint32_t* new_pop = pop1;
  uint8_t cycle_type = 0;
  uint32_t cycle = 0;
  curandState_t state;
  curand_init(100 + id, 0, 0, &state);
  gen_population_nums(pop, settings.pop_size, genes, &state);
  //gen_population_bits(pop, pop_size, genes, genes / 4);
  uint32_t *ranking = (uint32_t *) malloc(settings.pop_size * sizeof(uint32_t));
  uint32_t rank_sum = 0, max = 0;
  knapsack_node_t score = {0, 0};
  uint32_t parent1_score, parent2_score;
  int64_t parent1_idx, parent2_idx;

  for (uint32_t c = 0; c < settings.cycles; c++) {
    for (uint32_t x = 0; x < settings.pop_size; x +=2) {
      // Note: parents are chosen WITH replacement
      parent1_score = 0;
      parent2_score = 0;
      parent1_idx = -1;
      parent2_idx = -1;
      rank_sum = 0;

      for (uint32_t i = 0; i < settings.pop_size; i++) {
        score = test_fitness_kernel(&pop[i * genes], genes, knapsack, settings.max_weight);

        // Not all generations may be meet the constraints or have the highest
        // value, so store the best-scoring cell seen so far.
        if (score.weight < settings.max_weight && score.value > max) {
            max = score.value;
            memcpy(best_cell, &pop[i * genes], genes * sizeof(uint32_t));
        }

        ranking[i] = score.value;
        rank_sum += score.value;
      }

      if (rank_sum == 0) {
          return;
      }

      while (parent1_score == 0) {
        parent1_score = curand(&state) % rank_sum;
      }
      while (parent2_score == 0 || parent2_score == parent1_score) {
        parent2_score = curand(&state) % rank_sum;
      }
      rank_sum = 0;

      for (int64_t i = 0; i < settings.pop_size; i++) {
        rank_sum += ranking[i];
        if (parent1_score < rank_sum && parent1_idx == -1) {
          parent1_idx = i;
        }

        if (parent2_score < rank_sum && parent2_idx == -1) {
          parent2_idx = i;
        }

        if (parent1_idx != -1 && parent2_idx != -1) {
          break;
        }
      }

      breed(&new_pop[x * genes], &new_pop[(x + 1) * genes],
            &pop[parent1_idx * genes], &pop[parent2_idx * genes],
            genes / 2); // (2 * 8 * sizeof(uint32_t))

        if (curand(&state) % settings.mutation_prob == 0) {
          mutate(&new_pop[x * genes], genes, &state);
        }
        
        if (curand(&state) % settings.mutation_prob == 0) {
          mutate(&new_pop[(x + 1) * genes], genes, &state);
        }

    }

    if (cycle % settings.migration_freq == 0) {
      for (uint32_t i = 0; i < settings.migration_size; i++) {
        migr_indices[i] = curand(&state) % settings.pop_size;
        memcpy(&dMigration[i * ((id + 1) % islands) * genes + i],
               &new_pop[i * genes + migr_indices[i]],
               genes * sizeof(uint32_t));
      }
      __syncthreads();
      
      for (uint32_t i = 0; i < settings.migration_size; i++) {
        memcpy(&new_pop[i * genes + migr_indices[i]],
               &dMigration[i * id * genes + i],
               genes * sizeof(uint32_t));
      }
    }

    if (cycle_type == 0) {
      pop = pop1;
      new_pop = pop0;
      cycle_type = 1;
    } else {
      pop = pop0;
      new_pop = pop1;
      cycle_type = 0;
    }

    cycle++;
  }

  // Check the last generation.
  for (uint32_t i = 0; i < settings.pop_size; i++) {
      score = test_fitness_kernel(&pop[i * genes], genes, knapsack, settings.max_weight);
      if (score.weight < settings.max_weight && score.value > max) {
        max = score.value;
        memcpy(best_cell, &pop[i * genes], genes * sizeof(uint32_t));
      }
  }

  memcpy(&dResults[id * genes], best_cell, genes * sizeof(uint32_t));

  free(pop0);
  free(pop1);
  free(ranking);
  free(best_cell);
  free(migr_indices);
}

int main(int argc, char **argv)
{
  //cudaStream_t stream;
  parameters_t settings = {
    16, // pop_size
    16, // item_count
    100, // max_value
    100, // max_weight
    100, // cycles
    1, // mutation_prob
    16 / 10, // migration_size,
    100 / 10, // migration_freq,
  };
  uint32_t genes = settings.item_count;
  knapsack_node_t *dNodes;
  uint32_t islands = 10;
  uint32_t max = 0;
  knapsack_node_t score = {0, 0};
  uint32_t* best_cell = (uint32_t *) malloc(genes * sizeof(uint32_t));
  uint32_t *results = (uint32_t *) calloc(islands * genes, sizeof(uint32_t));
  uint32_t *dResults;
  uint32_t *dMigration;

  knapsack_node_t *nodes = (knapsack_node_t*) calloc(genes, sizeof(knapsack_node_t));
  gen_nodes(nodes, genes, settings.max_value, settings.max_weight);

  CHECK(cudaMalloc((void **)&dNodes, sizeof(knapsack_node_t) * genes));
  CHECK(cudaMemcpy(dNodes, nodes, genes * sizeof(knapsack_node_t), cudaMemcpyHostToDevice));

  CHECK(cudaMalloc((void **)&dResults, islands * genes * sizeof(uint32_t)));
  CHECK(cudaMalloc((void **)&dMigration, islands * genes * settings.migration_size * sizeof(uint32_t)));

  //CHECK(cudaStreamCreate(&stream));

  cycle<<<1, islands>>>(dNodes, settings, dMigration, dResults, islands);

  //CHECK(cudaStreamSynchronize(stream));
  cudaDeviceSynchronize();

  CHECK(cudaMemcpy(results, dResults, genes * islands * sizeof(uint32_t), cudaMemcpyDeviceToHost));

  for (uint16_t i = 0; i < islands; i++) {
      score = test_fitness(&results[i * genes], genes, nodes, settings.max_weight);
      if (score.weight < settings.max_weight && score.value > max) {
          max = score.value;
          memcpy(best_cell, &results[i * genes], genes * sizeof(uint32_t));
      }
  }

  if (max == 0) {
    printf("Nothing found.\n");
  } else {
    print_cell_nodes(best_cell, genes, nodes);
  }

  free(best_cell);
  free(nodes);
  CHECK(cudaFree(dNodes));
  free(results);
  CHECK(cudaFree(dResults));
  return 0;
}
