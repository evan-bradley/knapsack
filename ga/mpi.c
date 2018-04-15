#include "../common/kp.h"
#include <memory.h>
#include <mpi.h>

#define get_bit(arr, i) (arr[i / 32] & (1 << (i % 32)))
#define set_bit(arr, i) (arr[i / 32] |= (1 << (i % 32)))
#define unset_bit(arr, i) (arr[i / 32] &= ~(1 << (i % 32)))

// Generate roughly 50/50 distribution of enabled and disabled gene-bits.
void gen_population_nums(uint32_t *pop, uint32_t len, uint32_t genes, int seed_adjust)
{
  srand(100 + seed_adjust);
  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < genes; j++) {
      pop[i * genes + j] = rand() % 2;
    }
  }
}

// Generate genes with a specified number of bits enabled.
void gen_population_bits(uint32_t *pop, uint32_t len, uint32_t genes, uint32_t bits, int seed_adjust)
{
  srand(100 + seed_adjust);
  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < bits; j++) {
      pop[i * genes + (rand() % genes)] = 1;
    }
  }
}

void print_cell(uint32_t *cell, uint32_t len)
{
    for (uint32_t i = 0; i < len; i++) {
        printf("%u", cell[i]);
    }
    printf("\n");
}

void print_cell_nodes(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
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

void print_population(uint32_t *pop, uint32_t len, uint32_t genes)
{
  for (uint32_t i = 0; i < len; i++) {
      print_cell(&pop[i * genes], genes);
  }
}

knapsack_node_t test_fitness(uint32_t *cell, uint32_t len,
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

uint32_t get_weight(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
{
    uint32_t weight = 0;
    for (uint32_t i = 0; i < len; i++) {
        if (cell[i] == 1) {
            weight += knapsack[i].weight;
        }
    }

    return weight;
}

uint32_t get_value(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack)
{
    uint32_t value = 0;
    for (uint32_t i = 0; i < len; i++) {
        if (cell[i] == 1) {
            value += knapsack[i].value;
        }
    }

    return value;
}

void breed(uint32_t *child1, uint32_t *child2, uint32_t *parent1, uint32_t *parent2, uint32_t half)
{
  for (uint32_t i = 0; i < half; i++) {
    child1[i] = parent1[i];
    child1[i + half] = parent2[half - i];
    child2[i] = parent1[half - i];
    child2[i + half] = parent2[i];
  }
}

void mutate(uint32_t *cell, uint32_t len)
{
    cell[rand() % len] = !cell[rand() % len];
}

void copy_arr(uint32_t **arr0, uint32_t **arr1, uint32_t m, uint32_t n)
{
    for (uint32_t i = 0; i < m; i++) {
        for (uint32_t j = 0; j < n; j++) {
            arr1[i][j] = arr0[i][j];
        }
    }
}

void cycle(int id, int tasklen,
           knapsack_node_t *knapsack, parameters_t settings)
{
  uint32_t swap_size = settings.pop_size / 10;
  uint8_t swap_cycle = 10;
  uint32_t genes = settings.item_count;

  uint32_t* pop0 = (uint32_t *) calloc(settings.pop_size * genes, sizeof(uint32_t));
  uint32_t* pop1 = (uint32_t *) calloc(settings.pop_size * genes, sizeof(uint32_t));
  uint32_t* best_cell = (uint32_t *) calloc(genes, sizeof(uint32_t));
  uint32_t* send_buf = (uint32_t *) calloc(swap_size * genes, sizeof(uint32_t));
  uint32_t* recv_buf = (uint32_t *) calloc(swap_size * genes, sizeof(uint32_t));
  uint32_t* swap_indices = (uint32_t *) calloc(swap_size, sizeof(uint32_t));
  uint32_t* results = (uint32_t *) calloc(tasklen * genes, sizeof(uint32_t));

  MPI_Request r_send, r_recv;
  uint32_t* pop = pop0;
  uint32_t* new_pop = pop1;
  uint8_t cycle_type = 0, cycle = 0;
  //gen_population_nums(pop, pop_size, genes, id);
  gen_population_bits(pop, settings.pop_size, genes, settings.max_weight / genes, id);
  //gen_population_bits(pop, pop_size, genes, genes / 4);
  uint32_t *ranking = calloc(settings.pop_size, sizeof(uint32_t));
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
        score = test_fitness(&pop[i * genes], genes, knapsack, settings.max_weight);

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
        parent1_score = rand() % rank_sum;
      }
      while (parent2_score == 0 || parent2_score == parent1_score) {
        parent2_score = rand() % rank_sum;
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

        if (rand() % settings.mutation_prob == 0) {
          mutate(&new_pop[x * genes], genes);
        }

        if (rand() % settings.mutation_prob == 0) {
          mutate(&new_pop[(x + 1) * genes], genes);
        }

    }

    // NOTE: There can be duplicate indicies here.
    for (uint32_t i = 0; i < swap_size; i++) {
      swap_indices[i] = rand() % settings.pop_size;
      memcpy(&send_buf[i], &new_pop[i * genes + swap_indices[i]], genes * sizeof(uint32_t));
    }

    if (cycle % swap_cycle == 0) {
      MPI_Isend(send_buf, swap_size * genes, MPI_UNSIGNED,
                (id + 1) % tasklen, 0, MPI_COMM_WORLD, &r_send);
      MPI_Irecv(recv_buf, swap_size * genes, MPI_UNSIGNED,
                (id - 1) % tasklen, 0, MPI_COMM_WORLD, &r_recv);
    }
    MPI_Wait(&r_send, MPI_STATUS_IGNORE);
    MPI_Wait(&r_recv, MPI_STATUS_IGNORE);

    for (uint32_t i = 0; i < swap_size; i++) {
      memcpy(&new_pop[i * genes + swap_indices[i]], &recv_buf[i], genes * sizeof(uint32_t));
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
  }

  // Check the last generation.
  for (uint32_t i = 0; i < settings.pop_size; i++) {
      score = test_fitness(&pop[i * genes], genes, knapsack, settings.max_weight);
      if (score.weight < settings.max_weight && score.value > max) {
        max = score.value;
        memcpy(best_cell, &pop[i * genes], genes * sizeof(uint32_t));
      }
  }

  MPI_Gather(best_cell, genes, MPI_UNSIGNED,
             results,   genes, MPI_UNSIGNED,
             0, MPI_COMM_WORLD);

  if (id == 0) {
      for (uint32_t i = 0; i < settings.pop_size; i++) {
          score = test_fitness(&results[i * genes], genes, knapsack, settings.max_weight);
          if (score.weight < settings.max_weight && score.value > max) {
            max = score.value;
            memcpy(best_cell, &results[i * genes], genes * sizeof(uint32_t));
          }
      }


    if (max == 0) {
        printf("No solutions found\n");
    } else {
        print_cell_nodes(best_cell, genes, knapsack);
    }
  }

  free(pop0);
  free(pop1);
  free(send_buf);
  free(recv_buf);
  free(swap_indices);
  free(results);
  free(ranking);
  free(best_cell);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int taskid, tasklen;
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Comm_size(MPI_COMM_WORLD, &tasklen);

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

  knapsack_node_t *nodes = (knapsack_node_t*) calloc(settings.item_count, sizeof(knapsack_node_t));
  gen_nodes(nodes, settings.item_count, settings.max_value, settings.max_weight);

  cycle(taskid, tasklen, nodes, settings);

  free(nodes);

  MPI_Finalize();

  return 0;
}
