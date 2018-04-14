#include "kp.h"
#include <memory.h>

#define get_bit(arr, i) (arr[i / 32] & (1 << (i % 32)))
#define set_bit(arr, i) (arr[i / 32] |= (1 << (i % 32)))
#define unset_bit(arr, i) (arr[i / 32] &= ~(1 << (i % 32)))

// Generate roughly 50/50 distribution of enabled and disabled gene-bits.
void gen_population_nums(uint32_t *pop, uint32_t len, uint32_t genes)
{
  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < genes; j++) {
      pop[i * genes + j] = rand() % 2;
    }
  }
}

// Generate genes with a specified number of bits enabled.
void gen_population_bits(uint32_t *pop, uint32_t len, uint32_t genes,
                         uint32_t bits, int seed_adjust)
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
      print_cell(&pop[i], genes);
  }
}

uint32_t cell_dist(uint32_t *cell0, uint32_t *cell1, uint32_t genes) {
  uint32_t sum = 0;

  for (uint32_t i = 0; i < genes; i++) {
    sum += abs(cell0[i] - cell1[i]);
  }

  return floor(sqrt(sum));
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

void breed(uint32_t *child1, uint32_t *child2,
           uint32_t *parent1, uint32_t *parent2,
           uint32_t half)
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

void cycle(knapsack_node_t *knapsack, parameters_t settings)
{
  uint32_t* pop0 = (uint32_t *) calloc(settings.pop_size * settings.item_count,
                                       sizeof(uint32_t));
  uint32_t* pop1 = (uint32_t *) calloc(settings.pop_size * settings.item_count,
                                       sizeof(uint32_t));
  uint32_t* best_cell = (uint32_t *) calloc(settings.item_count, sizeof(uint32_t));

  uint32_t* pop = pop0;
  uint32_t* new_pop = pop1;
  uint8_t cycle = 0;
  //gen_population_nums(pop, pop_size, item_count);
  gen_population_bits(pop, settings.pop_size, settings.item_count, settings.item_count / 2, 0);
  uint32_t *ranking = calloc(settings.pop_size, sizeof(uint32_t));
  uint32_t rank_sum = 0, max = 0;
  knapsack_node_t score = {0, 0};
  uint32_t parent1_score, parent2_score;
  int64_t parent1_idx, parent2_idx;
  uint32_t genes = settings.item_count;

  for (uint32_t c = 0; c < settings.cycles; c++) {
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

    for (uint32_t x = 0; x < settings.pop_size; x +=2) {
      // Note: parents are chosen WITH replacement
      parent1_score = 0;
      parent2_score = 0;
      parent1_idx = -1;
      parent2_idx = -1;

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

    if (cycle == 0) {
      pop = pop1;
      new_pop = pop0;
      cycle = 1;
    } else {
      pop = pop0;
      new_pop = pop1;
      cycle = 0;
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

  if (max == 0) {
      printf("No solutions found\n");
  } else {
      print_cell_nodes(best_cell, genes, knapsack);
  }

  free(pop0);
  free(pop1);
  free(ranking);
  free(best_cell);
}

int main(int argc, char **argv)
{
  parameters_t settings = {
    16, // pop_size
    16, // item_count
    100, // max_value
    100, // max_weight
    100, // cycles
    1, // mutation_prob
  };

  knapsack_node_t *nodes = (knapsack_node_t*) calloc(settings.item_count, sizeof(knapsack_node_t));
  gen_nodes(nodes, settings.item_count, settings.max_value, settings.max_weight);
  //print_nodes(nodes, item_count);

  cycle(nodes, settings);

  free(nodes);
  return 0;
}
