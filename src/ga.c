#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define get_bit(arr, i) (arr[i / 32] & (1 << (i % 32)))
#define set_bit(arr, i) (arr[i / 32] |= (1 << (i % 32)))
#define unset_bit(arr, i) (arr[i / 32] &= ~(1 << (i % 32)))

typedef struct knapnode {
  uint32_t value;
  uint32_t weight;
} knapsack_node_t;

void gen_nodes(knapsack_node_t *knapsack, uint32_t nodes, uint32_t max_weight)
{
  for (uint32_t i = 0; i < nodes; i++) {
    knapsack[i].value = rand() % (nodes - 1) + 1;
    knapsack[i].weight = rand() % (max_weight - 1) + 1;
    }
}

void gen_population_nums(uint32_t **pop, uint32_t len, uint32_t genes)
{
  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < genes; j++) {
      pop[i][j] = rand() % (genes - 1);
    }
  }
}

void gen_population_bits(uint32_t **pop, uint32_t len, uint32_t genes, uint32_t bits)
{
  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < bits; j++) {
      set_bit(pop[i], rand() % (genes - 1));
    }
  }
}

uint64_t test_fitness(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack, uint64_t max_weight)
{
  uint64_t value_total = 0;
  uint64_t weight_total = 0;

  for (uint32_t i = 0; i < len; i++) {
    if (get_bit(cell, i) == 1) {
      weight_total += knapsack[i].weight;
      value_total += knapsack[i].value;

      // Note that breaking once we hit the max weight only counts the values up
      // to that point, and so might induce bias.
      if (weight_total > max_weight) {
        value_total -= knapsack[i].value;
        break;
      }
    }
  }

  if (weight_total == 0) {
      return 0;
  }

  // Prioritize the value, but give a bump for a better value-to-weight ratio.
  return value_total + (uint64_t) (value_total / weight_total);
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

void copy_arr(uint32_t **arr0, uint32_t **arr1, uint32_t m, uint32_t n)
{
    for (uint32_t i = 0; i < m; i++) {
        for (uint32_t j = 0; j < n; j++) {
            arr1[i][j] = arr0[i][j];
        }
    }
}

void cycle(knapsack_node_t *knapsack, uint32_t pop_size,
           uint32_t item_count, uint32_t max_weight, uint32_t cycles)
{
  uint32_t** pop0 = (uint32_t **) malloc(pop_size * sizeof(uint32_t *));
  for (uint32_t i = 0; i < pop_size; i++) {
    pop0[i] = (uint32_t *) calloc(item_count / 32, sizeof(uint32_t));
  }

  uint32_t** pop1 = (uint32_t **) malloc(pop_size * sizeof(uint32_t *));
  for (uint32_t i = 0; i < pop_size; i++) {
    pop1[i] = (uint32_t *) calloc(item_count / 32, sizeof(uint32_t));
  }

  uint32_t** pop = pop0;
  uint32_t** new_pop = pop1;
  uint8_t cycle = 0;
  gen_population_nums(pop, pop_size, item_count);
  //gen_population_bits(pop, pop_size, item_count, item_count / 4);
  uint64_t *ranking = calloc(pop_size, sizeof(uint64_t));
  uint64_t rank_sum, max = 0, score;
  uint64_t parent1_score, parent2_score;
  int64_t parent1_idx, parent2_idx;

  for (uint32_t c = 0; c < cycles; c++) {
    for (uint32_t x = 0; x < pop_size; x +=2) {
      // Note: parents are chosen WITH replacement
      parent1_score = 0;
      parent2_score = 0;
      parent1_idx = -1;
      parent2_idx = -1;
      rank_sum = 0;

      for (uint32_t i = 0; i < pop_size; i++) {
        score = test_fitness(pop[i], item_count, knapsack, max_weight);
        ranking[i] = score;
        rank_sum += score;
      }

      while (parent1_score == 0) {
      	parent1_score = rand() % rank_sum;
      }
      while (parent2_score == 0) {
      	parent2_score = rand() % rank_sum;
      }
      rank_sum = 0;

      for (int64_t i = 0; i < pop_size; i++) {
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

      breed(new_pop[x], new_pop[x + 1], pop[parent1_idx], pop[parent2_idx], pop_size / 64);
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

  for (uint16_t i = 0; i < pop_size; i++) {
    score = test_fitness(pop[i], item_count, knapsack, max_weight);
    max = score > max ? score : max;
  }

  printf("max weight: %llu\n", max);

  for (uint32_t i = 0; i < pop_size; i++) {
    free(pop[i]);
  }
  free(pop);

  for (uint32_t i = 0; i < pop_size; i++) {
    free(new_pop[i]);
  }
  free(new_pop);

  free(ranking);
}

int main(int argc, char **argv)
{
  srand(100);
  uint32_t pop_size = 100;
  uint32_t item_count = 1024; // Must be a power of two.
  uint32_t max_weight = 10000;
  uint32_t cycles = 100;

  knapsack_node_t *nodes = (knapsack_node_t*) calloc(item_count, sizeof(knapsack_node_t));
  gen_nodes(nodes, item_count, max_weight);

  cycle(nodes, pop_size, item_count, max_weight, cycles);

  free(nodes);
  return 0;
}
