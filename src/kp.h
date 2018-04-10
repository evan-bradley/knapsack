#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

typedef struct knapnode {
  uint32_t value;
  uint32_t weight;
} knapsack_node_t;

void gen_nodes(knapsack_node_t *knapsack, uint32_t nodes, uint32_t max_value, uint32_t max_weight)
{
  srand(100);
  for (uint32_t i = 0; i < nodes; i++) {
    knapsack[i].value = rand() % (max_value - 1) + 1;
    knapsack[i].weight = rand() % (max_weight - 1) + 1;
  }
}

void print_nodes(knapsack_node_t *knapsack, uint32_t nodes)
{
  srand(100);
  for (uint32_t i = 0; i < nodes; i++) {
    printf("{ value: %u\tweight: %u }\n", knapsack[i].value, knapsack[i].weight);
  }
}
