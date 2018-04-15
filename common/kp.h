#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

typedef struct knapsack_node_s {
  uint32_t value;
  uint32_t weight;
} knapsack_node_t;

typedef struct parameters_s {
  uint32_t pop_size;
  uint32_t item_count; // Must be a power of two.
  uint32_t max_value;
  uint32_t max_weight;
  uint32_t cycles;
  uint32_t mutation_prob;
  uint32_t migration_size;
  uint32_t migration_freq;
} parameters_t;

typedef struct file_header_s {
  uint32_t item_count;
  uint32_t max_value;
  uint32_t max_weight;
  uint32_t best_value;
  uint32_t best_weight;
} file_header_t;

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

void write_data(file_header_t *header, knapsack_node_t *nodes, uint32_t *best, const char *filename)
{
    FILE *data = fopen(filename, "wb");

    fwrite(header, sizeof(file_header_t), 1, data);
    fwrite(best, sizeof(uint32_t), header->item_count, data);
    fwrite(nodes, sizeof(knapsack_node_t), header->item_count, data);

    fclose(data);
}

void read_data(file_header_t *header, knapsack_node_t **nodes, uint32_t **best, const char *filename)
{
    FILE *data = fopen(filename, "rb");

    fread(header, sizeof(file_header_t), 1, data);

    *best = (uint32_t *) malloc(header->item_count * sizeof(uint32_t));
    fread(*best, sizeof(uint32_t), header->item_count, data);

    *nodes = (knapsack_node_t *) malloc(header->item_count * sizeof(knapsack_node_t));
    fread(*nodes, sizeof(knapsack_node_t), header->item_count, data);

    fclose(data);
}
