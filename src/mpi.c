#include "kp.h"
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

uint64_t test_fitness(uint32_t *cell, uint32_t len, knapsack_node_t *knapsack, uint64_t max_weight)
{
  uint32_t value_total = 0;
  uint32_t weight_total = 0;

  for (uint32_t i = 0; i < len; i++) {
    if (cell[i] == 1) {
      weight_total += knapsack[i].weight;
      value_total += knapsack[i].value;

      // Note that breaking once we hit the max weight only counts the values up
      // to that point, and so might induce bias.
      if (weight_total > max_weight) {
        value_total = value_total / 2;
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

void cycle(int id, int tasklen, knapsack_node_t *knapsack, uint32_t pop_size,
           uint32_t item_count, uint32_t max_weight, uint32_t cycles)
{
  uint32_t swap_size = pop_size / 10;
  uint8_t swap_cycle = 10;

  uint32_t* pop0 = (uint32_t *) calloc(pop_size * item_count, sizeof(uint32_t));
  uint32_t* pop1 = (uint32_t *) calloc(pop_size * item_count, sizeof(uint32_t));
  uint32_t* best_cell = (uint32_t *) calloc(item_count, sizeof(uint32_t));
  uint32_t* send_buf = (uint32_t *) calloc(swap_size * item_count, sizeof(uint32_t));
  uint32_t* recv_buf = (uint32_t *) calloc(swap_size * item_count, sizeof(uint32_t));
  uint32_t* swap_indices = (uint32_t *) calloc(swap_size, sizeof(uint32_t));
  uint32_t* results = (uint32_t *) calloc(tasklen * item_count, sizeof(uint32_t));

  MPI_Request r_send, r_recv;
  uint32_t* pop = pop0;
  uint32_t* new_pop = pop1;
  uint8_t cycle_type = 0, cycle = 0;
  //gen_population_nums(pop, pop_size, item_count, id);
  gen_population_bits(pop, pop_size, item_count, max_weight / item_count, id);
  //gen_population_bits(pop, pop_size, item_count, item_count / 4);
  uint32_t *ranking = calloc(pop_size, sizeof(uint32_t));
  uint32_t rank_sum = 0, max = 0, score = 0;
  uint32_t parent1_score, parent2_score;
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
        score = test_fitness(&pop[i * item_count], item_count, knapsack, max_weight);

        ranking[i] = score;
        rank_sum += score;
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

      breed(&new_pop[x * item_count], &new_pop[(x + 1) * item_count],
            &pop[parent1_idx * item_count], &pop[parent2_idx * item_count],
            pop_size / 2); // (2 * 8 * sizeof(uint32_t))

        mutate(&new_pop[x * item_count], item_count);
        mutate(&new_pop[(x + 1) * item_count], item_count);

    }

    // The population will not necessarily converge to the optimal solution,
    // so record the best seen cell out of all previously seen cells at each
    // generation.
    for (uint16_t i = 0; i < pop_size; i++) {
        uint32_t cell_weight = get_weight(&new_pop[i * item_count], item_count, knapsack);
        if (cell_weight < max_weight) {
            score = get_value(&new_pop[i * item_count], item_count, knapsack);

            if (score > max) {
                max = score;
                memcpy(best_cell, &new_pop[i * item_count], item_count * sizeof(uint32_t));
            }
        }
    }

    // NOTE: There can be duplicate indicies here.
    for (uint32_t i = 0; i < swap_size; i++) {
      swap_indices[i] = rand() % pop_size;
      memcpy(&send_buf[i], &new_pop[i * item_count], item_count * sizeof(uint32_t));
    }

    if (cycle % swap_cycle == 0) {
      MPI_Isend(send_buf, swap_size * item_count, MPI_UNSIGNED,
                (id + 1) % tasklen, 0, MPI_COMM_WORLD, &r_send);
      MPI_Irecv(recv_buf, swap_size * item_count, MPI_UNSIGNED,
                (id - 1) % tasklen, 0, MPI_COMM_WORLD, &r_recv);
    }
    MPI_Wait(&r_send, MPI_STATUS_IGNORE);
    MPI_Wait(&r_recv, MPI_STATUS_IGNORE);

    for (uint32_t i = 0; i < swap_size; i++) {
      swap_indices[i] = rand() % pop_size;
      memcpy(&new_pop[i * item_count], &recv_buf[i], item_count * sizeof(uint32_t));
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

  MPI_Gather(best_cell, item_count, MPI_UNSIGNED,
             results,   item_count, MPI_UNSIGNED,
             0, MPI_COMM_WORLD);

  if (id == 0) {
    for (uint16_t i = 0; i < tasklen; i++) {
        uint32_t cell_weight = get_weight(&results[i * item_count], item_count, knapsack);
        if (cell_weight < max_weight) {
            score = get_value(&results[i * item_count], item_count, knapsack);

            if (score > max) {
                max = score;
                memcpy(best_cell, &results[i * item_count], item_count * sizeof(uint32_t));
            }
        }
    }

    if (max == 0) {
        printf("No solutions found\n");
    } else {
        print_cell_nodes(best_cell, item_count, knapsack);
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

  uint32_t pop_size = 16;
  uint32_t item_count = 16; // Must be a power of two.
  uint32_t max_value = 100;
  uint32_t max_weight = 100;
  uint32_t cycles = 1000;

  knapsack_node_t *nodes = (knapsack_node_t*) calloc(item_count, sizeof(knapsack_node_t));
  gen_nodes(nodes, item_count, max_value, max_weight);

  cycle(taskid, tasklen, nodes, pop_size, item_count, max_weight, cycles);

  free(nodes);

  MPI_Finalize();

  return 0;
}
