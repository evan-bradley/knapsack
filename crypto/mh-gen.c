#include "mh.c"
#include <stdlib.h>
#include <stdio.h>
#include "../common/kp.h"
#include <time.h>

// Makes x a superincreasing set, randomly generated.
void gen_superinc(int x[], int n) {
  
  x[0] = rand()%1000;
  
  for (int i=1; i<n; i++) {
    x[i] = 2*x[i-1] + (rand()%x[i-1]) + 1;
  }
}

int sum_set(int x[], int n) {
  int total = 0;
  for (int i=0; i<n; i++)
    total += x[i];
  return total;
}

int isprime(int v) {
    int i;
    if (v < 0) v = -v;
    if (v < 2 || !((unsigned)v & 1))
      return 0;
    if (v == 2) return 1;
    for (i = 3; i * i <= v; i+=2)
      if (v % i == 0)
	return 0;
    return 1;
}

int find_prime(int a) {
  int done = 0;
  int result = a;
  while (!done) {
    result++;
    if (isprime(result))
      done = 1;
  }
  return result;
}

void make_nodes(knapsack_node_t *knapsack, int n, int x[]) {
  for (int i=0; i<n; i++) {
    knapsack[i].value = x[i];
    knapsack[i].weight = x[i];
  }
}

// Does the work to spit public and private keys into a file
// writes, in order
// header
// message (1010101010101010)
// nodes representing public set
// encrypted message (subset sum on public set)
void make_keyinfo(file_header_t *header, const char *filename, int n) {
  srand(time(NULL));
  int *x = malloc(n*sizeof(int));
  gen_superinc(x, n);
  int total = sum_set(x, n);
  int prime = find_prime(total);
  int b = (rand()%prime) + 3;
  int *pub_set = malloc(N*sizeof(int));
  memcpy(x, pub_set, n*sizeof(int));
  sis_to_pub(pub_set, n, prime, b);
  uint32_t *msg = (uint32_t *) calloc(n, sizeof(uint32_t));
  for (int i=0; i<n; i++)
    msg[i] = i%2;
  knapsack_node_t *nodes = (knapsack_node_t*) calloc(nm sizeof(knapsack_node_t));
  make_nodes(nodes, n, pub_set);
  int encrypted = encrypt(pub_set, msg, n);
  
  FILE *data = fopen(filename, "wb");
  fwrite(header, sizeof(file_header_t), 1, data);
  fwrite(msg, sizeof(uint32_t), n, data);
  fwrite(nodes, sizeof(knapsack_node_t), n, data);
  fwrite(encrypted, sizeof(int), 1, data);
  fclose(data);
}
