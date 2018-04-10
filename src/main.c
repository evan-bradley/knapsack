#include "kp.h"
#include <time.h>

// Keep[N][W]
uint32_t **keep;
// V[n + 1][W + 1]
uint32_t **V;

/**
 *
 * @param node The knapsack potential entries
 * @param n The total number of files
 * @param W The total size possible
 * @return The optimal solution
 */
int Knapsack(knapsack_node_t *nodes, uint32_t n, uint32_t W)
{

    for (uint32_t i = 0; i < W; i++)
    {
        V[0][i] = 0;
    }

    for (uint32_t i = 1; i <= n; i++)
    {
        knapsack_node_t node = nodes[i];
        for (uint32_t j = 0; j <= W; j++)
        {
            if ((node.weight <= j) && (node.value + V[i - 1][j - node.weight] > V[i - 1][j]))
            {
                V[i][j] = node.value + V[i - 1][j - node.weight];
                keep[i][j] = 1;
            }
            else
            {
                V[i][j] = V[i - 1][j];
                keep[i][j] = 0;
            }
        }
    }
    uint32_t K = W;
    for (int i = n; i > 0; i--)
    {
        knapsack_node_t node = nodes[i];
        if (keep[i][K] == 1)
        {
            // output i
          printf("Good i: %d; %u, %u\n", i, nodes[i].value, nodes[i].weight);
            K = K - node.weight;
        }
    }
    return V[n][W];
}

uint32_t* allocate_mem(uint32_t*** arr, uint32_t n, uint32_t m)
{
    *arr = (uint32_t**)malloc(n * sizeof(uint32_t*));
    uint32_t *arr_data = malloc( n * m * sizeof(uint32_t));
    for(uint32_t i=0; i<n; i++)
        (*arr)[i] = arr_data + i * m ;
    return arr_data; //free pouint32_t
}

void deallocate_mem(uint32_t*** arr, uint32_t* arr_data){
    free(arr_data);
    free(*arr);
}

int main()
{
    uint32_t node_count = 16;
    uint32_t max_value = 100;
    uint32_t max_weight = 100;

    // Generate keep and V
    uint32_t *keep_point = allocate_mem(&keep, node_count + 1, max_weight + 1);
    uint32_t *V_point = allocate_mem(&V, node_count + 1, max_weight + 1);

    // Generate the nodes
    knapsack_node_t *nodes = (knapsack_node_t*) calloc(node_count, sizeof(knapsack_node_t));
    gen_nodes(nodes, node_count, max_value, max_weight);
    print_nodes(nodes, node_count);

    // Evaluate
    uint32_t maxVal = Knapsack(nodes, node_count, max_weight);
    printf("Max values: %d\n", maxVal);

    free(nodes);
    deallocate_mem(&keep, keep_point);
    deallocate_mem(&V, V_point);

    return 0;
}
