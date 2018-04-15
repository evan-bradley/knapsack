#include "../common/kp.h"
#include <time.h>

 #define max(a,b) \
      ({ __typeof__ (a) _a = (a); \
         __typeof__ (b) _b = (b); \
        _a > _b ? _a : _b; })

/**
 *
 * @param node The knapsack potential entries
 * @param n The total number of files
 * @param W The total size possible
 * @return The optimal solution
 */
knapsack_node_t Knapsack(knapsack_node_t *nodes, uint32_t **V,
                         uint32_t **keep, file_header_t opt)
{
    knapsack_node_t best = {0, 0};
    uint32_t n = opt.item_count;
    uint32_t W = opt.max_weight;

    for (uint32_t i = 0; i <= n; i++)
    {
        // NOTE: if i == 0, node is garbage; just set to satisfy the compiler.
        knapsack_node_t node = i != 0 ? nodes[i - 1] : nodes[i];
        for (uint32_t j = 0; j <= W; j++)
        {
            if (i == 0 || j == 0)
            {
              V[i][j] = 0;
            }
            else if ((node.weight <= j))
            {
                //V[i][j] = max(V[i - 1][j], V[i - 1][j - node.weight] + node.value);
                //keep[i][j] = 1;
                if (node.value + V[i - 1][j - node.weight] > V[i - 1][j]) {
                  V[i][j] = V[i - 1][j - node.weight] + node.value;
                  keep[i][j] = 1;
                } else {
                  V[i][j] = V[i - 1][j];
                }
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
        knapsack_node_t node = nodes[i - 1];
        if (keep[i][K] == 1)
        {
            best.weight += nodes[i - 1].weight;
            K = K - node.weight;
        }
    }

    best.value = V[n][W];
    return best;
}

uint32_t* allocate_mem(uint32_t*** arr, uint32_t n, uint32_t m)
{
    *arr = (uint32_t**)calloc(n, sizeof(uint32_t*));
    uint32_t *arr_data = calloc(n * m, sizeof(uint32_t));
    for(uint32_t i=0; i<n; i++)
        (*arr)[i] = arr_data + i * m ;
    return arr_data; //free point
}

void deallocate_mem(uint32_t*** arr, uint32_t* arr_data)
{
    free(arr_data);
    free(*arr);
}

int main()
{
    file_header_t options = {
      1024,  // item_count
      100, // max_value
      2048, // max_weight
      0, // best_value
      0, // best_weight
    };

    // Keep[N][W]
    uint32_t **keep;
    // V[n + 1][W + 1]
    uint32_t **V;

    // Generate keep and V
    uint32_t *keep_point = allocate_mem(&keep, options.item_count + 1, options.max_weight + 1);
    uint32_t *V_point = allocate_mem(&V, options.item_count + 1, options.max_weight + 1);
    uint32_t *best = (uint32_t *) calloc(options.item_count, sizeof(uint32_t));

    // Generate the nodes
    knapsack_node_t *nodes = (knapsack_node_t*) calloc(options.item_count, sizeof(knapsack_node_t));
    gen_nodes(nodes, options.item_count, options.max_value, options.max_weight);
    //print_nodes(nodes, options.item_count);

    // Evaluate
    knapsack_node_t maxVal = Knapsack(nodes, V, keep, options);
    printf("Maximum: value %u, weight %u\n", maxVal.value, maxVal.weight);
    options.best_value = maxVal.value;
    options.best_weight = maxVal.weight;

    uint32_t K = options.max_weight;
    for (int i = options.item_count; i > 0; i--)
    {
        knapsack_node_t node = nodes[i - 1];
        best[i - 1] = keep[i][K];
        if (keep[i][K] == 1)
        {
            K = K - node.weight;
        }
    }

    write_data(&options, nodes, best, "test.kp");

    free(nodes);
    deallocate_mem(&keep, keep_point);
    deallocate_mem(&V, V_point);

    return 0;
}
