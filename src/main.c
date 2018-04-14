#include "kp.h"
#include <time.h>

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

    for (uint32_t i = 0; i < opt.max_weight; i++)
    {
        V[0][i] = 0;
    }

    for (uint32_t i = 1; i <= opt.item_count; i++)
    {
        knapsack_node_t node = nodes[i];
        for (uint32_t j = 0; j <= opt.max_weight; j++)
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
    uint32_t K = opt.max_weight;
    for (int i = opt.item_count; i > 0; i--)
    {
        knapsack_node_t node = nodes[i];
        if (keep[i][K] == 1)
        {
            best.weight += nodes[i].weight;
            K = K - node.weight;
        }
    }

    best.value = V[opt.item_count][opt.max_weight];
    return best;
}

uint32_t* allocate_mem(uint32_t*** arr, uint32_t n, uint32_t m)
{
    *arr = (uint32_t**)malloc(n * sizeof(uint32_t*));
    uint32_t *arr_data = malloc( n * m * sizeof(uint32_t));
    for(uint32_t i=0; i<n; i++)
        (*arr)[i] = arr_data + i * m ;
    return arr_data; //free pouint32_t
}

void deallocate_mem(uint32_t*** arr, uint32_t* arr_data)
{
    free(arr_data);
    free(*arr);
}

int main()
{
    file_header_t options = {
      16,  // item_count
      100, // max_value
      100, // max_weight
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

    for (uint32_t i = 0; i < options.item_count; i++) {
        best[i] = keep[i][options.max_weight];
    }

    write_data(&options, nodes, best, "test.kp");

    free(nodes);
    deallocate_mem(&keep, keep_point);
    deallocate_mem(&V, V_point);

    return 0;
}
