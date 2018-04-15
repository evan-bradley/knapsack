#include "../common/kp.h"
#include "../common/common.h"
#include <time.h>
#include <cuda.h>

 #define max(a,b) \
      ({ __typeof__ (a) _a = (a); \
         __typeof__ (b) _b = (b); \
        _a > _b ? _a : _b; })

__global__ void
solve_weight(uint32_t *V, uint32_t *keep, uint32_t i,
            knapsack_node_t *nodes, file_header_t *opt)
{
    uint32_t j = threadIdx.x;
    uint32_t W = opt->max_weight;

    if (j > W) {
        return;
    }

    // NOTE: if i == 0, node is garbage; just set to satisfy the compiler.
    knapsack_node_t node = i != 0 ? nodes[i - 1] : nodes[i];

    if (i == 0 || j == 0) {
        V[i + W * j] = 0;
    }
    else if ((node.weight <= j)) {
        if (V[(i - 1) + W * (j - node.weight)] + node.value > V[(i - 1) + W * j]) {
            V[i + W * j] = V[(i - 1) + W * (j - node.weight)] + node.value;
            keep[i + W * j] = 1;
        }
        else {
            V[i + W * j] = V[(i - 1) + W * j];
        }
    }
    else {
        V[i + W * j] = V[(i - 1) + W * j];
        keep[i + W * j] = 0;
    }
}

/**
 *
 * @param node The knapsack potential entries
 * @param n The total number of files
 * @param W The total size possible
 * @return The optimal solution
 */
__global__ void
kp01(knapsack_node_t *nodes, knapsack_node_t *best,
     uint32_t *V, uint32_t *keep, file_header_t *opt)
{
    uint32_t n = opt->item_count;
    uint32_t W = opt->max_weight;

    dim3 block_dim(1024, 1, 1);
    dim3 grid_dim(ceil((float)W / (float)block_dim.x), 1, 1);
    printf("grid: %u, block: %u\n", grid_dim.x, block_dim.x);

    best->value = 0;
    best->weight = 0;

    for (uint32_t i = 1; i <= n; i++) {
        solve_weight<<<grid_dim, block_dim>>>(V, keep, i, nodes, opt);
        cudaDeviceSynchronize();
    }

    uint32_t K = W;
    for (int i = n; i > 0; i--) {
        knapsack_node_t node = nodes[i - 1];
        if (keep[i + W * K] == 1) {
            best->weight += nodes[i - 1].weight;
            K = K - node.weight;
        }
    }

    best->value = V[n + W * W];
}

int
main(int argc, char **argv)
{
    file_header_t options = {
      1024,  // item_count
      100, // max_value
      2048, // max_weight
      0, // best_value
      0, // best_weight
    };

    uint32_t n = options.item_count;
    uint32_t W = options.max_weight;

    // The first rows are filled with 0s for simplicity, so the arrays
    // are N+1 x W+1
    // Keep[N + 1][W + 1]
    uint32_t *keep = (uint32_t *) calloc((n + 1) * (W + 1), sizeof(uint32_t));
    // V[n + 1][W + 1]
    uint32_t *V = (uint32_t *) calloc((n + 1) * (W + 1), sizeof(uint32_t));
    uint32_t *best = (uint32_t *) calloc(n, sizeof(uint32_t));
    knapsack_node_t kp_max;

    // Generate the nodes
    knapsack_node_t *nodes = (knapsack_node_t*) calloc(n, sizeof(knapsack_node_t));
    gen_nodes(nodes, n, options.max_value, W);
    //print_nodes(nodes, n);

    uint32_t *dKeep, *dV;
    knapsack_node_t *dNodes, *dKP_max;
    file_header_t *dOpt;

    CHECK(cudaMalloc((void **)&dKeep, (n + 1) * (W + 1) * sizeof(uint32_t)));
    CHECK(cudaMalloc((void **)&dV, (n + 1) * (W + 1) * sizeof(uint32_t)));
    CHECK(cudaMalloc((void **)&dKP_max, sizeof(knapsack_node_t)));
    CHECK(cudaMalloc((void **)&dNodes, n * sizeof(knapsack_node_t)));
    CHECK(cudaMalloc((void **)&dOpt, sizeof(file_header_t)));

    CHECK(cudaMemset(dKeep, 0, (n + 1) * (W + 1) * sizeof(uint32_t)));
    CHECK(cudaMemset(dV, 0, (n + 1) * (W + 1) * sizeof(uint32_t)));
    CHECK(cudaMemset(dKP_max, 0, sizeof(knapsack_node_t)));
    CHECK(cudaMemcpy(dNodes, nodes, n * sizeof(knapsack_node_t), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dOpt, &options, sizeof(file_header_t), cudaMemcpyHostToDevice));

    kp01<<<1, 1>>>(dNodes, dKP_max, dV, dKeep, dOpt);
    CHECK(cudaDeviceSynchronize());

    CHECK(cudaMemcpy(&kp_max, dKP_max, sizeof(knapsack_node_t), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(keep, dKeep, (n + 1) * (W + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));

    options.best_value = kp_max.value;
    options.best_weight = kp_max.weight;
    printf("Value: %u\tWeight: %u\n", kp_max.value, kp_max.weight);

    uint32_t K = W;
    for (int i = n; i > 0; i--) {
        knapsack_node_t node = nodes[i - 1];
        best[i - 1] = keep[i + W * K];
        if (keep[i + W * K] == 1) {
            K = K - node.weight;
        }
    }

    write_data(&options, nodes, best, "test.kp");

    CHECK(cudaFree(dKeep));
    CHECK(cudaFree(dV));
    CHECK(cudaFree(dKP_max));
    CHECK(cudaFree(dNodes));
    CHECK(cudaFree(dOpt));

    free(best);
    free(nodes);
    free(keep);
    free(V);

    return 0;
}
