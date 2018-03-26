#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>


// Keep[N][W]
int **keep;
// V[n + 1][W + 1]
int **V;


typedef struct knapnode
{
    int value;
    int weight;
} knapsack_node;

void generate_nodes(int N, int W, knapsack_node *nodes)
{
    srand(time(NULL));

    for (int i = 0; i < N; i++)
    {
        nodes[i].value = (rand() % (N - 1)) + 1;
        nodes[i].weight = rand() % W;
        printf("Node: v:%d, w:%d\n", nodes[i].value, nodes[i].weight);
    }
}

/**
 *
 * @param node The knapsack potential entries
 * @param n The total number of files
 * @param W The total size possible
 * @return The optimal solution
 */
int Knapsack(knapsack_node *nodes, int n, int W)
{
    
    for (int i = 0; i < W; i++)
    {
        V[0][i] = 0;
    }
    
    for (int i = 1; i < n; i++)
    {
        knapsack_node node = nodes[i];
        for (int j = 0; j < W; j++)
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
    int K = W;
    for (int i = n; i > 0; i--)
    {
        knapsack_node node = nodes[i];
        if (keep[i][K] == 1)
        {
            // output i
            printf("Good i: %d\n", i);
            K = K - node.weight;
        }
    }
    return V[n][W];
}

int* allocate_mem(int*** arr, int n, int m)
{
    *arr = (int**)malloc(n * sizeof(int*));
    int *arr_data = malloc( n * m * sizeof(int));
    for(int i=0; i<n; i++)
        (*arr)[i] = arr_data + i * m ;
    return arr_data; //free point
}

void deallocate_mem(int*** arr, int* arr_data){
    free(arr_data);
    free(*arr);
}

int main()
{
    int N = 100; // The number of things
    int W = 200; // The total possible weight

    // Generate keep and V
    int *keep_point = allocate_mem(&keep, N + 1, W + 1);
    int *V_point = allocate_mem(&V, N + 1, W + 1);

    // Generate the nodes
    knapsack_node *nodes = (knapsack_node*) calloc(N, sizeof(knapsack_node));
    generate_nodes(N, W, nodes);

    // Evaluate
    int maxVal = Knapsack(nodes, N, W);
    printf("Max values: %d\n", maxVal);

    free(nodes);
    deallocate_mem(&keep, keep_point);
    deallocate_mem(&V, V_point);

    return 0;
}