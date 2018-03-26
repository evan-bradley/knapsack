#include <stdio.h>



int keep[100][100];
// V[n + 1][W + 1]
int V[100][100];


int Knapsack(int *v, int*w, int n, int W)
{
    
    for (int i = 0; i < W; i++)
    {
        V[0][i] = 0;
    }
    
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < W; j++)
        {
            if ((w[i] < j) && (v[i] + V[i - 1][j - w[i]] > V[i - 1][j]))
            {
                V[i][j] = v[i] + V[i - 1][j - w[i]];
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
        if (keep[i][K] == 1)
        {
            // output i
            K = K - w[i];
        }
    }
    return V[n][W];
}

int main()
{
    
    
    
    return 0;
}