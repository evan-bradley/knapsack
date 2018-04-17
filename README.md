# knapsack
This repository contains code experimenting with parallel algorithms that solve
the 0/1 knapsack problem.

# Compilation
| Algorithm name | command |
| ---            | ---     |
| Serial DP | gcc -g dp.c -o dp |
| CUDA DP | nvcc -arch=compute_61 -rdc=true -g dp-cuda.cu -o dpcu -lcudadevrt |
| Serial GA | gcc -lm -g ga.c -o ga |
| MPI GA | mpicc -g mpi.c -o mpi |
| CUDA GA | nvcc -lcurand ga.cu -o ./gacu |

#Running
| Algorithm wame | command    |
| ---           | ---        |
| Serial DP |./dp|
| CUDA DP |./dpcu|
| Serial GA |./ga| 
| MPI GA |mpirun -np 4 ./mpi|
| CUDA GA |./gacu| 
