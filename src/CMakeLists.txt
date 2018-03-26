cmake_minimum_required(VERSION 3.9)
project(src C)

find_package(CUDA QUIET REQUIRED)

set(CMAKE_C_STANDARD 11)

add_executable(src main.c)