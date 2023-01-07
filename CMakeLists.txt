cmake_minimum_required(VERSION 3.16)
project(iSONG CUDA)

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
set(CUDA_INCLUDE /usr/local/cuda/include)


#查找cuda环境
# find_package(CUDA REQUIRED)
#头文件
include_directories(${CUDA_INCLUDE})
#库文件lib
#link_directories()
add_executable(iSONG main.cu bin_heap.cuh)

set_target_properties(iSONG PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
