#include <iostream>
#include <cuda_runtime.h>

int main() {
    std::cout << "Hello, World!" << std::endl;
    cudaDeviceReset();
    return 0;
}
