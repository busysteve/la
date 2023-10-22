
//  nvcc -o matrix_multiply matrix_multiply.cu



//Create a C++ project and include the necessary headers:
#include <iostream>
#include <vector>
#include <cuda_runtime.h>


//Define a function to allocate and initialize matrices:
void initializeMatrix(std::vector<float>& mat, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            mat[i * cols + j] = static_cast<float>(rand() % 100);
        }
    }
}


//Define the matrix multiplication kernel:
__global__ void matrixMultiplyKernel(const float* A, const float* B, float* C, int m, int n, int k) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < m && col < k) {
        float sum = 0.0f;
        for (int i = 0; i < n; ++i) {
            sum += A[row * n + i] * B[i * k + col];
        }
        C[row * k + col] = sum;
    }
}


//Define the main function:
int main() {
    const int m = 4;  // Number of rows in A
    const int n = 3;  // Number of columns in A and rows in B
    const int k = 5;  // Number of columns in B

    // Allocate and initialize matrices A and B
    std::vector<float> hostA(m * n);
    std::vector<float> hostB(n * k);
    std::vector<float> hostC(m * k);
    initializeMatrix(hostA, m, n);
    initializeMatrix(hostB, n, k);

    // Allocate device memory for A, B, and C
    float *deviceA, *deviceB, *deviceC;
    cudaMalloc(&deviceA, m * n * sizeof(float));
    cudaMalloc(&deviceB, n * k * sizeof(float));
    cudaMalloc(&deviceC, m * k * sizeof(float));

    // Copy matrices A and B from host to device
    cudaMemcpy(deviceA, hostA.data(), m * n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceB, hostB.data(), n * k * sizeof(float), cudaMemcpyHostToDevice);

    // Define grid and block dimensions for the kernel
    dim3 gridDim(2, 2); // You can adjust these dimensions based on your GPU and matrix size
    dim3 blockDim(2, 2);

    // Launch the matrix multiplication kernel
    matrixMultiplyKernel<<<gridDim, blockDim>>>(deviceA, deviceB, deviceC, m, n, k);
    cudaDeviceSynchronize();

    // Copy the result matrix C from the device to host
    cudaMemcpy(hostC.data(), deviceC, m * k * sizeof(float), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(deviceA);
    cudaFree(deviceB);
    cudaFree(deviceC);

    // Output the result
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < k; ++j) {
            std::cout << hostC[i * k + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}







