#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>

/*
Instead of calling getmax, call a
kernel getmaxcu() that outsources the job of finding the maximum number to the GPU. Let the main function allocate and populate the array in the host. Then you need to allocate memory in the device for that array and transfer the array from the host to the device, calculate the maximum in the device, then transfer that number back to the host.
*/

unsigned int getmax(unsigned int *, unsigned int);

__global__ void getmaxgpu(unsigned int *nums, unsigned int *global_max, unsigned int *size)
{
  __shared__ int block_max; // max of all threads in a block
  if (threadIdx.x == 0) { 
    block_max = 0;
  }
  __syncthreads();
  int index = threadIdx.x + (blockIdx.x * blockDim.x); // get thread index
  if(index < *size)
	  atomicMax(&block_max, nums[index]); // get max across all threads in the block

  __syncthreads();
  atomicMax(global_max, block_max); // get max across all blocks
}


int main(int argc, char *argv[])
{
    unsigned int size = 0;  // The size of the array
    unsigned int i;  // loop index
    unsigned int * numbers; //pointer to the array
    unsigned int * dev_numbers; // device copy of numbers
    unsigned int max = 0;
    unsigned int *d_max;
	unsigned int *d_size;
    int num_blocks;
    int num_threads = 256;
    
    if(argc !=2)
    {
       printf("usage: maxseq num\n");
       printf("num = size of the array\n");
       exit(1);
    }

    size = atol(argv[1]);

    if (size < num_threads) {
      num_threads = size;
      num_blocks = size / num_threads;
	}
	else {
      num_blocks = size / num_threads;
	  if (size%num_threads) {
		  num_blocks += 1;
	  }
    }

    numbers = (unsigned int *)malloc(size * sizeof(unsigned int));
    if( !numbers )
    {
       printf("Unable to allocate mem for an array of size %u\n", size);
       exit(1);
    }    

    srand(time(NULL)); // setting a seed for the random number generator
    // Fill-up the array with random numbers from 0 to size-1 
    for( i = 0; i < size; i++)
       numbers[i] = rand()  % size;

	//printf("The maximum number in the array (seq) is: %u\n", getmax(numbers, size));

    cudaMalloc((void **)&dev_numbers, size * sizeof(unsigned int));
    cudaMalloc((void **)&d_max, sizeof(unsigned int));
	cudaMalloc((void **)&d_size, sizeof(unsigned int));
	cudaMemcpy(d_max, &max, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_size, &size, sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_numbers, numbers, size * sizeof(unsigned int), cudaMemcpyHostToDevice);
	free(numbers);
	//cudaDeviceSynchronize();
    getmaxgpu<<<num_blocks, num_threads>>>(dev_numbers, d_max, d_size);

    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
      fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
    }

    cudaMemcpy(&max, d_max, sizeof(unsigned int), cudaMemcpyDeviceToHost);
    
    cudaFree(dev_numbers);
    cudaFree(&d_max);
	cudaFree(&d_size);

    cudaDeviceSynchronize();
    
    printf("The maximum number in the array (gpu) is: %u\n", max);

    exit(0);
}

/*
   input: pointer to an array of long int
          number of elements in the array
   output: the maximum number of the array
*/

unsigned int getmax(unsigned int num[], unsigned int size)
{
  unsigned int i;
  unsigned int max = num[0];

  for(i = 1; i < size; i++)
	if(num[i] > max)
	   max = num[i];

  return( max );
}


