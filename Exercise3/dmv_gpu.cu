/*
 *  dmv_gpu.cu -- Template for DMV GPU kernels
 *
 *  Copyright (C) 2010-2013, Computing Systems Laboratory (CSLab)
 *  Copyright (C) 2010-2013, Vasileios Karakasis
 */ 
#include <stdio.h>
#include "dmv.h"

/*
 *  Utility function to get the thread ID within the
 *  global working space.
 */ 
__device__ int get_global_tid()
{
    return (gridDim.x*blockIdx.y + blockIdx.x)*blockDim.x*blockDim.y +
        blockDim.x*threadIdx.y + threadIdx.x;
}

/*
 *  Utility function to get the thread ID within the
 *  local/block working space.
 */ 
__device__ int get_local_tid()
{
    return blockDim.x*threadIdx.y + threadIdx.x;
}

/*
 *  Naive kernel
 */ 
__global__ void dmv_gpu_naive(const value_t *a, const value_t *x, value_t *y,
                              size_t n)
{
	int i;
	int idx = get_global_tid();

    if (idx >= n) {
        return;
    }
	
	y[idx] = 0;
    for (i = 0; i < n; ++i) {
        y[idx] += a[idx*n+i]*x[i];
	}
        
}

/*
 *  Coalesced memory acceses
 */
__global__ void dmv_gpu_coalesced(const value_t *a, const value_t *x,
                                  value_t *y, size_t n)
{
    int i, j;
	int idx = get_global_tid();

    if (idx >= n) {
        return;
    }
	
	y[idx] = 0;
    for (i = 0; i < n; ++i) {
        y[idx] += a[idx*n+i]*x[i];
	} 
}

/*
 *  Use of shared memory
 */
__global__ void dmv_gpu_shmem(const value_t *a, const value_t *x, value_t *y,
                              size_t n)
{
	int i,j;
	int idx = get_global_tid();
	int idx2 = get_local_tid();
	extern __shared__ value_t shmem[];
	
	if (idx >= n) {
        return;
    }
	
	y[idx] = 0;
	for(j = idx2; j < n; j += blockDim.y) {
        shmem[j] = x[j];
        __syncthreads();
        for(i = j-idx2; i < j-idx2+blockDim.y; i++) {
                y[idx] += a[i*n+idx]*shmem[i];
        }
    }	
}
