#include "complex.h"
#include"kernel.cuh"
#include"kernel.h"
#include"whu_constants_and_definitions.h"
#include <assert.h>
#include<iostream>
#include<cuda_runtime.h>
#include<stdio.h>
#include <cuda.h>
using namespace std;

#define CUDA_CALL(x) {const cudaError_t a=(x); if(a != cudaSuccess) {printf("\nerror in line:%d CUDAError:%s(err_num=%d)\n",__LINE__,cudaGetErrorString(a),a); cudaDeviceReset(); assert(0); }}




