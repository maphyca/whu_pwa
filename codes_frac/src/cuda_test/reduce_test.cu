#include<iostream>
#include <unistd.h>
#include<cuda.h>
#include<stdio.h>
#include<time.h>
#include <sys/time.h>
const int threads_per_block=64;
clock_t start,end;
__global__ void reduce(int *arrays,int numbers,int *result)
{
  int id = threadIdx.x + blockDim.x*blockIdx.x;
  __shared__ int s_arrays[threads_per_block];
  if (id< numbers)
  {
    s_arrays[threadIdx.x]=arrays[id];
  }
  else
  {
    s_arrays[threadIdx.x]=0;
  }
      __syncthreads();
  for(int i = blockDim.x/2;i>=1;i/=2)
    {
      if(threadIdx.x<i)
      {
      s_arrays[threadIdx.x]+=s_arrays[threadIdx.x+i];
      }
      __syncthreads();
    }
  result[blockIdx.x]=s_arrays[0];
}
int main(int argc, char *argv[])
{
  int numbers=atoi(argv[1]);
  int *a = new int[numbers];
  for(int i=0;i<numbers;i++)
    a[i]=1;
  int *b = new int[numbers];
  int *c = new int[numbers];
  int *da,*db;
  cudaMalloc((void**)&da,numbers*sizeof(int));
  cudaMalloc((void**)&db,numbers*sizeof(int));
  cudaMemcpyAsync(da,a,numbers*sizeof(int),cudaMemcpyHostToDevice);
  int count = numbers;
  while(1)
    {
      reduce<<<(count+63)/64,threads_per_block>>>(da,count,db);
      cudaDeviceSynchronize();
      count=(count+63)/64;
      std::cout<<count<<std::endl;
      
      if(count==1)
        {
          cudaMemcpyAsync(b,db,numbers*sizeof(int),cudaMemcpyDeviceToHost);
          break;
        }
      reduce<<<(count+63)/64,threads_per_block>>>(db,count,da);
      cudaDeviceSynchronize();
      std::cout<<count<<std::endl;
      count=(count+63)/64;
      if(count==1)
        {
          cudaMemcpyAsync(b,da,numbers*sizeof(int),cudaMemcpyDeviceToHost);
          break;
        }

    }

  for (int i=0;i<numbers;i++)
    {
      c[0]+=a[i];
    }
  start=clock();
  std::cout<<"b0 "<<b[0]<<" c0 "<<c[0]<<std::endl;
  return 0;
}
