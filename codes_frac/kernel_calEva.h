/*************************************************************************
	> File Name: kernel_calEva.h
	> Author:
	> Mail:
	> Created Time: 2017年03月21日 星期二 18时20分53秒
 ************************************************************************/

#ifndef _KERNEL_CALEVA_H
#define _KERNEL_CALEVA_H
#include "/usr/local/cuda/include/cuda_runtime.h"
#include "cuComplex.h"
#include <iostream>
#include "cu_PWA_PARAS.h"
#include <vector>
#include <fstream>
#include <math.h>
#include "cu_DPFPropogator.h"
#include <assert.h>
#include <vector>
#include "MultDevice.h"
//double *d_float_pp=NULL;
//#include "cu_PWA_PARAS.h"
    //int initialize_data(std::vector<PWA_PARAS>&, DataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //int data_distribution(DataPointers&, CudaDataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //double calEva(const PWA_PARAS &pp, int idp);
    //double kernel_calEva(const PWA_PARAS &pp,int idp);
    //void func(DataPointers& cpu_data_pointers);

class cuda_kernel {
    public:

    cuda_kernel() {};
    double *d_fx[DEVICE_NUM];
    int *d_parameter[DEVICE_NUM];
    double *d_paraList[DEVICE_NUM];
    double2 * d_complex_para[DEVICE_NUM];
    double *d_mlk[DEVICE_NUM];

    void cu_malloc_h_pp(double *,double *&,int,int);
    int host_store_fx(std::vector<double *>,int *,double *,int , double *,double * ,int ,int );
    //int malloc_mem(int end, int begin, int para_size, int *h_parameter)
    //
    int warp_malloc_mem(int, int, int, int *);
    int malloc_mem(int, int, int, int*);
};
    __device__ double calEva(const cu_PWA_PARAS *pp, const int * parameter , double2 * complex_para ,const double * d_paraList,double *d_mlk,int idp) ;
    __global__ void kernel_store_fx(const double * float_pp,const int *parameter,double2 * d_complex_para ,const double *d_paraList,int para_size,double * d_fx,double *d_mlk,int end,int begin);
#endif

