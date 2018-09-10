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

//const double rk=0.493677;

/*//const double rp=0.139556995;
#ifndef CAFFE_COMMON_CUH_
#define CAFFE_COMMON_CUH_


#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600

#else
static __inline__ __device__ double atomicAdd(double *address, double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  if (val==0.0)
    return __longlong_as_double(old);
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}


#endif
#endif
*/
__device__ complex cro(
                       double sx,
                       double am1,
                       double am2) 
{
  double t1=(am1+am2) * (am1 + am2); // double t1=pow((am1+am2),2);
  double t2=(am1-am2) * (am1 - am2); // double t2=pow((am1-am2),2);
  double st=(sx-t1)*(sx-t2);
  double cro=sqrt(fabs(st))/sx;
  complex result;
  if (st<0.) result.y=cro;
  else result.x=cro;
  return  result;
}

__device__ complex propogator980(
                                 double mass,
                                 double g11,
                                 double g22,
                                 double sx)
{
  complex ci(0,1);
  double rm=mass*mass;
  complex propogator980=1.0/(rm-sx-ci*(g11*cro(sx,rp,rp)+g22*cro(sx,rk,rk)));
  return propogator980;
}

__device__ complex pip(
                       double sx)
{
  //?    complex ci(0,1);
  double xk2=sx-0.3116676;     //0.3116676=16.*0.139568*0.139568
  if(xk2<=0.)xk2=0.0;
  double r4pip=sqrt(xk2/sx)/(1.0+exp(9.8-3.5*sx));    //9.8=3.5*2.8
  return  make_complex(r4pip,0);
}

__device__ complex propogator(
                              double mass,
                              double width,
                              double sx) 
{
  complex ci(0,1);
  double am=mass;
  double g1=mass*width;
  complex prop=g1/(sx-am * am +ci*g1); // complex prop=g1/(sx-pow(am,2)+ci*g1);
  return prop;
}

__device__ complex propogator1270(
                                  double mass,
                                  double width,
                                  double sx) 
{
  complex ci(0,1);
  double rm=mass*mass;
  double gr=mass*width;
  double q2r=0.25*rm-0.0194792;
  double b2r=q2r*(q2r+0.1825)+0.033306;
  double g11270=gr*b2r/pow(q2r,2.5);
  double q2=0.25*sx-0.0194792;
  double b2=q2*(q2+0.1825)+0.033306;
  double g1=g11270*pow(q2,2.5)/b2;
  complex prop=gr/(sx-rm+ci*g1);
  return prop;
}
__global__ void cal_fCP(double *par, double *fCP_real, double *fCP_imag, int numbers)
{
  int id = threadIdx.x+blockDim.x*blockIdx.x;
  if(id<numbers){
    double rho0 = par[end_category * id + rho_category];
    double frac0 = par[end_category * id + frac_category];
    double phi0 = par[end_category * id + phi_category];
    rho0 *= exp(frac0);
    fCP_real[id]=rho0*cos(phi0);
    fCP_imag[id]=rho0*sin(phi0);
  }
}
__global__ void propogator1(
                            double mass,
                            double width,
                            double *sx,
                            double *b2qjvf2,
                            double *wu,
                            double *w0p22,
                            double *fCF_real,
                            double *fCF_imag,
                            int numbers)
{
  int id = threadIdx.x+blockDim.x*blockIdx.x;
  if(id<numbers){
    complex crp1 = propogator(mass, width, sx[id]);
    complex cr0p11 = crp1 / b2qjvf2[id];
    complex result;
    //01 contribution
    result = wu[id] * crp1;
    fCF_real[id] =result.x; 
    fCF_imag[id] =result.y; 
    result = wu[id+numbers] * crp1;
    fCF_real[id+numbers] =result.x; 
    fCF_imag[id+numbers] =result.y; 

    //02 contribution
    result = w0p22[id] * cr0p11;
    fCF_real[id+numbers*2] =result.x; 
    fCF_imag[id+numbers*2] =result.y; 
    result = w0p22[id+numbers] * cr0p11;
    fCF_real[id+numbers*3] =result.x; 
    fCF_imag[id+numbers*3] =result.y; 
  }

}
__global__ void propogator2(
                            double mass,
                            double g11,
                            double g22,
                            double *sx,
                            double *b2qjvf2,
                            double *wu,
                            double *w0p22,
                            double *fCF_real,
                            double *fCF_imag,
                            int numbers)
{
  int id = threadIdx.x+blockDim.x*blockIdx.x;
  if(id<numbers){
    id=threadIdx.x+blockDim.x*blockIdx.x;
    complex crp1 = propogator980(mass, g11, g22, sx[id]);
    complex cr0p11 = crp1 / b2qjvf2[id];
    complex result;
    //01 contribution
    result = wu[id] * crp1;
    fCF_real[id] =result.x; 
    fCF_imag[id] =result.y; 
    result = wu[id+numbers] * crp1;
    fCF_real[id+numbers] =result.x; 
    fCF_imag[id+numbers] =result.y; 

    //02 contribution
    result = w0p22[id] * cr0p11;
    fCF_real[id+numbers*2] =result.x; 
    fCF_imag[id+numbers*2] =result.y; 
    result = w0p22[id+numbers] * cr0p11;
    fCF_real[id+numbers*3] =result.x; 
    fCF_imag[id+numbers*3] =result.y;
  }

}
__global__ void propogator7(
                            double mass,
                            double width,
                            double *sv2,
                            double *sv3,
                            double *b1qjv2,
                            double *b1qbv2,
                            double *b1qjv3,
                            double *b1qbv3,
                            double *w1m12,
                            double *w1m13,
                            double *fCF_real,
                            double *fCF_imag,
                            int numbers)
{
  int id = threadIdx.x+blockDim.x*blockIdx.x;
  if(id<numbers){
    complex crp1 = propogator(mass, width, sv2[id]);
    complex crp11 = propogator(mass, width, sv3[id]);
    complex cr1m12_1 = crp1 / b1qjv2[id] / b1qbv2[id];
    complex cr1m13_1 = crp11 / b1qjv3[id] / b1qbv3[id];
    complex result;
    //1-__1 contribution

    result = w1m12[id] * cr1m12_1 + w1m13[id] * cr1m13_1;
    fCF_real[id] =result.x; 
    fCF_imag[id] =result.y; 
    result = w1m12[id+numbers] * cr1m12_1 + w1m13[id+numbers] * cr1m13_1;
    fCF_real[id+numbers] =result.x; 
    fCF_imag[id+numbers] =result.y; 
  }
}
__global__ void propogator8(
                            double mass,
                            double width,
                            double *sv2,
                            double *sv3,
                            double *b2qbv2,
                            double *b2qbv3,
                            double *b2qjv2,
                            double *b2qjv3,
                            double *w1p12_1,
                            double *w1p13_1,
                            double *w1p12_2,
                            double *w1p13_2,
                            double *w1p12_3,
                            double *w1p13_3,
                            double *w1p12_4,
                            double *w1p13_4,
                            double *fCF_real,
                            double *fCF_imag,
                            int numbers )
{
  int id = threadIdx.x+blockDim.x*blockIdx.x;
  if(id<numbers){
    complex crp1 = propogator(mass, width, sv2[id]);
    complex crp11 = propogator(mass, width, sv3[id]);
    complex c1p12_12 = crp1 / b2qbv2[id];
    complex c1p13_12 = crp11 / b2qbv3[id];
    complex c1p12_13 = crp1 / b2qjv2[id];
    complex c1p13_13 = crp11 / b2qjv3[id];
    complex c1p12_14 = c1p12_12 / b2qjv2[id];
    complex c1p13_14 = c1p13_12 / b2qjv3[id];
    complex result;
    // z 1+ 1
    result = w1p12_1[id] * crp1 + w1p13_1[id] * crp11;
    fCF_real[id] =result.x; 
    fCF_imag[id] =result.y; 
    result = w1p12_1[id+numbers] * crp1 + w1p13_1[id+numbers] * crp11;
    fCF_real[id+numbers] =result.x; 
    fCF_imag[id+numbers] =result.y; 

    // z 1+ 2
    result = w1p12_2[id] * c1p12_12 + w1p13_2[id] * c1p13_12;
    fCF_real[id+numbers*2] =result.x; 
    fCF_imag[id+numbers*2] =result.y; 
    result = w1p12_2[id+numbers] * c1p12_12 + w1p13_2[id+numbers] * c1p13_12;
    fCF_real[id+numbers*3] =result.x; 
    fCF_imag[id+numbers*3] =result.y; 

    // z 1+ 3
    result = w1p12_3[id] * c1p12_13 + w1p13_3[id] * c1p13_13;
    fCF_real[id+numbers*4] =result.x; 
    fCF_imag[id+numbers*4] =result.y; 
    result = w1p12_3[id+numbers] * c1p12_13 + w1p13_3[id+numbers] * c1p13_13;
    fCF_real[id+numbers*5] =result.x; 
    fCF_imag[id+numbers*5] =result.y; 

    // z 1+ 4
    result = w1p12_4[id] * c1p12_14 + w1p13_4[id] * c1p13_14;
    fCF_real[id+numbers*6] =result.x; 
    fCF_imag[id+numbers*6] =result.y; 
    result = w1p12_4[id+numbers] * c1p12_14 + w1p13_4[id+numbers] * c1p13_14;
    fCF_real[id+numbers*7] =result.x; 
    fCF_imag[id+numbers*7] =result.y; 
  }
}
__global__ void propogator6(
                            double mass,
                            double width,
                            double *sx,
                            double *b2qf2xx,
                            double *b2qjvf2,
                            double *b4qjvf2,
                            double *w2p1,
                            double *w2p2,
                            double *w2p3,
                            double *w2p4,
                            double *w2p5,
                            double *fCF_real,
                            double *fCF_imag,
                            int numbers)
{
  int id = threadIdx.x+blockDim.x*blockIdx.x;
  if(id<numbers){
    complex crp1 = propogator1270(mass, width, sx[id]);
    complex cw2p11 = crp1 / b2qf2xx[id];
    complex cw2p12 = cw2p11 / b2qjvf2[id];
    complex cw2p15 = cw2p11 / b4qjvf2[id];
    complex result;
    //21 contribution
    result = w2p1[id] * cw2p11;
    fCF_real[id] =result.x; 
    fCF_imag[id] =result.y; 
    result = w2p1[id+numbers] * cw2p11;
    fCF_real[id+numbers] =result.x; 
    fCF_imag[id+numbers] =result.y; 

    //22 contribution
    result=w2p2[id] * cw2p12;
    fCF_real[id+numbers*2] = result.x;
    fCF_imag[id+numbers*2] = result.y;
    result= w2p2[id+numbers] * cw2p12;
    fCF_real[id+numbers*3]=result.x;
    fCF_imag[id+numbers*3]=result.y;

    //23 contribution
    result = w2p3[id] * cw2p12;
    fCF_real[id+numbers*4] =result.x; 
    fCF_imag[id+numbers*4] =result.y; 
    result = w2p3[id+numbers] * cw2p12;
    fCF_real[id+numbers*5] =result.x; 
    fCF_imag[id+numbers*5] =result.y; 

    //24 contribution
    result = w2p4[id] * cw2p12;
    fCF_real[id+numbers*6] =result.x; 
    fCF_imag[id+numbers*6] =result.y; 
    result = w2p4[id+numbers] * cw2p12;
    fCF_real[id+numbers*7] =result.x; 
    fCF_imag[id+numbers*7] =result.y; 

    //25 contribution
    result = w2p5[id] * cw2p15;
    fCF_real[id+numbers*8] =result.x; 
    fCF_imag[id+numbers*8] =result.y; 
    result = w2p5[id+numbers] * cw2p15;
    fCF_real[id+numbers*9] =result.x; 
    fCF_imag[id+numbers*9] =result.y; 
  }
}
void kernel::calEva()
{
  //cout << "number_of_amplitudes = " << number_of_amplitudes << endl;
  //cout << "number_of_events_ = " << number_of_events_ << endl;
  int i = 0;
  while (i < number_of_amplitudes)
    {
      int propType_now = h_par[end_category * i + propType_category];
      cout << "amplitude i = " << i << endl;
      cout << "propType_now = " << propType_now << endl;
      switch(propType_now)
        {
        case 1: // f0
          {
            double mass0 = h_par[end_category * i + mass_category];
            double width0 = h_par[end_category * i + width_category];
            bool _not_changed =
              ((mass0 == h_par_back[end_category * i + mass_category])
               && (width0 == h_par_back[end_category * i + width_category]));
            if (!_not_changed) {
              cout << "prop = " << propType_now << " : work one time!!!" << endl;
              cudaDeviceSynchronize();
              propogator1<<<Blocks,Threads>>>(
                                              mass0,
                                              width0,
                                              d_s23,
                                              d_b2qjvf2,
                                              d_wu,
                                              d_w0p22,
                                              fCF_real+2*i*number_of_data,
                                              fCF_imag+2*i*number_of_data,
                                              number_of_data);
            }
            i = i + 2;
          }
          break;
          //	Flatte   Propagator Contribution
        case 2: // f0 980
          {
            double mass980 = h_par[end_category * i + mass_category];
            double g10 = h_par[end_category * i + g1_category];
            double g20 = h_par[end_category * i + g2_category];
            bool _not_changed =
              ((mass980 == h_par_back[end_category * i + mass_category])
               && (g10 == h_par_back[end_category * i + g1_category])
               && (g20 == h_par_back[end_category * i + g2_category]));
            if (!_not_changed) {
              cout << "prop = " << propType_now << " : work one time!!!" << endl;
              cudaDeviceSynchronize();
              propogator2<<<Blocks,Threads>>>(
                                              mass980,
                                              g10,
                                              g20,
                                              d_s23,
                                              d_b2qjvf2,
                                              d_w0p22,
                                              d_wu,
                                              fCF_real+2*i*number_of_data,
                                              fCF_imag+2*i*number_of_data,
                                              number_of_data);
            }
            i = i + 2;
          }
          break;
        case 7: //1m1800
          {
            double mass0 = h_par[end_category * i + mass_category];
            double width0 = h_par[end_category * i + width_category];
            bool _not_changed =
              ((mass0 == h_par_back[end_category * i + mass_category])
               && (width0 == h_par_back[end_category * i + width_category]));
            if (!_not_changed)
              {
                cout << "prop = " << propType_now << " : work one time!!!" << endl;
                cudaDeviceSynchronize();
                propogator7<<<Blocks,Threads>>>(
                                                mass0,
                                                width0,
                                                d_sv2,
                                                d_sv3,
                                                d_b1qjv2,
                                                d_b1qbv2,
                                                d_b1qjv3,
                                                d_b1qbv3,
                                                d_w1m12,
                                                d_w1m13,
                                                fCF_real+2*number_of_data*i,
                                                fCF_imag+2*number_of_data*i,
                                                number_of_data);
              }
            i = i + 1;
          }
          break;
        case 8: //1p1800
          {
            double mass0 = h_par[end_category * i + mass_category];
            double width0 = h_par[end_category * i + width_category];
            bool _not_changed =
              ((mass0 == h_par_back[end_category * i + mass_category])
               && (width0 == h_par_back[end_category * i + width_category]));
            if (!_not_changed)
              {
                cout << "prop = " << propType_now << " : work one time!!!" << endl;
                cudaDeviceSynchronize();
                propogator8<<<Blocks,Threads>>>(
                                                mass0,
                                                width0,
                                                d_sv2,
                                                d_sv3,
                                                d_b2qbv2,
                                                d_b2qbv3,
                                                d_b2qjv2,
                                                d_b2qjv3,
                                                d_w1p12_1,
                                                d_w1p13_1,
                                                d_w1p12_2,
                                                d_w1p13_2,
                                                d_w1p12_3,
                                                d_w1p13_3,
                                                d_w1p12_4,
                                                d_w1p13_4,
                                                fCF_real+2*i*number_of_data,
                                                fCF_imag+2*i*number_of_data,
                                                number_of_data);
              }
            i = i + 4;
          }
          break;
        case 6: //f2
          {
            double mass0 = h_par[end_category * i + mass_category];
            double width0 = h_par[end_category * i + width_category];
            bool _not_changed =
              ((mass0 == h_par_back[end_category * i + mass_category])
               && (width0 == h_par_back[end_category * i + width_category]));
            if (!_not_changed) {
              cout << "prop = " << propType_now << " : work one time!!!" << endl;
              cudaDeviceSynchronize();
              propogator6<<<Blocks,Threads>>>(
                                              mass0,
                                              width0,
                                              d_s23,
                                              d_b2qf2xx,
                                              d_b2qjvf2,
                                              d_b4qjvf2,
                                              d_w2p1,
                                              d_w2p2,
                                              d_w2p3,
                                              d_w2p4,
                                              d_w2p5,
                                              fCF_real+2*i*number_of_data,
                                              fCF_imag+2*i*number_of_data,
                                              number_of_data);
            }
            i = i + 5;
          }
          break;
        default :
          cout << "Do not know how to deal with prop type " << propType_now << endl;
          exit(1);
          ;
        }
    }

  cudaDeviceSynchronize();
  cal_fCP<<<1,number_of_amplitudes>>>(d_par,fCP_real,fCP_imag,number_of_amplitudes);
  cudaDeviceSynchronize();

}

__global__ void reduce(double *arrays,int numbers,double *result)
{
  int id = threadIdx.x + blockDim.x*blockIdx.x;
  extern __shared__ double s_arrays[];
  if (id< numbers) s_arrays[threadIdx.x]=arrays[id];
  else s_arrays[threadIdx.x]=0;
  __syncthreads();
  for(int i = blockDim.x/2;i>=1;i/=2)
    {
      if(threadIdx.x<i)
      s_arrays[threadIdx.x]+=s_arrays[threadIdx.x+i];
      __syncthreads();
    }
  result[blockIdx.x]=s_arrays[0];
}


__global__ void cal_phsp(double *fCP_real, double *fCP_imag, double *fCF_real, double *fCF_imag, double *result,int number_of_amplitudes,int numbers)
{

  int id=threadIdx.x+blockIdx.x*blockDim.x;
  if(id<numbers)
    {
      complex cw1,cw2,fCP,fCF;
      cw1=cw2=complex(0,0);
      for(int i = 0; i < number_of_amplitudes; i+=1)
        {
          fCP = make_complex(fCP_real[i],fCP_imag[i]);
          fCF = make_complex(fCF_real[id+2*i*numbers],fCF_imag[id+2*i*numbers]);
          cw1 = cw1 + fCP * fCF;

          fCF = make_complex(fCF_real[id+(i*2+1)*numbers],fCF_imag[id+(i*2+1)*numbers]);
  
   
  
          cw2 = cw2 + fCP * fCF;
        }
      
      result[id]=(real(cw1) * real(cw1) + imag(cw1) * imag(cw1) + real(cw2) * real(cw2) + imag(cw2) * imag(cw2)) / 2.0;
    }
}



__global__ void cal_likelihood(double *fCP_real, double *fCP_imag, double *fCF_real, double *fCF_imag, double *fx, int number_of_amplitudes,int numbers)
{
  int id=threadIdx.x+blockIdx.x*blockDim.x;
  if(id<numbers)
    {
      complex cw1,cw2,fCP,fCF;
      cw1=cw2=complex(0,0);
      for(int i = 0; i < number_of_amplitudes; i+=1)
        {
          fCP = make_complex(fCP_real[i],fCP_imag[i]);
          fCF = make_complex(fCF_real[id+2*i*numbers],fCF_imag[id+2*i*numbers]);
          cw1 = cw1 + fCP * fCF;

          fCF = make_complex(fCF_real[id+(i*2+1)*numbers],fCF_imag[id+(i*2+1)*numbers]);
          cw2 = cw2 + fCP * fCF;
        }
      fx[id] = -log((real(cw1) * real(cw1) + imag(cw1) * imag(cw1) + real(cw2) * real(cw2) + imag(cw2) * imag(cw2)) / 2.0);
    }
}

__global__ void cal_penalty(double *fCP_real, double *fCP_imag, double *fCF_real, double *fCF_imag,double *result,  int number_of_amplitudes,int numbers)
{

   int id=threadIdx.x+blockIdx.x*blockDim.x;
   complex fCP,fCF,cw1,cw2;
  double temp;
  cw1=cw2=complex(0,0);
  temp=0;
  if(id<numbers)
    {
      for(int i = 0; i < number_of_amplitudes; i+=1)
        {
          fCP = make_complex(fCP_real[i],fCP_imag[i]);
          fCF = make_complex(fCF_real[id+2*i*numbers],fCF_imag[id+2*i*numbers]);
          cw1 = fCP * conj(fCP);
          cw2 = fCF * conj(fCF) / 2.0;
          fCF = make_complex(fCF_real[id+(i*2+1)*numbers],fCF_imag[id+(i*2+1)*numbers]);
          cw2 = cw2 + fCF * conj(fCF) / 2.0;
          temp+= real(cw1) * real(cw2);
        }
      result[id]=temp;
  
    }
}


kernel::kernel(std::vector<double *> Data, int Device_id, int start, int end,int nAmps,int numbers)
{
  h_par_back=new double[nAmps*end_category];
  h_par=new double[nAmps*end_category];
  h_phsp_container=new double[numbers];
  Threads = threads_per_block;
  Blocks = (end - start + Threads -1)/Threads;
  number_of_data = end-start;
  number_of_amplitudes = nAmps;
  CUDA_CALL(cudaSetDevice(Device_id));
  CUDA_CALL(cudaMalloc((void **)&d_par,number_of_amplitudes*sizeof(double)*end_category));
  CUDA_CALL(cudaMalloc((void **)&d_container,(end-start)*sizeof(double)));

  //prop1
  CUDA_CALL(cudaMalloc((void **)&d_w0p22, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_wu, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b2qjvf2, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_s23, (end-start)*sizeof(double)));
  //prop6
  CUDA_CALL(cudaMalloc((void **)&d_b2qf2xx, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b4qjvf2, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w2p1, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w2p2, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w2p3, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w2p4, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w2p5, 2*(end-start)*sizeof(double)));
  //prop7
  CUDA_CALL(cudaMalloc((void **)&d_sv2, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_sv3, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b1qjv2, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b1qbv2, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b1qjv3, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b1qbv3, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1m12, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1m13, 2*(end-start)*sizeof(double)));
  //prop8
  CUDA_CALL(cudaMalloc((void **)&d_b2qbv2, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b2qbv3, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b2qjv2, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_b2qjv3, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p12_1, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p13_1, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p12_2, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p13_2, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p12_3, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p13_3, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p12_4, 2*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_w1p13_4, 2*(end-start)*sizeof(double)));



  CUDA_CALL(cudaMalloc((void **)&fCP_real, nAmps*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&fCP_imag, nAmps*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&fCF_real, 2*nAmps*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&fCF_imag, 2*nAmps*(end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_phsp, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_likelihood, (end-start)*sizeof(double)));
  CUDA_CALL(cudaMalloc((void **)&d_penalty, (end-start)*sizeof(double)));


  //prop1
  CUDA_CALL(cudaMemcpyAsync(d_b2qjvf2, Data[b2qjvf2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));

  CUDA_CALL(cudaMemcpyAsync(d_wu, Data[wu]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_wu+(end-start), Data[wu]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w0p22, Data[w0p22]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w0p22+(end-start), Data[w0p22]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_s23, Data[s23]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));

  //prop6
  CUDA_CALL(cudaMemcpyAsync(d_b2qf2xx, Data[b2qf2xx]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b4qjvf2, Data[b4qjvf2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p1, Data[w2p1]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p1+(end-start), Data[w2p1]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p2, Data[w2p2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p2+(end-start), Data[w2p2]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p3, Data[w2p3]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p3+(end-start), Data[w2p3]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p4, Data[w2p4]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p4+(end-start), Data[w2p4]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p5, Data[w2p5]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w2p5+(end-start), Data[w2p5]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));

  //prop7
  CUDA_CALL(cudaMemcpyAsync(d_sv3, Data[sv3]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_sv2, Data[sv2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b1qjv2, Data[b1qjv2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b1qjv3, Data[b1qjv3]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b1qbv2, Data[b1qbv2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b1qbv3, Data[b1qbv3]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));

  CUDA_CALL(cudaMemcpyAsync(d_w1m12, Data[w1m12]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1m12+(end-start), Data[w1m12]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1m13, Data[w1m13]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1m13+(end-start), Data[w1m13]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));


  //prop8
  CUDA_CALL(cudaMemcpyAsync(d_b2qbv2, Data[b2qbv2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b2qbv3, Data[b2qbv3]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b2qjv2, Data[b2qjv2]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_b2qjv3, Data[b2qjv3]+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));


  CUDA_CALL(cudaMemcpyAsync(d_w1p12_1, Data[w1p12_1]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p12_1+(end-start), Data[w1p12_1]+numbers+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_1, Data[w1p13_1]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_1+(end-start), Data[w1p13_1]+numbers+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p12_2, Data[w1p12_2]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p12_2+(end-start), Data[w1p12_2]+numbers+start, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_2, Data[w1p13_2]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_2+(end-start), Data[w1p13_2]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p12_3, Data[w1p12_3]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p12_3+(end-start), Data[w1p12_3]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_3, Data[w1p13_3]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_3+(end-start), Data[w1p13_3]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p12_4, Data[w1p12_4]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p12_4+(end-start), Data[w1p12_4]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_4, Data[w1p13_4]+start, 2*(end-start)*sizeof(double),cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpyAsync(d_w1p13_4+(end-start), Data[w1p13_4]+start+numbers, (end-start)*sizeof(double),cudaMemcpyHostToDevice));
  cudaDeviceSynchronize();


 
  

}

void kernel::par_trans(const std::vector<double> par) const
{
  cout<<"number_of_amplitudes: "<<number_of_amplitudes<<endl;
  CUDA_CALL(cudaMemcpyAsync(d_par, &par[0],end_category*number_of_amplitudes*sizeof(double),cudaMemcpyHostToDevice));
  for(int i=0;i<number_of_amplitudes*end_category;i++)
    {
      h_par[i]=par[i];
    }
}


kernel::kernel()
{
}

double kernel::sum_penalty()
{
  cal_penalty<<<Blocks,Threads>>>(fCP_real, fCP_imag, fCF_real,fCF_imag,d_penalty,  number_of_amplitudes, number_of_data);
  int count=number_of_data;
    while(1)
    {
      reduce<<<(count+Threads-1)/Threads,Threads,Threads*sizeof(double)>>>(d_penalty,count,d_container);
      count=(count+Threads-1)/Threads;
      
      if(count==1)
        {
          cudaMemcpyAsync(&h_penalty,d_container,sizeof(double),cudaMemcpyDeviceToHost);
          break;
        }
      reduce<<<(count+Threads-1)/Threads,Threads,Threads*sizeof(double)>>>(d_container,count,d_penalty);
      count=(count+Threads-1)/Threads;
      if(count==1)
        {
          cudaMemcpyAsync(&h_penalty,d_penalty,sizeof(double),cudaMemcpyDeviceToHost);
          break;
        }

        }
    return h_penalty;
}

double kernel::sum_phsp()
{
    cal_phsp<<<Blocks,Threads>>>(fCP_real, fCP_imag, fCF_real,fCF_imag,d_phsp,number_of_amplitudes,number_of_data);
  int count=number_of_data;
     while(1)
    {
      reduce<<<(count+Threads-1)/Threads,Threads,Threads*sizeof(double)>>>(d_phsp,count,d_container);
      count=(count+Threads-1)/Threads;
      
      if(count==1)
        {
          cudaMemcpyAsync(&h_phsp,d_container,sizeof(double),cudaMemcpyDeviceToHost);
          break;
        }
      reduce<<<(count+Threads-1)/Threads,Threads,Threads*sizeof(double)>>>(d_container,count,d_phsp);
      count=(count+Threads-1)/Threads;
      if(count==1)
        {
          cudaMemcpyAsync(&h_phsp,d_phsp,sizeof(double),cudaMemcpyDeviceToHost);
          break;
        }
      
        }
     return h_phsp;
}

double  kernel::sum_likelihood()
{
  cal_likelihood<<<Blocks,Threads>>>(fCP_real, fCP_imag, fCF_real,fCF_imag,d_likelihood, number_of_amplitudes,number_of_data);
   int count=number_of_data;
     while(1)
    {
      reduce<<<(count+Threads-1)/Threads,Threads,Threads*sizeof(double)>>>(d_likelihood,count,d_container);
      count=(count+Threads-1)/Threads;
      
      if(count==1)
        {
          cudaMemcpyAsync(&h_likelihood,d_container,sizeof(double),cudaMemcpyDeviceToHost);
          break;
        }
      reduce<<<(count+Threads-1)/Threads,Threads,Threads*sizeof(double)>>>(d_container,count,d_likelihood);
      count=(count+Threads-1)/Threads;
      if(count==1)
        {
          cudaMemcpyAsync(&h_likelihood,d_likelihood,sizeof(double),cudaMemcpyDeviceToHost);
          break;
        }

        }
     return h_likelihood;
}
 
void kernel::trans_phsp()
{
    cal_phsp<<<Blocks,Threads>>>(fCP_real, fCP_imag, fCF_real,fCF_imag,d_phsp,number_of_amplitudes,number_of_data);
    cudaMemcpyAsync(h_phsp_container,d_phsp,number_of_data*sizeof(double),cudaMemcpyDeviceToHost);

}

