#include <iostream> 
#include <stdlib.h> 
#include <time.h>
#include "complex.h"
#include "whu_constants_and_definitions.h"
#include "TComplex.h"
using namespace std;
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

__global__ void propogator980(
                                 double mass,
                                 double g11,
                                 double g22,
                                 double sx,
                                 double *test_r,
                                 double *test_i)
{
  complex ci(0,1);
  double rm=mass*mass;
  complex propogator980=1.0/(rm-sx-ci*(g11*cro(sx,rp,rp)+g22*cro(sx,rk,rk)));
  test_r= propogator980.x;
  test_i= propogator980.y;
}
/*__global__ void propogator2(
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
TComplex cpu_cro(double sx, double am1, double am2)  {
    TComplex ci(0, 1);
    double t1 = (am1 + am2) * (am1 + am2);  // double t1=pow((am1+am2),2);
    double t2 = (am1 - am2) * (am1 - am2);  // double t2=pow((am1-am2),2);
    double st = (sx - t1) * (sx - t2);
    double cro = sqrt(fabs(st)) / sx;
    TComplex result = cro;
    if (st < 0.) result = cro * ci;
    return result;
}
*/
TComplex cpu_propogator980(double mass, double g11, double g22,
                                    double sx)  {
    TComplex ci(0, 1);
    double rm = mass * mass;
    TComplex propogator980 =
        1.0 / (rm - sx - ci * (g11 * cpu_cro(sx, rp, rp) + g22 * cpu_cro(sx, rk, rk)));
    return propogator980;
}
void cpu_propogator2(double mass, double g11, double g22,
                                  double *sx, double *b2qjvf2, double *wu,
                                  double *w0p22, TComplex *fCF0, TComplex *fCF1,
                                  int vec_size) {
    for (int i = 0; i < vec_size; i++) {
        TComplex crp1 = cpu_propogator980(mass, g11, g22, sx[i]);
        TComplex cr0p11 = crp1 / b2qjvf2[i];

        // 01 contribution
        fCF0[i + vec_size * 0] = wu[i + vec_size * 0] * crp1;
        fCF0[i + vec_size * 1] = wu[i + vec_size * 1] * crp1;

        // 02 contribution
        fCF1[i + vec_size * 0] = w0p22[i + vec_size * 0] * cr0p11;
        fCF1[i + vec_size * 1] = w0p22[i + vec_size * 1] * cr0p11;
    }
}

int main()
{
  int numbers = 1000;
  double *pa,*pb,*pc,*pd;
  double *pda,*pdb,*pdc,*pdd,*pde,*pdf,*pfr,*pfi;
  TComplex *pe;
  pa = new double[numbers];
  pb = new double[numbers];
  pc = new double[numbers*2];
  pd = new double[numbers*2];
  pfr = new double[numbers*4];
  pfi = new double[numbers*4];

  pe = new TComplex[numbers*4];
  for(int i=0;i<numbers;i++)
    {
      pa[i]=rand()/(double)RAND_MAX;
      pb[i]=rand()/(double)RAND_MAX;
      pc[i]=rand()/(double)RAND_MAX;
      pd[i]=rand()/(double)RAND_MAX;
      pc[i+numbers]=rand()/(double)RAND_MAX;
      pd[i+numbers]=rand()/(double)RAND_MAX;
    }

  cudaMalloc((void**)&pda,numbers*sizeof(double));
  cudaMalloc((void**)&pdb,numbers*sizeof(double));
  cudaMalloc((void**)&pdc,2*numbers*sizeof(double));
  cudaMalloc((void**)&pdd,2*numbers*sizeof(double));
  cudaMalloc((void**)&pde,4*numbers*sizeof(double));
  cudaMalloc((void**)&pdf,4*numbers*sizeof(double));

  cudaMemcpyAsync(pda,pa,numbers*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpyAsync(pdb,pb,numbers*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpyAsync(pdc,pc,2*numbers*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpyAsync(pdd,pd,2*numbers*sizeof(double),cudaMemcpyHostToDevice);

  propogator2<<<(numbers+63)/64,64>>>(1,1,1,pda,pdb,pdc,pdd,pde,pdf,numbers);
  cudaMemcpyAsync(pfr,pde,4*numbers*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpyAsync(pfi,pdf,4*numbers*sizeof(double),cudaMemcpyDeviceToHost);
  cpu_propogator980(1,1,1,pa,pb,pc,pd,pe,,numbers);
  for(int i=0;i<20;i++)
    {
      cout<<"cpu p2 "<<" : "<<pe[0]<<"  gpu :"<<pfr[0]<<"+ "<<pfi[0]<<endl;
    }
   return 0;
}


