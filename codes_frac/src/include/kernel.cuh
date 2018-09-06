
__device__ complex cro(
                       double sx,
                       double am1,
                       double am2);

__device__ complex propogator980(
                                 double mass,
                                 double g11,
                                 double g22,
                                 double sx);
__device__ complex propogator(
                              double mass,
                              double width,
                              double sx) ;
__device__ complex propogator1270(
                                  double mass,
                                  double width,
                                  double sx) ;
__global__ void cal_fCP(double *par, double *fCP_real, double *fCP_imag, int numbers);

__global__ void propogator1(
                            double mass,
                            double width,
                            double *sx,
                            double *b2qjvf2,
                            double *wu,
                            double *w0p22,
                            double *fCF_real,
                            double *fCF_imag,
                            int numbers);

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
                            int numbers);

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
                            int numbers);

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
                            int numbers );

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
                            int numbers);

__global__ void reduce(double *arrays,int numbers,double *result);
__global__ void cal_phsp(double *fCP_real, double *fCP_imag, double *fCF_real, double *fCF_imag, double *result,int number_of_amplitudes,int numbers);
__global__ void cal_likelihood(double *fCP_real, double *fCP_imag, double *fCF_real, double *fCF_imag,double *fx,  int number_of_amplitudes, int numbers);
__global__ void cal_penalty(double *fCP_real, double *fCP_imag, double *fCF_real, double *fCF_imag,double *result, int number_of_amplitudes, int number);

