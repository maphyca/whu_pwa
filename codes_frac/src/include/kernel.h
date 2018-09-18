#include<vector>
class kernel{
  //pp 指针
public:
  double *d_b2qbv2;
  double *d_b2qbv3;
  double *d_b2qjv2;
  double *d_b2qjv3;
  double *d_b2qjvf2;
  double *d_b2qf2xx;
  double *d_b1qbv2;
  double *d_b1qbv3;
  double *d_b1q2r23;
  double *d_b4qjvf2;
  double *d_b1qjv2;
  double *d_b1qjv3;
  double *d_sv2;
  double *d_sv3;
  double *d_s23;
  double *d_wu;
  double *d_wpf22;
  double *d_w0p22;
  double *d_w2p1;
  double *d_w2p2;
  double *d_w2p3;
  double *d_w2p4;
  double *d_w2p5;
  double *d_w1p12_1;
  double *d_w1p13_1;
  double *d_w1p12_2;
  double *d_w1p13_2;
  double *d_w1p12_3;
  double *d_w1p13_3;
  double *d_w1p12_4;
  double *d_w1p13_4;
  double *d_w1m12;
  double *d_w1m13;
  //fCP,fCF指针
  double *fCP_real;
  double *fCP_imag;
  double *fCF_real;
  double *fCF_imag;

  double *h_par;
  double *h_par_back;
  double *d_par;
  double *d_phsp;
  double *d_penalty;
  double *d_likelihood;
  double *d_container;
  double *d_weight;

  double h_penalty;
  double h_phsp;
  double h_likelihood;

  double *h_phsp_container;

  int number_of_amplitudes;
  int Threads,Blocks;
  int number_of_data;
  kernel(std::vector<double*>, int, int, int, int, int);
  kernel();
  void par_trans(const std::vector<double>  ) const;
  void calEva();
  double sum_phsp();
  double sum_likelihood();
  double sum_penalty();
  void trans_phsp();


};
