
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "MultDevice.h"
#include "data_obj.h"
#include "DPFAngular.h"
#include <iostream>

using namespace std;

DataObject::~DataObject()
{
  for(int i = 0; i < number_of_events_; i++) {
    delete[] mcp1[i];
    delete[] mcp2[i];
    delete[] mcp3[i];
    delete[] mcp4[i];
    delete[] mcp5[i];
  }
  delete[] mcp1;
  delete[] mcp2;
  delete[] mcp3;
  delete[] mcp4;
  delete[] mcp5;

  pwa_paras_.resize(0);
}
void DataObject::parameters_vector_resize()
{
    wu.resize(number_of_events_, vector<double>(4));
    w0p22.resize(number_of_events_, vector<double>(4));
    ak23w.resize(number_of_events_, vector<double>(4));
    w2p2.resize(number_of_events_, vector<double>(4));
    w2p1.resize(number_of_events_, vector<double>(4));
    w2p3.resize(number_of_events_, vector<double>(4));
    w2p4.resize(number_of_events_, vector<double>(4));
    w2p5.resize(number_of_events_, vector<double>(4));
    b2qf2xx.resize(number_of_events_, 0);
    b4qjvf2.resize(number_of_events_, 0);
    b2qjv2.resize(number_of_events_, 0);
    b2qjv3.resize(number_of_events_, 0);
    b2qbv2.resize(number_of_events_, 0);
    b2qbv3.resize(number_of_events_, 0);
    b1qjv2.resize(number_of_events_, 0);
    b1qjv3.resize(number_of_events_, 0);
    b1qbv2.resize(number_of_events_, 0);
    b1qbv3.resize(number_of_events_, 0);
    b1q2r23.resize(number_of_events_, 0);
    sv.resize(number_of_events_, 0);
    s23.resize(number_of_events_, 0);
    sv2.resize(number_of_events_, 0);
    sv3.resize(number_of_events_, 0);
    wpf22.resize(number_of_events_, vector<double>(2));
    b2qjvf2.resize(number_of_events_, 0);
    w1p12_1.resize(number_of_events_, vector<double>(4));
    w1p13_1.resize(number_of_events_, vector<double>(4));
    w1p12_2.resize(number_of_events_, vector<double>(4));
    w1p13_2.resize(number_of_events_, vector<double>(4));
    w1p12_3.resize(number_of_events_, vector<double>(2));
    w1p13_3.resize(number_of_events_, vector<double>(2));
    w1p12_4.resize(number_of_events_, vector<double>(2));
    w1p13_4.resize(number_of_events_, vector<double>(2));
    w1m12.resize(number_of_events_, vector<double>(2));
    w1m13.resize(number_of_events_, vector<double>(2));
}
void DataObject::read_events() {
    initialize_mcp();
    load_mcp_from_dat_file();
    cout << "There are " << number_of_events_ << " events in " << dat_file_name_ << std::endl;
    //convert_mcp_to_pwa_paras();
    //cout << "convert data in " << dat_file_name_ << " to PWA_PARAS format." << endl;

}
void DataObject::read_weight_file() {

}
void DataObject::initialize_mcp()
{
  mcp1=new double*[number_of_events_];
  mcp2=new double*[number_of_events_];
  mcp3=new double*[number_of_events_];
  mcp4=new double*[number_of_events_];
  mcp5=new double*[number_of_events_];
  for(int i=0; i<number_of_events_; i++){
    mcp1[i]=new double[4];
    mcp2[i]=new double[4];
    mcp3[i]=new double[4];
    mcp4[i]=new double[4];
    mcp5[i]=new double[4];
  }

}

void DataObject::load_mcp_from_dat_file()
{
  double fx1,fy1,fz1,ft1,fx2,fy2,fz2,ft2,fx3,fy3,fz3,ft3,fx4,fy4,fz4,ft4,fx5,fy5,fz5,ft5;
  FILE *fp;
  if((fp=fopen(dat_file_name_,"r"))==NULL)
  {printf("can't open input file");
    return;
  }
  //cout << "------->start input mcp(pshp)" << _dp->_phspfile << endl;
  int i=0;
  while(fscanf(fp,"%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n",&fx1,&fy1,&fz1,&ft1,&fx2,&fy2,&fz2,&ft2,&fx3,&fy3,&fz3,&ft3,&fx4,&fy4,&fz4,&ft4,&fx5,&fy5,&fz5,&ft5)!=EOF)
  {
    //  //cout<<"haha: "<< __LINE__ << endl;
    mcp1[i][0]=fx1;mcp1[i][1]=fy1;mcp1[i][2]=fz1;mcp1[i][3]=ft1;
    mcp2[i][0]=fx2;mcp2[i][1]=fy2;mcp2[i][2]=fz2;mcp2[i][3]=ft2;
    mcp3[i][0]=fx3;mcp3[i][1]=fy3;mcp3[i][2]=fz3;mcp3[i][3]=ft3;
    mcp4[i][0]=fx4;mcp4[i][1]=fy4;mcp4[i][2]=fz4;mcp4[i][3]=ft4;
    mcp5[i][0]=fx5;mcp5[i][1]=fy5;mcp5[i][2]=fz5;mcp5[i][3]=ft5;
    i++;
  }
  fclose(fp);
  //    Nmc = count_lines() / 5;
  if (i != number_of_events_) {
      std::cout << "Line:" << __LINE__ << "  There is memory allocated error!" << endl;
    exit(1);
  }

}
int DataObject::count_lines() {
  ifstream in(dat_file_name_);
  string line;
  int n = 0;
  while (getline(in, line)) {
    n++;
  }
  //cout << df << " have " << n << " lines." << endl;
  return n;
}
void DataObject::convert_mcp_to_pwa_paras() {
  DPFAngular _amp;
  _amp.setdp(_dp);
  pwa_paras_.resize(0);
  for(int i = 0; i < number_of_events_; i++) {
    PWA_PARAS _the_pwa_paras;
    _amp.calculate0p(
        mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],
        mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],
        mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],
        mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],
        mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3],
        _the_pwa_paras);
    pwa_paras_.push_back(_the_pwa_paras);
  }
}
void DataObject::cpu_convert_mcp_to_pwa_paras() {
  for(int i = 0; i < number_of_events_; i++) {
    cpu_calculate0p(i);
  }
}

double DataObject::scalar(vector<double> &a1, vector<double> &a2) const
{
  double (*fDel)[4];
  fDel = _dp->fDel;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      if(i==j){
        if(i<3)
          fDel[i][j]=-1;
        else
          fDel[i][j]=1;
      } else {
        fDel[i][j]=0;
      }
    }
  }
  double scal=0;

  for(Int_t i=0;i<4;i++){
    scal+=a1[i]*a2[i]*_dp->fDel[i][i];
  }
  return scal;
}
double DataObject::scalar(double *a1, double *a2) const
{
  double (*fDel)[4];
  fDel = _dp->fDel;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      if(i==j){
        if(i<3)
          fDel[i][j]=-1;
        else
          fDel[i][j]=1;
      } else {
        fDel[i][j]=0;
      }
    }
  }
  double scal=0;

  for(Int_t i=0;i<4;i++){
    scal+=a1[i]*a2[i]*_dp->fDel[i][i];
  }
  return scal;
}
double DataObject::cpu_calculate0p(
    //double _p11, double _p12, double _p13, double _p14,
    //double _p21, double _p22, double _p23, double _p24,
    //double _p31, double _p32, double _p33, double _p34,
    //double _p41, double _p42, double _p43, double _p44,
    //double _p51, double _p52, double _p53, double _p54,
    int id)
{
  double ap23[4],apv2[4],apv3[4];
  double ak1[4],ak2[4],ak3[4],ak4[4],ak5[4];
  ak1[0]=mcp1[id][0]; ak1[1]=mcp1[id][1]; ak1[2]=mcp1[id][2]; ak1[3]=mcp1[id][3];
  ak2[0]=mcp2[id][0]; ak2[1]=mcp2[id][1]; ak2[2]=mcp2[id][2]; ak2[3]=mcp2[id][3];
  ak3[0]=mcp3[id][0]; ak3[1]=mcp3[id][1]; ak3[2]=mcp3[id][2]; ak3[3]=mcp3[id][3];
  ak4[0]=mcp4[id][0]; ak4[1]=mcp4[id][1]; ak4[2]=mcp4[id][2]; ak4[3]=mcp4[id][3];
  ak5[0]=mcp5[id][0]; ak5[1]=mcp5[id][1]; ak5[2]=mcp5[id][2]; ak5[3]=mcp5[id][3];
  //ak1[0]=_p11; ak1[1]=_p12; ak1[2]=_p13; ak1[3]=_p14;
  //ak2[0]=_p21; ak2[1]=_p22; ak2[2]=_p23; ak2[3]=_p24;
  //ak3[0]=_p31; ak3[1]=_p32; ak3[2]=_p33; ak3[3]=_p34;
  //ak4[0]=_p41; ak4[1]=_p42; ak4[2]=_p43; ak4[3]=_p44;
  //ak5[0]=_p51; ak5[1]=_p52; ak5[2]=_p53; ak5[3]=_p54;

  double ak23[4],ak45[4],ar[4],ak23u[4],ak45u[4],aru[4];
  double ak12[4],ak13[4],akv2m3[4],akv3m2[4],ak12u[4],ak13u[4],akv2m3u[4],akv3m2u[4];
  double delv2w[4][4],delv3w[4][4],ak12w[4],ak13w[4],akv2m3w[4],akv3m2w[4];
  double del45w[4][4],t2wvf[4][4],wt[4];
  double arw[4];
  double ak23wu[4];
  double fPCMS[4];
  double arwdarw,tmp0;
  fPCMS[0]=fPCMS[1]=fPCMS[2]=0; fPCMS[3]=psi_mass;
  double fFUD=0.22/3.0;
  //double fMK2=kaon_mass*kaon_mass,fMP2=pion_mass*pion_mass
  double fMPsi2=psi_mass*psi_mass;
  double q2r45,b1q2r45,qjvf2;
  double qf2xx,tmpf2,qjv2,qjv3,qbv2,qbv3;
  double q2r23;
  double tmp12w,tmp13w,tmpv2m3,tmpv3m2;
  //double temp;
  double tmp1,tmp2,tmp;
  double (*fDel)[4],(*fGel)[4],(*E)[4][4][4],(*G1)[4][4][4],(*G3)[4][4][4][4][4],t2wf[4][4],del23w[4][4],t2wfu[4][4],w2p1u[4],ttfw[4][4][4],t2p3[4][4],t4wvf[4][4][4][4];
  double t2v2[4][4],t2v3[4][4],t2b3[4][4],t2b2[4][4];
  double w1p12_1u[4],w1p13_1u[4],w1p12_2u[4],w1p13_2u[4];

  //vector<double> ap23(4),apv2(4),apv3(4);
  //vector<double> ak1(4),ak2(4),ak3(4),ak4(4),ak5(4);
  //ak1[0]=_p11; ak1[1]=_p12; ak1[2]=_p13; ak1[3]=_p14;
  //ak2[0]=_p21; ak2[1]=_p22; ak2[2]=_p23; ak2[3]=_p24;
  //ak3[0]=_p31; ak3[1]=_p32; ak3[2]=_p33; ak3[3]=_p34;
  //ak4[0]=_p41; ak4[1]=_p42; ak4[2]=_p43; ak4[3]=_p44;
  //ak5[0]=_p51; ak5[1]=_p52; ak5[2]=_p53; ak5[3]=_p54;

  //vector<double> ak23(4),ak45(4),ar(4),ak23u(4),ak45u(4),aru(4);
  //vector<double> ak12(4),ak13(4),akv2m3(4),akv3m2(4),ak12u(4),ak13u(4),akv2m3u(4),akv3m2u(4);
  //vector<vector<double> > delv2w(4, vector<double>(4));
  //vector<vector<double> > delv3w(4, vector<double>(4));

  //vector<double> ak12w(4);
  //vector<double> ak13w(4);
  //vector<double> akv2m3w(4);
  //vector<double> akv3m2w(4);
  //vector<vector<double> > del45w(4, vector<double>(4));
  //vector<vector<double> > t2wvf(4, vector<double>(4));
  //vector<double> wt(4);
  //vector<double> arw(4);
  //vector<double> ak23wu(4);
  //vector<double> fPCMS(4);
  //double arwdarw,tmp0;
  //fPCMS[0]=fPCMS[1]=fPCMS[2]=0; fPCMS[3]=psi_mass;
  //double fFUD=0.22/3.0;
  ////double fMK2=kaon_mass*kaon_mass,fMP2=pion_mass*pion_mass
  //double fMPsi2=psi_mass*psi_mass;
  //double q2r45,b1q2r45,qjvf2;
  //double qf2xx,tmpf2,qjv2,qjv3,qbv2,qbv3;
  //double q2r23;
  //double tmp12w,tmp13w,tmpv2m3,tmpv3m2;
  ////double temp;
  //double tmp1,tmp2,tmp;
  //double (*fDel)[4],(*fGel)[4],(*E)[4][4][4],(*G1)[4][4][4],(*G3)[4][4][4][4][4],t2wf[4][4],del23w[4][4],t2wfu[4][4],w2p1u[4],ttfw[4][4][4],t2p3[4][4],t4wvf[4][4][4][4];
  //double t2v2[4][4],t2v3[4][4],t2b3[4][4],t2b2[4][4];
  //vector<double> w1p12_1u(4),w1p13_1u(4),w1p12_2u(4),w1p13_2u(4);

  fDel = _dp->fDel;
  fGel = _dp->fGel;
  E = _dp->E;
  G1 = _dp->G1;
  G3 = _dp->G3;

  //#include dpfangular_cc_commet.cc 影响程序的格式，因此这部分被注释的内容放到了dpfangular_cc_commet中

  //cout<<"haha: "<< __LINE__ << endl;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      if(i==j){
        if(i<3)
          fDel[i][j]=-1;
        else
          fDel[i][j]=1;
      } else {
        fDel[i][j]=0;
      }
    }
  }
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      if(i==j){
        if(i<3)
          fGel[i][j]=-1;
        else
          fGel[i][j]=0;
      } else {
        fGel[i][j]=0;
      }
    }
  }
  //cout<<"haha: "<< __LINE__ << endl;
  E[0][1][2][3]= 1.0;
  E[0][1][3][2]=-1.0;
  E[0][2][1][3]=-1.0;
  E[0][2][3][1]= 1.0;
  E[0][3][1][2]= 1.0;
  E[0][3][2][1]=-1.0;
  E[1][0][2][3]=-1.0;
  E[1][0][3][2]= 1.0;
  E[1][2][0][3]= 1.0;
  E[1][2][3][0]=-1.0;
  E[1][3][0][2]=-1.0;
  E[1][3][2][0]= 1.0;
  E[2][0][1][3]= 1.0;
  E[2][0][3][1]=-1.0;
  E[2][1][0][3]=-1.0;
  E[2][1][3][0]= 1.0;
  E[2][3][0][1]= 1.0;
  E[2][3][1][0]=-1.0;
  E[3][0][1][2]=-1.0;
  E[3][0][2][1]= 1.0;
  E[3][1][0][2]= 1.0;
  E[3][1][2][0]=-1.0;
  E[3][2][0][1]=-1.0;
  E[3][2][1][0]= 1.0;
  for(Int_t I =0;I<4;I++){
    for(Int_t J=0;J<4;J++){
      for(Int_t K=0;K<4;K++){
        for(Int_t L=0;L<4;L++){
          // G1[I][J][K][L] = (fGel[J][K]*fGel[I][L] + fGel[K][I]*fGel[J][L])/2 -(fGel[I][J]*fGel[K][L])/3 ;
          G1[I][J][K][L] = fGel[I][J]*fGel[K][L] + fGel[J][K]*fGel[I][L] +fGel[K][I]*fGel[J][L]; // last (...), i.e. combination of \tilde{g} in Eq.(37) of PRD48,1225
        }
      }
    }
  }// Eq.(20)
  for(Int_t I1=0;I1<4;I1++){
    for(Int_t I2=0;I2<4;I2++){
      for(Int_t I3=0;I3<4;I3++){
        for(Int_t I4=0;I4<4;I4++){
          for(Int_t I5=0;I5<4;I5++){
            for(Int_t I6=0;I6<4;I6++){
              G3[I1][I2][I3][I4][I5][I6] =
                ( fGel[I1][I2]*fGel[I4][I5]*fGel[I3][I6] +
                  fGel[I1][I2]*fGel[I5][I6]*fGel[I3][I4] +
                  fGel[I1][I2]*fGel[I4][I6]*fGel[I3][I5] +
                  fGel[I1][I3]*fGel[I4][I6]*fGel[I2][I5] +
                  fGel[I1][I3]*fGel[I4][I5]*fGel[I2][I6] +
                  fGel[I1][I3]*fGel[I5][I6]*fGel[I2][I4] +
                  fGel[I2][I3]*fGel[I5][I6]*fGel[I1][I4] +
                  fGel[I2][I3]*fGel[I4][I5]*fGel[I1][I6] +
                  fGel[I2][I3]*fGel[I4][I6]*fGel[I1][I5] )/15.0 -
                ( fGel[I1][I4]*fGel[I2][I5]*fGel[I3][I6] +
                  fGel[I1][I4]*fGel[I2][I6]*fGel[I3][I5] +
                  fGel[I1][I5]*fGel[I2][I4]*fGel[I3][I6] +
                  fGel[I1][I5]*fGel[I2][I6]*fGel[I3][I4] +
                  fGel[I1][I6]*fGel[I2][I5]*fGel[I3][I4] +
                  fGel[I1][I6]*fGel[I2][I4]*fGel[I3][I5] )/6.0;
            }
          }
        }
      }
    }
  }// Eq. (21)

  // 1 phi 2 pi+ 3 pi- 4 K+ 5 K-
  // ak1 ~ ak5, four momentum from input, ak1 is the sum of 3 and 4, etc. the phi's momentum
  for(Int_t i=0;i<4;i++){
    ap23[i]=ak2[i]+ak3[i]; // pi+ + pi-
    apv2[i]=ak1[i]+ak2[i]; // phi + pi+
    apv3[i]=ak1[i]+ak3[i]; // phi + pi-
    ak23[i]=ak2[i]-ak3[i]; // r_34 (relative difference between pi+ and pi-)
    ak45[i]=ak4[i]-ak5[i]; // r_12 (relative difference between K+ and K-)
    ar[i]=ak1[i]-ap23[i];  // r^\mu (relative difference between phi and two pions)
    //
    ak23u[i]=ak23[i]*fDel[i][i];
    ak45u[i]=ak45[i]*fDel[i][i];
    aru[i]=ar[i]*fDel[i][i];// r_\mu(phi f)
  }
  sv[id]=scalar(ak1,ak1);      // p^2(phi)
  s23[id]=scalar(ap23,ap23);   // p^2(f)
  sv2[id]=scalar(apv2,apv2);  //  p^2(123)
  sv3[id]=scalar(apv3,apv3); //   p^2(124)

  // 0+ contribution
  // prepare tensors for Eq.(36)
  q2r45=sv[id]*0.25 - (_dp->_m[0])*(_dp->_m[0]); // q2r45 is Q^2_abc, Eq.(13), phi(a)-> K+(b) K-(c)
  b1q2r45=sqrt(2/(q2r45+fFUD)); // b1q2r45 is B1(Q_abc), Eq.(14), it can be ignored for phi is narrow, i.e. it should be always close to a constant
  // when final states' momenta change
  //cout << "b1q2r45 " << b1q2r45 << endl;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      del45w[i][j]=fDel[i][j]-ak1[i]*ak1[j]/sv[id];
      //cout << "del45w " << del45w[i][j] << endl;
    }
  }// del45w is \tilde{g}^{\mu\nu}(\phi)

  for(Int_t i=0;i<4;i++){
    wt[i]=0.0;
    for(Int_t j=0;j<4;j++){
      wt[i]=wt[i]+del45w[i][j]*ak45u[j];
    }
    wt[i]=wt[i]*b1q2r45; // Eq. (10), \tilde{r}^\mu = t^{(1)\mu}(12)
    wu[id][i]=wt[i]*fDel[i][i]; // \tilde{r}^\mu = t^{(1)}_\mu(12)
  }
  // prepare tensors for Eq.(37)
  for(Int_t i=0;i<4;i++){
    arw[i]=0.0;
    for(Int_t j=0;j<4;j++){
      arw[i]=arw[i]+fGel[i][j]*aru[j];  // \tilde(r)^\mu(phi f0)
    }
  }
  arwdarw=scalar(arw,arw);
  tmp0=(-arwdarw)/3;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      t2wvf[i][j]=arw[i]*arw[j]+tmp0*fGel[i][j]; // \tilde{T}^{(2)\mu\nu}(phi f0), from Eq.(11)
      //cout << "t2wvf " << t2wvfuu[i][j] << endl;
    }
  }
  qjvf2=0.25*(pow((psi_mass*psi_mass+sv[id]-s23[id]),2))/(psi_mass*psi_mass)-sv[id]; // Eq.(13), a is psi, b is phi, c is resonance
  b2qjvf2[id]=sqrt(13./(pow(qjvf2,2)+3.*qjvf2*fFUD+9.*pow(fFUD,2))); // B2(Qabc) in Eq.(15), can be neglected since psi is very narrow
  //cout << "b2qjvf2 " << b2qjvf2 << endl;

  //0+ contribution
  for(Int_t i=0;i<2;i++){
    w0p22[id][i]=0.0;
    for(Int_t j=0;j<4;j++){
      w0p22[id][i]=w0p22[id][i]+t2wvf[i][j]*wu[id][j]; // the \tilde{T}^{(2)\mu\nv}(phi f0) \tilde{t}^{(1)}_{\nu}(12)  in Eq.(37), only prepare \mu = 1 and 2
    }
  }
  // for(int i=0; i<4; i++) cout << "w0pp2 " << w0p22[i] << endl;
  // Till now, the tensors in Eq.(36) [wt] and (37) [w0p22] have been prepared.
  // The missing parts are propagators, etc. BW.

  //   2+ contribution
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      del23w[i][j]=fDel[i][j]-ap23[i]*ap23[j]/s23[id]; // \tilde{g}^{\mu\nu}(f2) in Eq.(11)
    }
  }
  for(Int_t i=0;i<4;i++){
    ak23w[id][i]=0.0;
    for(Int_t j=0;j<4;j++){
      ak23w[id][i]=ak23w[id][i]+del23w[i][j]*ak23u[j]; // \tilde{r}^{\mu}(f2) in Eq.(11), r(a)-> pi+(b) pi-(c)
    }
  }
  qf2xx=0.25*s23[id]-(_dp->_m[2])*(_dp->_m[2]);
  b2qf2xx[id]=sqrt(13./(pow(qf2xx,2)+3.*qf2xx*fFUD+9.*pow(fFUD,2))); // B2 in Eq.(15) of pi+ pi-
  b4qjvf2[id]=sqrt(12746./(pow(qjvf2,4)+10*pow(qjvf2,3)*fFUD+135*pow(qjvf2,2)*pow(fFUD,2)+1575*qjvf2*pow(fFUD,3)+11025*pow(fFUD,4))); // B4 in Eq.(17) of psi
  tmpf2=-scalar(ak23w[id],ak23w[id])/3;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      t2wf[i][j]=(ak23w[id][i]*ak23w[id][j]+tmpf2*del23w[i][j])*b2qf2xx[id]; // \tilde{t}^{(2)\mu\nu}(f2) in Eq.(11), \tilde{t}^{(2)\mu\nu}(34) in Eq.(41)
      t2wfu[i][j]=t2wf[i][j]*fDel[j][j];  // \tilde{t}^{(2)\mu}_{\nu}(34)
      //cout << "t2wfu " << t2wfu[i][j] << endl;
    }
  }
  for(Int_t i=0;i<4;i++){
    w2p1[id][i]=0.0;
    for(Int_t j=0;j<4;j++){
      w2p1[id][i]=w2p1[id][i]+t2wf[i][j]*wu[id][j]; // \tilde{t}^{(2)\mu\nu}(34) \tilde{t}^{(1)}_{\nu}(12) the second and third terms in Eq.(42)
    }
    w2p1u[i]=w2p1[id][i]*fDel[i][i]; // \tilde{t}^{(2)}_{\mu\nu}(34) \tilde{t}^{(1)\nu}(12)
  }
  for(Int_t i=0;i<2;i++){
    w2p2[id][i]=0.0;
    for(Int_t j=0;j<4;j++){
      w2p2[id][i]=w2p2[id][i]+t2wvf[i][j]*w2p1u[j]; // \tilde{T}^{(2)\mu\alpha}(phi f2)\tilde{t}^{(2)}_{\alpha\nu}(34) \tilde{t}^{(1)\nu}(12)， the first three terms in Eq.(42)
    }
  }
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      for(Int_t k=0;k<4;k++){
        ttfw[i][j][k]=t2wfu[i][j]*wt[k]; // \tilde{t}^{(2)\lambda}_\delta(34) \tilde{t}^{(1)\nu}(12)
      }
    }
  }
  //  cout<<"haha: "<< __LINE__ << endl;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<4;j++){
      t2p3[i][j]=0.0;
      for(Int_t k1=0;k1<3;k1++){
        for(Int_t k2=0;k2<3;k2++){
          t2p3[i][j]=t2p3[i][j]+(E[i][k1][3][k2]*ttfw[k1][j][k2]+E[j][k1][3][k2]*ttfw[k1][i][k2]); // [...]p^\sigma_psi \tilde{t}^{(1)\nu}(12)
        }
      }
    }
  }
  for(Int_t i=0;i<2;i++){
    w2p3[id][i]=0.0;
    for(Int_t j1=0;j1<3;j1++){
      for(Int_t j2=0;j2<3;j2++){
        for(Int_t j3=0;j3<4;j3++){
          w2p3[id][i]=w2p3[id][i]+E[i][j1][j2][3]*t2wvf[j1][j3]*t2p3[j2][j3]; // epsilon^{\mu\alpha\beta\gamma}p_{\psi\alpha}\tilde^{(2)\lambda}_{\beta}(\phi f2)[...] in Eq. (43)
        }
      }
    }
  }
  double test[4][4][4];
  for(int i=0;i<4;i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<4; k++){
        test[i][j][k]=0;
        for(int m=0; m<4; m++){
          for(int n=0; n<4; n++){
            test[i][j][k] = test[i][j][k] + ttfw[i][m][n]*fDel[m][j]*fDel[n][k];
          }
        }
      }
    }
  }
  for(Int_t i=0;i<2;i++){
    w2p4[id][i]=0.0;
    for(Int_t j1=0;j1<4;j1++){
      for(Int_t j2=0;j2<4;j2++){
        for(Int_t j3=0;j3<4;j3++){
          for(Int_t j4=0;j4<4;j4++){
            for(Int_t j5=0;j5<4;j5++){
              w2p4[id][i]=w2p4[id][i]+G3[i][j1][j2][j3][j4][j5]*t2wvf[j1][j2]*ttfw[j3][j4][j5]; //  P^{(3)}\tilde{T}^{(2)}(\phi f2)\tilde{t}^{(2)}(34)\tilde{t}^{(1)}(12) Eq.(44)
              // w2p4[i]=w2p4[i]+G3[i][j1][j2][j3][j4][j5]*t2wvf[j1][j2]*test[j3][j4][j5];
            }
          }
        }
      }
    }
  }
  tmp1=(-arwdarw)/7;
  tmp2=(pow(arwdarw,2))/35;
  for(Int_t i=0;i<2;i++){
    for(Int_t j=0;j<4;j++){
      for(Int_t k=0;k<4;k++){
        for(Int_t l=0;l<4;l++){
          t4wvf[i][j][k][l]=arw[i]*arw[j]*arw[k]*arw[l]+tmp1*(fGel[i][j]*arw[k]*arw[l]+fGel[j][k]*arw[i]*arw[l]+fGel[k][i]*arw[j]*arw[l]+fGel[i][l]*arw[j]*arw[k]+fGel[j][l]*arw[k]*arw[i]+fGel[k][l]*arw[i]*arw[j])+tmp2*G1[i][j][k][l]; // defined in Eq(37) of PRD48,1225
        }
      }
    }
  }
  for(Int_t i=0;i<2;i++){
    w2p5[id][i]=0.0;
    for(Int_t j1=0;j1<4;j1++){
      for(Int_t j2=0;j2<4;j2++){
        for(Int_t j3=0;j3<4;j3++){
          w2p5[id][i]=w2p5[id][i]+t4wvf[i][j1][j2][j3]*ttfw[j2][j3][j1]; // \tilde{T}^{(4)}(\phi f2) \tilde{t}^{(1)}(12) \tilde{t}^{(2)}(34), in Eq.(45)
        }
      }
    }
  }
  // 1+ and 1- contribution
  qjv2=0.25*pow((fMPsi2+sv2[id]-(_dp->_m[2])*(_dp->_m[2])),2)/fMPsi2-sv2[id]; // Q^2(abc) = Q^2(\psi' 123 \pi4)
  qjv3=0.25*pow((fMPsi2+sv3[id]-(_dp->_m[2])*(_dp->_m[2])),2)/fMPsi2-sv3[id]; // Q^2(\psi 124 \pi3)
  qbv2=0.25*pow((sv2[id]+sv[id]-(_dp->_m[2])*(_dp->_m[2])),2)/sv2[id]-sv[id];   // Q^2(123 12 \pi3)
  qbv3=0.25*pow((sv3[id]+sv[id]-(_dp->_m[2])*(_dp->_m[2])),2)/sv3[id]-sv[id];   // Q^2(124 12 \pi4)
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      delv2w[i][j]=fDel[i][j]-apv2[i]*apv2[j]/sv2[id];   // \tilde{g}^{\mu\nu}(123)
      delv3w[i][j]=fDel[i][j]-apv3[i]*apv3[j]/sv3[id];   // \tilde{g}^{\mu\nu}(124)
    }
  }
  for(int i=0;i<4;i++){
    ak12[i]=ak1[i]-ak2[i];        // r^\mu(phi pi3) i.e. p_{12}-p_3
    ak13[i]=ak1[i]-ak3[i];        // r^\mu(phi pi4) i.e. p_{12}-p_4
    akv2m3[i]=apv2[i]-ak3[i];     // r^\mu(rho' pi4) i.e. p_{123}-p_4
    akv3m2[i]=apv3[i]-ak2[i];     // r^\mu(rho' pi3) i.e. p_{124}-p_3
    ak12u[i]=ak12[i]*fDel[i][i];  // r_\mu(phi pi3)
    ak13u[i]=ak13[i]*fDel[i][i];  // r_\mu(phi pi4)
    akv2m3u[i]=akv2m3[i]*fDel[i][i]; // r_\mu(rho' pi4)
    akv3m2u[i]=akv3m2[i]*fDel[i][i]; // r_\mu(rho' pi3)
  }
  for(int i=0;i<4;i++){
    ak12w[i]=0.0;
    ak13w[i]=0.0;
    akv2m3w[i]=0.0;
    akv3m2w[i]=0.0;
    for(int j=0;j<4;j++){
      ak12w[i]=ak12w[i]+delv2w[i][j]*ak12u[j];      // \tilde{r}^\mu (phi pi3)
      ak13w[i]=ak13w[i]+delv3w[i][j]*ak13u[j];      // \tilde{r}^\mu (phi pi4)
      akv2m3w[i]=akv2m3w[i]+fGel[i][j]*akv2m3u[j];  // \tilde{r}^\mu (b pi4)
      akv3m2w[i]=akv3m2w[i]+fGel[i][j]*akv3m2u[j];  // \tilde{r}^\mu (b pi3)
    }
  }
  // 1+ contribution
  b2qjv2[id]=sqrt(13./(pow(qjv2,2)+3.*qjv2*fFUD+9.*pow(fFUD,2)));  // B2(Q psi b pi4)
  b2qjv3[id]=sqrt(13./(pow(qjv3,2)+3.*qjv3*fFUD+9.*pow(fFUD,2)));  // B2(Q psi b pi3)
  b2qbv2[id]=sqrt(13./(pow(qbv2,2)+3.*qbv2*fFUD+9.*pow(fFUD,2)));  // B2(Q b phi pi3)
  b2qbv3[id]=sqrt(13./(pow(qbv3,2)+3.*qbv3*fFUD+9.*pow(fFUD,2)));  // B2(Q b phi pi4)
  tmp12w=-scalar(ak12w,ak12w)/3.0;
  tmp13w=-scalar(ak13w,ak13w)/3.0;
  tmpv2m3=-scalar(akv2m3w,akv2m3w)/3.0;
  tmpv3m2=-scalar(akv3m2w,akv3m2w)/3.0;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      t2v2[i][j]=(ak12w[i]*ak12w[j]+tmp12w*delv2w[i][j])*b2qbv2[id];    // \tilde{t}^{(2)\mu\nu}(phi pi3)
      t2v3[i][j]=(ak13w[i]*ak13w[j]+tmp13w*delv3w[i][j])*b2qbv3[id];    // \tilde{t}^{(2)\mu\nu}(phi pi4)
      t2b3[i][j]=(akv2m3w[i]*akv2m3w[j]+tmpv2m3*fGel[i][j])*b2qjv2[id]; // \tilde{T}^{(2)\mu\nu}(b pi4)
      t2b2[i][j]=(akv3m2w[i]*akv3m2w[j]+tmpv3m2*fGel[i][j])*b2qjv3[id]; // \tilde{T}^{(2)\mu\nu}(b pi3)
    }
  }
  for(int i=0;i<4;i++){
    w1p12_1[id][i]=0.0;
    w1p13_1[id][i]=0.0;
    w1p12_2[id][i]=0.0;
    w1p13_2[id][i]=0.0;
    for(int j=0;j<4;j++){
      w1p12_1[id][i]=w1p12_1[id][i]+delv2w[i][j]*wu[id][j]; // \tilde{g}^{\mu\nu}(123) \tilde{t}^{(1)}_\nu(12) Eq.(47)
      w1p13_1[id][i]=w1p13_1[id][i]+delv3w[i][j]*wu[id][j]; // \tilde{g}^{\mu\nu}(124) \tilde{t}^{(1)}_\nu(12) Eq.(47)
      w1p12_2[id][i]=w1p12_2[id][i]+t2v2[i][j]*wu[id][j];   // \tilde{t}_{(2)\mu\nu}(phi 3) \tilde{t}^{(1\nu}(12) Eq.(48)
      w1p13_2[id][i]=w1p13_2[id][i]+t2v3[i][j]*wu[id][j];   // \tilde{t}^{(2)\mu\nu}(phi 4) \tilde{t}^{(1)\nu}(12) Eq.(48)
    }
    // cout << "w1p12_1 " << w1p12_1[i] << "  w1p13_1 " << w1p13_1[i] << "  w1p12_2 " << w1p12_2[i] << "  w1p13_2 " << w1p13_2[i] << endl;
    w1p12_1u[i]=w1p12_1[id][i]*fDel[i][i];
    w1p13_1u[i]=w1p13_1[id][i]*fDel[i][i];
    w1p12_2u[i]=w1p12_2[id][i]*fDel[i][i];
    w1p13_2u[i]=w1p13_2[id][i]*fDel[i][i];
  }
  for(int i=0;i<2;i++){
    w1p12_3[id][i]=0.0;
    w1p13_3[id][i]=0.0;
    w1p12_4[id][i]=0.0;
    w1p13_4[id][i]=0.0;
    for(int j=0;j<4;j++){
      w1p12_3[id][i]=w1p12_3[id][i]+t2b3[i][j]*w1p12_1u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 4} \tilde{g}_{(123)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(49)
      w1p13_3[id][i]=w1p13_3[id][i]+t2b2[i][j]*w1p13_1u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 3} \tilde{g}_{(124)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(49)
      w1p12_4[id][i]=w1p12_4[id][i]+t2b3[i][j]*w1p12_2u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 4} \tilde{t}_{(phi 3)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(50)
      w1p13_4[id][i]=w1p13_4[id][i]+t2b2[i][j]*w1p13_2u[j]; // \tilde{T}^{(2)\mu\lambda}_{b1 3} \tilde{g}_{(phi 4)\lambda \nu} \tilde{t}^{(1)\nu}_{(12)} in Eq.(50)
    }
    // cout << "w1p12_3 " << w1p12_3[i] << "  w1p13_3 " << w1p13_3[i] << "  w1p12_4 " << w1p12_4[i] << "  w1p13_4 " << w1p13_4[i] << endl;
  }
  //  1- contribution
  b1qjv2[id]=sqrt(2./(qjv2+fFUD)); // B1(Q psi rho' pi4), rho'(123)
  b1qjv3[id]=sqrt(2./(qjv3+fFUD)); // B1(Q psi rho' pi3), rho'(124)
  b1qbv2[id]=sqrt(2./(qbv2+fFUD)); // B1(Q rho' phi pi3), rho'(123)
  b1qbv3[id]=sqrt(2./(qbv3+fFUD)); // B1(Q rho' phi pi4), rho'(124)
  for(int i=0;i<2;i++){
    w1m12[id][i]=0.0;
    w1m13[id][i]=0.0;
    for(int j1=0;j1<3;j1++){
      for(int j2=0;j2<3;j2++){
        for(int j3=0;j3<3;j3++){
          for(int j4=0;j4<3;j4++){
            tmp=E[i][j1][j2][4]*E[j2][j3][j4][4];
            w1m12[id][i]=w1m12[id][i]+tmp*akv2m3w[j1]*ak12w[j3]*wt[j4]*b1qbv2[id]; // The second term in Eq.(46)
            w1m13[id][i]=w1m13[id][i]+tmp*akv3m2w[j1]*ak13w[j3]*wt[j4]*b1qbv3[id]; // The first term in Eq.(46)
          }
        }
      }
    }
  }
  // phi(1650)f0(980)  contribution // mass 1680 +- 20 MeV, width 150 +- 50 MeV, dominant decay KK*(892)
  q2r23=0.25*s23[id]-(_dp->_m[2])*(_dp->_m[2]); // Q^2_{abc}(f->pi+ pi-)
  b1q2r23[id]=sqrt(2./(q2r23+fFUD));  // B1(Q_abc) = B1(Q f pi+ pi-)
  for(int i=0;i<4;i++){
    ak23wu[i]=ak23w[id][i]*fDel[i][i]*b1q2r23[id]; // \tilde{t}_mu(34)
  }
  for(int i=0;i<2;i++){
    wpf22[id][i]=0.0;
    for(int j=0;j<4;j++){
      wpf22[id][i]=wpf22[id][i]-t2wvf[i][j]*ak23wu[j]; // -\tilde{T}^{(2)\mu\nu}(phi f0) \tilde{t}_nu(34)
    }
  }
  return 0;
}


bool DataObject::check_integraty()
{
    for(int i = 0; i < number_of_events_; i++) {
        if (pwa_paras_[i].b2qjvf2 != b2qjvf2[i]) return false;
        if (pwa_paras_[i].w0p22[0] != w0p22[i][0]) return false;
        if (pwa_paras_[i].w0p22[1] != w0p22[i][1]) return false;
        if (pwa_paras_[i].ak23w[1] != ak23w[i][1]) return false;
        if (pwa_paras_[i].w1m12[1] != w1m12[i][1]) return false;
    }
    return true;

}

