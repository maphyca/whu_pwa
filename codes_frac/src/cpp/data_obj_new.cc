#ifndef DATA_H
#define DATA_H
#include"data_obj_new.h"
#endif

#ifndef CONFIG_H 
#define CONFIG_H
#include"pwa_conf.h"
#endif
#ifndef TCOMPLEX_H
#define TCOMPLEX_H
#include "TComplex.h"
#endif

#include "whu_constants_and_definitions.h"
using namespace std;

#define ss(x) (hp[x])
#define vv(x, y) (hp[x + y])
#define vp(x) (&hp[x])
#define mm(x, y, z) (hp[x + y * 4 + z])

DataObject::DataObject()
{
  cerr<<"warning : a new empty class is being created ";
}

DataObject::~DataObject()
{
  delete [] _mcp;
  delete [] _weight;

}

DataObject::DataObject(string d_name, DPFPWAPoint* pwa_point)
{
  data_file_name = d_name;
  read_data();
  _weight = new double[_number_of_events];
  fill(_weight,_weight+_number_of_events,1.0);//权重统一初始化为1，避免出错
  _dp = pwa_point;
}

DataObject::DataObject(string d_name, string w_name,DPFPWAPoint* pwa_point)
{
  data_file_name = d_name;
  weight_file_name = w_name;
  read_data();
  _weight = new double[_number_of_events];
  read_weight();
  _dp = pwa_point;
}

void DataObject::read_weight()
{
  ifstream iwf;
  iwf.open(weight_file_name,ios::in);
  if(iwf.is_open())
    {
      cout<<"weight file "<<weight_file_name<<endl;
      for(int i=0;i<_number_of_events;i++)
        {
          iwf>>_weight[i];
        }
    }
  else
   {
      throw runtime_error("can not open weight file, please check your file and filename!");
    }

}

void DataObject::read_data()
{
  ifstream idf;
  idf.open(data_file_name,ios::in);
  cout<<"file"<<data_file_name<<endl;
  if (idf.is_open())
    {
      string line;
      int n = 0;
      while (getline(idf, line))
        {
          n++;
        }
      _number_of_events = n/5;
      cout<<" n evenet "<<_number_of_events<<endl;
      if(n%5 !=0)
        cerr<<"warning : the number of lines in data file is not a multiple of 5, please check your data file "<<endl;
      _mcp = new double[20*_number_of_events];
      idf.clear();
      idf.seekg(0,ios::beg);
      for(int i=0;i<5*_number_of_events;i++)
        {
          idf>>_mcp[i*4]>>_mcp[i*4+1]>>_mcp[i*4+2]>>_mcp[i*4+3];
        }
      _hpv = new double[_number_of_events*parameter_vector_index_end];
    }
  else
    {
      throw runtime_error("can not open data file, please check your file and filename!");
    }
}
  //cout << df << " have " << n << " lines." << endl;

void DataObject::read_pwa()
{
  for(int i=0;i<_number_of_events;i++)
    {
      double* hp = new double[parameter_vector_index_end];
      calculate0p(i,hp);
      store_parameters(i, hp);
    }
}
void DataObject::calculate0p(int id, double* hp)
{
    double ap23[4], apv2[4], apv3[4];
    double ak1[4], ak2[4], ak3[4], ak4[4], ak5[4];
    ak1[0] = _mcp[id*20+0];
    ak1[1] = _mcp[id*20+1];
    ak1[2] = _mcp[id*20+2];
    ak1[3] = _mcp[id*20+3];
    ak2[0] = _mcp[id*20+4];
    ak2[1] = _mcp[id*20+5];
    ak2[2] = _mcp[id*20+6];
    ak2[3] = _mcp[id*20+7];
    ak3[0] = _mcp[id*20+8];
    ak3[1] = _mcp[id*20+9];
    ak3[2] = _mcp[id*20+10];
    ak3[3] = _mcp[id*20+11];
    ak4[0] = _mcp[id*20+12];
    ak4[1] = _mcp[id*20+13];
    ak4[2] = _mcp[id*20+14];
    ak4[3] = _mcp[id*20+15];
    ak5[0] = _mcp[id*20+16];
    ak5[1] = _mcp[id*20+17];
    ak5[2] = _mcp[id*20+18];
    ak5[3] = _mcp[id*20+19];
    double(*fDel)[4], (*fGel)[4], (*E)[4][4][4], (*G1)[4][4][4],
        (*G3)[4][4][4][4][4];
    double t4wvf[4][4][4][4], ttfw[4][4][4];

    fDel = _dp->fDel;
    fGel = _dp->fGel;
    E = _dp->E;
    G1 = _dp->G1;
    G3 = _dp->G3;

//    // cout<<"haha: "<< __LINE__ << endl;
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            if (i == j) {
//                if (i < 3)
//                    fDel[i][j] = -1;
//                else
//                    fDel[i][j] = 1;
//            } else {
//                fDel[i][j] = 0;
//            }
//        }
//    }
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            if (i == j) {
//                if (i < 3)
//                    fGel[i][j] = -1;
//                else
//                    fGel[i][j] = 0;
//            } else {
//                fGel[i][j] = 0;
//            }
//        }
//    }
//    // cout<<"haha: "<< __LINE__ << endl;
//    E[0][1][2][3] = 1.0;
//    E[0][1][3][2] = -1.0;
//    E[0][2][1][3] = -1.0;
//    E[0][2][3][1] = 1.0;
//    E[0][3][1][2] = 1.0;
//    E[0][3][2][1] = -1.0;
//    E[1][0][2][3] = -1.0;
//    E[1][0][3][2] = 1.0;
//    E[1][2][0][3] = 1.0;
//    E[1][2][3][0] = -1.0;
//    E[1][3][0][2] = -1.0;
//    E[1][3][2][0] = 1.0;
//    E[2][0][1][3] = 1.0;
//    E[2][0][3][1] = -1.0;
//    E[2][1][0][3] = -1.0;
//    E[2][1][3][0] = 1.0;
//    E[2][3][0][1] = 1.0;
//    E[2][3][1][0] = -1.0;
//    E[3][0][1][2] = -1.0;
//    E[3][0][2][1] = 1.0;
//    E[3][1][0][2] = 1.0;
//    E[3][1][2][0] = -1.0;
//    E[3][2][0][1] = -1.0;
//    E[3][2][1][0] = 1.0;
//    for (int I = 0; I < 4; I++) {
//        for (int J = 0; J < 4; J++) {
//            for (int K = 0; K < 4; K++) {
//                for (int L = 0; L < 4; L++) {
//                    // G1[I][J][K][L] = (fGel[J][K]*fGel[I][L] +
//                    // fGel[K][I]*fGel[J][L])/2 -(fGel[I][J]*fGel[K][L])/3 ;
//                    G1[I][J][K][L] = fGel[I][J] * fGel[K][L] +
//                                     fGel[J][K] * fGel[I][L] +
//                                     fGel[K][I] * fGel[J][L];
//                    // last (...), i.e.
//                    // combination of \tilde{g} in
//                    // Eq.(37) of PRD48,1225
//                }
//            }
//        }
//    }  // Eq.(20)
//    for (int I1 = 0; I1 < 4; I1++) {
//        for (int I2 = 0; I2 < 4; I2++) {
//            for (int I3 = 0; I3 < 4; I3++) {
//                for (int I4 = 0; I4 < 4; I4++) {
//                    for (int I5 = 0; I5 < 4; I5++) {
//                        for (int I6 = 0; I6 < 4; I6++) {
//                            G3[I1][I2][I3][I4][I5][I6] =
//                                (fGel[I1][I2] * fGel[I4][I5] * fGel[I3][I6] +
//                                 fGel[I1][I2] * fGel[I5][I6] * fGel[I3][I4] +
//                                 fGel[I1][I2] * fGel[I4][I6] * fGel[I3][I5] +
//                                 fGel[I1][I3] * fGel[I4][I6] * fGel[I2][I5] +
//                                 fGel[I1][I3] * fGel[I4][I5] * fGel[I2][I6] +
//                                 fGel[I1][I3] * fGel[I5][I6] * fGel[I2][I4] +
//                                 fGel[I2][I3] * fGel[I5][I6] * fGel[I1][I4] +
//                                 fGel[I2][I3] * fGel[I4][I5] * fGel[I1][I6] +
//                                 fGel[I2][I3] * fGel[I4][I6] * fGel[I1][I5]) /
//                                    15.0 -
//                                (fGel[I1][I4] * fGel[I2][I5] * fGel[I3][I6] +
//                                 fGel[I1][I4] * fGel[I2][I6] * fGel[I3][I5] +
//                                 fGel[I1][I5] * fGel[I2][I4] * fGel[I3][I6] +
//                                 fGel[I1][I5] * fGel[I2][I6] * fGel[I3][I4] +
//                                 fGel[I1][I6] * fGel[I2][I5] * fGel[I3][I4] +
//                                 fGel[I1][I6] * fGel[I2][I4] * fGel[I3][I5]) /
//                                    6.0;
//                        }
//                    }
//                }
//            }
//        }
//    }  // Eq. (21)

    // 1 phi 2 pi+ 3 pi- 4 K+ 5 K-
    // ak1 ~ ak5, four momentum from input, ak1 is the sum of 3 and 4, etc. the
    // phi's momentum
    for (int i = 0; i < 4; i++) {
        ap23[i] = ak2[i] + ak3[i];  // pi+ + pi-
        apv2[i] = ak1[i] + ak2[i];  // phi + pi+
        apv3[i] = ak1[i] + ak3[i];  // phi + pi-
        vv(ak23, i) = ak2[i] - ak3[i];
        // r_34 (relative difference between pi+ and pi-)
        vv(ak45, i) = ak4[i] - ak5[i];
        // r_12 (relative difference between K+ and K-)
        vv(ar, i) = ak1[i] - ap23[i];
        // r^\mu (relative difference between phi and two pions)
        vv(ak23u, i) = vv(ak23, i) * fDel[i][i];
        vv(ak45u, i) = vv(ak45, i) * fDel[i][i];
        vv(aru, i) = vv(ar, i) * fDel[i][i];  // r_\mu(phi f)
    }
    ss(sv) = scalar(ak1, ak1);     // p^2(phi)
    ss(s23) = scalar(ap23, ap23);  // p^2(f)
    ss(sv2) = scalar(apv2, apv2);  //  p^2(123)
    ss(sv3) = scalar(apv3, apv3);  //   p^2(124)

    // 0+ contribution
    // prepare tensors for Eq.(36)
    ss(q2r45) = ss(sv) * 0.25 - (_dp->_m[0]) * (_dp->_m[0]);
    // q2r45 is Q^2_abc, Eq.(13), phi(a)-> K+(b) K-(c)
    ss(b1q2r45) = sqrt(2 / (ss(q2r45) + fFUD));
    // b1q2r45 is B1(Q_abc), Eq.(14), it can
    // be ignored for phi is narrow, i.e. it
    // should be always close to a constant
    // when final states' momenta change
    // cout << "b1q2r45 " << b1q2r45 << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mm(del45w, i, j) = fDel[i][j] - ak1[i] * ak1[j] / ss(sv);
            // cout << "del45w " << del45w[i][j] << endl;
        }
    }  // del45w is \tilde{g}^{\mu\nu}(\phi)

    for (int i = 0; i < 4; i++) {
        vv(wt, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(wt, i) = vv(wt, i) + mm(del45w, i, j) * vv(ak45u, j);
        }
        vv(wt, i) = vv(wt, i) * ss(b1q2r45);
        // Eq. (10), \tilde{r}^\mu = t^{(1)\mu}(12)
        vv(wu, i) = vv(wt, i) * fDel[i][i];
        // \tilde{r}^\mu = t^{(1)}_\mu(12)
    }
    // prepare tensors for Eq.(37)
    for (int i = 0; i < 4; i++) {
        vv(arw, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(arw, i) =
                vv(arw, i) + fGel[i][j] * vv(aru, j);  // \tilde(r)^\mu(phi f0)
        }
    }
    ss(arwdarw) = scalar(vp(arw), vp(arw));
    double tmp0 = (-ss(arwdarw)) / 3;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mm(t2wvf, i, j) = vv(arw, i) * vv(arw, j) + tmp0 * fGel[i][j];
            // \tilde{T}^{(2)\mu\nu}(phi f0), from Eq.(11)
            // cout << "t2wvf " << t2wvfuu[i][j] << endl;
        }
    }
    ss(qjvf2) =
        0.25 * (pow((mpsip * mpsip + ss(sv) - ss(s23)), 2)) / (mpsip * mpsip) -
        ss(sv);
    // Eq.(13), a is psi, b is phi, c is resonance
    ss(b2qjvf2) = sqrt(
        13. / (pow(ss(qjvf2), 2) + 3. * ss(qjvf2) * fFUD + 9. * pow(fFUD, 2)));
    // B2(Qabc) in Eq.(15), can
    // be neglected since psi is
    // very narrow
    // cout << "b2qjvf2 " << b2qjvf2 << endl;

    // 0+ contribution
    for (int i = 0; i < 2; i++) {
        vv(w0p22, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(w0p22, i) = vv(w0p22, i) + mm(t2wvf, i, j) * vv(wu, j);
            // the \tilde{T}^{(2)\mu\nv}(phi
            // f0) \tilde{t}^{(1)}_{\nu}(12)
            // in Eq.(37), only prepare \mu =
            // 1 and 2
        }
    }
    // for(int i=0; i<4; i++) cout << "w0pp2 " << w0p22[i] << endl;
    // Till now, the tensors in Eq.(36) [wt] and (37) [w0p22] have been
    // prepared.
    // The missing parts are propagators, etc. BW.

    //   2+ contribution
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mm(del23w, i, j) = fDel[i][j] - ap23[i] * ap23[j] / ss(s23);
            // \tilde{g}^{\mu\nu}(f2) in Eq.(11)
        }
    }
    for (int i = 0; i < 4; i++) {
        vv(ak23w, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(ak23w, i) = vv(ak23w, i) + mm(del23w, i, j) * vv(ak23u, j);
            // \tilde{r}^{\mu}(f2) in Eq.(11),
            // r(a)-> pi+(b) pi-(c)
        }
    }
    ss(qf2xx) = 0.25 * ss(s23) - (_dp->_m[2]) * (_dp->_m[2]);
    ss(b2qf2xx) = sqrt(
        13. / (pow(ss(qf2xx), 2) + 3. * ss(qf2xx) * fFUD + 9. * pow(fFUD, 2)));
    // B2 in Eq.(15) of pi+ pi-
    ss(b4qjvf2) =
        sqrt(12746. / (pow(ss(qjvf2), 4) + 10 * pow(ss(qjvf2), 3) * fFUD +
                       135 * pow(ss(qjvf2), 2) * pow(fFUD, 2) +
                       1575 * ss(qjvf2) * pow(fFUD, 3) + 11025 * pow(fFUD, 4)));
    // B4 in Eq.(17) of psi
    double tmpf2 = -scalar(vp(ak23w), vp(ak23w)) / 3;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mm(t2wf, i, j) =
                (vv(ak23w, i) * vv(ak23w, j) + tmpf2 * mm(del23w, i, j)) *
                ss(b2qf2xx);
            // \tilde{t}^{(2)\mu\nu}(f2) in Eq.(11),
            // \tilde{t}^{(2)\mu\nu}(34) in Eq.(41)
            mm(t2wfu, i, j) = mm(t2wf, i, j) * fDel[j][j];
            // \tilde{t}^{(2)\mu}_{\nu}(34)
            // cout << "t2wfu " << t2wfu[i][j] << endl;
        }
    }
    for (int i = 0; i < 4; i++) {
        vv(w2p1, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(w2p1, i) = vv(w2p1, i) + mm(t2wf, i, j) * vv(wu, j);
            // \tilde{t}^{(2)\mu\nu}(34)
            // \tilde{t}^{(1)}_{\nu}(12) the
            // second and third terms in Eq.(42)
        }
        vv(w2p1u, i) = vv(w2p1, i) * fDel[i][i];
        // \tilde{t}^{(2)}_{\mu\nu}(34) \tilde{t}^{(1)\nu}(12)
    }
    for (int i = 0; i < 2; i++) {
        vv(w2p2, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(w2p2, i) = vv(w2p2, i) + mm(t2wvf, i, j) * vv(w2p1u, j);
            // \tilde{T}^{(2)\mu\alpha}(phi
            // f2)\tilde{t}^{(2)}_{\alpha\nu}(34)
            // \tilde{t}^{(1)\nu}(12)， the first three
            // terms in Eq.(42)
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                ttfw[i][j][k] = mm(t2wfu, i, j) * vv(wt, k);
                // \tilde{t}^{(2)\lambda}_\delta(34)
                // \tilde{t}^{(1)\nu}(12)
            }
        }
    }
    //  cout<<"haha: "<< __LINE__ << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mm(t2p3, i, j) = 0.0;
            for (int k1 = 0; k1 < 3; k1++) {
                for (int k2 = 0; k2 < 3; k2++) {
                    mm(t2p3, i, j) =
                        mm(t2p3, i, j) + (E[i][k1][3][k2] * ttfw[k1][j][k2] +
                                          E[j][k1][3][k2] * ttfw[k1][i][k2]);
                    // [...]p^\sigma_psi
                    // \tilde{t}^{(1)\nu}(12)
                }
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        vv(w2p3, i) = 0.0;
        for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 3; j2++) {
                for (int j3 = 0; j3 < 4; j3++) {
                    vv(w2p3, i) =
                        vv(w2p3, i) +
                        E[i][j1][j2][3] * mm(t2wvf, j1, j3) * mm(t2p3, j2, j3);
                    // epsilon^{\mu\alpha\beta\gamma}p_{\psi\alpha}\tilde^{(2)\lambda}_{\beta}(\phi
                    // f2)[...] in Eq. (43)
                }
            }
        }
    }
    double test[4][4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                test[i][j][k] = 0;
                for (int m = 0; m < 4; m++) {
                    for (int n = 0; n < 4; n++) {
                        test[i][j][k] = test[i][j][k] +
                                        ttfw[i][m][n] * fDel[m][j] * fDel[n][k];
                    }
                }
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        vv(w2p4, i) = 0.0;
        for (int j1 = 0; j1 < 4; j1++) {
            for (int j2 = 0; j2 < 4; j2++) {
                for (int j3 = 0; j3 < 4; j3++) {
                    for (int j4 = 0; j4 < 4; j4++) {
                        for (int j5 = 0; j5 < 4; j5++) {
                            vv(w2p4, i) = vv(w2p4, i) +
                                          G3[i][j1][j2][j3][j4][j5] *
                                              mm(t2wvf, j1, j2) *
                                              ttfw[j3][j4][j5];
                            //  P^{(3)}\tilde{T}^{(2)}(\phi
                            //  f2)\tilde{t}^{(2)}(34)\tilde{t}^{(1)}(12)
                            //  Eq.(44)
                            // w2p4[i]=w2p4[i]+G3[i][j1][j2][j3][j4][j5]*t2wvf[j1][j2]*test[j3][j4][j5];
                        }
                    }
                }
            }
        }
    }
    double tmp1 = (-ss(arwdarw)) / 7;
    double tmp2 = (pow(ss(arwdarw), 2)) / 35;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                for (int l = 0; l < 4; l++) {
                    t4wvf[i][j][k][l] =
                        vv(arw, i) * vv(arw, j) * vv(arw, k) * vv(arw, l) +
                        tmp1 * (fGel[i][j] * vv(arw, k) * vv(arw, l) +
                                fGel[j][k] * vv(arw, i) * vv(arw, l) +
                                fGel[k][i] * vv(arw, j) * vv(arw, l) +
                                fGel[i][l] * vv(arw, j) * vv(arw, k) +
                                fGel[j][l] * vv(arw, k) * vv(arw, i) +
                                fGel[k][l] * vv(arw, i) * vv(arw, j)) +
                        tmp2 * G1[i][j][k][l];
                    // defined in Eq(37) of PRD48,1225
                }
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        vv(w2p5, i) = 0.0;
        for (int j1 = 0; j1 < 4; j1++) {
            for (int j2 = 0; j2 < 4; j2++) {
                for (int j3 = 0; j3 < 4; j3++) {
                    vv(w2p5, i) =
                        vv(w2p5, i) + t4wvf[i][j1][j2][j3] * ttfw[j2][j3][j1];
                    // \tilde{T}^{(4)}(\phi
                    // f2)
                    // \tilde{t}^{(1)}(12)
                    // \tilde{t}^{(2)}(34),
                    // in Eq.(45)
                }
            }
        }
    }
    // 1+ and 1- contribution
    ss(qjv2) = 0.25 * pow((_dp->_M2 + ss(sv2) - (_dp->_m[2]) * (_dp->_m[2])), 2) /
                   _dp->_M2 -
               ss(sv2);
    // Q^2(abc) = Q^2(\psi' 123 \pi4)
    ss(qjv3) = 0.25 * pow((_dp->_M2 + ss(sv3) - (_dp->_m[2]) * (_dp->_m[2])), 2) /
                   _dp->_M2 -
               ss(sv3);
    // Q^2(\psi 124 \pi3)
    ss(qbv2) = 0.25 * pow((ss(sv2) + ss(sv) - (_dp->_m[2]) * (_dp->_m[2])), 2) /
                   ss(sv2) -
               ss(sv);
    // Q^2(123 12 \pi3)
    ss(qbv3) = 0.25 * pow((ss(sv3) + ss(sv) - (_dp->_m[2]) * (_dp->_m[2])), 2) /
                   ss(sv3) -
               ss(sv);
    // Q^2(124 12 \pi4)
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mm(delv2w, i, j) = fDel[i][j] - apv2[i] * apv2[j] / ss(sv2);
            // \tilde{g}^{\mu\nu}(123)
            mm(delv3w, i, j) = fDel[i][j] - apv3[i] * apv3[j] / ss(sv3);
            // \tilde{g}^{\mu\nu}(124)
        }
    }
    for (int i = 0; i < 4; i++) {
        vv(ak12, i) = ak1[i] - ak2[i];     // r^\mu(phi pi3) i.e. p_{12}-p_3
        vv(ak13, i) = ak1[i] - ak3[i];     // r^\mu(phi pi4) i.e. p_{12}-p_4
        vv(akv2m3, i) = apv2[i] - ak3[i];  // r^\mu(rho' pi4) i.e. p_{123}-p_4
        vv(akv3m2, i) = apv3[i] - ak2[i];  // r^\mu(rho' pi3) i.e. p_{124}-p_3
        vv(ak12u, i) = vv(ak12, i) * fDel[i][i];      // r_\mu(phi pi3)
        vv(ak13u, i) = vv(ak13, i) * fDel[i][i];      // r_\mu(phi pi4)
        vv(akv2m3u, i) = vv(akv2m3, i) * fDel[i][i];  // r_\mu(rho' pi4)
        vv(akv3m2u, i) = vv(akv3m2, i) * fDel[i][i];  // r_\mu(rho' pi3)
    }
    for (int i = 0; i < 4; i++) {
        vv(ak12w, i) = 0.0;
        vv(ak13w, i) = 0.0;
        vv(akv2m3w, i) = 0.0;
        vv(akv3m2w, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(ak12w, i) = vv(ak12w, i) + mm(delv2w, i, j) * vv(ak12u, j);
            // \tilde{r}^\mu (phi pi3)
            vv(ak13w, i) = vv(ak13w, i) + mm(delv3w, i, j) * vv(ak13u, j);
            // \tilde{r}^\mu (phi pi4)
            vv(akv2m3w, i) = vv(akv2m3w, i) + fGel[i][j] * vv(akv2m3u, j);
            // \tilde{r}^\mu (b pi4)
            vv(akv3m2w, i) = vv(akv3m2w, i) + fGel[i][j] * vv(akv3m2u, j);
            // \tilde{r}^\mu (b pi3)
        }
    }
    // 1+ contribution
    ss(b2qjv2) = sqrt(
        13. / (pow(ss(qjv2), 2) + 3. * ss(qjv2) * fFUD + 9. * pow(fFUD, 2)));
    // B2(Q psi b pi4)
    ss(b2qjv3) = sqrt(
        13. / (pow(ss(qjv3), 2) + 3. * ss(qjv3) * fFUD + 9. * pow(fFUD, 2)));
    // B2(Q psi b pi3)
    ss(b2qbv2) = sqrt(
        13. / (pow(ss(qbv2), 2) + 3. * ss(qbv2) * fFUD + 9. * pow(fFUD, 2)));
    // B2(Q b phi pi3)
    ss(b2qbv3) = sqrt(
        13. / (pow(ss(qbv3), 2) + 3. * ss(qbv3) * fFUD + 9. * pow(fFUD, 2)));
    // B2(Q b phi pi4)
    double tmp12w = -scalar(vp(ak12w), vp(ak12w)) / 3.0;
    double tmp13w = -scalar(vp(ak13w), vp(ak13w)) / 3.0;
    double tmpv2m3 = -scalar(vp(akv2m3w), vp(akv2m3w)) / 3.0;
    double tmpv3m2 = -scalar(vp(akv3m2w), vp(akv3m2w)) / 3.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mm(t2v2, i, j) =
                (vv(ak12w, i) * vv(ak12w, j) + tmp12w * mm(delv2w, i, j)) *
                ss(b2qbv2);
            // \tilde{t}^{(2)\mu\nu}(phi pi3)
            mm(t2v3, i, j) =
                (vv(ak13w, i) * vv(ak13w, j) + tmp13w * mm(delv3w, i, j)) *
                ss(b2qbv3);
            // \tilde{t}^{(2)\mu\nu}(phi pi4)
            mm(t2b3, i, j) =
                (vv(akv2m3w, i) * vv(akv2m3w, j) + tmpv2m3 * fGel[i][j]) *
                ss(b2qjv2);
            // \tilde{T}^{(2)\mu\nu}(b pi4)
            mm(t2b2, i, j) =
                (vv(akv3m2w, i) * vv(akv3m2w, j) + tmpv3m2 * fGel[i][j]) *
                ss(b2qjv3);
            // \tilde{T}^{(2)\mu\nu}(b pi3)
        }
    }
    for (int i = 0; i < 4; i++) {
        vv(w1p12_1, i) = 0.0;
        vv(w1p13_1, i) = 0.0;
        vv(w1p12_2, i) = 0.0;
        vv(w1p13_2, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(w1p12_1, i) = vv(w1p12_1, i) + mm(delv2w, i, j) * vv(wu, j);
            // \tilde{g}^{\mu\nu}(123)
            // \tilde{t}^{(1)}_\nu(12) Eq.(47)
            vv(w1p13_1, i) = vv(w1p13_1, i) + mm(delv3w, i, j) * vv(wu, j);
            // \tilde{g}^{\mu\nu}(124)
            // \tilde{t}^{(1)}_\nu(12) Eq.(47)
            vv(w1p12_2, i) = vv(w1p12_2, i) + mm(t2v2, i, j) * vv(wu, j);
            // \tilde{t}_{(2)\mu\nu}(phi 3)
            // \tilde{t}^{(1\nu}(12) Eq.(48)
            vv(w1p13_2, i) = vv(w1p13_2, i) + mm(t2v3, i, j) * vv(wu, j);
            // \tilde{t}^{(2)\mu\nu}(phi 4)
            // \tilde{t}^{(1)\nu}(12) Eq.(48)
        }
        // cout << "w1p12_1 " << w1p12_1[i] << "  w1p13_1 " << w1p13_1[i] << "
        // w1p12_2 " << w1p12_2[i] << "  w1p13_2 " << w1p13_2[i] << endl;
        vv(w1p12_1u, i) = vv(w1p12_1, i) * fDel[i][i];
        vv(w1p13_1u, i) = vv(w1p13_1, i) * fDel[i][i];
        vv(w1p12_2u, i) = vv(w1p12_2, i) * fDel[i][i];
        vv(w1p13_2u, i) = vv(w1p13_2, i) * fDel[i][i];
    }
    for (int i = 0; i < 2; i++) {
        vv(w1p12_3, i) = 0.0;
        vv(w1p13_3, i) = 0.0;
        vv(w1p12_4, i) = 0.0;
        vv(w1p13_4, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(w1p12_3, i) = vv(w1p12_3, i) + mm(t2b3, i, j) * vv(w1p12_1u, j);
            // \tilde{T}^{(2)\mu\lambda}_{b1 4}
            // \tilde{g}_{(123)\lambda \nu}
            // \tilde{t}^{(1)\nu}_{(12)} in Eq.(49)
            vv(w1p13_3, i) = vv(w1p13_3, i) + mm(t2b2, i, j) * vv(w1p13_1u, j);
            // \tilde{T}^{(2)\mu\lambda}_{b1 3}
            // \tilde{g}_{(124)\lambda \nu}
            // \tilde{t}^{(1)\nu}_{(12)} in Eq.(49)
            vv(w1p12_4, i) = vv(w1p12_4, i) + mm(t2b3, i, j) * vv(w1p12_2u, j);
            // \tilde{T}^{(2)\mu\lambda}_{b1 4}
            // \tilde{t}_{(phi 3)\lambda \nu}
            // \tilde{t}^{(1)\nu}_{(12)} in Eq.(50)
            vv(w1p13_4, i) = vv(w1p13_4, i) + mm(t2b2, i, j) * vv(w1p13_2u, j);
            // \tilde{T}^{(2)\mu\lambda}_{b1 3}
            // \tilde{g}_{(phi 4)\lambda \nu}
            // \tilde{t}^{(1)\nu}_{(12)} in Eq.(50)
        }
        // cout << "w1p12_3 " << w1p12_3[i] << "  w1p13_3 " << w1p13_3[i] << "
        // w1p12_4 " << w1p12_4[i] << "  w1p13_4 " << w1p13_4[i] << endl;
    }
    //  1- contribution
    ss(b1qjv2) = sqrt(2. / (ss(qjv2) + fFUD));  // B1(Q psi rho' pi4), rho'(123)
    ss(b1qjv3) = sqrt(2. / (ss(qjv3) + fFUD));  // B1(Q psi rho' pi3), rho'(124)
    ss(b1qbv2) = sqrt(2. / (ss(qbv2) + fFUD));  // B1(Q rho' phi pi3), rho'(123)
    ss(b1qbv3) = sqrt(2. / (ss(qbv3) + fFUD));  // B1(Q rho' phi pi4), rho'(124)
    for (int i = 0; i < 2; i++) {
        vv(w1m12, i) = 0.0;
        vv(w1m13, i) = 0.0;
        for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 3; j2++) {
                for (int j3 = 0; j3 < 3; j3++) {
                    for (int j4 = 0; j4 < 3; j4++) {
                        double tmp = E[i][j1][j2][3] * E[j2][j3][j4][3];
                        vv(w1m12, i) = vv(w1m12, i) +
                                       tmp * vv(akv2m3w, j1) * vv(ak12w, j3) *
                                           vv(wt, j4) * ss(b1qbv2);
                        // The second term in Eq.(46)
                        vv(w1m13, i) = vv(w1m13, i) +
                                       tmp * vv(akv3m2w, j1) * vv(ak13w, j3) *
                                           vv(wt, j4) * ss(b1qbv3);
                        // The first term in Eq.(46)
                    }
                }
            }
        }
    }
    // phi(1650)f0(980)  contribution // mass 1680 +- 20 MeV, width 150 +- 50
    // MeV, dominant decay KK*(892)
    ss(q2r23) = 0.25 * ss(s23) - (_dp->_m[2]) * (_dp->_m[2]);
    // Q^2_{abc}(f->pi+ pi-)
    ss(b1q2r23) = sqrt(2. / (ss(q2r23) + fFUD));
    // B1(Q_abc) = B1(Q f pi+ pi-)
    for (int i = 0; i < 4; i++) {
        vv(ak23wu, i) =
            vv(ak23w, i) * fDel[i][i] * ss(b1q2r23);  // \tilde{t}_mu(34)
    }
    for (int i = 0; i < 2; i++) {
        vv(wpf22, i) = 0.0;
        for (int j = 0; j < 4; j++) {
            vv(wpf22, i) = vv(wpf22, i) - mm(t2wvf, i, j) * vv(ak23wu, j);
            // -\tilde{T}^{(2)\mu\nu}(phi f0) \tilde{t}_nu(34)
        }
    }
}
void DataObject::store_parameters(int id, const double* hp) {
    for (int i = parameter_vector_index_start; i < parameter_vector_index_end;i++)
      {
        _hpv[i * _number_of_events + id] = hp[i];
      }
}

double DataObject::scalar(double* a1, double* a2) const
{
  double scal=0;

  for(int i=0;i<4;i++)
    {
      scal+=a1[i]*a2[i]*_dp->fDel[i][i];
    }
  return scal;
}
double* DataObject::Get_mcp()
{
  return _mcp;
}

int DataObject::Get_number_of_events()
{
  return _number_of_events;
}

double* DataObject::Get_weight()
{
  return _weight;
}

double* DataObject::Get_hpv()
{
  return _hpv;
}
ostream &operator <<(ostream &out, DataObject &d)
{
  for(int i=0;i<5*d._number_of_events;i++)
    {
      out<<d._mcp[i*4]<<" "<<d._mcp[i*4+1]<<" "<<d._mcp[i*4+2]<<" "<<d._mcp[i*4+3]<<endl;
    }
  return out;
}
