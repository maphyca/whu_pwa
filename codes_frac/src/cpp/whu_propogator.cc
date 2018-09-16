#include "whu_propogator.h"
//#include "DPFAngular.h"
#include "TComplex.h"
#include "whu_constants_and_definitions.h"
// const double rk=0.493677;
// const double rp=0.139556995;

TComplex CPUWaveFunc::cro(double sx, double am1, double am2) const {
    TComplex ci(0, 1);
    double t1 = (am1 + am2) * (am1 + am2);  // double t1=pow((am1+am2),2);
    double t2 = (am1 - am2) * (am1 - am2);  // double t2=pow((am1-am2),2);
    double st = (sx - t1) * (sx - t2);
    double cro = sqrt(fabs(st)) / sx;
    TComplex result = cro;
    if (st < 0.) result = cro * ci;
    return result;
}
TComplex CPUWaveFunc::propogator980(double mass, double g11, double g22,
                                    double sx) const {
    TComplex ci(0, 1);
    double rm = mass * mass;
    TComplex propogator980 =
        1.0 / (rm - sx - ci * (g11 * cro(sx, rp, rp) + g22 * cro(sx, rk, rk)));
    return propogator980;
}
TComplex CPUWaveFunc::pip(double sx) const {
    TComplex ci(0, 1);
    double xk2 = sx - 0.3116676;  // 0.3116676=16.*0.139568*0.139568
    if (xk2 <= 0.) xk2 = 0.0;
    double r4pip = sqrt(xk2 / sx) / (1.0 + exp(9.8 - 3.5 * sx));  // 9.8=3.5*2.8
    return r4pip;
}
TComplex CPUWaveFunc::propogator(double mass, double width, double sx) const {
    TComplex ci(0, 1);
    double am = mass;
    double g1 = mass * width;
    TComplex prop = g1 / (sx - am * am + ci * g1);
    // TComplex prop=g1/(sx-pow(am,2)+ci*g1);
    return prop;
}
TComplex CPUWaveFunc::propogator1270(double mass, double width,
                                     double sx) const {
    TComplex ci(0, 1);
    double rm = mass * mass;
    double gr = mass * width;
    double q2r = 0.25 * rm - 0.0194792;
    double b2r = q2r * (q2r + 0.1825) + 0.033306;
    double g11270 = gr * b2r / pow(q2r, 2.5);
    double q2 = 0.25 * sx - 0.0194792;
    double b2 = q2 * (q2 + 0.1825) + 0.033306;
    double g1 = g11270 * pow(q2, 2.5) / b2;
    TComplex prop = gr / (sx - rm + ci * g1);
    return prop;
}
void CPUWaveFunc::cpu_propogator1(double mass, double width, double *sx,
                                  double *b2qjvf2, double *wu, double *w0p22,
                                  TComplex *fCF0, TComplex *fCF1,
                                  int vec_size) {
    for (int i = 0; i < vec_size; i++) {
        TComplex crp1 = propogator(mass, width, sx[i]);
        TComplex cr0p11 = crp1 / b2qjvf2[i];

        // 01 contribution
        fCF0[i + vec_size * 0] = wu[i + vec_size * 0] * crp1;
        fCF0[i + vec_size * 1] = wu[i + vec_size * 1] * crp1;

        // 02 contribution
        fCF1[i + vec_size * 0] = w0p22[i + vec_size * 0] * cr0p11;
        fCF1[i + vec_size * 1] = w0p22[i + vec_size * 1] * cr0p11;
    }
}
void CPUWaveFunc::cpu_propogator2(double mass, double g11, double g22,
                                  double *sx, double *b2qjvf2, double *wu,
                                  double *w0p22, TComplex *fCF0, TComplex *fCF1,
                                  int vec_size) {
    for (int i = 0; i < vec_size; i++) {
        TComplex crp1 = propogator980(mass, g11, g22, sx[i]);
        TComplex cr0p11 = crp1 / b2qjvf2[i];

        // 01 contribution
        fCF0[i + vec_size * 0] = wu[i + vec_size * 0] * crp1;
        fCF0[i + vec_size * 1] = wu[i + vec_size * 1] * crp1;

        // 02 contribution
        fCF1[i + vec_size * 0] = w0p22[i + vec_size * 0] * cr0p11;
        fCF1[i + vec_size * 1] = w0p22[i + vec_size * 1] * cr0p11;
    }
}
void CPUWaveFunc::cpu_propogator7(double mass, double width, double *sv2,
                                  double *sv3, double *b1qjv2, double *b1qbv2,
                                  double *b1qjv3, double *b1qbv3, double *w1m12,
                                  double *w1m13, TComplex *fCF, int vec_size) {
    for (int i = 0; i < vec_size; i++) {
        TComplex crp1 = propogator(mass, width, sv2[i]);
        TComplex crp11 = propogator(mass, width, sv3[i]);
        TComplex cr1m12_1 = crp1 / b1qjv2[i] / b1qbv2[i];
        TComplex cr1m13_1 = crp11 / b1qjv3[i] / b1qbv3[i];

        // 1-__1 contribution
        fCF[i + vec_size * 0] = w1m12[i + vec_size * 0] * cr1m12_1 +
                                w1m13[i + vec_size * 0] * cr1m13_1;
        fCF[i + vec_size * 1] = w1m12[i + vec_size * 1] * cr1m12_1 +
                                w1m13[i + vec_size * 1] * cr1m13_1;
    }
}
void CPUWaveFunc::cpu_propogator8(
    double mass, double width, double *sv2, double *sv3, double *b2qbv2,
    double *b2qbv3, double *b2qjv2, double *b2qjv3, double *w1p12_1,
    double *w1p13_1, double *w1p12_2, double *w1p13_2, double *w1p12_3,
    double *w1p13_3, double *w1p12_4, double *w1p13_4, TComplex *fCF0,
    TComplex *fCF1, TComplex *fCF2, TComplex *fCF3, int vec_size) {
    for (int i = 0; i < vec_size; i++) {
        TComplex crp1 = propogator(mass, width, sv2[i]);
        TComplex crp11 = propogator(mass, width, sv3[i]);
        TComplex c1p12_12 = crp1 / b2qbv2[i];
        TComplex c1p13_12 = crp11 / b2qbv3[i];
        TComplex c1p12_13 = crp1 / b2qjv2[i];
        TComplex c1p13_13 = crp11 / b2qjv3[i];
        TComplex c1p12_14 = c1p12_12 / b2qjv2[i];
        TComplex c1p13_14 = c1p13_12 / b2qjv3[i];

        // z 1+ 1
        fCF0[i + vec_size * 0] = w1p12_1[i + vec_size * 0] * crp1 +
                                 w1p13_1[i + vec_size * 0] * crp11;
        fCF0[i + vec_size * 1] = w1p12_1[i + vec_size * 1] * crp1 +
                                 w1p13_1[i + vec_size * 1] * crp11;

        // z 1+ 2
        fCF1[i + vec_size * 0] = w1p12_2[i + vec_size * 0] * c1p12_12 +
                                 w1p13_2[i + vec_size * 0] * c1p13_12;
        fCF1[i + vec_size * 1] = w1p12_2[i + vec_size * 1] * c1p12_12 +
                                 w1p13_2[i + vec_size * 1] * c1p13_12;

        // z 1+ 3
        fCF2[i + vec_size * 0] = w1p12_3[i + vec_size * 0] * c1p12_13 +
                                 w1p13_3[i + vec_size * 0] * c1p13_13;
        fCF2[i + vec_size * 1] = w1p12_3[i + vec_size * 1] * c1p12_13 +
                                 w1p13_3[i + vec_size * 1] * c1p13_13;

        // z 1+ 4
        fCF3[i + vec_size * 0] = w1p12_4[i + vec_size * 0] * c1p12_14 +
                                 w1p13_4[i + vec_size * 0] * c1p13_14;
        fCF3[i + vec_size * 1] = w1p12_4[i + vec_size * 1] * c1p12_14 +
                                 w1p13_4[i + vec_size * 1] * c1p13_14;
    }
}
void CPUWaveFunc::cpu_propogator6(double mass, double width, double *sx,
                                  double *b2qf2xx, double *b2qjvf2,
                                  double *b4qjvf2, double *w2p1, double *w2p2,
                                  double *w2p3, double *w2p4, double *w2p5,
                                  TComplex *fCF0, TComplex *fCF1,
                                  TComplex *fCF2, TComplex *fCF3,
                                  TComplex *fCF4, int vec_size) {
    for (int i = 0; i < vec_size; i++) {
        TComplex crp1 = propogator1270(mass, width, sx[i]);
        TComplex cw2p11 = crp1 / b2qf2xx[i];
        TComplex cw2p12 = cw2p11 / b2qjvf2[i];
        TComplex cw2p15 = cw2p11 / b4qjvf2[i];

        // 21 contribution
        fCF0[i + vec_size * 0] = w2p1[i + vec_size * 0] * cw2p11;
        fCF0[i + vec_size * 1] = w2p1[i + vec_size * 1] * cw2p11;

        // 22 contribution
        fCF1[i + vec_size * 0] = w2p2[i + vec_size * 0] * cw2p12;
        fCF1[i + vec_size * 1] = w2p2[i + vec_size * 1] * cw2p12;

        // 23 contribution
        fCF2[i + vec_size * 0] = w2p3[i + vec_size * 0] * cw2p12;
        fCF2[i + vec_size * 1] = w2p3[i + vec_size * 1] * cw2p12;

        // 24 contribution
        fCF3[i + vec_size * 0] = w2p4[i + vec_size * 0] * cw2p12;
        fCF3[i + vec_size * 1] = w2p4[i + vec_size * 1] * cw2p12;

        // 25 contribution
        fCF4[i + vec_size * 0] = w2p5[i + vec_size * 0] * cw2p15;
        fCF4[i + vec_size * 1] = w2p5[i + vec_size * 1] * cw2p15;
    }
}

void CPUWaveFunc::parameters_vector_resize() {
    whether_hold_flag.resize(parameter_vector_index_end, false);
    hpv.resize(parameter_vector_index_end, NULL);
    //hp.resize(parameter_vector_index_end, 0);
    for (int i = parameter_vector_index_start; i < parameter_vector_index_end;
         i++) {
        if (i > one_dimensional_start && i < one_dimensional_end) {
            hpv[i] = new double[number_of_events_];
        }
        if (i > two_dimensional_start && i < two_dimensional_end &&
            (i - two_dimensional_start) % 2 == 1) {
            hpv[i] = new double[number_of_events_ * 2];
            hpv[i + 1] = hpv[i] + number_of_events_;
        }
        if (i > four_dimensional_start && i < four_dimensional_end &&
            (i - four_dimensional_start) % 4 == 1) {
            hpv[i] = new double[number_of_events_ * 4];
            hpv[i + 1] = hpv[i] + number_of_events_ * 1;
            hpv[i + 2] = hpv[i] + number_of_events_ * 2;
            hpv[i + 3] = hpv[i] + number_of_events_ * 3;
        }
    }

    // wu.resize(number_of_events_, vector<double>(4));
    // w0p22.resize(number_of_events_, vector<double>(4));
    // ak23w.resize(number_of_events_, vector<double>(4));
    // w2p2.resize(number_of_events_, vector<double>(4));
    // w2p1.resize(number_of_events_, vector<double>(4));
    // w2p3.resize(number_of_events_, vector<double>(4));
    // w2p4.resize(number_of_events_, vector<double>(4));
    // w2p5.resize(number_of_events_, vector<double>(4));
    // b2qf2xx.resize(number_of_events_, 0);
    // b4qjvf2.resize(number_of_events_, 0);
    // b2qjv2.resize(number_of_events_, 0);
    // b2qjv3.resize(number_of_events_, 0);
    // b2qbv2.resize(number_of_events_, 0);
    // b2qbv3.resize(number_of_events_, 0);
    // b1qjv2.resize(number_of_events_, 0);
    // b1qjv3.resize(number_of_events_, 0);
    // b1qbv2.resize(number_of_events_, 0);
    // b1qbv3.resize(number_of_events_, 0);
    // b1q2r23.resize(number_of_events_, 0);
    // sv.resize(number_of_events_, 0);
    // s23.resize(number_of_events_, 0);
    // sv2.resize(number_of_events_, 0);
    // sv3.resize(number_of_events_, 0);
    // wpf22.resize(number_of_events_, vector<double>(2));
    // b2qjvf2.resize(number_of_events_, 0);
    // w1p12_1.resize(number_of_events_, vector<double>(4));
    // w1p13_1.resize(number_of_events_, vector<double>(4));
    // w1p12_2.resize(number_of_events_, vector<double>(4));
    // w1p13_2.resize(number_of_events_, vector<double>(4));
    // w1p12_3.resize(number_of_events_, vector<double>(2));
    // w1p13_3.resize(number_of_events_, vector<double>(2));
    // w1p12_4.resize(number_of_events_, vector<double>(2));
    // w1p13_4.resize(number_of_events_, vector<double>(2));
    // w1m12.resize(number_of_events_, vector<double>(2));
    // w1m13.resize(number_of_events_, vector<double>(2));
}
void CPUWaveFunc::cpu_convert_mcp_to_pwa_paras() {
    pwa_paras_.resize(0);
    for (int i = 0; i < number_of_events_; i++) {
        vector<double> hp(parameter_vector_index_end);
        cpu_calculate0p(i, hp);
        store_parameters(i, hp);
    }
}
// void CPUWaveFunc::convert_mcp_to_pwa_paras() {
//    DPFAngular _amp;
//    _amp.setdp(_dp);
//    pwa_paras_.resize(0);
//    for (int i = 0; i < number_of_events_; i++) {
//        PWA_PARAS _the_pwa_paras;
//        _amp.calculate0p(mcp1[i][0], mcp1[i][1], mcp1[i][2], mcp1[i][3],
//                         mcp2[i][0], mcp2[i][1], mcp2[i][2], mcp2[i][3],
//                         mcp3[i][0], mcp3[i][1], mcp3[i][2], mcp3[i][3],
//                         mcp4[i][0], mcp4[i][1], mcp4[i][2], mcp4[i][3],
//                         mcp5[i][0], mcp5[i][1], mcp5[i][2], mcp5[i][3],
//                         _the_pwa_paras);
//        pwa_paras_.push_back(_the_pwa_paras);
//    }
//}
bool CPUWaveFunc::check_integraty() {
    for (int i = 0; i < number_of_events_; i++) {
        // if (pwa_paras_[i].b2qjvf2 != b2qjvf2[i]) return false;
        // if (pwa_paras_[i].w0p22[0] != w0p22[i][0]) return false;
        // if (pwa_paras_[i].w0p22[1] != w0p22[i][1]) return false;
        // if (pwa_paras_[i].ak23w[1] != ak23w[i][1]) return false;
        // if (pwa_paras_[i].w1m12[1] != w1m12[i][1]) return false;
        // if (pwa_paras_[i].w1m13[1] != w1m13[i][1]) return false;
        // if (pwa_paras_[i].w1m13[0] != w1m13[i][0]) return false;
    }
    return true;
}
void CPUWaveFunc::cpu_resize_intermediate_variables(int number_of_amplitudes) {
    fCP.resize(number_of_amplitudes);
    // crp1.resize(number_of_amplitudes);
    // crp11.resize(number_of_amplitudes);
    fCF.resize(number_of_amplitudes);
    for (int i = 0; i < number_of_amplitudes; i++) {
        // crp1[i].resize(number_of_events_);
        // crp11[i].resize(number_of_events_);
        // fCF[i].resize(number_of_events_, vector<TComplex>(4));
        fCF[i] = new TComplex[number_of_events_ * 4];
    }
}
double CPUWaveFunc::cpu_calEva(const vector<double> &par,
                               vector<double> &par_back,
                               int number_of_amplitudes) {
    // cout << "number_of_amplitudes = " << number_of_amplitudes << endl;
    // cout << "number_of_events_ = " << number_of_events_ << endl;
    int i = 0;
    while (i < number_of_amplitudes) {
        int propType_now = par[end_category * i + propType_category];
        cout << "amplitude i = " << i << endl;
        cout << "propType_now = " << propType_now << endl;
        switch (propType_now) {
            case 1:  // f0
            {
                double mass0 = par[end_category * i + mass_category];
                double width0 = par[end_category * i + width_category];
                bool _not_changed =
                    ((mass0 == par_back[end_category * i + mass_category]) &&
                     (width0 == par_back[end_category * i + width_category]));
                if (!_not_changed) {
                    cout << "prop = " << propType_now << " : work one time!!!"
                         << endl;
                    cpu_propogator1(mass0, width0, hpv[s23], hpv[b2qjvf2],
                                    hpv[wu], hpv[w0p22], fCF[i], fCF[i + 1],
                                    number_of_events_);
                }
                i = i + 2;
            } break;
            //	Flatte   Propagator Contribution
            case 2:  // f0 980
            {
                double mass980 = par[end_category * i + mass_category];
                double g10 = par[end_category * i + g1_category];
                double g20 = par[end_category * i + g2_category];
                bool _not_changed =
                    ((mass980 == par_back[end_category * i + mass_category]) &&
                     (g10 == par_back[end_category * i + g1_category]) &&
                     (g20 == par_back[end_category * i + g2_category]));
                if (!_not_changed) {
                    cout << "prop = " << propType_now << " : work one time!!!"
                         << endl;
                    cpu_propogator2(mass980, g10, g20, hpv[s23], hpv[b2qjvf2],
                                    hpv[w0p22], hpv[wu], fCF[i], fCF[i + 1],
                                    number_of_events_);
                }
                i = i + 2;
            } break;
            case 7:  // 1m1800
            {
                double mass0 = par[end_category * i + mass_category];
                double width0 = par[end_category * i + width_category];
                bool _not_changed =
                    ((mass0 == par_back[end_category * i + mass_category]) &&
                     (width0 == par_back[end_category * i + width_category]));
                if (!_not_changed) {
                    cout << "prop = " << propType_now << " : work one time!!!"
                         << endl;
                    cpu_propogator7(mass0, width0, hpv[sv2], hpv[sv3],
                                    hpv[b1qjv2], hpv[b1qbv2], hpv[b1qjv3],
                                    hpv[b1qbv3], hpv[w1m12], hpv[w1m13], fCF[i],
                                    number_of_events_);
                }
                i = i + 1;
            } break;
            case 8:  // 1p1800
            {
                double mass0 = par[end_category * i + mass_category];
                double width0 = par[end_category * i + width_category];
                bool _not_changed =
                    ((mass0 == par_back[end_category * i + mass_category]) &&
                     (width0 == par_back[end_category * i + width_category]));
                if (!_not_changed) {
                    cout << "prop = " << propType_now << " : work one time!!!"
                         << endl;
                    cpu_propogator8(
                        mass0, width0, hpv[sv2], hpv[sv3], hpv[b2qbv2],
                        hpv[b2qbv3], hpv[b2qjv2], hpv[b2qjv3], hpv[w1p12_1],
                        hpv[w1p13_1], hpv[w1p12_2], hpv[w1p13_2], hpv[w1p12_3],
                        hpv[w1p13_3], hpv[w1p12_4], hpv[w1p13_4], fCF[i],
                        fCF[i + 1], fCF[i + 2], fCF[i + 3], number_of_events_);
                }
                i = i + 4;
            } break;
            case 6:  // f2
            {
                double mass0 = par[end_category * i + mass_category];
                double width0 = par[end_category * i + width_category];
                bool _not_changed =
                    ((mass0 == par_back[end_category * i + mass_category]) &&
                     (width0 == par_back[end_category * i + width_category]));
                if (!_not_changed) {
                    cout << "prop = " << propType_now << " : work one time!!!"
                         << endl;
                    cpu_propogator6(mass0, width0, hpv[s23], hpv[b2qf2xx],
                                    hpv[b2qjvf2], hpv[b4qjvf2], hpv[w2p1],
                                    hpv[w2p2], hpv[w2p3], hpv[w2p4], hpv[w2p5],
                                    fCF[i], fCF[i + 1], fCF[i + 2], fCF[i + 3],
                                    fCF[i + 4], number_of_events_);
                }
                i = i + 5;
            } break;
            default:
                cout << "Do not know how to deal with prop type "
                     << propType_now << endl;
                exit(1);
                ;
        }
    }

    for (int i = 0; i < number_of_amplitudes; i++) {
        double rho0 = par[end_category * i + rho_category];
        double frac0 = par[end_category * i + frac_category];
        double phi0 = par[end_category * i + phi_category];
        rho0 *= TMath::Exp(frac0);
        fCP[i] = TComplex(rho0 * TMath::Cos(phi0), rho0 * TMath::Sin(phi0));
    }

    return 0;
}

double CPUWaveFunc::sum_phsp(int number_of_amplitudes) {
    double _sum(0);
    for (int _event = 0; _event < number_of_events_; _event++) {
        TComplex cw1(0, 0), cw2(0, 0);
        for (int i = 0; i < number_of_amplitudes; i++) {
            cw1 = cw1 + fCP[i] * fCF[i][_event + number_of_events_ * 0];
            cw2 = cw2 + fCP[i] * fCF[i][_event + number_of_events_ * 1];
        }
        _sum += (cw1.Re() * cw1.Re() + cw1.Im() * cw1.Im() +
                 cw2.Re() * cw2.Re() + cw2.Im() * cw2.Im()) /
                2.0;
    }
    return _sum;
}
double CPUWaveFunc::sum_likelihood(int number_of_amplitudes) {
    double _sum(0);
    for (int _event = 0; _event < number_of_events_; _event++) {
        TComplex cw1(0, 0), cw2(0, 0);
        for (int i = 0; i < number_of_amplitudes; i++) {
            cw1 = cw1 + fCP[i] * fCF[i][_event + number_of_events_ * 0];
            cw2 = cw2 + fCP[i] * fCF[i][_event + number_of_events_ * 1];
        }
        _sum += -log((cw1.Re() * cw1.Re() + cw1.Im() * cw1.Im() +
                      cw2.Re() * cw2.Re() + cw2.Im() * cw2.Im()) /
                     2.0);
    }
    return _sum;
}
double CPUWaveFunc::sum_penalty(int number_of_amplitudes) {
    double _sum(0);
    for (int _event = 0; _event < number_of_events_; _event++) {
        for (int i = 0; i < number_of_amplitudes; i++) {
            TComplex cw1 = fCP[i] * TComplex::Conjugate(fCP[i]);
            TComplex cw2(0, 0);
            cw2 += fCF[i][_event + number_of_events_ * 0] *
                   TComplex::Conjugate(fCF[i][_event + number_of_events_ * 0]) /
                   2.0;
            cw2 += fCF[i][_event + number_of_events_ * 1] *
                   TComplex::Conjugate(fCF[i][_event + number_of_events_ * 1]) /
                   2.0;
            _sum += cw1.Re() * cw2.Re();
        }
    }
    return _sum;
}

void CPUWaveFunc::test_generate_root_file_kk(kernel *ker,string res)
{
  //FILE *idp = fopen("/home/arthur/test/old/whu_pwa/newbase/idp_kk_all.dat","r");
  TString kk_weight_file_name = "../test/"+res+"/kk_weight_"+res+".root";
  double idp_list[number_of_events_];
  int idp_num=0;
  double tmp=0;
  /*  while(fscanf(idp,"%d%lf\n",&idp_num,&tmp)!=EOF){
    idp_list[idp_num] = tmp;
    cout<<"test idp "<<tmp<<"test  number "<<idp_num<<endl;
  }
  fclose(idp);*/
    DATA_PWA_KK ss;
    DATA_ANA_KK aa;
    TFile *fout = new TFile(kk_weight_file_name, "RECREATE");
    TTree *ana_tr = new TTree("ana_tr", "ana information");
    
    ana_tr->Branch("kk", &aa, MEM_ANA_KK);
    cout<<"test number_of_events : "<<number_of_events_<<endl;
    for(Int_t i = 0; i < number_of_events_; i++) {
        ss.Kp2X = mcp2[i][0]; ss.Kp2Y = mcp2[i][1]; ss.Kp2Z = mcp2[i][2]; ss.Kp2E = mcp2[i][3];
        ss.Km2X = mcp3[i][0]; ss.Km2Y = mcp3[i][1]; ss.Km2Z = mcp3[i][2]; ss.Km2E = mcp3[i][3];
        ss.Kp1X = mcp4[i][0]; ss.Kp1Y = mcp4[i][1]; ss.Kp1Z = mcp4[i][2]; ss.Kp1E = mcp4[i][3];
        ss.Km1X = mcp5[i][0]; ss.Km1Y = mcp5[i][1]; ss.Km1Z = mcp5[i][2]; ss.Km1E = mcp5[i][3];
        if(i==0) cout<<"test ss : "<<endl<<"Kp2x : "<<ss.Kp2X<<endl;
    //    PWA_PARAS pp;
    //    _amp.calculate0p(
    //            mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],
    //            mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],
    //            mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],
    //            mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],
    //            mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3],
    //            pp);

    //    ss.weight = calEva(pp, i);
        ss.weight = ker->h_phsp_container[i];
        
        //pwa_tr->Fill();
        DATA_ORIG_KK tt;
        pwa_to_orig(ss, tt);
        orig_to_ana(tt, aa);
        ana_tr->Fill();
    }
    
        cout << "haha: " << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
    //pwa_tr->Write();
    ana_tr->Write();
    //fout->Write();
    fout->Close();
    cout << kk_weight_file_name << " is created!!!!" << endl;
  
}
void CPUWaveFunc::test_generate_root_file_pp(kernel *ker,string res)
{
 TString phsp_weight_file_name = "../test/"+res+"/pipi_weight_"+res+".root";
    DATA_PWA_PIPI ss;
    DATA_ANA_PIPI aa;
    TFile *fout = new TFile(phsp_weight_file_name, "RECREATE");
    //TTree *pwa_tr = new TTree("pwa_tr", "pwa information");
    TTree *ana_tr = new TTree("ana_tr", "ana information");
    //pwa_tr->Branch("pp", &ss, MEM_PWA_PIPI);
    ana_tr->Branch("pp", &aa, MEM_ANA_PIPI);
    for(Int_t i = 0; i < number_of_events_; i++) {
        ss.pipX = mcp2[i][0]; ss.pipY = mcp2[i][1]; ss.pipZ = mcp2[i][2]; ss.pipE = mcp2[i][3];
        ss.pimX = mcp3[i][0]; ss.pimY = mcp3[i][1]; ss.pimZ = mcp3[i][2]; ss.pimE = mcp3[i][3];
        ss.KpX = mcp4[i][0]; ss.KpY = mcp4[i][1]; ss.KpZ = mcp4[i][2]; ss.KpE = mcp4[i][3];
        ss.KmX = mcp5[i][0];  ss.KmY = mcp5[i][1]; ss.KmZ = mcp5[i][2]; ss.KmE = mcp5[i][3];
    //    PWA_PARAS pp;
    //    _amp.calculate0p(
    //            mcp1[i][0],mcp1[i][1],mcp1[i][2],mcp1[i][3],
    //            mcp2[i][0],mcp2[i][1],mcp2[i][2],mcp2[i][3],
    //            mcp3[i][0],mcp3[i][1],mcp3[i][2],mcp3[i][3],
    //            mcp4[i][0],mcp4[i][1],mcp4[i][2],mcp4[i][3],
    //            mcp5[i][0],mcp5[i][1],mcp5[i][2],mcp5[i][3],
    //            pp);

    //    ss.weight = calEva(pp, i);
        ss.weight = ker->h_phsp_container[i];
        //pwa_tr->Fill();
        DATA_ORIG_PIPI tt;
        pwa_to_orig(ss, tt);
        orig_to_ana(tt, aa);
        ana_tr->Fill();
    }
        cout << "haha: " << __FILE__ << " " << __FUNCTION__ << "   Line:" << __LINE__ << endl;
    //pwa_tr->Write();
    ana_tr->Write();
    //fout->Write();
    fout->Close();
    cout << phsp_weight_file_name << " is created!!!!" << endl;
}


void test_record(double *val,int number,string fname)
{
  
}
