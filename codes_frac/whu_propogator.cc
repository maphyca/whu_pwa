#include "whu_propogator.h"
#include "TComplex.h"
#include "whu_constants_and_definitions.h"
#include "DPFAngular.h"

//const double rk=0.493677;
//const double rp=0.139556995;

TComplex CPUWaveFunc::cro(
        double sx,
        double am1,
        double am2) const
{
    TComplex ci(0,1);
    double t1=(am1+am2) * (am1 + am2); // double t1=pow((am1+am2),2);
    double t2=(am1-am2) * (am1 - am2); // double t2=pow((am1-am2),2);
    double st=(sx-t1)*(sx-t2);
    double cro=sqrt(fabs(st))/sx;
    TComplex result = cro;
    if (st<0.) result=cro*ci;
    return  result;
}
TComplex CPUWaveFunc::propogator980(
        double mass,
        double g11,
        double g22,
        double sx) const
{
    TComplex ci(0,1);
    double rm=mass*mass;
    TComplex propogator980=1.0/(rm-sx-ci*(g11*cro(sx,rp,rp)+g22*cro(sx,rk,rk)));
    return propogator980;
}
TComplex CPUWaveFunc::pip(
        double sx) const
{
    TComplex ci(0,1);
    double xk2=sx-0.3116676;     //0.3116676=16.*0.139568*0.139568
    if(xk2<=0.)xk2=0.0;
    double r4pip=sqrt(xk2/sx)/(1.0+exp(9.8-3.5*sx));    //9.8=3.5*2.8
    return  r4pip;
}
//TComplex CPUWaveFunc::propogator600(
//        double mass,
//        double b1,
//        double b2,
//        double b3,
//        double b4,
//        double b5,
//        double sx) const
//{
//    TComplex ci(0,1);
//    double am1=mass;
//    double as=am1*am1;
//    double cgam1=am1*(b1+b2*sx)*cro(sx,rp,rp)/cro(as,rp,rp)*(sx-0.0097)/(as-0.0097)*exp(-(sx-as)/b3);
//    double cgam2=am1*b4*pip(sx)/pip(as);
//    TComplex propogator600=1.0/(as-sx-ci*b5*(cgam1+cgam2));
//    return propogator600;
//}
TComplex CPUWaveFunc::propogator(
        double mass,
        double width,
        double sx) const
{
    TComplex ci(0,1);
    double am=mass;
    double g1=mass*width;
    TComplex prop=g1/(sx-am * am +ci*g1); // TComplex prop=g1/(sx-pow(am,2)+ci*g1);
    return prop;
}
TComplex CPUWaveFunc::propogator1270(
        double mass,
        double width,
        double sx) const
{
    TComplex ci(0,1);
    double rm=mass*mass;
    double gr=mass*width;
    double q2r=0.25*rm-0.0194792;
    double b2r=q2r*(q2r+0.1825)+0.033306;
    double g11270=gr*b2r/pow(q2r,2.5);
    double q2=0.25*sx-0.0194792;
    double b2=q2*(q2+0.1825)+0.033306;
    double g1=g11270*pow(q2,2.5)/b2;
    TComplex prop=gr/(sx-rm+ci*g1);
    return prop;
}
void CPUWaveFunc::cpu_propogator1(
        double mass,
        double width,
        const vector<double> &sx,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &wu,
        const vector<vector<double> > &w0p22,
        vector<vector<TComplex> > &fCF0,
        vector<vector<TComplex> > &fCF1,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double crp1 = propogator(mass, width, sx[i]);
        double cr0p11 = crp1 / b2qjvf2[i];

        //01 contribution
        fCF0[i][0] = wu[i][0] * crp1;
        fCF0[i][1] = wu[i][1] * crp1;

        //02 contribution
        fCF1[i][0] = w0p22[i][0] * cr0p11;
        fCF1[i][1] = w0p22[i][1] * cr0p11;
    }

}
void CPUWaveFunc::cpu_propogator2(
        double mass,
        double g11,
        double g22,
        const vector<double> &sx,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &wu,
        const vector<vector<double> > &w0p22,
        vector<vector<TComplex> > &fCF0,
        vector<vector<TComplex> > &fCF1,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double crp1 = propogator980(mass, g11, g22, sx[i]);
        double cr0p11 = crp1 / b2qjvf2[i];

        //01 contribution
        fCF0[i][0] = wu[i][0] * crp1;
        fCF0[i][1] = wu[i][1] * crp1;

        //02 contribution
        fCF1[i][0] = w0p22[i][0] * cr0p11;
        fCF1[i][1] = w0p22[i][1] * cr0p11;
    }

}
void CPUWaveFunc::cpu_propogator7(
        double mass,
        double width,
        const vector<double> &sv2,
        const vector<double> &sv3,
        const vector<double> &b1qjv2,
        const vector<double> &b1qbv2,
        const vector<double> &b1qjv3,
        const vector<double> &b1qbv3,
        const vector<vector<double> > &w1m12,
        const vector<vector<double> > &w1m13,
        vector<vector<TComplex> > &fCF,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double crp1 = propogator(mass, width, sv2[i]);
        double crp11 = propogator(mass, width, sv3[i]);
        double cr1m12_1 = crp1 / b1qjv2[i] / b1qbv2[i];
        double cr1m13_1 = crp11 / b1qjv3[i] / b1qbv3[i];

        //1-__1 contribution
        fCF[i][0] = w1m12[i][0] * cr1m12_1 + w1m13[i][0] * cr1m13_1;
        fCF[i][1] = w1m12[i][1] * cr1m12_1 + w1m13[i][1] * cr1m13_1;
    }
}
void CPUWaveFunc::cpu_propogator8(
        double mass,
        double width,
        const vector<double> &sv2,
        const vector<double> &sv3,
        const vector<double> &b2qbv2,
        const vector<double> &b2qbv3,
        const vector<double> &b2qjv2,
        const vector<double> &b2qjv3,
        const vector<vector<double> > &w1p12_1,
        const vector<vector<double> > &w1p13_1,
        const vector<vector<double> > &w1p12_2,
        const vector<vector<double> > &w1p13_2,
        const vector<vector<double> > &w1p12_3,
        const vector<vector<double> > &w1p13_3,
        const vector<vector<double> > &w1p12_4,
        const vector<vector<double> > &w1p13_4,
        vector<vector<TComplex> > &fCF0,
        vector<vector<TComplex> > &fCF1,
        vector<vector<TComplex> > &fCF2,
        vector<vector<TComplex> > &fCF3,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double crp1 = propogator(mass, width, sv2[i]);
        double crp11 = propogator(mass, width, sv3[i]);
        double c1p12_12 = crp1 / b2qbv2[i];
        double c1p13_12 = crp11 / b2qbv3[i];
        double c1p12_13 = crp1 / b2qjv2[i];
        double c1p13_13 = crp11 / b2qjv3[i];
        double c1p12_14 = c1p12_12 / b2qjv2[i];
        double c1p13_14 = c1p13_12 / b2qjv3[i];

        // z 1+ 1
        fCF0[i][0] = w1p12_1[i][0] * crp1 + w1p13_1[i][0] * crp11;
        fCF0[i][1] = w1p12_1[i][1] * crp1 + w1p13_1[i][1] * crp11;

        // z 1+ 2
        fCF1[i][0] = w1p12_2[i][0] * c1p12_12 + w1p13_2[i][0] * c1p13_12;
        fCF1[i][1] = w1p12_2[i][1] * c1p12_12 + w1p13_2[i][1] * c1p13_12;

        // z 1+ 3
        fCF2[i][0] = w1p12_3[i][0] * c1p12_13 + w1p13_3[i][0] * c1p13_13;
        fCF2[i][1] = w1p12_3[i][1] * c1p12_13 + w1p13_3[i][1] * c1p13_13;

        // z 1+ 4
        fCF3[i][0] = w1p12_4[i][0] * c1p12_14 + w1p13_4[i][0] * c1p13_14;
        fCF3[i][1] = w1p12_4[i][1] * c1p12_14 + w1p13_4[i][1] * c1p13_14;
    }
}
void CPUWaveFunc::cpu_propogator6(
        double mass,
        double width,
        const vector<double> &sx,
        const vector<double> &b2qf2xx,
        const vector<double> &b2qjvf2,
        const vector<double> &b4qjvf2,
        const vector<vector<double> > &w2p1,
        const vector<vector<double> > &w2p2,
        const vector<vector<double> > &w2p3,
        const vector<vector<double> > &w2p4,
        const vector<vector<double> > &w2p5,
        vector<vector<TComplex> > &fCF0,
        vector<vector<TComplex> > &fCF1,
        vector<vector<TComplex> > &fCF2,
        vector<vector<TComplex> > &fCF3,
        vector<vector<TComplex> > &fCF4,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double crp1 = propogator1270(mass, width, sx[i]);
        double cw2p11 = crp1 / b2qf2xx[i];
        double cw2p12 = cw2p11 / b2qjvf2[i];
        double cw2p15 = cw2p11 / b4qjvf2[i];

        //21 contribution
        fCF0[i][0] = w2p1[i][0] * cw2p11;
        fCF0[i][1] = w2p1[i][1] * cw2p11;

        //22 contribution
        fCF1[i][0]=w2p2[i][0] * cw2p12;
        fCF1[i][1]=w2p2[i][1] * cw2p12;

        //23 contribution
        fCF2[i][0] = w2p3[i][0] * cw2p12;
        fCF2[i][1] = w2p3[i][1] * cw2p12;

        //24 contribution
        fCF3[i][0] = w2p4[i][0] * cw2p12;
        fCF3[i][1] = w2p4[i][1] * cw2p12;

        //25 contribution
        fCF4[i][0] = w2p5[i][0] * cw2p15;
        fCF4[i][1] = w2p5[i][1] * cw2p15;
    }
}

void CPUWaveFunc::cpu_cast_spin801(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<TComplex> &crp11,
        const vector<vector<double> > &w1p12_1,
        const vector<vector<double> > &w1p13_1,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        fCF[i][0] = w1p12_1[i][0] * crp1[i] + w1p13_1[i][0] * crp11[i];
        fCF[i][1] = w1p12_1[i][1] * crp1[i] + w1p13_1[i][1] * crp11[i];
    }
}
void CPUWaveFunc::cpu_cast_spin802(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<TComplex> &crp11,
        const vector<double> &b2qbv2,
        const vector<double> &b2qbv3,
        const vector<vector<double> > &w1p12_2,
        const vector<vector<double> > &w1p13_2,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double c1p12_12 = crp1[i] / b2qbv2[i];
        double c1p13_12 = crp11[i] / b2qbv3[i];
        fCF[i][0] = w1p12_2[i][0] * c1p12_12 + w1p13_2[i][0] * c1p13_12;
        fCF[i][1] = w1p12_2[i][1] * c1p12_12 + w1p13_2[i][1] * c1p13_12;
    }
}
void CPUWaveFunc::cpu_cast_spin803(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<TComplex> &crp11,
        const vector<double> &b2qjv2,
        const vector<double> &b2qjv3,
        const vector<vector<double> > &w1p12_3,
        const vector<vector<double> > &w1p13_3,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double c1p12_13 = crp1[i] / b2qjv2[i];
        double c1p13_13 = crp11[i] / b2qjv3[i];
        fCF[i][0] = w1p12_3[i][0] * c1p12_13 + w1p13_3[i][0] * c1p13_13;
        fCF[i][1] = w1p12_3[i][1] * c1p12_13 + w1p13_3[i][1] * c1p13_13;
    }
}
void CPUWaveFunc::cpu_cast_spin804(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<TComplex> &crp11,
        const vector<double> &b2qbv2,
        const vector<double> &b2qjv2,
        const vector<double> &b2qbv3,
        const vector<double> &b2qjv3,
        const vector<vector<double> > &w1p12_4,
        const vector<vector<double> > &w1p13_4,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double c1p12_12 = crp1[i] / b2qbv2[i];
        double c1p13_12 = crp11[i] / b2qbv3[i];
        double c1p12_14 = c1p12_12 / b2qjv2[i];
        double c1p13_14 = c1p13_12 / b2qjv3[i];
        fCF[i][0] = w1p12_4[i][0] * c1p12_14 + w1p13_4[i][0] * c1p13_14;
        fCF[i][1] = w1p12_4[i][1] * c1p12_14 + w1p13_4[i][1] * c1p13_14;
    }
}
void CPUWaveFunc::cpu_cast_spin701(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<TComplex> &crp11,
        const vector<double> &b1qjv2,
        const vector<double> &b1qbv2,
        const vector<double> &b1qjv3,
        const vector<double> &b1qbv3,
        const vector<vector<double> > &w1m12,
        const vector<vector<double> > &w1m13,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cr1m12_1 = crp1[i] / b1qjv2[i] / b1qbv2[i];
        double cr1m13_1 = crp11[i] / b1qjv3[i] / b1qbv3[i];
        fCF[i][0] = w1m12[i][0] * cr1m12_1 + w1m13[i][0] * cr1m13_1;
        fCF[i][1] = w1m12[i][1] * cr1m12_1 + w1m13[i][1] * cr1m13_1;
    }
}
void CPUWaveFunc::cpu_cast_spin101(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<vector<double> > &wu,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        fCF[i][0] = wu[i][0] * crp1[i];
        fCF[i][1] = wu[i][1] * crp1[i];
    }
}
void CPUWaveFunc::cpu_cast_spin102(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &w0p22,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cr0p11 = crp1[i] / b2qjvf2[i];
        fCF[i][0] = w0p22[i][0] * cr0p11;
        fCF[i][1] = w0p22[i][1] * cr0p11;
    }
}
void CPUWaveFunc::cpu_cast_spin201(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<vector<double> > &wu,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        fCF[i][0] = wu[i][0] * crp1[i];
        fCF[i][1] = wu[i][1] * crp1[i];
    }
}
void CPUWaveFunc::cpu_cast_spin202(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &w0p22,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cr0p11 = crp1[i] / b2qjvf2[i];
        fCF[i][0] = w0p22[i][0] * cr0p11;
        fCF[i][1] = w0p22[i][1] * cr0p11;
    }
}
void CPUWaveFunc::cpu_cast_spin601(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<double> &b2qf2xx,
        const vector<vector<double> > &w2p1,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cw2p11 = crp1[i] / b2qf2xx[i];
        fCF[i][0] = w2p1[i][0] * cw2p11;
        fCF[i][1] = w2p1[i][1] * cw2p11;
    }
}
void CPUWaveFunc::cpu_cast_spin602(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<double> &b2qf2xx,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &w2p2,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cw2p11 = crp1[i] / b2qf2xx[i];
        double cw2p12 = cw2p11 / b2qjvf2[i];
        fCF[i][0]=w2p2[i][0] * cw2p12;
        fCF[i][1]=w2p2[i][1] * cw2p12;
    }
}
void CPUWaveFunc::cpu_cast_spin603(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<double> &b2qf2xx,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &w2p3,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cw2p11 = crp1[i] / b2qf2xx[i];
        double cw2p12 = cw2p11 / b2qjvf2[i];
        fCF[i][0] = w2p3[i][0] * cw2p12;
        fCF[i][1] = w2p3[i][1] * cw2p12;
    }

}
void CPUWaveFunc::cpu_cast_spin604(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<double> &b2qf2xx,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &w2p4,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cw2p11 = crp1[i] / b2qf2xx[i];
        double cw2p12 = cw2p11 / b2qjvf2[i];
        fCF[i][0] = w2p4[i][0] * cw2p12;
        fCF[i][1] = w2p4[i][1] * cw2p12;
    }
}
void CPUWaveFunc::cpu_cast_spin605(
        vector<vector<TComplex> > &fCF,
        const vector<TComplex> &crp1,
        const vector<double> &b2qf2xx,
        const vector<double> &b4qjvf2,
        const vector<vector<double> > &w2p5,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double cw2p11 = crp1[i] / b2qf2xx[i];
        double cw2p15 = cw2p11 / b4qjvf2[i];
        fCF[i][0] = w2p5[i][0] * cw2p15;
        fCF[i][1] = w2p5[i][1] * cw2p15;
    }
}
void CPUWaveFunc::parameters_vector_resize()
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
void CPUWaveFunc::cpu_convert_mcp_to_pwa_paras() {
    pwa_paras_.resize(0);
    for(int i = 0; i < number_of_events_; i++) {
        cpu_calculate0p(i);
    }
}
void CPUWaveFunc::convert_mcp_to_pwa_paras() {
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
bool CPUWaveFunc::check_integraty()
{
    for(int i = 0; i < number_of_events_; i++) {
        if (pwa_paras_[i].b2qjvf2 != b2qjvf2[i]) return false;
        if (pwa_paras_[i].w0p22[0] != w0p22[i][0]) return false;
        if (pwa_paras_[i].w0p22[1] != w0p22[i][1]) return false;
        if (pwa_paras_[i].ak23w[1] != ak23w[i][1]) return false;
        if (pwa_paras_[i].w1m12[1] != w1m12[i][1]) return false;
        if (pwa_paras_[i].w1m13[1] != w1m13[i][1]) return false;
        if (pwa_paras_[i].w1m13[0] != w1m13[i][0]) return false;
    }
    return true;

}
void CPUWaveFunc::cpu_resize_intermediate_variables(int number_of_amplitudes)
{
    fCP.resize(number_of_amplitudes);
    crp1.resize(number_of_amplitudes);
    crp11.resize(number_of_amplitudes);
    fCF.resize(number_of_amplitudes);
    for(unsigned i = 0; i < crp1.size(); i++)
    {
        crp1[i].resize(number_of_events_);
        crp11[i].resize(number_of_events_);
        fCF[i].resize(number_of_events_, vector<TComplex>(4));
    }
}
double CPUWaveFunc::cpu_calEva(const vector<double> &par, vector<double> &par_back, int number_of_amplitudes)
{
    vector<bool> not_changed(number_of_amplitudes);
    int i = 0;
    while (i < number_of_amplitudes)
    {
        int propType_now = par[end_category * i + propType_category];
        switch(propType_now)
        {
            case 1: // f0
                {
                    double mass0 = par[end_category * i + mass_category];
                    double width0 = par[end_category * i + width_category];
                    not_changed[i] =
                        ((mass0 == par_back[end_category * i + mass_category])
                         && (width0 == par_back[end_category * i + width_category]));
                    if (!not_changed[i]) {
                        cpu_propogator1(
                                mass0,
                                width0,
                                s23,
                                b2qjvf2,
                                wu,
                                w0p22,
                                fCF[i],
                                fCF[i + 1],
                                number_of_events_);
                        i = i + 2;
                    }
                }
                break;
                //	Flatte   Propagator Contribution
            case 2: // f0 980
                {
                    double mass980 = par[end_category * i + mass_category];
                    double g10 = par[end_category * i + g1_category];
                    double g20 = par[end_category * i + g2_category];
                    not_changed[i] =
                        ((mass980 == par_back[end_category * i + mass_category])
                         && (g10 == par_back[end_category * i + g1_category])
                         && (g20 == par_back[end_category * i + g2_category]));
                    if (!not_changed[i]) {
                        cpu_propogator2(
                                mass980,
                                g10,
                                g20,
                                s23,
                                b2qjvf2,
                                w0p22,
                                wu,
                                fCF[i],
                                fCF[i + 1],
                                number_of_events_);
                        i = i + 2;
                    }
                }
                break;
            case 7: //1m1800
                {
                    double mass0 = par[end_category * i + mass_category];
                    double width0 = par[end_category * i + width_category];
                    not_changed[i] =
                        ((mass0 == par_back[end_category * i + mass_category])
                         && (width0 == par_back[end_category * i + width_category]));
                    if (!not_changed[i])
                    {
                        cpu_propogator7(
                                mass0,
                                width0,
                                sv2,
                                sv3,
                                b1qjv2,
                                b1qbv2,
                                b1qjv3,
                                b1qbv3,
                                w1m12,
                                w1m13,
                                fCF[i],
                                number_of_events_);
                        i = i + 1;
                    }
                }
                break;
            case 8: //1p1800
                {
                    double mass0 = par[end_category * i + mass_category];
                    double width0 = par[end_category * i + width_category];
                    not_changed[i] =
                        ((mass0 == par_back[end_category * i + mass_category])
                         && (width0 == par_back[end_category * i + width_category]));
                    if (!not_changed[i])
                    {
                        cpu_propogator8(
                                mass0,
                                width0,
                                sv2,
                                sv3,
                                b2qbv2,
                                b2qbv3,
                                b2qjv2,
                                b2qjv3,
                                w1p12_1,
                                w1p13_1,
                                w1p12_2,
                                w1p13_2,
                                w1p12_3,
                                w1p13_3,
                                w1p12_4,
                                w1p13_4,
                                fCF[i],
                                fCF[i + 1],
                                fCF[i + 2],
                                fCF[i + 3],
                                number_of_events_);
                        i = i + 4;
                    }
                }
                break;
            case 6: //f2
                {
                    double mass0 = par[end_category * i + mass_category];
                    double width0 = par[end_category * i + width_category];
                    not_changed[i] =
                        ((mass0 == par_back[end_category * i + mass_category])
                         && (width0 == par_back[end_category * i + width_category]));
                    if (!not_changed[i]) {
                        cpu_propogator6(
                                mass0,
                                width0,
                                s23,
                                b2qf2xx,
                                b2qjvf2,
                                b4qjvf2,
                                w2p1,
                                w2p2,
                                w2p3,
                                w2p4,
                                w2p5,
                                fCF[i],
                                fCF[i + 1],
                                fCF[i + 2],
                                fCF[i + 3],
                                fCF[i + 4],
                                number_of_events_);
                        i = i + 5;
                    }
                }
            default :
                cout << "Do not know how to deal with prop type " << propType_now << endl;
                exit(1);
                ;
        }
    }

    par_back = par;

    for(int i = 0; i < number_of_amplitudes; i++) {
        double rho0 = par[end_category * i + rho_category];
        double frac0 = par[end_category * i + frac_category];
        double phi0 = par[end_category * i + phi_category];
        int propType_now = par[end_category * i + propType_category];
        rho0 *= TMath::Exp(frac0);
        fCP[i]=TComplex(rho0*TMath::Cos(phi0),rho0*TMath::Sin(phi0));
    }
    vector<double> cw1(number_of_events_, 0), cw2(number_of_events_, 0);
    for(int _event = 0; _event < number_of_events_; _event++)
    {
        for(int i = 0; i < number_of_amplitudes; i++)
        {
            cw1[_event] = cw1[_event] + fCP[i] * fCF[i][_event][0];
            cw2[_event] = cw2[_event] + fCP[i] * fCF[i][_event][1];
        }
        //value[_event] = (cw1[_event] * TComplex::Conjugate(cw1[_event]) + cw2[_event] * TComplex::Conjugate(cw2[_event]) ) / 2.0;
    }
    for(int _event = 0; _event < number_of_events_; _event++)
    {
        for(int i = 0; i < number_of_amplitudes; i++)
        {
            TComplex cw = fCP[i] * TComplex::Conjugate(fCP[i]);
            double pa = cw.Re();
            cw = TComplex(0, 0);
            cw += fCF[i][_event][0] * TComplex::Conjugate(fCF[i][_event][0]) / 2.0;
            cw += fCF[i][_event][1] * TComplex::Conjugate(fCF[i][_event][1]) / 2.0;
            double fu = cw.Re();
        //mlk[_event][i] = pa * fu;
        }
    }


    for(int _event = 0; _event < number_of_events_; _event++)
    {
        double carry(0);
        for(int i = 0; i < number_of_amplitudes; i++)
        {
            for(int j = 0; j < number_of_amplitudes; j++)
            {
                double pa, fu;
                TComplex cw = fCP[i] * TComplex::Conjugate(fCP[j]);
                if (i==j) pa = cw.Re();
                else if(i<j) pa = 2*cw.Re();
                else pa = 2 * cw.Im();

                cw = TComplex(0,0);
                for(int k=0; k < 2; k++)
                {
                    cw += fCF[i][_event][k] * TComplex::Conjugate(fCF[j][_event][k]) / 2.0;
                }
                if (i <= j) fu = cw.Re();
                if (i > j) fu = -cw.Im();

                carry += pa * fu;
            }
        }
        //value[_event] = (carry <= 0) ? 1e-30 : carry;
    }
    for(int _event = 0; _event < number_of_events_; _event++)
    {
        for(int i = 0; i < number_of_amplitudes; i++)
        {
                TComplex cw = fCP[i] * TComplex::Conjugate(fCP[i]);
                double pa = cw.Re();

                cw = TComplex(0,0);
                for(int k=0; k < 2; k++)
                {
                    cw += fCF[i][_event][k] * TComplex::Conjugate(fCF[i][_event][k]) / 2.0;
                }
                double fu = cw.Re();
                //mlk[_event][i] = pa * fu;
        }
    }
    return 0;

}
