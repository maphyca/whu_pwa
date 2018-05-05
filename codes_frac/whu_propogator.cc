#include "whu_propogator.h"
#include "TComplex.h"
#include "whu_constants_and_definitions.h"

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
TComplex CPUWaveFunc::propogator600(
        double mass,
        double b1,
        double b2,
        double b3,
        double b4,
        double b5,
        double sx) const
{
    TComplex ci(0,1);
    double am1=mass;
    double as=am1*am1;
    double cgam1=am1*(b1+b2*sx)*cro(sx,rp,rp)/cro(as,rp,rp)*(sx-0.0097)/(as-0.0097)*exp(-(sx-as)/b3);
    double cgam2=am1*b4*pip(sx)/pip(as);
    TComplex propogator600=1.0/(as-sx-ci*b5*(cgam1+cgam2));
    return propogator600;
}
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
void CPUWaveFunc::cpu_propogator(
        double mass,
        double width,
        vector<TComplex> &crp1,
        const vector<double> &sx,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        crp1[i] = propogator(mass, width, sx[i]);
    }
}
void CPUWaveFunc::cpu_propogator980(
        double mass,
        double g11,
        double g22,
        vector<TComplex> &crp1,
        const vector<double> &sx,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        crp1[i] = propogator980(mass, g11, g22, sx[i]);
    }
}
void CPUWaveFunc::cpu_propogator600(
        double mass,
        double b1,
        double b2,
        double b3,
        double b4,
        double b5,
        vector<TComplex> &crp1,
        const vector<double> &sx,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        crp1[i] = propogator600(mass, b1, b2, b3, b4, b5, sx[i]);
    }
}
void CPUWaveFunc::cpu_propogator1270(
        double mass,
        double width,
        vector<TComplex> &crp1,
        const vector<double> &sx,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        crp1[i] = propogator1270(mass, width, sx[i]);
    }
}
void CPUWaveFunc::cpu_cast_spin11(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<double> &crp11,
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
void CPUWaveFunc::cpu_cast_spin12(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<double> &crp11,
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
void CPUWaveFunc::cpu_cast_spin13(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<double> &crp11,
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
void CPUWaveFunc::cpu_cast_spin14(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<double> &crp11,
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
void CPUWaveFunc::cpu_cast_spin111(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<double> &crp11,
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
void CPUWaveFunc::cpu_cast_spin191(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<double> &crp11,
        const vector<double> &b1q2r23,
        const vector<vector<double> > &ak23w,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double crpf1 = crp1[i]*crp11[i]/b1q2r23[i];
        fCF[i][0] = ak23w[i][0]*crpf1;
        fCF[i][1] = ak23w[i][1]*crpf1;
    }
}
void CPUWaveFunc::cpu_cast_spin192(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<double> &crp11,
        const vector<double> &b1q2r23,
        const vector<double> &b2qjvf2,
        const vector<vector<double> > &wpf22,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        double crpf1 = crp1[i] * crp11[i] / b1q2r23[i];
        double crpf2 = crpf1 / b2qjvf2[i];
        fCF[i][0] = wpf22[i][0] * crpf2;
        fCF[i][1] = wpf22[i][1] * crpf2;
    }
}
void CPUWaveFunc::cpu_cast_spin1(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
        const vector<vector<double> > &wu,
        int vec_size)
{
    for(int i = 0; i < vec_size; i++)
    {
        fCF[i][0] = wu[i][0] * crp1[i];
        fCF[i][1] = wu[i][1] * crp1[i];
    }
}
void CPUWaveFunc::cpu_cast_spin2(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
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
void CPUWaveFunc::cpu_cast_spin21(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
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
void CPUWaveFunc::cpu_cast_spin22(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
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
void CPUWaveFunc::cpu_cast_spin23(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
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
void CPUWaveFunc::cpu_cast_spin24(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
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
void CPUWaveFunc::cpu_cast_spin25(
        vector<vector<TComplex> > &fCF,
        const vector<double> &crp1,
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
