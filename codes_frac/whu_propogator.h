#ifndef PROPOGATOR_HH
#define PROPOGATOR_HH

#include <iostream>
#include <fstream>
#include <math.h>
#include "TNamed.h"

//#if defined(USEROOT) || defined(__CINT__)
//#include "RooComplex.h"
//#else
//#include "RooFitCore/RooComplex.hh"
//#endif

#include "DPFPWAPoint.h"
#include "TComplex.h"

using namespace std;

class Propogator {

public:
  Propogator() {};

  virtual ~Propogator() {}


  TComplex cro(Double_t sx, Double_t am1, Double_t am2)const;
  TComplex propogator980(Double_t mass, Double_t g11, Double_t g22,Double_t sx)const;
  TComplex pip(Double_t sx)const;
  TComplex propogator600(Double_t mass, Double_t b1, Double_t b2, Double_t b3, Double_t b4, Double_t b5, Double_t sx)const;
  TComplex propogator(Double_t mass,Double_t width,Double_t sx) const;
  TComplex propogator1270(Double_t mass,Double_t width,Double_t sx) const;

  void cpu_propogator(double, double, vector<TComplex> &, const vector<double>&, int);
  void cpu_propogator980(double, double, double, vector<TComplex> &, const vector<double>&, int);
  void cpu_propogator600(double, double, double, double, double, double, vector<TComplex> &, const vector<double>&, int);
  void cpu_propogator1270(double, double, vector<TComplex> &, const vector<double>&, int);

void cpu_cast_spin11(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin12(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin13(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin14(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin111(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin191(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin192(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin1(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin2(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin21(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin22(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin23(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin24(
        vector<vector<TComplex> > &fCF,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);
void cpu_cast_spin25(
        vector<vector<TComplex> > &,
        const vector<double> &,
        const vector<double> &,
        const vector<double> &,
        const vector<vector<double> > &,
        int);


private:

//  DPFPWAPoint* _dp;

//  ClassDef(Propogator,1)
};

#endif
