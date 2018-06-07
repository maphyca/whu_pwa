#ifndef PROPOGATOR_HH
#define PROPOGATOR_HH

#include <math.h>
#include "TComplex.h"
#include "data_obj.h"

using namespace std;

class CPUWaveFunc : public DataObject {
   public:
    // CPUWaveFunc() {};
    CPUWaveFunc(TString dat_file_name, DPFPWAPoint *pwa_point)
        : DataObject(dat_file_name, pwa_point) {
        parameters_vector_resize();
        cpu_convert_mcp_to_pwa_paras();
        //            convert_mcp_to_pwa_paras();
        //        if (!check_integraty()) {
        //            cout << "test integraty fail!!!!" << endl;
        //            exit(1);
        //
        //        }
    };
    CPUWaveFunc(TString dat_file_name, TString weight_file_name,
                DPFPWAPoint *pwa_point)
        : DataObject(dat_file_name, weight_file_name, pwa_point) {
        parameters_vector_resize();
        cpu_convert_mcp_to_pwa_paras();
    };

    virtual ~CPUWaveFunc() {}
    void parameters_vector_resize();
    void cpu_convert_mcp_to_pwa_paras();
    double cpu_calculate0p(int, std::vector<double>&);
    void store_parameters(int, const std::vector<double>&);
    // void convert_mcp_to_pwa_paras();
    bool check_integraty();
    double cpu_calEva(const vector<double> &, vector<double> &, int);
    void cpu_resize_intermediate_variables(int);

    TComplex cro(double, double, double) const;
    TComplex propogator980(double, double, double, double) const;
    TComplex pip(double) const;
    TComplex propogator(double, double, double) const;
    TComplex propogator1270(double, double, double) const;

    std::vector<PWA_PARAS> pwa_paras() const { return pwa_paras_; }

    void cpu_propogator1(double mass, double width, double *sx, double *b2qjvf2,
                         double *wu, double *w0p22, TComplex *fCF0,
                         TComplex *fCF1, int vec_size);
    void cpu_propogator2(double mass, double g11, double g22, double *sx,
                         double *b2qjvf2, double *wu, double *w0p22,
                         TComplex *fCF0, TComplex *fCF1, int vec_size);
    void cpu_propogator7(double mass, double width, double *sv2, double *sv3,
                         double *b1qjv2, double *b1qbv2, double *b1qjv3,
                         double *b1qbv3, double *w1m12, double *w1m13,
                         TComplex *fCF, int vec_size);
    void cpu_propogator8(double mass, double width, double *sv2, double *sv3,
                         double *b2qbv2, double *b2qbv3, double *b2qjv2,
                         double *b2qjv3, double *w1p12_1, double *w1p13_1,
                         double *w1p12_2, double *w1p13_2, double *w1p12_3,
                         double *w1p13_3, double *w1p12_4, double *w1p13_4,
                         TComplex *fCF0, TComplex *fCF1, TComplex *fCF2,
                         TComplex *fCF3, int vec_size);
    void cpu_propogator6(double mass, double width, double *sx, double *b2qf2xx,
                         double *b2qjvf2, double *b4qjvf2, double *w2p1,
                         double *w2p2, double *w2p3, double *w2p4, double *w2p5,
                         TComplex *fCF0, TComplex *fCF1, TComplex *fCF2,
                         TComplex *fCF3, TComplex *fCF4, int vec_size);

    double sum_likelihood(int);
    double sum_phsp(int);
    double sum_penalty(int);

   private:
    std::vector<PWA_PARAS> pwa_paras_;  // 利用四动量转换而来的中间变量
    std::vector<bool> whether_hold_flag;
    std::vector<double *> hpv;  // hold parameter vector for speedup
    //std::vector<double> hp;     // used parameter in computation


    vector<TComplex> fCP;
    vector<TComplex *> fCF;
};

#endif
