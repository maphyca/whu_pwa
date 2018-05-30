#ifndef FITPROXY_H
#define FITPROXY_H

#include "DPFPWAPdf.h"
#include "DPFPWAPoint.h"
#include "TString.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "PWA_CTRL.H"


//void Info(RooRealVar *bb);

class fitproxy {
private:
    TString outf_phipp_file_name_, outf_phikk_file_name_;
    TString data_phipp_file_name_, data_phikk_file_name_;
    TString proj_phipp_file_name_, proj_phikk_file_name_;
    TString phsp_phipp_file_name_, phsp_phikk_file_name_;
    TString idx_pp_file_name_, idx_kk_file_name_;
    RooArgSet theSet;

    vector<RooRealVar> vv;
    vector< vector<RooRealVar> > vpars;
    vector<TString> vn;

    RooRealVar idp;

    RooArgSet allparas, fitparas;
    RooDataSet *datapp, *datakk;

    //DPFPWAPoint *dphipp, *dphikk;
    DPFPWAPdf *pdfphipp, *pdfphikk;

    DataObject *phikk_phsp_events_, *phipp_phsp_events_;
    DataObject *phikk_data_events_, *phikk_data_events_;

    FitParametersInterface *fit_parameters_interface_;

public:
    fitproxy(const PWA_CTRL & pwa_ctrl) {
        outf_phipp_file_name_ = pwa_ctrl.outf_phipp;
        outf_phikk_file_name_ = pwa_ctrl.outf_phikk;
        proj_phipp_file_name_ = pwa_ctrl.proj_phipp;
        proj_phikk_file_name_ = pwa_ctrl.proj_phikk;

        phsp_phipp_file_name_ = pwa_ctrl.phsp_phipp;
        phsp_phikk_file_name_ = pwa_ctrl.phsp_phikk;
        data_phipp_file_name_ = pwa_ctrl.data_phipp;
        data_phikk_file_name_ = pwa_ctrl.data_phikk;

        idx_pp_file_name_ = pwa_ctrl.idx_pp;
        idx_kk_file_name_ = pwa_ctrl.idx_kk;

        phikk_phsp_events_ = new DataObject(phsp_phikk_file_name_, pwa_points_phikk);
        phikk_data_events_ = new DataObject(data_phikk_file_name_, idx_kk_file_name_, pwa_points_phikk);
        phipp_phsp_events_ = new DataObject(data_phipp_file_name_, pwa_points_phipp);
        phipp_data_events_ = new DataObject(data_phipp_file_name_, idx_pp_file_name_, pwa_points_phipp);

    };
    ~fitproxy() {
        cout << "*****" << endl;
        cout << "delete dphipp" << endl;
        delete dphipp;
        delete dphikk;
        cout << "delete pdfphipp" << endl;
        delete pdfphipp;
        delete pdfphikk;
        cout << "delete datapp" << endl;
        datapp->Delete();
        delete datapp;
        datakk->Delete();
        delete datakk;
        cout << "End delete" << endl;
        vpars.resize(0);
        vv.resize(0);
        vn.resize(0);
        cout << "end vpars resize" << endl;
        theSet.Delete();
        allparas.Delete();
        fitparas.Delete();
        vv.resize(0);
        vpars.resize(0);
        vn.resize(0);
    };
    //void init_all_files_name(const PWA_CTRL&);
    //void init_input_argset();
    void read_data();
    void init_pdf(const PWA_CTRL&);
    void setup_resonances();
    void print_all_paras();
    void print_fit_paras();
    void add_res0_list(TString, double, double, double, double);
    void add_res2_list(TString, double, double, double, double);
    void add_res980_list(TString, double, double);
    void add_res1m_list(TString, double, double, double, double);
    void add_res1p_list(TString, double, double, double, double);
    void createlist_allparas();
    void add_fitparas(TString rn);

    void act_res0(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res2(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res980(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res1m(TString rn, bool act_pp = true, bool act_kk = true);
    void act_res1p(TString rn, bool act_pp = true, bool act_kk = true);
    RooRealVar* gp(TString tn);

    void set_all_constant(TString rn);
    void set_all_constant();
    void store_fit_paras(TString fname);
    void store_all_paras(TString fname);
    void reload_paras(TString fname);

    void FIT();
    void Prepare_Figs(const PWA_CTRL &);

};

#endif


