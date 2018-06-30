#ifndef DATA_OBJECT
#define DATA_OBJECT

#include <vector>
#include <iostream>
#include "TString.h"
#include "PWA_PARAS.h"
#include "DPFPWAPoint.h"

class DataObject {
    public:
        DataObject(TString dat_file_name, DPFPWAPoint *pwa_point) :
            dat_file_name_(dat_file_name),
            _dp(pwa_point)
    {
        number_of_events_ = count_lines() / 5;
        cout << "lines = " << number_of_events_ << endl;
        read_events();
    };
        DataObject(TString dat_file_name, TString weight_file_name, DPFPWAPoint *pwa_point) :
            dat_file_name_(dat_file_name),
            weight_file_name_(weight_file_name),
            _dp(pwa_point)
    {
        number_of_events_ = count_lines() / 5;
        cout << "lines = " << number_of_events_ << endl;
        read_events();
        read_weight_file();
    };
        ~DataObject();
        //std::vector<PWA_PARAS> pwa_paras() const { return pwa_paras_; }
        std::vector<double> data_weights() const { return data_weights_; }
        TString dat_file_name() const { return dat_file_name_; }
        int number_of_events() const { return number_of_events_; }
        //void parameters_vector_resize();
        double scalar(vector<double> &, vector<double> &) const;
        double scalar(double *, double *) const;

    protected:
        double **mcp1; // 用来读取格式化的四动量和权重数
        double **mcp2;
        double **mcp3;
        double **mcp4;
        double **mcp5;
        std::vector<double> data_weights_; // 每个事例的权重数据
        TString dat_file_name_;
        TString weight_file_name_;
        int number_of_events_;
        DPFPWAPoint *_dp;

        void load_mcp_from_dat_file();
        //void convert_mcp_to_pwa_paras();
        void initialize_mcp();
        int  count_lines();
        void read_weight_file();
        void read_events();

};

#endif
