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
        parameters_vector_resize();
        read_events();
        cout << "begin cpu convert mcp to pwa paras!!!" << endl;
        cpu_convert_mcp_to_pwa_paras();
        //if (!check_integraty()) {
        //    cout << "test integraty fail!!!!" << endl;
        //    exit(1);

        //}
    };
        DataObject(TString dat_file_name, TString weight_file_name, DPFPWAPoint *pwa_point) :
            dat_file_name_(dat_file_name),
            weight_file_name_(weight_file_name),
            _dp(pwa_point)
    {
        number_of_events_ = count_lines() / 5;
        parameters_vector_resize();
        read_events();
        convert_mcp_to_pwa_paras();
        read_weight_file();
    };
        ~DataObject();
        std::vector<PWA_PARAS> pwa_paras() const { return pwa_paras_; }
        std::vector<double> data_weights() const { return data_weights_; }
        TString dat_file_name() const { return dat_file_name_; }
        int number_of_events() const { return number_of_events_; }
        void parameters_vector_resize();
        double scalar(vector<double> &, vector<double> &) const;
        double scalar(double *, double *) const;
        double cpu_calculate0p(
                //double, double, double, double,
                //double, double, double, double,
                //double, double, double, double,
                //double, double, double, double,
                //double, double, double, double,
                int);
        void cpu_convert_mcp_to_pwa_paras();
        bool check_integraty();

    protected:
        double **mcp1; // 用来读取格式化的四动量和权重数
        double **mcp2;
        double **mcp3;
        double **mcp4;
        double **mcp5;
        std::vector<PWA_PARAS> pwa_paras_; // 利用四动量转换而来的中间变量
        std::vector<double> data_weights_; // 每个事例的权重数据
        TString dat_file_name_;
        TString weight_file_name_;
        int number_of_events_;
        DPFPWAPoint *_dp;

        void load_mcp_from_dat_file();
        void convert_mcp_to_pwa_paras();
        void initialize_mcp();
        int  count_lines();
        void read_weight_file();
        void read_events();

    vector<vector<double> > wu,w0p22,ak23w,w2p2,w2p1;
    vector<vector<double> > w2p3,w2p4,w2p5;
    vector<double> b2qf2xx,b4qjvf2,b2qjv2,b2qjv3,b2qbv2,b2qbv3,b1qjv2,b1qjv3,b1qbv2,b1qbv3,b1q2r23;
    vector<double> sv,s23,sv2,sv3;
    vector<vector<double> > wpf22;
    vector<double> b2qjvf2;
    vector<vector<double> > w1p12_1,w1p13_1,w1p12_2,w1p13_2;
    vector<vector<double> > w1p12_3,w1p13_3,w1p12_4,w1p13_4;
    vector<vector<double> > w1m12,w1m13;
};

#endif
