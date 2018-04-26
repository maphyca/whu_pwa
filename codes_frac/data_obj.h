#ifndef DATA_OBJECT
#define DATA_OBJECT

#include <vector>
#include <iostream>
#include "TString.h"
#include "PWA_PARAS.h"
#include "DPFPWAPoint.h"

#ifdef TIMING
#include <sys/time.h>
#endif

class DataObject {
    public:
        DataObject(TString dat_file_name, DPFPWAPoint *pwa_point) :
            dat_file_name_(dat_file_name),
            pwa_point_(pwa_point)
    {

#ifdef TIMING
        gettimeofday(&timer_, NULL);
        start_time_ = timer_.tv_sec + timer_.tv_usec / 1000000.0;
#endif

        read_events_and_convert_to_pwa_paras();

#ifdef TIMING
        gettimeofday(&timer_, NULL);
        end_time_ = timer_.tv_sec + timer_.tv_usec / 1000000.0;
        std::cout << "Load mcp from " << dat_file_name_ << " takting time: " << end_time - start_time << std::endl;
#endif

    };
        DataObject(TString dat_file_name, TString weight_file_name, DPFPWAPoint *pwa_point) :
            dat_file_name_(dat_file_name),
            weight_file_name_(weight_file_name),
            pwa_point_(pwa_point)
    {

#ifdef TIMING
        gettimeofday(&timer_, NULL);
        start_time_ = timer_.tv_sec + timer_.tv_usec / 1000000.0;
#endif

        read_events_and_convert_to_pwa_paras();
        read_weight_file();



#ifdef TIMING
        gettimeofday(&timer_, NULL);
        end_time_ = timer_.tv_sec + timer_.tv_usec / 1000000.0;
        std::cout << "Load mcp from " << dat_file_name_ << " takting time: " << end_time - start_time << std::endl;
#endif

    };
        ~DataObject();
        std::vector<PWA_PARAS> pwa_paras() const { return pwa_paras_; }
        std::vector<double> data_weights() const { return data_weights_; }
        TString dat_file_name() const { return dat_file_name_; }
        int number_of_events() const { return number_of_events_; }

    private:
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
        DPFPWAPoint *pwa_point_;

        void load_mcp_from_dat_file();
        void convert_mcp_to_pwa_paras();
        void initialize_mcp();
        int  count_lines();
        void read_weight_file();
        void read_events_and_convert_to_pwa_paras();

#ifdef TIMING
        struct timeval timer_;
        double start_time_, end_time_, step_time_;
#endif
};

#endif
