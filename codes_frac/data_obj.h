#ifndef DATA_OBJECT 
#define DATA_OBJECT

#include <vector>
#include <iostream>
#include "TString.h"
#include "PWA_PARAS.h"

#ifdef TIMING
#include <sys/time.h>
#endif

class DataObject {
  public:
    DataObject(TString dat_file_name) : dat_file_name_(dat_file_name)
  {
#ifdef TIMING
    gettimeofday(&timer_, NULL);
    start_time_ = timer_.tv_sec + timer_.tv_usec / 1000000.0;
#endif  
    std::cout << "TEST" << std::endl;
    number_of_events_ = count_lines() / 5;
    initialize_mcp();
    load_mcp_from_dat_file();
    std::cout << "lines: " << number_of_events_ << std::endl;
    convert_mcp_to_pwa_paras();
#ifdef TIMING
    gettimeofday(&timer_, NULL);
    end_time_ = timer_.tv_sec + timer_.tv_usec / 1000000.0;
    std::cout << "Load mcp from " << dat_file_name_ << " takting time: " << end_time - start_time << std::endl;
#endif  
  };
    ~DataObject();
    std::vector<PWA_PARAS> pwa_paras() const { return pwa_paras_; }
    TString dat_file_name() const { return dat_file_name_; }
    int number_of_events() const { return number_of_events_; }

  private:
    double **mcp1; // 用来读取格式化的四动量和权重数据 
    double **mcp2;  
    double **mcp3;  
    double **mcp4;  
    double **mcp5;  
    std::vector<PWA_PARAS> pwa_paras_; // 利用四动量转换而来的中间变量
    TString dat_file_name_;
    int number_of_events_;

    void load_mcp_from_dat_file();
    void convert_mcp_to_pwa_paras();
    void initialize_mcp();
    int count_lines();
#ifdef TIMING
    struct timeval timer_;
    double start_time_, end_time_, step_time_;
#endif
};

#endif
