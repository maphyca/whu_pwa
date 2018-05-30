#ifndef FIT_PARAMETERS_INTERFACE
#define FIT_PARAMETERS_INTERFACE

#include "PWA_CTRL.H"
#include <vector>
#include <string>
#include <map>
#include "whu_constants_and_definitions.h"
#include "whu_my_parameter.h"



class FitParametersInterface {
    public:
        FitParametersInterface() {};
        FitParametersInterface(const PWA_CTRL &);
        std::vector<double> paraList() { return paraList_; }
        std::vector<int> tagParaList() const { return tagParaList_; }
        std::vector<MyParameter> myParameterTable() { return myParameterTable_; }
        int number_of_amplitudes() const { return nAmps; }
//        void format_parameter(MyParameter &, std::string, double);
//        void format_parameter(MyParameter &, std::string, double, double);
//        void format_parameter(MyParameter &, std::string, double, double, double, double);
        void prepare_my_parameter_table(const std::string);
        int position_in_fit_parameter_list(const MyParameter &);
        void initialize_fit_parameter_mapping();
        void add_pp_amplitude_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_pp_amplitude600_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_pp_amplitude980_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_pp_amplitude1680_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_kk_amplitude_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_kk_amplitude600_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_kk_amplitude980_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_kk_amplitude1680_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
    MyParameter gp(std::string parameter_name) {
        for(int i = 0; i < myParameterTable_.size(); i++) {
            if (parameter_name == myParameterTable_[i].get_name()) {
                return myParameterTable_[i];
            }
        }
        std::cout << "Cannot find the parameter " << parameter_name << " in myParameterTable!!" << std::endl;
        exit(1);
    };
    void act_pp_resonance_1p(std::string);
    void act_kk_resonance_1p(std::string);
    void act_pp_resonance_980(std::string);
    void act_kk_resonance_980(std::string);
    void act_pp_resonance_f0(std::string);
    void act_kk_resonance_f0(std::string);
    void act_pp_resonance_f2(std::string);
    void act_kk_resonance_f2(std::string);

    double str2float(const std::string);
    void string_to_vector(std::string, std::vector<std::string>  &);
    void shape_of_mapping_kk();
    void shape_of_mapping_pp();

    private:
        std::vector<MyParameter> myParameterTable_;
        std::vector<MyParameter> fit_parameter_list_;
        std::vector<std::vector<int> > fit_parameter_mapping_kk_;
        std::vector<std::vector<int> > fit_parameter_mapping_pp_;

        std::vector<double> paraList_;
        std::vector<int> tagParaList_;
        int _CN_spinList;
        int _CN_massList;
        int _CN_mass2List;
        int _CN_widthList;
        int _CN_g1List;
        int _CN_g2List;
        int _CN_b1List;
        int _CN_b2List;
        int _CN_b3List;
        int _CN_b4List;
        int _CN_b5List;
        int _CN_rhoList;
        int _CN_fracList;
        int _CN_phiList;
        int _CN_propList;
        int _CN_end;

        int _spinIter;
        int _massIter;
        int _mass2Iter;
        int _widthIter;
        int _g1Iter;
        int _g2Iter;
        int _b1Iter;
        int _b2Iter;
        int _b3Iter;
        int _b4Iter;
        int _b5Iter;
        int _rhoIter;
        int _fracIter;
        int _phiIter;
        int _propIter;

        std::string *nameList;
        std::string *titleList;
        std::string *titleListT;
        int nAmps;
        int nStates;
        int nStatesb;
        int nStatesg1g2;
        int nStateswidth;
        int nStatesmass2;
};


#endif
