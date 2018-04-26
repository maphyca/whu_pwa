#ifndef FIT_PARAMETERS_INTERFACE
#define FIT_PARAMETERS_INTERFACE

#include "PWA_CTRL.H"
#include <vector>
#include <string>
#include <map>
#include "whu_constants_and_definitions.h"
#include "whu_my_parameter.h"

double str2float(const std::string);
void string_to_vector(std::string, std::vector<std::string>  &);


class FitParametersInterface {
    public:
        FitParametersInterface() {};
        static std::vector<MyParameter> myParameterTable() { return my_parameter_table_; }
        static void prepare_my_parameter_table(const std::string);
        static int position_in_fit_parameter_list(const MyParameter &);
        MyParameter gp(std::string);
        static std::vector<MyParameter> get_my_parameter_table() { return my_parameter_table_; };
        static std::vector<MyParameter> get_fit_parameter_list() { return fit_parameter_list_; };
        static void information_of_parameter_table();
        static void information_of_parameter_list_sent_to_minuit();
    protected:
        static std::vector<MyParameter> my_parameter_table_;
        static std::vector<MyParameter> fit_parameter_list_;
};

class AmplitudeMethodWithFitParametersInterface: public FitParametersInterface {
    public:
        AmplitudeMethodWithFitParametersInterface() {
            initialize_fit_parameter_mapping();
        };
        void add_amplitude_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_amplitude600_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_amplitude980_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);
        void add_amplitude1680_to_fit_parameter_list(const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&, const MyParameter&);

        void initialize_fit_parameter_mapping();
        void shape_of_mapping();
        std::vector<std::vector<int> > get_fit_parameter_mapping() { return fit_parameter_mapping_; };
        int number_of_amplitudes();

    protected:
        std::vector<std::vector<int> > fit_parameter_mapping_;

};
class FitParametersOfPhiPP : public AmplitudeMethodWithFitParametersInterface {
    public:
        FitParametersOfPhiPP():AmplitudeMethodWithFitParametersInterface() {};
        void act_resonance_1p(std::string);
        void act_resonance_980(std::string);
        void act_resonance_f0(std::string);
        void act_resonance_f2(std::string);
};
class FitParametersOfPhiKK : public AmplitudeMethodWithFitParametersInterface {
    public:
        FitParametersOfPhiKK():AmplitudeMethodWithFitParametersInterface() {};
        void act_resonance_1p(std::string);
        void act_resonance_980(std::string);
        void act_resonance_f0(std::string);
        void act_resonance_f2(std::string);
};

#endif

