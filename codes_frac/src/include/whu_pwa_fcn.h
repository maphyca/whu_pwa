#ifndef PWA_FCN
#define PWA_FCN

#include "Minuit2/FCNBase.h"
#include "data_obj.h"
#include <iostream>
#include "whu_constants_and_definitions.h"
#include "fit_parameter_interface.h"
#include "whu_propogator.h"

#include <vector>

#endif

namespace ROOT {

    namespace Minuit2 {


        class PWAFcn : public FCNBase {

            public:

                //PWAFcn(std::vector<DataObject*>& data_set, FitParametersInterface **parameter_list_set) :
        PWAFcn(std::vector<DataObject*>& data_set, std::vector<AmplitudeMethodWithFitParametersInterface*>& parameter_list_set,std::vector<kernel*>kernel_set ,int number_of_parameters) :
          data_set_(data_set),
            parameter_list_set_(parameter_list_set),
            kernel_set_(kernel_set),
            number_of_parameters_(number_of_parameters)
            {
                //for(int i = 0; i < end_data_object_index; i++) {
                //    std::cout << data_set_[i]->dat_file_name() << std::endl;
                //    std::cout << data_set_[i]->number_of_events() << std::endl;
                //}
                for(int i = 0; i < end_list_index; i++) {
                    if (parameter_list_set_[i] == NULL) continue;
                    parameter_list_set_[i]->shape_of_mapping();
                }
                ((CPUWaveFunc*)data_set_[phipp_phsp_index])
                    ->cpu_resize_intermediate_variables(
                            parameter_list_set_[phipp_list_index]->number_of_amplitudes());
                ((CPUWaveFunc*)data_set_[phipp_data_index])
                    ->cpu_resize_intermediate_variables(
                            parameter_list_set_[phipp_list_index]->number_of_amplitudes());
                ((CPUWaveFunc*)data_set_[phikk_phsp_index])
                    ->cpu_resize_intermediate_variables(
                            parameter_list_set_[phikk_list_index]->number_of_amplitudes());
                ((CPUWaveFunc*)data_set_[phikk_data_index])
                    ->cpu_resize_intermediate_variables(
                            parameter_list_set_[phikk_list_index]->number_of_amplitudes());
                number_of_events.resize(end_data_object_index);
                for(int i = 0; i < end_data_object_index; i++) {
                    number_of_events[i] = ((CPUWaveFunc*)data_set_[i])->number_of_events();
                }
            };

                ~PWAFcn() {}

                virtual double Up() const {return fErrorDef;}
                virtual double operator()(const std::vector<double>&) const;
                std::vector<DataObject*> data_set() const { return data_set_; }
                std::vector<kernel*> kernel_set() const { return kernel_set_; }
                void SetErrorDef(double def) {fErrorDef = def;}

            private:


                //FitParametersInterface **parameter_list_set_;
                std::vector<DataObject*> data_set_ ;
                std::vector<kernel*> kernel_set_ ;
                std::vector<AmplitudeMethodWithFitParametersInterface*> parameter_list_set_;
                unsigned number_of_parameters_;
                double fErrorDef;
                std::vector<int> number_of_events;

        };

    }  // namespace Minuit2

}  // namespace ROOT
