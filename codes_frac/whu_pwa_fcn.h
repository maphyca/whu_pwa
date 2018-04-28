#ifndef PWA_FCN
#define PWA_FCN

#include "Minuit2/FCNBase.h"
#include "data_obj.h"
#include <iostream>
#include "whu_constants_and_definitions.h"
#include "fit_parameter_interface.h"

#include <vector>

#endif

namespace ROOT {

    namespace Minuit2 {


        class PWAFcn : public FCNBase {

            public:

                //PWAFcn(std::vector<DataObject*>& data_set, FitParametersInterface **parameter_list_set) :
                PWAFcn(std::vector<DataObject*>& data_set, std::vector<AmplitudeMethodWithFitParametersInterface*>& parameter_list_set, int number_of_parameters) :
                    data_set_(data_set),
                    parameter_list_set_(parameter_list_set),
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

            };

                ~PWAFcn() {}

                virtual double Up() const {return fErrorDef;}
                virtual double operator()(const std::vector<double>&) const;
                std::vector<DataObject*> data_set() const { return data_set_; }
                void SetErrorDef(double def) {fErrorDef = def;}

            private:


                //FitParametersInterface **parameter_list_set_;
                std::vector<DataObject*> data_set_ ;
                std::vector<AmplitudeMethodWithFitParametersInterface*> parameter_list_set_;
                double fErrorDef;
                unsigned number_of_parameters_;
                std::vector<std::vector<int> > tag_list_;
                std::vector<std::vector<int> > par_list_;
        };

    }  // namespace Minuit2

}  // namespace ROOT
