#include "whu_pwa_fcn.h"
#include "whu_constants_and_definitions.h"


#include <cassert>

namespace ROOT {

    namespace Minuit2 {


        double PWAFcn::operator()(const std::vector<double>& par) const {
            assert(par.size() == number_of_parameters_);
            for(int i = 0; i < end_list_index; i++) {
                if (parameter_list_set_[i] == NULL) continue;
                parameter_list_set_[i]->assignment_of_minuit_parameters_from_par(par);
                //parameter_list_set_[i]->shape_of_mapping();
                //parameter_list_set_[i]->shape_of_minuit_parameters();
            }
            return 0;
        }


    }  // namespace Minuit2

}  // namespace ROOT
