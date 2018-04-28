#include "whu_pwa_fcn.h"
#include "whu_constants_and_definitions.h"


#include <cassert>

namespace ROOT {

    namespace Minuit2 {


        double PWAFcn::operator()(const std::vector<double>& par) const {

            assert(par.size() == number_of_parameters_);
//            for(int i = 0; i < end_list_index; i++)
//            {
//                int _id = 0;
//                for(int j = start_category; j < end_category; j++)
//                {
//                    for(unsigned k = 0; k < parameter_list_set_[i]->get_local_fit_parameter_mapping()[j]; k++)
//                    {
//                        par_list_[i][_id++] = parameter_list_set_[i]->get_local_fit_parameter_mapping()[j][k]);
//                    }
//                }
//            }

//            GaussFunction gauss(par[0], par[1], par[2]);
//
//            double chi2 = 0.;
//            for(unsigned int n = 0; n < fMeasurements.size(); n++) {
//                chi2 += ((gauss(fPositions[n]) - fMeasurements[n])*(gauss(fPositions[n]) - fMeasurements[n])/fMVariances[n]);
//            }
//
//            return chi2;
                return 0;
        }


    }  // namespace Minuit2

}  // namespace ROOT
