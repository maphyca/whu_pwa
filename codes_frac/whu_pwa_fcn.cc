#include "whu_pwa_fcn.h"
#include "whu_constants_and_definitions.h"


#include <cassert>

namespace ROOT {

    namespace Minuit2 {


        double PWAFcn::operator()(const std::vector<double>& par) const {
            //cout << par.size() << " === " << number_of_parameters_ << endl;
            //cout << par[36] << endl;
            assert(par.size() == number_of_parameters_);
            for(int i = 0; i < end_list_index; i++) {
                if (parameter_list_set_[i] == NULL) continue;
                parameter_list_set_[i]->assignment_of_minuit_parameters_from_par(par);
                //parameter_list_set_[i]->shape_of_mapping();
                //parameter_list_set_[i]->shape_of_minuit_parameters();
            }
            //cout << ">>>>>>>>>>>>" << endl;
            //for(unsigned i = 0; i < par.size(); i++) {
            //    cout << par[i] << endl;
            //}

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
