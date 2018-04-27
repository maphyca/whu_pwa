#include "whu_pwa_fcn.h"
#include "whu_constants_and_definitions.h"


#include <cassert>

namespace ROOT {

    namespace Minuit2 {


        double PWAFcn::operator()(const std::vector<double>& par) const {

//            assert(par.size() == 3);
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
