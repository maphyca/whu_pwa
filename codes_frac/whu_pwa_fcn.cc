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

            //parameter_list_set_[phipp_list_index]->assignment_of_minuit_parameters_from_par(par);
            ((CPUWaveFunc*)data_set_[phipp_phsp_index])->cpu_calEva(
                parameter_list_set_[phipp_list_index]->get_minuit_parameters(),
                parameter_list_set_[phipp_list_index]->get_minuit_parameters_back(),
                parameter_list_set_[phipp_list_index]->number_of_amplitudes());
            ((CPUWaveFunc*)data_set_[phipp_data_index])->cpu_calEva(
                parameter_list_set_[phipp_list_index]->get_minuit_parameters(),
                parameter_list_set_[phipp_list_index]->get_minuit_parameters_back(),
                parameter_list_set_[phipp_list_index]->number_of_amplitudes());
            ((CPUWaveFunc*)data_set_[phikk_phsp_index])->cpu_calEva(
                parameter_list_set_[phikk_list_index]->get_minuit_parameters(),
                parameter_list_set_[phikk_list_index]->get_minuit_parameters_back(),
                parameter_list_set_[phikk_list_index]->number_of_amplitudes());
            ((CPUWaveFunc*)data_set_[phikk_data_index])->cpu_calEva(
                parameter_list_set_[phikk_list_index]->get_minuit_parameters(),
                parameter_list_set_[phikk_list_index]->get_minuit_parameters_back(),
                parameter_list_set_[phikk_list_index]->number_of_amplitudes());

            double phsp_phipp, likelihood_phipp, penalty_phipp;
            double phsp_phikk, likelihood_phikk, penalty_phikk;

            phsp_phipp =
                ((CPUWaveFunc*)data_set_[phipp_phsp_index])->sum_phsp(
                    parameter_list_set_[phipp_list_index]->number_of_amplitudes());
            penalty_phipp =
                ((CPUWaveFunc*)data_set_[phipp_phsp_index])->sum_penalty(
                    parameter_list_set_[phipp_list_index]->number_of_amplitudes());
            cout << "phipp phsp integral = " << phsp_phipp << endl;
            cout << "phipp phsp penalty = " << penalty_phipp << endl;

            phsp_phikk =
                ((CPUWaveFunc*)data_set_[phikk_phsp_index])->sum_phsp(
                    parameter_list_set_[phikk_list_index]->number_of_amplitudes());
            penalty_phikk =
                ((CPUWaveFunc*)data_set_[phikk_phsp_index])->sum_penalty(
                    parameter_list_set_[phikk_list_index]->number_of_amplitudes());
            cout << "phikk phsp integral = " << phsp_phikk << endl;
            cout << "phikk phsp penalty = " << penalty_phikk << endl;

            likelihood_phipp =
                ((CPUWaveFunc*)data_set_[phipp_data_index])->sum_likelihood(
                    parameter_list_set_[phipp_list_index]->number_of_amplitudes());
            likelihood_phikk =
                ((CPUWaveFunc*)data_set_[phikk_data_index])->sum_likelihood(
                    parameter_list_set_[phikk_list_index]->number_of_amplitudes());
            cout << "phipp data likelihood = " << likelihood_phipp << endl;
            cout << "phikk data likelihood = " << likelihood_phikk << endl;

            for(int i = 0; i < end_list_index; i++) {
                if (parameter_list_set_[i] == NULL) continue;
                parameter_list_set_[i]->copy_minuit_parameter_to_back();
                //parameter_list_set_[i]->shape_of_mapping();
                //parameter_list_set_[i]->shape_of_minuit_parameters();
            }

            double _rec = phsp_phipp + likelihood_phipp + penalty_phipp + phsp_phikk + likelihood_phikk + penalty_phikk;
            cout << "rec = " << _rec << endl;
            return _rec;
        }


    }  // namespace Minuit2

}  // namespace ROOT
