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
            /*              //parameter_list_set_[phipp_list_index]->assignment_of_minuit_parameters_from_par(par);
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

            phsp_phikk =
                ((CPUWaveFunc*)data_set_[phikk_phsp_index])->sum_phsp(
                    parameter_list_set_[phikk_list_index]->number_of_amplitudes());
            penalty_phikk =
                ((CPUWaveFunc*)data_set_[phikk_phsp_index])->sum_penalty(
                    parameter_list_set_[phikk_list_index]->number_of_amplitudes());

            likelihood_phipp =
                ((CPUWaveFunc*)data_set_[phipp_data_index])->sum_likelihood(
                    parameter_list_set_[phipp_list_index]->number_of_amplitudes());
            likelihood_phikk =
                ((CPUWaveFunc*)data_set_[phikk_data_index])->sum_likelihood(
                    parameter_list_set_[phikk_list_index]->number_of_amplitudes());
            */
            //gpu part

              (kernel_set_[phipp_phsp_index])->par_trans(
                parameter_list_set_[phipp_list_index]->get_minuit_parameters()
                                                         );
            (kernel_set_[phipp_data_index])->par_trans(
                parameter_list_set_[phipp_list_index]->get_minuit_parameters()
                                                       );
            (kernel_set_[phikk_phsp_index])->par_trans(
                parameter_list_set_[phikk_list_index]->get_minuit_parameters()
                                                       );
            (kernel_set_[phikk_data_index])->par_trans(
                parameter_list_set_[phikk_list_index]->get_minuit_parameters()
                                                       );

            (kernel_set_[phipp_phsp_index])->calEva();
            (kernel_set_[phipp_data_index])->calEva();
            (kernel_set_[phikk_phsp_index])->calEva();
            (kernel_set_[phikk_data_index])->calEva();

            double phsp_phipp, likelihood_phipp, penalty_phipp;
            double phsp_phikk, likelihood_phikk, penalty_phikk;

            phsp_phipp =
              (kernel_set_[phipp_phsp_index])->sum_phsp();
            penalty_phipp =
                (kernel_set_[phipp_phsp_index])->sum_penalty();

            phsp_phikk =
                (kernel_set_[phikk_phsp_index])->sum_phsp();
            penalty_phikk =
                (kernel_set_[phikk_phsp_index])->sum_penalty();

            likelihood_phipp =
                (kernel_set_[phipp_data_index])->sum_likelihood();
            likelihood_phikk =
                (kernel_set_[phikk_data_index])->sum_likelihood();
            

            //message

              cout << " phipp phsp integral = " << phsp_phipp << endl;
              cout << " phipp phsp penalty = " << penalty_phipp << endl;


              cout << " phikk phsp integral = " << phsp_phikk << endl;
              cout << " phikk phsp penalty = " << penalty_phikk << endl;
              
              
              cout << " phipp data likelihood = " << likelihood_phipp << endl;
              cout << " phikk data likelihood = " << likelihood_phikk << endl;

                     for(int i = 0; i < end_list_index; i++) {
              if (parameter_list_set_[i] == NULL) continue;
              parameter_list_set_[i]->copy_minuit_parameter_to_back();
              //parameter_list_set_[i]->shape_of_mapping();
              //parameter_list_set_[i]->shape_of_minuit_parameters();
            }

            //double _rec = likelihood_phipp + likelihood_phikk;
            double _rec = + log(phsp_phipp * number_of_events[phipp_data_index] / number_of_events[phipp_phsp_index]) * number_of_events[phipp_data_index] + likelihood_phipp + log(phsp_phikk * number_of_events[phikk_data_index] / number_of_events[phikk_phsp_index]) * number_of_events[phikk_data_index] + likelihood_phikk;
            //double _rec = + log(phsp_phipp) * 1000 + likelihood_phipp + penalty_phipp + log(phsp_phikk) * 100 + likelihood_phikk + penalty_phikk;
            cout << "rec = " << _rec << endl;
            return _rec;
        }


    }  // namespace Minuit2

}  // namespace ROOT
