
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TString.h"
#include "data_obj.h"
#include "whu_constants_and_definitions.h"
#include "fit_parameter_interface.h"

int main()
{

    DPFPWAPoint *pwa_point_phipp = new DPFPWAPoint(mka, mka, mpi, mpi, mpsip);
    DataObject phsp_pp("../newbase/phsp_pwa_pp_1000.dat", pwa_point_phipp);
    //DataObject phsp_pp_200000("../newbase/phsp_pwa_pp_200000.dat", pwa_point_phipp);

    DPFPWAPoint *pwa_point_phikk = new DPFPWAPoint(mka, mka, mka, mka, mpsip);
    DataObject phsp_kk("../newbase/phsp_pwa_kk_1000.dat", pwa_point_phikk);
    //DataObject phsp_kk_200000("../newbase/phsp_pwa_kk_200000.dat", pwa_point_phikk);

    FitParametersInterface::prepare_my_parameter_table("whu_pwa_parameter_table.txt");
    FitParametersInterface::information_of_parameter_table();
    FitParametersInterface::information_of_parameter_list_sent_to_minuit();


    FitParametersOfPhiPP phipp_parameters;
    FitParametersOfPhiKK phikk_parameters;

    phipp_parameters.act_resonance_980("f00980");
    phikk_parameters.act_resonance_f0("f01000");
    phikk_parameters.act_resonance_980("f00980");
    phikk_parameters.act_resonance_f2("f21270");

    FitParametersInterface::information_of_parameter_list_sent_to_minuit();

    cout << "There are " << phipp_parameters.number_of_amplitudes() << " amplitudes for phipp" << endl;
    cout << "There are " << phikk_parameters.number_of_amplitudes() << " amplitudes for phikk" << endl;

    phipp_parameters.shape_of_mapping();
    phikk_parameters.shape_of_mapping();


    return 0;
}
