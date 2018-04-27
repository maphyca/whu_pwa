
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TString.h"
#include "data_obj.h"
#include "whu_constants_and_definitions.h"
#include "fit_parameter_interface.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include "whu_pwa_fcn.h"

using namespace ROOT::Minuit2;

int main()
{

    DPFPWAPoint *pwa_point_phipp = new DPFPWAPoint(mka, mka, mpi, mpi, mpsip);
    DPFPWAPoint *pwa_point_phikk = new DPFPWAPoint(mka, mka, mka, mka, mpsip);

    vector<DataObject*> data_set(end_data_object_index);
    data_set[phipp_phsp_index] = new DataObject("../newbase/phsp_pwa_pp_1000.dat", pwa_point_phipp);
    data_set[phipp_data_index] = new DataObject("../newbase/data_pwa_pipi_weight_superlength_all.dat", pwa_point_phipp);
    //DataObject phsp_pp_200000("../newbase/phsp_pwa_pp_200000.dat", pwa_point_phipp);

    data_set[phikk_phsp_index] = new DataObject("../newbase/phsp_pwa_kk_1000.dat", pwa_point_phikk);
    data_set[phikk_data_index] = new DataObject("../newbase/data_pwa_kk_weight_superlength_all.dat", pwa_point_phikk);
    //DataObject phsp_kk_200000("../newbase/phsp_pwa_kk_200000.dat", pwa_point_phikk);


    FitParametersInterface::prepare_my_parameter_table("whu_pwa_parameter_table.txt");
    FitParametersInterface::information_of_parameter_table();
    //FitParametersInterface::information_of_parameter_list_sent_to_minuit();


    vector<AmplitudeMethodWithFitParametersInterface*> parameter_list_set(end_list_index, NULL);
    parameter_list_set[phipp_list_index]=new FitParametersOfPhiPP();
    parameter_list_set[phikk_list_index]=new FitParametersOfPhiKK();

    ((FitParametersOfPhiPP*)parameter_list_set[phipp_list_index])->act_resonance_980("f00980");
    ((FitParametersOfPhiPP*)parameter_list_set[phipp_list_index])->act_resonance_1p("1p1800");
    ((FitParametersOfPhiPP*)parameter_list_set[phipp_list_index])->act_resonance_1m("1m1800");
    ((FitParametersOfPhiKK*)parameter_list_set[phikk_list_index])->act_resonance_f0("f01000");
    ((FitParametersOfPhiKK*)parameter_list_set[phikk_list_index])->act_resonance_980("f00980");
    ((FitParametersOfPhiKK*)parameter_list_set[phikk_list_index])->act_resonance_f2("f21270");

    vector<int> parameter_list_for_minuit(end_list_index);
    for(unsigned i = 0; i < parameter_list_set.size(); i++) {
        if (parameter_list_set[i] == NULL) continue;
        parameter_list_set[i]->remap_local_mapping_for_minuit(parameter_list_for_minuit);
    }
    //FitParametersInterface::information_of_parameter_list_sent_to_minuit();
    //
    parameter_list_set[phipp_list_index]->shape_of_mapping();
    parameter_list_set[phipp_list_index]->shape_of_local_mapping(parameter_list_for_minuit);
    parameter_list_set[phikk_list_index]->shape_of_mapping();
    parameter_list_set[phikk_list_index]->shape_of_local_mapping(parameter_list_for_minuit);
    exit(1);




    PWAFcn my_pwa_fcn(data_set, parameter_list_set);

    //FitParametersInterface::information_of_parameter_list_sent_to_minuit();
    return 0;
}
