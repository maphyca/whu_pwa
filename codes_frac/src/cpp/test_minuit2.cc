#include <iostream>
#include <fstream>
#include <string>

#include <sstream>
#include "TString.h"
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
#include "data_obj.h"
#include "whu_constants_and_definitions.h"
#include "fit_parameter_interface.h"
#include "whu_propogator.h"

using namespace ROOT::Minuit2;

void setup_mnuserparameters(MnUserParameters& upar, vector<int> & parameter_list_for_minuit)
{
    cout << parameter_list_for_minuit.size() << "*****" << endl;
    for(unsigned i = 0; i < parameter_list_for_minuit.size(); i++) {
        MyParameter the_parameter(FitParametersInterface::myParameterTable()[parameter_list_for_minuit[i]]);
        upar.Add(the_parameter.get_name(), the_parameter.get_value(), 0.1);
        cout << i << "***" << the_parameter.get_name() << endl;
    }
}
void update_mnuserparameters(MnUserParameters& upar, vector<int> & parameter_list_for_minuit)
{
    for(unsigned i = 0; i < parameter_list_for_minuit.size(); i++) {
        MyParameter the_parameter(FitParametersInterface::myParameterTable()[parameter_list_for_minuit[i]]);
        if (the_parameter.get_type() == both_limits) upar.SetLimits(the_parameter.get_name(), the_parameter.get_lowlimit(), the_parameter.get_uplimit());
        else if (the_parameter.get_type() == only_lowlimit) upar.SetLowerLimit(the_parameter.get_name(), the_parameter.get_lowlimit());
        else if (the_parameter.get_type() == only_uplimit) upar.SetUpperLimit(the_parameter.get_name(), the_parameter.get_uplimit());
    }
}
void update_mnmigrad(MnMigrad & migrad, vector<int> & parameter_list_for_minuit)
{
    for(unsigned i = 0; i < parameter_list_for_minuit.size(); i++) {
        MyParameter the_parameter(FitParametersInterface::myParameterTable()[parameter_list_for_minuit[i]]);
        if (the_parameter.get_type() == fixxd) migrad.Fix(i);
        //else if (the_parameter.get_type() == freed) migrad.Release(i);
    }
}

int main()
{

    DPFPWAPoint *pwa_point_phipp = new DPFPWAPoint(mka, mka, mpi, mpi, mpsip);
    DPFPWAPoint *pwa_point_phikk = new DPFPWAPoint(mka, mka, mka, mka, mpsip);

    vector<DataObject*> data_set(end_data_object_index);
    data_set[phipp_phsp_index] = new CPUWaveFunc("../../newbase/phsp_pwa_pp_200000.dat", pwa_point_phipp);
    data_set[phipp_data_index] = new CPUWaveFunc("../../newbase/data_pwa_pipi_weight_superlength_all.dat", pwa_point_phipp);
    //DataObject phsp_pp_200000("../newbase/phsp_pwa_pp_200000.dat", pwa_point_phipp);

    data_set[phikk_phsp_index] = new CPUWaveFunc("../../newbase/phsp_pwa_kk_200000.dat", pwa_point_phikk);
    data_set[phikk_data_index] = new CPUWaveFunc("../../newbase/data_pwa_kk_weight_superlength_all.dat", pwa_point_phikk);
    //DataObject phsp_kk_200000("../newbase/phsp_pwa_kk_200000.dat", pwa_point_phikk);



    FitParametersInterface::prepare_my_parameter_table("../resource/whu_pwa_parameter_table.txt");
    FitParametersInterface::information_of_parameter_table();
    //FitParametersInterface::information_of_parameter_list_sent_to_minuit();


    vector<AmplitudeMethodWithFitParametersInterface*> parameter_list_set(end_list_index, NULL);
    parameter_list_set[phipp_list_index]=new FitParametersOfPhiPP();
    parameter_list_set[phikk_list_index]=new FitParametersOfPhiKK();

    vector<string> resonances = {
        "f00980",
        "1p1800",
        "1m1800",
        "f01000",
        "f01370",
        "f01500",
        "f01710",
        "f02020",
        "f02100",
        "f02200",
        "f02330",
        "f21000",
        "f21270",
        "f21525",
        "f21501",
        "f21810",
        "f21910",
        "f21950",
        "f22150",
        "f22300"};
    vector<vector<string> > resonance_list_set(end_list_index);
    resonance_list_set[phipp_list_index] = resonances;
    resonance_list_set[phikk_list_index] = resonances;

    ((FitParametersOfPhiPP*)parameter_list_set[phipp_list_index])->act_resonances(resonance_list_set[phipp_list_index]);
    ((FitParametersOfPhiKK*)parameter_list_set[phikk_list_index])->act_resonances(resonance_list_set[phikk_list_index]);

    vector<int> parameter_list_for_minuit;
    for(unsigned i = 0; i < parameter_list_set.size(); i++) {
        if (parameter_list_set[i] == NULL) continue;
        parameter_list_set[i]->remap_local_mapping_for_minuit(parameter_list_for_minuit);
        //parameter_list_set[i]->create_category_tags();
        parameter_list_set[i]->create_minuit_mapping();
    }
    for(unsigned i = 0; i < parameter_list_set.size(); i++) {
        if (parameter_list_set[i] == NULL) continue;
        int _id = 0;
            for(unsigned k = 0; k < parameter_list_set[i]->number_of_amplitudes(); k++)
        {
        for(int j = 0; j < end_category; j++)
            {
              cout << FitParametersInterface::get_my_parameter_table()[parameter_list_for_minuit[parameter_list_set[i]->get_minuit_mapping()[_id++]]].get_name() << " == " << FitParametersInterface::get_my_parameter_table()[parameter_list_for_minuit[parameter_list_set[i]->get_local_fit_parameter_mapping()[j][k]]].get_name() << endl;
            }
        }
    }
    for(unsigned i = 0; i < parameter_list_set.size(); i++) {
        if (parameter_list_set[i] == NULL) continue;
        parameter_list_set[i]->shape_of_local_mapping(parameter_list_for_minuit);
    }
    //for(unsigned i = 0; i < parameter_list_set.size(); i++) {
    //    if (parameter_list_set[i] == NULL) continue;
    //    cout << "****" << parameter_list_set[i]->get_minuit_mapping().size() << endl;
    //    for(unsigned j = 0; j < parameter_list_set[i]->get_minuit_mapping().size(); j++) {
    //        cout << parameter_list_set[i]->get_minuit_mapping()[j] << endl;
    //    }
    //}


    vector<kernel*> kernel_set(end_data_object_index);

    kernel_set[phipp_phsp_index] = new kernel(((CPUWaveFunc*)data_set[phipp_phsp_index])->hpv,0,0,((CPUWaveFunc*)data_set[phipp_phsp_index])->number_of_events(),parameter_list_set[phipp_list_index]->number_of_amplitudes(),((CPUWaveFunc*)data_set[phipp_phsp_index])->number_of_events());
    kernel_set[phipp_data_index] = new kernel(((CPUWaveFunc*)data_set[phipp_data_index])->hpv,0,0,((CPUWaveFunc*)data_set[phipp_data_index])->number_of_events(),parameter_list_set[phipp_list_index]->number_of_amplitudes(),((CPUWaveFunc*)data_set[phipp_data_index])->number_of_events());
    kernel_set[phikk_phsp_index] = new kernel(((CPUWaveFunc*)data_set[phikk_phsp_index])->hpv,0,0,((CPUWaveFunc*)data_set[phikk_phsp_index])->number_of_events(),parameter_list_set[phikk_list_index]->number_of_amplitudes(),((CPUWaveFunc*)data_set[phikk_phsp_index])->number_of_events());
    kernel_set[phikk_data_index] = new kernel(((CPUWaveFunc*)data_set[phikk_data_index])->hpv,0,0,((CPUWaveFunc*)data_set[phikk_data_index])->number_of_events(),parameter_list_set[phikk_list_index]->number_of_amplitudes(),((CPUWaveFunc*)data_set[phikk_data_index])->number_of_events());

     /*   
    std::vector<double> test_par(68*end_category),test_par_back(68*end_category);
    for(int i=0;i<68*end_category;i++)
      {
        if(i%end_category==propType_category){
          test_par[i]=2;
          test_par_back[i]=2;
        }
        else{
          test_par[i]=rand()/(double)RAND_MAX;
          test_par_back[i]=rand()/(double)RAND_MAX;
        }
      }
    ((CPUWaveFunc*) data_set[phikk_phsp_index])->cpu_resize_intermediate_variables(68);
    ((CPUWaveFunc*) data_set[phikk_phsp_index])->cpu_calEva(test_par,test_par_back,68);
    kernel_set[phikk_phsp_index]->par_trans(test_par);
    kernel_set[phikk_phsp_index]->calEva();
    double gpu_phsp=kernel_set[phikk_phsp_index]->sum_phsp();

    double test_phsp=((CPUWaveFunc*) data_set[phikk_phsp_index])->sum_phsp(68);
    cout<<"cpu phsp : "<<test_phsp<<"   gpu phsp : "<<gpu_phsp<<endl;

    */

    {
    MnUserParameters upar;
    setup_mnuserparameters(upar, parameter_list_for_minuit);
    update_mnuserparameters(upar, parameter_list_for_minuit);

    PWAFcn my_pwa_fcn(data_set, parameter_list_set,kernel_set, parameter_list_for_minuit.size());

    MnMigrad migrad(my_pwa_fcn, upar);

    update_mnmigrad(migrad, parameter_list_for_minuit);

    FunctionMinimum min = migrad();
    }
            for(int i = 0; i < end_list_index; i++) {
                if (parameter_list_set[i] == NULL) continue;
                parameter_list_set[i]->shape_of_mapping();
                parameter_list_set[i]->shape_of_minuit_parameters();
            }

            //FitParametersInterface::information_of_parameter_list_sent_to_minuit();*/
    return 0;
}
