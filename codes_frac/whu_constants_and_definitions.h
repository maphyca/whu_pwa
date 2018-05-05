#ifndef CONSTANTS_AND_DEFINITIONS
#define CONSTANTS_AND_DEFINITIONS

#include <string>

const double mpsip=3.686,mka=0.493677,mpi=0.13957;

const bool binding_phi = true;

const int kBlank = -1;
const std::string kBlankParameter = "xxxxxxxxxxxx";
const double kParameterLimit = 1e30;

const double rk=0.493677;
const double rp=0.139556995;


enum ParameterCategory {
    //start_category,
    spin_category, mass_category, mass2_category, width_category,
    g1_category, g2_category, b1_category, b2_category, b3_category,
    b4_category, b5_category, rho_category,
    frac_category, phi_category, propType_category,
    end_category
};

enum DataObjectIndex {
    phipp_data_index,
    phipp_phsp_index,
    phikk_data_index,
    phikk_phsp_index,
    end_data_object_index
};

enum ParameterListIndex {
    phipp_list_index,
    phikk_list_index,
    end_list_index
};
enum ParameterType {
    fixxd,
    freed,
    only_uplimit,
    only_lowlimit,
    both_limits
};
#endif

