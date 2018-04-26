#ifndef CONSTANTS_AND_DEFINITIONS
#define CONSTANTS_AND_DEFINITIONS

const double mpsip=3.686,mka=0.493677,mpi=0.13957;

const bool binding_phi = true;

//struct MyParameter {
//    TString name;
//    double value;
//    double error;
//    double uplimit, lowlimit;
//    int type;
//};

enum ParameterCategory {
    start_category,
    spin_category, mass_category, mass2_category, width_category,
    g1_category, g2_category, b1_category, b2_category, b3_category,
    b4_category, b5_category, rho_category,
    frac_category, phi_category, propType_category,
    end_category
};

#endif

