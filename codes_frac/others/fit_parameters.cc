#include "fit_parameters.h"
#include <string>
#include <sstream>
#include "whu_my_parameter.h"

using namespace std;

//void FitParametersInterface::format_parameter(MyParameter & the_parameter, string name, double value)
//{
//    cout << "in my_parameter" << endl;
//    cout << name << endl;
//    cout << value << endl;
//    cout << "-----------------" << endl;
//    the_parameter.name = name;
//    the_parameter.value = value;
//    the_parameter.error = 0.;
//    the_parameter.lowlimit = value;
//    the_parameter.uplimit = value;
//    the_parameter.type = 0; // 固定参数
//}
//void FitParametersInterface::format_parameter(MyParameter & the_parameter, string name, double value, double error)
//{
//    cout << "in my_parameter" << endl;
//    cout << name << endl;
//    cout << value << endl;
//    cout << error << endl;
//    cout << "-----------------" << endl;
//    the_parameter.name = name;
//    the_parameter.value = value;
//    the_parameter.error = error;
//    the_parameter.lowlimit = -1e30;
//    the_parameter.uplimit = +1e30;
//    the_parameter.type = 1; // 自由参数
//}
//void FitParametersInterface::format_parameter(MyParameter & the_parameter, string name, double value, double error, double lowlimit, double uplimit)
//{
//    cout << "in my_parameter" << endl;
//    cout << name << endl;
//    cout << value << endl;
//    cout << error << endl;
//    cout << lowlimit << endl;
//    cout << uplimit << endl;
//    cout << "-----------------" << endl;
//    the_parameter.name = name;
//    the_parameter.value = value;
//    the_parameter.error = error;
//    the_parameter.lowlimit = lowlimit;
//    the_parameter.uplimit = uplimit;
//    if (lowlimit > -1e29 && uplimit < 1e29) {
//        the_parameter.type = 2; // 设置了上下界
//    } else if (lowlimit > -1e29) {
//        the_parameter.type = 3; // 设置了下界
//    } else if (uplimit < 1e29) {
//        the_parameter.type = 4; // 设置了上界
//    } else {
//        cout << "Formatting parameter " << name << " failed as limits cannot be comparable with 1e20!!!" << endl;
//        exit(1);
//    }
//}
void FitParametersInterface::string_to_vector(string str, vector<string> &vec)
{
    while(str.find(',') != string::npos) {
        int i  = str.find(',');
        vec.push_back(str.substr(0, i));
        str = str.substr(i + 1);
    }
    vec.push_back(str);
}
double FitParametersInterface::str2float(const string str)
{
    double _tmp;
    stringstream _ss(str);
    _ss >> _tmp;
    return _tmp;
}
void FitParametersInterface::prepare_my_parameter_table(const string table_file_name)
{
    fstream parameter_table_stream;
    //cfgFile.open("dfa");//打开文件
    const char* fname = table_file_name.c_str();
    parameter_table_stream.open(fname);//打开文件
    if( ! parameter_table_stream.is_open())
    {
        cout<<"can not open parameter table file " << table_file_name << " !!!"<<endl;
        exit(1);
    }
    char _tmp[10000];
    cout << "Begin read parameter table file-->" << endl;
    while(!parameter_table_stream.eof())//循环读取每一行
    {
        parameter_table_stream.getline(_tmp,1000);//每行读取前1000个字符
        string line(_tmp);
        cout << line << endl;
        size_t pos = line.find('#');
        if(pos==string::npos) continue;
        string fit_parameter_name = line.substr(0, pos);
        cout << "fit parmaeter " << fit_parameter_name  << " is getted!" << endl;
        string fit_parameter_specifies = line.substr(pos + 1);
        vector<string> specifi_list;
        string_to_vector(fit_parameter_specifies, specifi_list);
        //MyParameter _the_parameter;
        if (specifi_list.size() == 1) {
            myParameterTable_.push_back(MyParameter(fit_parameter_name, str2float(specifi_list[0])));
        } else if (specifi_list.size() == 2) {
            myParameterTable_.push_back(MyParameter(fit_parameter_name, str2float(specifi_list[0]), str2float(specifi_list[1])));
        } else if (specifi_list.size() == 4) {
            myParameterTable_.push_back(MyParameter(fit_parameter_name, str2float(specifi_list[0]), str2float(specifi_list[1]), str2float(specifi_list[2]), str2float(specifi_list[3])));
        } else {
            cout << "Format of the parameter " << fit_parameter_name << " is not correct in the parameter table file " << table_file_name << endl;
            exit(0);
        }
        myParameterTable_.back().info();
    }
    cout << "There are total " << myParameterTable_.size() << " parameters in this pwa analysis!" << endl;
}
int FitParametersInterface::position_in_fit_parameter_list(const MyParameter & the_parameter)
{
    for(int i = 0; i < fit_parameter_list_.size(); i++) {
        if (fit_parameter_list_[i].get_name() == the_parameter.get_name()) {
            return i;
        }
    }
    fit_parameter_list_.push_back(the_parameter);
    return fit_parameter_list_.size() - 1;
}
void FitParametersInterface::initialize_fit_parameter_mapping()
{
    fit_parameter_mapping_pp_.resize(end_category);
    for(int i = start_category; i < end_category; i++) {
        fit_parameter_mapping_pp_[i].resize(0);
    }

    //fit_parameter_mapping_pp_[spin_category].resize(0);
    //fit_parameter_mapping_pp_[mass_category].resize(0);
    //fit_parameter_mapping_pp_[mass2_category].resize(0);
    //fit_parameter_mapping_pp_[width_category].resize(0);
    //fit_parameter_mapping_pp_[g1_category].resize(0);
    //fit_parameter_mapping_pp_[g2_category].resize(0);
    //fit_parameter_mapping_pp_[b1_category].resize(0);
    //fit_parameter_mapping_pp_[b2_category].resize(0);
    //fit_parameter_mapping_pp_[b3_category].resize(0);
    //fit_parameter_mapping_pp_[b4_category].resize(0);
    //fit_parameter_mapping_pp_[b5_category].resize(0);
    //fit_parameter_mapping_pp_[rho_category].resize(0);
    //fit_parameter_mapping_pp_[frac_category].resize(0);
    //fit_parameter_mapping_pp_[phi_category].resize(0);
    //fit_parameter_mapping_pp_[propType_category].resize(0);

    fit_parameter_mapping_kk_.resize(end_category);
    for(int i = start_category; i < end_category; i++) {
        fit_parameter_mapping_kk_[i].resize(0);
    }

    //fit_parameter_mapping_kk_[spin_category].resize(0);
    //fit_parameter_mapping_kk_[mass_category].resize(0);
    //fit_parameter_mapping_kk_[mass2_category].resize(0);
    //fit_parameter_mapping_kk_[width_category].resize(0);
    //fit_parameter_mapping_kk_[g1_category].resize(0);
    //fit_parameter_mapping_kk_[g2_category].resize(0);
    //fit_parameter_mapping_kk_[b1_category].resize(0);
    //fit_parameter_mapping_kk_[b2_category].resize(0);
    //fit_parameter_mapping_kk_[b3_category].resize(0);
    //fit_parameter_mapping_kk_[b4_category].resize(0);
    //fit_parameter_mapping_kk_[b5_category].resize(0);
    //fit_parameter_mapping_kk_[rho_category].resize(0);
    //fit_parameter_mapping_kk_[frac_category].resize(0);
    //fit_parameter_mapping_kk_[phi_category].resize(0);
    //fit_parameter_mapping_kk_[propType_category].resize(0);
}
void FitParametersInterface::add_pp_amplitude_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& width, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_pp_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_pp_[mass_category].push_back(position_in_fit_parameter_list(mass));
    fit_parameter_mapping_pp_[width_category].push_back(position_in_fit_parameter_list(width));
    fit_parameter_mapping_pp_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_pp_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_pp_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_pp_[propType_category].push_back(position_in_fit_parameter_list(propType));
}
void FitParametersInterface::add_pp_amplitude1680_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass1680, const MyParameter& mass2, const MyParameter& width, const MyParameter& g1, const MyParameter& g2, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_pp_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_pp_[mass_category].push_back(position_in_fit_parameter_list(mass1680));
    fit_parameter_mapping_pp_[mass2_category].push_back(position_in_fit_parameter_list(mass2));
    fit_parameter_mapping_pp_[width_category].push_back(position_in_fit_parameter_list(width));
    fit_parameter_mapping_pp_[g1_category].push_back(position_in_fit_parameter_list(g1));
    fit_parameter_mapping_pp_[g2_category].push_back(position_in_fit_parameter_list(g2));
    fit_parameter_mapping_pp_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_pp_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_pp_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_pp_[propType_category].push_back(position_in_fit_parameter_list(propType));

}
void FitParametersInterface::add_pp_amplitude600_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& b1, const MyParameter& b2, const MyParameter& b3, const MyParameter& b4, const MyParameter& b5, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_pp_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_pp_[mass_category].push_back(position_in_fit_parameter_list(mass));
    fit_parameter_mapping_pp_[b1_category].push_back(position_in_fit_parameter_list(b1));
    fit_parameter_mapping_pp_[b2_category].push_back(position_in_fit_parameter_list(b2));
    fit_parameter_mapping_pp_[b3_category].push_back(position_in_fit_parameter_list(b3));
    fit_parameter_mapping_pp_[b4_category].push_back(position_in_fit_parameter_list(b4));
    fit_parameter_mapping_pp_[b5_category].push_back(position_in_fit_parameter_list(b5));
    fit_parameter_mapping_pp_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_pp_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_pp_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_pp_[propType_category].push_back(position_in_fit_parameter_list(propType));
}
void FitParametersInterface::add_pp_amplitude980_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& g1, const MyParameter& g2, const MyParameter& rho, const MyParameter& frac, const MyParameter&phi, const MyParameter& propType)
{
    fit_parameter_mapping_pp_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_pp_[mass_category].push_back(position_in_fit_parameter_list(mass));
    fit_parameter_mapping_pp_[g1_category].push_back(position_in_fit_parameter_list(g1));
    fit_parameter_mapping_pp_[g2_category].push_back(position_in_fit_parameter_list(g2));
    fit_parameter_mapping_pp_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_pp_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_pp_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_pp_[propType_category].push_back(position_in_fit_parameter_list(propType));
}
void FitParametersInterface::add_kk_amplitude_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& width, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_kk_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_kk_[mass_category].push_back(position_in_fit_parameter_list(mass));
    fit_parameter_mapping_kk_[width_category].push_back(position_in_fit_parameter_list(width));
    fit_parameter_mapping_kk_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_kk_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_kk_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_kk_[propType_category].push_back(position_in_fit_parameter_list(propType));
}
void FitParametersInterface::add_kk_amplitude1680_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass1680, const MyParameter& mass2, const MyParameter& width, const MyParameter& g1, const MyParameter& g2, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_kk_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_kk_[mass_category].push_back(position_in_fit_parameter_list(mass1680));
    fit_parameter_mapping_kk_[mass2_category].push_back(position_in_fit_parameter_list(mass2));
    fit_parameter_mapping_kk_[width_category].push_back(position_in_fit_parameter_list(width));
    fit_parameter_mapping_kk_[g1_category].push_back(position_in_fit_parameter_list(g1));
    fit_parameter_mapping_kk_[g2_category].push_back(position_in_fit_parameter_list(g2));
    fit_parameter_mapping_kk_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_kk_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_kk_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_kk_[propType_category].push_back(position_in_fit_parameter_list(propType));

}
void FitParametersInterface::add_kk_amplitude600_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& b1, const MyParameter& b2, const MyParameter& b3, const MyParameter& b4, const MyParameter& b5, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_kk_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_kk_[mass_category].push_back(position_in_fit_parameter_list(mass));
    fit_parameter_mapping_kk_[b1_category].push_back(position_in_fit_parameter_list(b1));
    fit_parameter_mapping_kk_[b2_category].push_back(position_in_fit_parameter_list(b2));
    fit_parameter_mapping_kk_[b3_category].push_back(position_in_fit_parameter_list(b3));
    fit_parameter_mapping_kk_[b4_category].push_back(position_in_fit_parameter_list(b4));
    fit_parameter_mapping_kk_[b5_category].push_back(position_in_fit_parameter_list(b5));
    fit_parameter_mapping_kk_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_kk_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_kk_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_kk_[propType_category].push_back(position_in_fit_parameter_list(propType));
}
void FitParametersInterface::add_kk_amplitude980_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& g1, const MyParameter& g2, const MyParameter& rho, const MyParameter& frac, const MyParameter&phi, const MyParameter& propType)
{
    fit_parameter_mapping_kk_[spin_category].push_back(position_in_fit_parameter_list(spin));
    fit_parameter_mapping_kk_[mass_category].push_back(position_in_fit_parameter_list(mass));
    fit_parameter_mapping_kk_[g1_category].push_back(position_in_fit_parameter_list(g1));
    fit_parameter_mapping_kk_[g2_category].push_back(position_in_fit_parameter_list(g2));
    fit_parameter_mapping_kk_[rho_category].push_back(position_in_fit_parameter_list(rho));
    fit_parameter_mapping_kk_[frac_category].push_back(position_in_fit_parameter_list(frac));
    fit_parameter_mapping_kk_[phi_category].push_back(position_in_fit_parameter_list(phi));
    fit_parameter_mapping_kk_[propType_category].push_back(position_in_fit_parameter_list(propType));
}
void FitParametersInterface::act_pp_resonance_1p(string rn)
{
    if (binding_phi) {
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
    } else {
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "p_phi3"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "p_phi4"), gp(rn + "a_typ_"));
    }
}
void FitParametersInterface::act_kk_resonance_1p(string rn)
{
    if (binding_phi) {
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
    } else {
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "k_phi3"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "k_phi4"), gp(rn + "a_typ_"));
    }
}
void FitParametersInterface::act_pp_resonance_980(string rn) {
    if (binding_phi) {
        add_pp_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_pp_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
    } else {
        add_pp_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_pp_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersInterface::act_kk_resonance_980(string rn) {
    if (binding_phi) {
        add_kk_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));

    } else {
        add_kk_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersInterface::act_pp_resonance_f0(string rn) {
    if (binding_phi) {
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
    } else {
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersInterface::act_kk_resonance_f0(string rn) {
    if (binding_phi) {
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
    } else {
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersInterface::act_pp_resonance_f2(string rn) {
    if (binding_phi) {
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc5"), gp(rn + "a_phi5"), gp(rn + "a_typ_"));
    } else {
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "p_phi3"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "p_phi4"), gp(rn + "a_typ_"));
        add_pp_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc5"), gp(rn + "p_phi5"), gp(rn + "a_typ_"));
    }
}
void FitParametersInterface::act_kk_resonance_f2(string rn) {
    if (binding_phi) {
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc5"), gp(rn + "a_phi5"), gp(rn + "a_typ_"));
    } else {
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "k_phi3"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "k_phi4"), gp(rn + "a_typ_"));
        add_kk_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc5"), gp(rn + "k_phi5"), gp(rn + "a_typ_"));
    }
}

FitParametersInterface::FitParametersInterface(const PWA_CTRL & pwa_ctrl)
{

}
//void FitParameters::addResonance980(const Char_t* name, const Char_t* title, RooAbsReal& spin,
//        RooAbsReal& mass, RooAbsReal& g1,RooAbsReal& g2,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType)
//{       _spinList.add(spin);
//    //  //cout<<"haha: "<< __LINE__ << endl;
//    _massList.add(mass);
//    _g1List.add(g1);
//    _g2List.add(g2);
//    _rhoList.add(rho);
//    _fracList.add(frac);
//    _phiList.add(phi);
//    _propList.add(propType);
//    nameList[nAmps]= name;
//    titleListT[nAmps]= title;
//    titleList[nStates]= title;
//    nAmps++;
//    nStates++;
//    nStatesg1g2++;
//}
//void DPFPWAPdf::addResonance1680(const Char_t* name, const Char_t* title, RooAbsReal& spin,
//        RooAbsReal& mass1680, RooAbsReal& mass2,RooAbsReal& width,RooAbsReal& g1,RooAbsReal& g2,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,
//        RooAbsReal& propType)
//{       _spinList.add(spin);
//    //  //cout<<"haha: "<< __LINE__ << endl;
//    _massList.add(mass1680);
//    _mass2List.add(mass2);
//    _widthList.add(width);
//    _g1List.add(g1);
//    _g2List.add(g2);
//    _rhoList.add(rho);
//    _fracList.add(frac);
//    _phiList.add(phi);
//    _propList.add(propType);
//    nameList[nAmps]= name;
//    titleListT[nAmps]= title;
//    titleList[nStates]= title;
//    nAmps++;
//    nStates++;
//    nStatesg1g2++;
//    nStatesmass2++;
//    nStateswidth++;
//}
//void DPFPWAPdf::addResonance(const string name, const string title, RooAbsReal& spin,
//        RooAbsReal& mass, RooAbsReal& width,RooAbsReal& rho, RooAbsReal& frac, RooAbsReal& phi,RooAbsReal& propType)
//{       _spinList.add(spin);
//    //  //cout<<"haha: "<< __LINE__ << endl;
//    _massList.add(mass);
//    _widthList.add(width);
//    _rhoList.add(rho);
//    _fracList.add(frac);
//    _phiList.add(phi);
//    _propList.add(propType);
//    nameList[nAmps]= name;
//    titleListT[nAmps]= title;
//    titleList[nStates]= title;
//    nAmps++;
//    nStates++;
//    nStateswidth++;
//}
void FitParametersInterface::shape_of_mapping_kk()
{
    cout << "The transported parameters to minuit2: ";
    for(int i = 0; i < fit_parameter_list_.size(); i++)
    {
        cout << fit_parameter_list_[i].get_name() << " ";
    }
    cout << endl;
    for(int i = 0; i < end_category; i++)
    {
        cout << i << " ";
        for(int j = 0; j < fit_parameter_mapping_kk_[i].size(); j++)
        {
            cout << fit_parameter_list_[fit_parameter_mapping_kk_[i][j]].get_name() << " ";
        }
        cout << endl;
    }
}
void FitParametersInterface::shape_of_mapping_pp()
{
    for(int i = 0; i < end_category; i++)
    {
        cout << i << " ";
        for(int j = 0; j < fit_parameter_mapping_pp_[i].size(); j++)
        {
            cout << fit_parameter_list_[fit_parameter_mapping_pp_[i][j]].get_name() << " ";
        }
        cout << endl;
    }

}
