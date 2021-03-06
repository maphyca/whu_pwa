#include "fit_parameter_interface.h"

using namespace std;

void string_to_vector(string str, vector<string> &vec)
{
    while(str.find(',') != string::npos) {
        int i  = str.find(',');
        vec.push_back(str.substr(0, i));
        str = str.substr(i + 1);
    }
    vec.push_back(str);
}
double str2float(const string str)
{
    double _tmp;
    stringstream _ss(str);
    _ss >> _tmp;
    return _tmp;
}

vector<MyParameter> FitParametersInterface::my_parameter_table_(0);
//vector<MyParameter> FitParametersInterface::fit_parameter_list_(0);

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
            my_parameter_table_.push_back(MyParameter(fit_parameter_name, str2float(specifi_list[0])));
        } else if (specifi_list.size() == 2) {
            my_parameter_table_.push_back(MyParameter(fit_parameter_name, str2float(specifi_list[0]), str2float(specifi_list[1])));
        } else if (specifi_list.size() == 4) {
            my_parameter_table_.push_back(MyParameter(fit_parameter_name, str2float(specifi_list[0]), str2float(specifi_list[1]), str2float(specifi_list[2]), str2float(specifi_list[3])));
        } else {
            cout << "Format of the parameter " << fit_parameter_name << " is not correct in the parameter table file " << table_file_name << endl;
            exit(0);
        }
        //my_parameter_table_.back().info();
    }
    cout << "There are total " << my_parameter_table_.size() << " parameters in this pwa analysis!" << endl;
}
MyParameter FitParametersInterface::gp(std::string parameter_name) {
    for(unsigned i = 0; i < my_parameter_table_.size(); i++) {
        if (parameter_name == my_parameter_table_[i].get_name()) {
            return my_parameter_table_[i];
        }
    }
    std::cout << "Cannot find the parameter " << parameter_name << " in myParameterTable!!" << std::endl;
    exit(1);
};
int FitParametersInterface::position_in_my_parameter_table(const MyParameter & the_parameter)
{
    for(unsigned i = 0; i < my_parameter_table_.size(); i++) {
        if (my_parameter_table_[i].get_name() == the_parameter.get_name()) {
            return i;
        }
    }
    exit(1);
    return -1;
}
void FitParametersInterface::information_of_parameter_table()
{
    if (my_parameter_table_.size() == 0 ) {
        cout << "There is no parameter in the parameter table!" << endl;
    } else {
        cout << "All fit parameters in parameter table: " << endl;
        for(unsigned i = 0; i < my_parameter_table_.size(); i++)
        {
            cout << my_parameter_table_[i].get_name() << " ";
        }
        cout << endl;
    }
}
//void FitParametersInterface::information_of_parameter_list_sent_to_minuit()
//{
//    if (fit_parameter_list_.size() == 0) {
//        cout << "There is no parameter in the list sent to minuit!" << endl;
//    } else {
//        cout << "All fit parameters in the list sent to minuit: " << endl;
//        for(unsigned i = 0; i < fit_parameter_list_.size(); i++)
//        {
//            cout << fit_parameter_list_[i].get_name() << " ";
//        }
//        cout << endl;
//    }
//}
unsigned AmplitudeMethodWithFitParametersInterface::number_of_amplitudes() {
    return fit_parameter_mapping_[0].size();
    //int _num = 0;
    //for(int i = 0; i < end_category; i++)
    //{
    //    _num += fit_parameter_mapping_[i].size();
    //}
    //return _num;
}
void AmplitudeMethodWithFitParametersInterface::shape_of_minuit_parameters()
{
    for(int i = 0; i < end_category; i++) {
        cout << i << " ";
        for(unsigned j = 0; j < number_of_amplitudes(); j++) {
            cout << minuit_parameters_[i + end_category * j] << " ";

        }
        cout << endl;
    }
}
void AmplitudeMethodWithFitParametersInterface::assignment_of_minuit_parameters_from_par(const std::vector<double> &par)
{
    for(unsigned j = 0; j < minuit_parameters_.size(); j++)
    {
        if(minuit_mapping_[j] == kBlank)
            minuit_parameters_[j] = kParameterLimit;
        else
            minuit_parameters_[j] = par[minuit_mapping_[j]];
    }
    //minuit_parameters_back_.resize(minuit_parameters_.size(), -kParameterLimit);
}
void AmplitudeMethodWithFitParametersInterface::create_minuit_mapping()
{
    minuit_mapping_.resize(number_of_amplitudes() * end_category);
    int _id = 0;
    for(unsigned k = 0; k < number_of_amplitudes(); k++)
    {
        for(int j = 0; j < end_category; j++)
        {
            minuit_mapping_[_id++] = local_fit_parameter_mapping_[j][k];
        }
    }
    minuit_parameters_.resize(number_of_amplitudes() * end_category);
    minuit_parameters_back_.resize(number_of_amplitudes() * end_category);
}

void AmplitudeMethodWithFitParametersInterface::initialize_fit_parameter_mapping()
{
    fit_parameter_mapping_.resize(end_category);
    local_fit_parameter_mapping_.resize(end_category);
    for(int i = 0; i < end_category; i++) {
        fit_parameter_mapping_[i].resize(0);
        local_fit_parameter_mapping_[i].resize(0);
    }
    cout << "fit_parameter_mapping_ is initialized and its size is " << fit_parameter_mapping_.size() << endl;
}
int AmplitudeMethodWithFitParametersInterface::position_in_parameter_list_for_minuit(int the_mapping, std::vector<int>& parameter_list_for_minuit)
{
    for(unsigned i = 0; i < parameter_list_for_minuit.size(); i++) {
        if (the_mapping == parameter_list_for_minuit[i]) return i;
    }
    parameter_list_for_minuit.push_back(the_mapping);
    return parameter_list_for_minuit.size() - 1;
}
void AmplitudeMethodWithFitParametersInterface::remap_local_mapping_for_minuit(std::vector<int>& parameter_list_for_minuit)
{
    for(int i = 0; i < end_category; i++) {
        local_fit_parameter_mapping_[i].resize(fit_parameter_mapping_[i].size());
        for(unsigned j = 0; j < fit_parameter_mapping_[i].size(); j++) {
            if(fit_parameter_mapping_[i][j] == kBlank)
                local_fit_parameter_mapping_[i][j] = kBlank;
            else
                local_fit_parameter_mapping_[i][j] = position_in_parameter_list_for_minuit(fit_parameter_mapping_[i][j], parameter_list_for_minuit);
        }
    }
//    create_reduction_list_of_propogator_computation();
}
//void AmplitudeMethodWithFitParametersInterface::create_reduction_list_of_propogator_computation()
//{
//    int crp_id(0);
//    for(unsigned i = 0; i < number_of_amplitudes(); i++)
//    {
//        int _pt = my_parameter_table_(fit_parameter_mapping_[propType_category][i]).get_value();
//        switch(_pt)
//        {
//            case 1:
//            case 2:
//                crp_mapping_[i] = crp_id;
//                crp_mapping_[i + 1] = crp_id;
//                crp_id++;
//                ignore_switch_[i] = true;
//                ignore_switch_[i + 1] = false;
//                i = i + 1;
//                break;
//            case 4:
//                crp_mapping_[i] = crp_id;
//                crp_mapping_[i + 1] = crp_id;
//                crp_mapping_[i + 2] = crp_id;
//                crp_mapping_[i + 3] = crp_id;
//                crp_id++;
//                ignore_switch_[i] = true;
//                ignore_switch_[i + 1] = false;
//                ignore_switch_[i + 2] = false;
//                ignore_switch_[i + 3] = false;
//                i = i + 3;
//                break;
//            case 6:
//                crp_mapping_[i] = crp_id;
//                crp_mapping_[i + 1] = crp_id;
//                crp_mapping_[i + 2] = crp_id;
//                crp_mapping_[i + 3] = crp_id;
//                crp_mapping_[i + 4] = crp_id;
//                crp_id++;
//                ignore_switch_[i] = true;
//                ignore_switch_[i + 1] = false;
//                ignore_switch_[i + 2] = false;
//                ignore_switch_[i + 3] = false;
//                ignore_switch_[i + 4] = false;
//                i = i + 4;
//                break;
//            default:
//                cout << "Do not know how to set the reduction list of propogator computation with prop type " << _pt << endl;
//                exit(1);
//        }
//    }
//    cout << "The crp mapping is " << endl;
//    for(unsigned i = 0; i < number_of_amplitudes(); i++)
//    {
//        cout << crp_mapping_[i] << " ";
//    }
//    cout << endl;
//    cout << "The ignore switch mapping is " << endl;
//    for(unsigned i = 0; i < number_of_amplitudes(); i++)
//    {
//        cout << ignore_switch_[i] << " ";
//    }
//    cout << endl;
//}
//void AmplitudeMethodWithFitParametersInterface::create_category_tags()
//{
//    category_tags_.resize(end_category);
//    category_tags_[start_category] = 0;
//    for(int i = start_category + 1; i < end_category; i++)
//    {
//        category_tags_[i] = category_tags_[i - 1] + local_fit_parameter_mapping_[i - 1].size();
//    }
//}
void AmplitudeMethodWithFitParametersInterface::shape_of_mapping()
{
    cout << "The parameters transported to minuit2: " << endl;
    for(int i = 0; i < end_category; i++)
    {
        cout << i << " ";
        for(unsigned j = 0; j < fit_parameter_mapping_[i].size(); j++)
        {
            if(fit_parameter_mapping_[i][j] == kBlank)
                cout << kBlankParameter << " ";
            else
                cout << my_parameter_table_[fit_parameter_mapping_[i][j]].get_name() << " ";
        }
        cout << endl;
    }
}
void AmplitudeMethodWithFitParametersInterface::shape_of_local_mapping(vector<int> &parameter_list_for_minuit)
{
    cout << "The parameters for local mapping transported to minuit2: " << endl;
    for(int i = 0; i < end_category; i++)
    {
        cout << i << " ";
        for(unsigned j = 0; j < fit_parameter_mapping_[i].size(); j++)
        {
            if(local_fit_parameter_mapping_[i][j] == kBlank)
                cout << kBlankParameter << " ";
            else
                cout << my_parameter_table_[parameter_list_for_minuit[local_fit_parameter_mapping_[i][j]]].get_name() << " ";
        }
        cout << endl;
    }
}
bool AmplitudeMethodWithFitParametersInterface::already_active(string resonance_name)
{
    for(unsigned i = 0; i < resonance_name_list_.size(); i++) {
        cout << resonance_name << " has been already activated!" << endl;
        if (resonance_name == resonance_name_list_[i]) return true;
    }
    resonance_name_list_.push_back(resonance_name);
    return false;
}
void AmplitudeMethodWithFitParametersInterface::fill_blank_of_fit_parameter_mapping()
{
    unsigned _max_length(0);
    for(int i = 0; i < end_category; i++) {
        if(_max_length < fit_parameter_mapping_[i].size()) _max_length = fit_parameter_mapping_[i].size();
    }
    for(int i = 0; i < end_category; i++) {
        if(_max_length > fit_parameter_mapping_[i].size()) fit_parameter_mapping_[i].push_back(kBlank);
    }
    for(int i = 0; i < end_category; i++) {
        if(_max_length != fit_parameter_mapping_[i].size()) {
            cout << "Fails filling blannk of fit parameter mapping!!!!" << endl;
            exit(1);
        }
    }
}
void AmplitudeMethodWithFitParametersInterface::add_amplitude_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& width, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_[spin_category].push_back(position_in_my_parameter_table(spin));
    fit_parameter_mapping_[mass_category].push_back(position_in_my_parameter_table(mass));
    fit_parameter_mapping_[width_category].push_back(position_in_my_parameter_table(width));
    fit_parameter_mapping_[rho_category].push_back(position_in_my_parameter_table(rho));
    fit_parameter_mapping_[frac_category].push_back(position_in_my_parameter_table(frac));
    fit_parameter_mapping_[phi_category].push_back(position_in_my_parameter_table(phi));
    fit_parameter_mapping_[propType_category].push_back(position_in_my_parameter_table(propType));
    fill_blank_of_fit_parameter_mapping();
}
void AmplitudeMethodWithFitParametersInterface::add_amplitude1680_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass1680, const MyParameter& mass2, const MyParameter& width, const MyParameter& g1, const MyParameter& g2, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_[spin_category].push_back(position_in_my_parameter_table(spin));
    fit_parameter_mapping_[mass_category].push_back(position_in_my_parameter_table(mass1680));
    fit_parameter_mapping_[mass2_category].push_back(position_in_my_parameter_table(mass2));
    fit_parameter_mapping_[width_category].push_back(position_in_my_parameter_table(width));
    fit_parameter_mapping_[g1_category].push_back(position_in_my_parameter_table(g1));
    fit_parameter_mapping_[g2_category].push_back(position_in_my_parameter_table(g2));
    fit_parameter_mapping_[rho_category].push_back(position_in_my_parameter_table(rho));
    fit_parameter_mapping_[frac_category].push_back(position_in_my_parameter_table(frac));
    fit_parameter_mapping_[phi_category].push_back(position_in_my_parameter_table(phi));
    fit_parameter_mapping_[propType_category].push_back(position_in_my_parameter_table(propType));
    fill_blank_of_fit_parameter_mapping();
}
void AmplitudeMethodWithFitParametersInterface::add_amplitude600_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& b1, const MyParameter& b2, const MyParameter& b3, const MyParameter& b4, const MyParameter& b5, const MyParameter& rho, const MyParameter& frac, const MyParameter& phi, const MyParameter& propType)
{
    fit_parameter_mapping_[spin_category].push_back(position_in_my_parameter_table(spin));
    fit_parameter_mapping_[mass_category].push_back(position_in_my_parameter_table(mass));
    fit_parameter_mapping_[b1_category].push_back(position_in_my_parameter_table(b1));
    fit_parameter_mapping_[b2_category].push_back(position_in_my_parameter_table(b2));
    fit_parameter_mapping_[b3_category].push_back(position_in_my_parameter_table(b3));
    fit_parameter_mapping_[b4_category].push_back(position_in_my_parameter_table(b4));
    fit_parameter_mapping_[b5_category].push_back(position_in_my_parameter_table(b5));
    fit_parameter_mapping_[rho_category].push_back(position_in_my_parameter_table(rho));
    fit_parameter_mapping_[frac_category].push_back(position_in_my_parameter_table(frac));
    fit_parameter_mapping_[phi_category].push_back(position_in_my_parameter_table(phi));
    fit_parameter_mapping_[propType_category].push_back(position_in_my_parameter_table(propType));
    fill_blank_of_fit_parameter_mapping();
}
void AmplitudeMethodWithFitParametersInterface::add_amplitude980_to_fit_parameter_list(const MyParameter& spin, const MyParameter& mass, const MyParameter& g1, const MyParameter& g2, const MyParameter& rho, const MyParameter& frac, const MyParameter&phi, const MyParameter& propType)
{
    fit_parameter_mapping_[spin_category].push_back(position_in_my_parameter_table(spin));
    fit_parameter_mapping_[mass_category].push_back(position_in_my_parameter_table(mass));
    fit_parameter_mapping_[g1_category].push_back(position_in_my_parameter_table(g1));
    fit_parameter_mapping_[g2_category].push_back(position_in_my_parameter_table(g2));
    fit_parameter_mapping_[rho_category].push_back(position_in_my_parameter_table(rho));
    fit_parameter_mapping_[frac_category].push_back(position_in_my_parameter_table(frac));
    fit_parameter_mapping_[phi_category].push_back(position_in_my_parameter_table(phi));
    fit_parameter_mapping_[propType_category].push_back(position_in_my_parameter_table(propType));
    fill_blank_of_fit_parameter_mapping();
}
void FitParametersOfPhiPP::act_resonance(string rn)
{
    if (rn.find("f00980") != string::npos) act_resonance_980(rn);
    else if (rn.find("f0") != string::npos) act_resonance_f0(rn);
    else if (rn.find("f2") != string::npos) act_resonance_f2(rn);
    else if (rn.find("1p") != string::npos) act_resonance_1p(rn);
    else if (rn.find("1m") != string::npos) act_resonance_1m(rn);
    else {
        cout << "Not know how to act resonance " << rn << endl;
        exit(1);
    }
    //if (rn.find("1m") != string::npos) {
    //    act_resonance_1m(rn);
    //    continue;
    //}

}
void FitParametersOfPhiKK::act_resonance(string rn)
{
    if (rn.find("f00980") != string::npos) act_resonance_980(rn);
    else if (rn.find("f0") != string::npos) act_resonance_f0(rn);
    else if (rn.find("f2") != string::npos) act_resonance_f2(rn);
    else if (rn.find("1p") != string::npos) act_resonance_1p(rn);
    else if (rn.find("1m") != string::npos) act_resonance_1m(rn);
    else {
        cout << "Not know how to act resonance " << rn << endl;
        exit(1);
    }
    //if (rn.find("1m") != string::npos) {
    //    act_resonance_1m(rn);
    //    continue;
    //}

}
void FitParametersOfPhiPP::act_resonance_1m(string rn)
{
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiKK::act_resonance_1m(string rn)
{
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiPP::act_resonance_1p(string rn)
{
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "p_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "p_phi4"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiKK::act_resonance_1p(string rn)
{
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "k_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "k_phi4"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiPP::act_resonance_980(string rn) {
    if (binding_phi) {
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
    } else {
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiKK::act_resonance_980(string rn) {
    if (binding_phi) {
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));

    } else {
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_amplitude980_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_g10_"), gp(rn + "a_g20_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiPP::act_resonance_f0(string rn) {
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiKK::act_resonance_f0(string rn) {
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiPP::act_resonance_f2(string rn) {
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc5"), gp(rn + "a_phi5"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc1"), gp(rn + "p_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc2"), gp(rn + "p_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc3"), gp(rn + "p_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc4"), gp(rn + "p_phi4"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "p_rho_"), gp(rn + "a_frc5"), gp(rn + "p_phi5"), gp(rn + "a_typ_"));
    }
}
void FitParametersOfPhiKK::act_resonance_f2(string rn) {
    if (binding_phi) {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "a_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "a_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "a_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "a_phi4"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc5"), gp(rn + "a_phi5"), gp(rn + "a_typ_"));
    } else {
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn1"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc1"), gp(rn + "k_phi1"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn2"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc2"), gp(rn + "k_phi2"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn3"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc3"), gp(rn + "k_phi3"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn4"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc4"), gp(rn + "k_phi4"), gp(rn + "a_typ_"));
        add_amplitude_to_fit_parameter_list(gp(rn + "a_spn5"), gp(rn + "a_mss_"), gp(rn + "a_wdt_"), gp(rn + "k_rho_"), gp(rn + "a_frc5"), gp(rn + "k_phi5"), gp(rn + "a_typ_"));
    }
}
void AmplitudeMethodWithFitParametersInterface::act_resonances(vector<string> & resonance_list)
{
    for(unsigned i = 0; i < resonance_list.size(); i++) {
        if (!already_active(resonance_list[i]))
            act_resonance(resonance_list[i]);
    }
}
void AmplitudeMethodWithFitParametersInterface::copy_minuit_parameter_to_back()
{
    minuit_parameters_back_ = minuit_parameters_;
}
