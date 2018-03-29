#include "../codes_frac/fitproxy.h"
#include "../codes_frac/PWA_CTRL.H"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TString.h"

void pwa_ana(PWA_CTRL & pwa_ctrl) {
    fitproxy *A = new fitproxy();
    A->init_all_files_name(pwa_ctrl);
    A->init_input_argset();
    A->setup_resonances();
    // A->add_res0_list("f01750", 1.7, 2.2 0.15, 0.7) here!!!
    A->createlist_allparas();

    for(vector<string>::iterator it = pwa_ctrl.paraConstList.begin(); it != pwa_ctrl.paraConstList.end(); it++) {
        A->gp(*it)->setConstant();
    }
    A->print_all_paras();

    A->read_data();
    A->init_pdf(pwa_ctrl);
    cout << "YYYYY" <<endl;

    for(vector<string>::iterator it = pwa_ctrl.actResList.begin(); it != pwa_ctrl.actResList.end(); it++) {
        cout << "test = " << *it << endl;
        if (*it == "f00980") {
            A->act_res980("f00980");
            continue;
        }
        if (it->find("1m") != string::npos) {
            A->act_res1m(*it);
            continue;
        }
        if (it->find("1p") != string::npos) {
            A->act_res1p(*it);
            continue;
        }
        if (it->find("f0") != string::npos) {
            A->act_res0(*it);
            continue;
        }
        if (it->find("f2") != string::npos) {
            A->act_res2(*it);
            continue;
        }
    }
    cout << "Active Res PP" << endl;
    if (pwa_ctrl.actResListPP.size() > 0) {
        for(vector<string>::iterator it = pwa_ctrl.actResListPP.begin(); it != pwa_ctrl.actResListPP.end(); it++) {
            cout << "Res PP " << *it << endl;
            if (*it == "f00980") {
                A->act_res980("f00980", true, false);
                continue;
            }
            if (it->find("1m") != string::npos) {
                A->act_res1m(*it, true, false);
                continue;
            }
            if (it->find("1p") != string::npos) {
                A->act_res1p(*it, true, false);
                continue;
            }
            if (it->find("f0") != string::npos) {
                A->act_res0(*it, true, false);
                continue;
            }
            if (it->find("f2") != string::npos) {
                A->act_res2(*it, true, false);
                continue;
            }
        }
    }
    cout << "Active Res KK " << endl;
    if (pwa_ctrl.actResListKK.size() > 0) {
        for(vector<string>::iterator it = pwa_ctrl.actResListKK.begin(); it != pwa_ctrl.actResListKK.end(); it++) {
            cout << "Res KK " << *it << endl;
            if (*it == "f00980") {
                A->act_res980("f00980", false, true);
                continue;
            }
            if (it->find("1m") != string::npos) {
                A->act_res1m(*it, false, true);
                continue;
            }
            if (it->find("1p") != string::npos) {
                A->act_res1p(*it, false, true);
                continue;
            }
            if (it->find("f0") != string::npos) {
                A->act_res0(*it, false, true);
                continue;
            }
            if (it->find("f2") != string::npos) {
                A->act_res2(*it, false, true);
                continue;
            }
        }
    }

    if (pwa_ctrl.reloadFitParaFile != "NONE") {
        cout << "Start reload fit results from " << pwa_ctrl.reloadFitParaFile << endl;
        A->reload_paras(pwa_ctrl.reloadFitParaFile);
    }
    cout << "finished reload_paras!!!!" << endl;
    A->print_fit_paras();
    cout << "print fit_paras!!!" << endl;

    //A->init_pdf(pwa_ctrl);
    A->FIT();
    A->print_fit_paras();
    A->store_fit_paras(pwa_ctrl.storeFitParaFile);
    cout << pwa_ctrl.storeFitParaFile << " is created!!!" << endl;

//    A->pdfphipp->projectpipi(pwa_ctrl);
//    A->pdfphikk->projectkk(pwa_ctrl);
    delete A;

}

void pwa_project(fitproxy *A, PWA_CTRL &pwa_ctrl, string res) {
        A->print_fit_paras();
        cout << "Begin " << res << " ^^^^^^^^" << endl;
        cout << "pwa_ctrl.reloadFitParaFile = " << pwa_ctrl.reloadFitParaFile << endl;
        A->reload_paras(pwa_ctrl.reloadFitParaFile);
        cout << "reload line: " << __LINE__ << endl;
        A->print_fit_paras();
        cout << "line: " << __LINE__ << endl;
        cout << "Start project " << res << "&&&&&&&&&" << endl;
        A->pdfphipp->showNumAll();
        A->pdfphikk->showNumAll();
        pwa_ctrl.actResList.resize(0);
        pwa_ctrl.actResList.push_back(res);
        A->pdfphikk->setup_iter_vec();
        A->pdfphipp->setup_iter_vec();
        A->pdfphikk->paras_getval();
        A->pdfphipp->paras_getval();
        A->pdfphikk->projectkk(pwa_ctrl);
        A->pdfphipp->projectpipi(pwa_ctrl);
}

void pwa_fraction(PWA_CTRL & pwa_ctrl) {
    cout << "Begin pwa_fraction" << endl;
    fitproxy *A = new fitproxy();
    A->init_all_files_name(pwa_ctrl);
    A->init_input_argset();
    A->setup_resonances();
    // A->add_res0_list("f01750", 1.7, 2.2 0.15, 0.7) here!!!
    A->createlist_allparas();


    //string paraConstStr;
    //readConfigFile(config_file, "const_paras", paraConstStr);
    ////string paraConstStr = "f01000a_mss_,f01000a_wdt_,f21000a_mss_,f21000a_wdt_";
    //vector<string> paraConstList;
    //string_to_vector(paraConstStr, paraConstList);
    for(vector<string>::iterator it = pwa_ctrl.paraConstList.begin(); it != pwa_ctrl.paraConstList.end(); it++) {
        A->gp(*it)->setConstant();
    }
    //A->gp("f01000a_mss_")->setConstant();
    //A->gp("f01000a_wdt_")->setConstant();
    //A->gp("f21000a_mss_")->setConstant();
    //A->gp("f21000a_wdt_")->setConstant();
    //A->print_all_paras();

    A->print_all_paras();
    A->read_data();
    A->init_pdf(pwa_ctrl);


    vector<string> rname(pwa_ctrl.actResList.begin(), pwa_ctrl.actResList.end());
    for(vector<string>::iterator it = rname.begin(); it != rname.end(); it++) {
        if (*it == "f00980") {
            A->act_res980("f00980");
            continue;
        }
        if (it->find("1m") != string::npos) {
            A->act_res1m(*it);
            continue;
        }
        if (it->find("1p") != string::npos) {
            A->act_res1p(*it);
            continue;
        }
        if (it->find("f0") != string::npos) {
            A->act_res0(*it);
            continue;
        }
        if (it->find("f2") != string::npos) {
            A->act_res2(*it);
        }
    }
    cout << "Active Res PP" << endl;
    if (pwa_ctrl.actResListPP.size() > 0) {
        for(vector<string>::iterator it = pwa_ctrl.actResListPP.begin(); it != pwa_ctrl.actResListPP.end(); it++) {
            cout << "Res PP " << *it << endl;
            if (*it == "f00980") {
                A->act_res980("f00980", true, false);
                continue;
            }
            if (it->find("1m") != string::npos) {
                A->act_res1m(*it, true, false);
                continue;
            }
            if (it->find("1p") != string::npos) {
                A->act_res1p(*it, true, false);
                continue;
            }
            if (it->find("f0") != string::npos) {
                A->act_res0(*it, true, false);
                continue;
            }
            if (it->find("f2") != string::npos) {
                A->act_res2(*it, true, false);
                continue;
            }
        }
    }
    cout << "Active Res KK " << endl;
    if (pwa_ctrl.actResListKK.size() > 0) {
        for(vector<string>::iterator it = pwa_ctrl.actResListKK.begin(); it != pwa_ctrl.actResListKK.end(); it++) {
            cout << "Res KK " << *it << endl;
            if (*it == "f00980") {
                A->act_res980("f00980", false, true);
                continue;
            }
            if (it->find("1m") != string::npos) {
                A->act_res1m(*it, false, true);
                continue;
            }
            if (it->find("1p") != string::npos) {
                A->act_res1p(*it, false, true);
                continue;
            }
            if (it->find("f0") != string::npos) {
                A->act_res0(*it, false, true);
                continue;
            }
            if (it->find("f2") != string::npos) {
                A->act_res2(*it, false, true);
                continue;
            }
        }
    }


    A->reload_paras(pwa_ctrl.reloadFitParaFile);


    A->pdfphikk->projectkk(pwa_ctrl);
    cout << "finish projectkk23" << endl;
    A->print_fit_paras();
    exit(1);
    A->pdfphipp->projectpipi(pwa_ctrl);
    cout << "finish projectpipi23" << endl;
    cout << "finish projectpipi and projectkk" << endl;


    for(vector<string>::iterator it = rname.begin(); it != rname.end(); it++) {
        cout << "line: " << __LINE__ << endl;
        A->pdfphipp->removeAll();
        A->pdfphikk->removeAll();
        cout << "line: " << __LINE__ << endl;
        if (*it == "f00980") {
        cout << "line: " << __LINE__ << endl;
            A->act_res980("f00980");
        cout << "line: " << __LINE__ << endl;
            pwa_project(A, pwa_ctrl, *it);
        cout << "line: " << __LINE__ << endl;
            continue;
        }
        if (it->find("1m") != string::npos) {
            A->act_res1m(*it);
            pwa_project(A, pwa_ctrl, *it);
            continue;
        }
        if (it->find("1p") != string::npos) {
            A->act_res1p(*it);
            pwa_project(A, pwa_ctrl, *it);
            continue;
        }
        if (it->find("f0") != string::npos) {
        cout << "line: " << __LINE__ << endl;
            A->act_res0(*it);
        cout << "line: " << __LINE__ << endl;
    A->print_fit_paras();
    exit(1);
        cout << "line: " << __LINE__ << endl;
            pwa_project(A, pwa_ctrl, *it);
        cout << "line: " << __LINE__ << endl;
            continue;
        }
        if (it->find("f2") != string::npos) {
            A->act_res2(*it);
            pwa_project(A, pwa_ctrl, *it);
            continue;
        }
    }
    if (pwa_ctrl.actResListPP.size() > 0) {
        for(vector<string>::iterator it = pwa_ctrl.actResListPP.begin(); it != pwa_ctrl.actResListPP.end(); it++) {
            A->pdfphipp->removeAll();
            A->pdfphikk->removeAll();
        cout << "line: " << __LINE__ << endl;
            cout << "Res PP " << *it << endl;
            if (*it == "f00980") {
                A->act_res980("f00980", true, false);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("1m") != string::npos) {
                A->act_res1m(*it, true, false);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("1p") != string::npos) {
                A->act_res1p(*it, true, false);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("f0") != string::npos) {
                A->act_res0(*it, true, false);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("f2") != string::npos) {
                A->act_res2(*it, true, false);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
        }
    }
    if (pwa_ctrl.actResListKK.size() > 0) {
        for(vector<string>::iterator it = pwa_ctrl.actResListKK.begin(); it != pwa_ctrl.actResListKK.end(); it++) {
            A->pdfphipp->removeAll();
            A->pdfphikk->removeAll();
        cout << "line: " << __LINE__ << endl;
            cout << "Res KK " << *it << endl;
            if (*it == "f00980") {
                A->act_res980("f00980", false, true);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("1m") != string::npos) {
                A->act_res1m(*it, false, true);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("1p") != string::npos) {
                A->act_res1p(*it, false, true);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("f0") != string::npos) {
                A->act_res0(*it, false, true);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
            if (it->find("f2") != string::npos) {
                A->act_res2(*it, false, true);
                pwa_project(A, pwa_ctrl, *it);
                continue;
            }
        }
    }
    cout << "finished fracton_pwa" << endl;

    //A->act_res980("f00980");
    //A->act_res0("f01000");

    //A->act_res1m("1m1800");
    //A->act_res1p("1p1800");

    //A->act_res0("f01500");
    //    A->act_res2("f21270");
    //    A->act_res2("f21525");
    // A->act_res0("f01750");  here!!!

    //A->FIT();
    //A->print_fit_paras();
    //A->store_fit_paras(pwa_ctrl.storeFitParaFile);
    //cout << pwa_ctrl.storeFitParaFile << " is created!!!" << endl;

    delete A;

}

int main(int argc, char *argv[])
{
    stringstream ss;
    ss << argv[1];
    string config_file = ss.str(); // config_file 配置文件
    PWA_CTRL pwa_ctrl; // PWA_CTRL 读取配置文件以后参数存放的结构体
    cout << "pwa_config file is " << config_file << endl;
    refresh_pwa_ctrl(config_file, pwa_ctrl);

    cout << pwa_ctrl.reloadFitParaFile << endl;
    if (argc == 2) {
        pwa_ana(pwa_ctrl);
    } else {
        stringstream tt;
        tt << argv[2];
        string fitResult = tt.str();
        pwa_ctrl.reloadFitParaFile = fitResult;

        cout << "pwa_ctrl.reloadFitParaFile = " << pwa_ctrl.reloadFitParaFile << endl;
        cout << "TESTTEST" << endl;
        pwa_fraction(pwa_ctrl);
    //    PWA_CTRL fr_pwa_ctrl;
    //    for(vector<string>::iterator it = pwa_ctrl.actResList.begin(); it != pwa_ctrl.actResList.end(); it++) {
    //        fr_pwa_ctrl = pwa_ctrl;
    //        fr_pwa_ctrl.actResList.resize(0);
    //        fr_pwa_ctrl.actResList.push_back(*it);
    //        pwa_fraction(fr_pwa_ctrl);
    //    }
    }
}


