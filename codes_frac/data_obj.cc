
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "MultDevice.h"
#include "data_obj.h"
#include "DPFAngular.h"
#include <iostream>

using namespace std;

DataObject::~DataObject()
{
    for(int i = 0; i < number_of_events_; i++) {
        delete[] mcp1[i];
        delete[] mcp2[i];
        delete[] mcp3[i];
        delete[] mcp4[i];
        delete[] mcp5[i];
    }
    delete[] mcp1;
    delete[] mcp2;
    delete[] mcp3;
    delete[] mcp4;
    delete[] mcp5;

}
void DataObject::read_events() {
    initialize_mcp();
    load_mcp_from_dat_file();
    cout << "There are " << number_of_events_ << " events in " << dat_file_name_ << std::endl;
}
void DataObject::read_weight_file() {

}
void DataObject::initialize_mcp()
{
    mcp1=new double*[number_of_events_];
    mcp2=new double*[number_of_events_];
    mcp3=new double*[number_of_events_];
    mcp4=new double*[number_of_events_];
    mcp5=new double*[number_of_events_];
    for(int i=0; i<number_of_events_; i++){
        mcp1[i]=new double[4];
        mcp2[i]=new double[4];
        mcp3[i]=new double[4];
        mcp4[i]=new double[4];
        mcp5[i]=new double[4];
    }

}

void DataObject::load_mcp_from_dat_file()
{
    double fx1,fy1,fz1,ft1,fx2,fy2,fz2,ft2,fx3,fy3,fz3,ft3,fx4,fy4,fz4,ft4,fx5,fy5,fz5,ft5;
    FILE *fp;
    if((fp=fopen(dat_file_name_,"r"))==NULL)
    {printf("can't open input file");
        return;
    }
    //cout << "------->start input mcp(pshp)" << _dp->_phspfile << endl;
    int i=0;
    while(fscanf(fp,"%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n%lf%lf%lf%lf\n",&fx1,&fy1,&fz1,&ft1,&fx2,&fy2,&fz2,&ft2,&fx3,&fy3,&fz3,&ft3,&fx4,&fy4,&fz4,&ft4,&fx5,&fy5,&fz5,&ft5)!=EOF)
    {
        //  //cout<<"haha: "<< __LINE__ << endl;
        mcp1[i][0]=fx1;mcp1[i][1]=fy1;mcp1[i][2]=fz1;mcp1[i][3]=ft1;
        mcp2[i][0]=fx2;mcp2[i][1]=fy2;mcp2[i][2]=fz2;mcp2[i][3]=ft2;
        mcp3[i][0]=fx3;mcp3[i][1]=fy3;mcp3[i][2]=fz3;mcp3[i][3]=ft3;
        mcp4[i][0]=fx4;mcp4[i][1]=fy4;mcp4[i][2]=fz4;mcp4[i][3]=ft4;
        mcp5[i][0]=fx5;mcp5[i][1]=fy5;mcp5[i][2]=fz5;mcp5[i][3]=ft5;
        i++;
    }
    fclose(fp);
    //    Nmc = count_lines() / 5;
    if (i != number_of_events_) {
        std::cout << "Line:" << __LINE__ << "  There is memory allocated error!" << endl;
        exit(1);
    }

}
int DataObject::count_lines() {
    ifstream in(dat_file_name_);
    string line;
    int n = 0;
    while (getline(in, line)) {
        n++;
    }
    //cout << df << " have " << n << " lines." << endl;
    return n;
}
double DataObject::scalar(vector<double> &a1, vector<double> &a2) const
{
    double (*fDel)[4];
    fDel = _dp->fDel;
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            if(i==j){
                if(i<3)
                    fDel[i][j]=-1;
                else
                    fDel[i][j]=1;
            } else {
                fDel[i][j]=0;
            }
        }
    }
    double scal=0;

    for(Int_t i=0;i<4;i++){
        scal+=a1[i]*a2[i]*_dp->fDel[i][i];
    }
    return scal;
}
double DataObject::scalar(double *a1, double *a2) const
{
    double (*fDel)[4];
    fDel = _dp->fDel;
    for(Int_t i=0;i<4;i++){
        for(Int_t j=0;j<4;j++){
            if(i==j){
                if(i<3)
                    fDel[i][j]=-1;
                else
                    fDel[i][j]=1;
            } else {
                fDel[i][j]=0;
            }
        }
    }
    double scal=0;

    for(Int_t i=0;i<4;i++){
        scal+=a1[i]*a2[i]*_dp->fDel[i][i];
    }
    return scal;
}
