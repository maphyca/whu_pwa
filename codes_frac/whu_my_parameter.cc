#include "whu_my_parameter.h"
#include <iostream>

using namespace std;

MyParameter::MyParameter(const MyParameter & the_parameter)
{
    name = the_parameter.name;
    value = the_parameter.value;
    error = the_parameter.error;
    uplimit = the_parameter.uplimit;
    lowlimit = the_parameter.lowlimit;
    type = the_parameter.type;
}
MyParameter MyParameter::operator=(const MyParameter & the_parameter)
{
    MyParameter this_parameter;
    this_parameter.name = the_parameter.name;
    this_parameter.value = the_parameter.value;
    this_parameter.error = the_parameter.error;
    this_parameter.uplimit = the_parameter.uplimit;
    this_parameter.lowlimit = the_parameter.lowlimit;
    this_parameter.type = the_parameter.type;
    return this_parameter;
}
MyParameter::MyParameter(string the_name, double the_value)
{
    name = the_name;
    value = the_value;
    error = 0;
    uplimit = the_value;
    lowlimit = the_value;
    type = 0;
}
MyParameter::MyParameter(string the_name, double the_value, double the_error)
{
    name = the_name;
    value = the_value;
    error = the_error;
    uplimit = +parameter_limit;
    lowlimit = -parameter_limit;
    type = 1;
}
MyParameter::MyParameter(string the_name, double the_value, double the_error, double the_lowlimit, double the_uplimit)
{
    name = the_name;
    value = the_value;
    error = the_error;
    lowlimit = the_lowlimit;
    uplimit = the_uplimit;
    if (lowlimit > -parameter_limit && uplimit < +parameter_limit) {
        type = 2; // 设置了上下界
    } else if (lowlimit > -parameter_limit) {
        type = 3; // 设置了下界
    } else if (uplimit < +parameter_limit) {
        type = 4; // 设置了上界
    } else {
        cout << "Formatting parameter " << name << " failed as limits cannot be comparable with 1e20!!!" << endl;
        exit(1);
    }
}
void MyParameter::info()
{
    cout << "The parameter is named as " << name << endl;
    cout << "Its value is " << value << endl;
    cout << "Its error is " << error << endl;
    cout << "Its lowlimit is " << lowlimit << " and its uplimit is " << uplimit << endl;
    cout << "Its parameter type is " << type << endl;
}
