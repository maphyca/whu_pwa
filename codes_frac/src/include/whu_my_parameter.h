#ifndef MY_PARAMETER
#define MY_PARAMETER
#include <string>

const double parameter_limit = 1e30;

class MyParameter {
    public:
        MyParameter(const MyParameter &);
        MyParameter operator=(const MyParameter &);
        ~MyParameter() {};
        MyParameter() {};
        MyParameter(std::string, double);
        MyParameter(std::string, double, double);
        MyParameter(std::string, double, double, double, double);

        std::string get_name() const { return name; };
        double get_value() const { return value; };
        double get_error() const { return error; };
        double get_uplimit() const { return uplimit; };
        double get_lowlimit() const { return lowlimit; };
        int get_type() const { return type; };

        void set_name(std::string the_name) { name = the_name; };
        void set_value(double the_value) { value = the_value; };
        void set_error(double the_error) { error = the_error; };
        void set_uplimit(double the_uplimit) { uplimit = the_uplimit; };
        void set_lowlimit(double the_lowlimit) { lowlimit = the_lowlimit; };
        void set_type(int the_type) { type = the_type; };
        void info();

    private:
        std::string name;
        double value;
        double error;
        double uplimit, lowlimit;
        int type;
};

#endif
