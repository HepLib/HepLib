#include "SD.h"

using namespace HepLib;

#include <iomanip>
#include <iostream>
#include <vector>

#include <nlopt.hpp>



double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    return sqrt(x[1]+x[0]);
}


int main() {

nlopt::opt opt(nlopt::LN_COBYLA, 2);
std::vector<double> lb(2);
lb[0] = 0; lb[1] = 0;
opt.set_lower_bounds(lb);
std::vector<double> ub(2);
ub[0] = 10; ub[1] = 10;
opt.set_upper_bounds(ub);

opt.set_max_objective(myvfunc, NULL);


opt.set_xtol_rel(1e-4);
std::vector<double> x(2);
x[0] = 1.234;
x[1] = 5.678;
double minf;

try{
    nlopt::result result = opt.optimize(x, minf);
    std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
        << std::setprecision(10) << minf << std::endl;
}
catch(std::exception &e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
}

}
