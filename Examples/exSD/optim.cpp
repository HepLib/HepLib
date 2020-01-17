#include "optim.hpp"

//
// Ackley function

double ackley_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data) {
    const double x = vals_inp(0);
    const double y = vals_inp(1);
    const double pi = arma::datum::pi;
    
    if(x>1 || x<0 || y>1 || y<0) return 1000;

    double obj_val = -20*std::exp( -0.2*std::sqrt(0.5*(x*x + y*y)) ) - std::exp( 0.5*(std::cos(2*pi*x) + std::cos(2*pi*y)) ) + 22.718282L;

    return -obj_val;
}

int main() {
    // initial values:
    arma::vec x(2);
    x(0) = 0.5;
    x(1) = 0.5;
    
    arma::cout << x.size() << std::endl;

    bool success = optim::de_prmm(x,ackley_fn,nullptr);

    arma::cout << "\nde: solution to Ackley test:\n" << x << arma::endl;
    arma::cout << ackley_fn(x, NULL, NULL) << arma::endl;

    return 0;
}
