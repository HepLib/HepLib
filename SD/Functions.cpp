/**
 * @file
 * @brief Functions for SD
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "SD.h"
#include <cmath>

namespace HepLib::SD {

    namespace {
        /*-----------------------------------------------------*/
        // Functions used in GiNaC
        /*-----------------------------------------------------*/
        static void print_VEO(const ex & ex1_in, const ex & ex2_in, const print_context & c) {
            auto oDigits = Digits;
            Digits = 50;
            static ex log10 = log(numeric("10"));
            if(is_zero(ex1_in) && is_zero(ex2_in)) {
                c.s << "0(0)";
            } else if(is_zero(ex1_in)) {
                numeric ex2 = ex_to<numeric>(Rationalize(ex2_in));
                int n2 = floor(ex_to<numeric>(log(ex2)/log10).to_double());
                double d2 = round(ex_to<numeric>(ex2/power(10,n2)).to_double());
                if(d2==10) {d2=1; n2++;}
                c.s << "0(" << d2 << ")";
                if(n2!=0) c.s << "E" << n2;
            } else if(is_zero(ex2_in)) {
                numeric ex1 = ex_to<numeric>(Rationalize(ex1_in));
                bool neg = bool(ex1<0);
                if(neg) ex1 = 0-ex1;
                int n1 = floor(ex_to<numeric>(log(ex1)/log10).to_double());
                int n12 = VEO_Digits;
                ostringstream oss;
                oss << evalf(ex1/power(10,n1));
                string sd1 = oss.str();
                if(sd1.length()<2+n12) {
                    oss.clear();
                    oss.str("");
                    oss << evalf(ex1/power(10,n1)+power(10,-3-n12));
                    sd1 = oss.str();
                } else {
                    int ci = stoi(sd1.substr(2+n12,1));
                    if(ci>4) {
                        oss.clear();
                        oss.str("");
                        oss << evalf(ex1/power(10,n1)+power(10,-n12));
                        sd1 = oss.str();
                    }
                }
                sd1 = sd1.substr(0, 2+n12);
                
                if(neg) c.s << "(-";
                c.s << sd1 << "(0)";
                if(n1!=0) c.s << "E" << n1;
                if(neg) c.s << ")";
            } else {
                numeric ex1 = ex_to<numeric>(Rationalize(ex1_in));
                numeric ex2 = ex_to<numeric>(Rationalize(ex2_in));
                bool neg = bool(ex1<0);
                if(neg) ex1 = 0-ex1;
                
                int n1 = floor(ex_to<numeric>(log(ex1)/log10).to_double());
                int n2 = floor(ex_to<numeric>(log(ex2)/log10).to_double());
                int d2 = round(ex_to<numeric>(ex2/power(10,n2)).to_double());
                if(d2==10) {d2=1; n2++;}
                
                if(n2>n1) {
                    if(n2==n1+1 && (ex1>5*power(10,n1)))  c.s << "1";
                    else c.s << "0";
                    c.s << "(" << d2 << ")";
                    if(n2!=0) c.s << "E" << n2;
                } else {
                    int n12 = n1-n2;
                    if(n12>VEO_Digits) {
                        n12 = VEO_Digits;
                        d2 = 0;
                    }
                    ostringstream oss;
                    oss << evalf(ex1/power(10,n1));
                    string sd1 = oss.str();
                    if(sd1.length()<3+n12) {
                        oss.clear();
                        oss.str("");
                        oss << evalf(ex1/power(10,n1)+power(10,-4-n12));
                        sd1 = oss.str();
                    } else {
                        int ci = stoi(sd1.substr(2+n12,1));
                        if(ci>4) {
                            oss.clear();
                            oss.str("");
                            oss << evalf(ex1/power(10,n1)+power(10,-n12));
                            sd1 = oss.str();
                        }
                    }
                    sd1 = sd1.substr(0, 2+n12);
                    if(sd1[sd1.length()-1]=='.') sd1 = sd1.substr(0, sd1.length()-1);
                    if(neg) c.s << "(-";
                    c.s << sd1 << "(" << d2 << ")";
                    if(n1!=0) c.s << "E" << n1;
                    if(neg) c.s << ")";
                }
            }
            Digits = oDigits;
            return;
        }
                
        static ex NoDiff_1P(const ex & x, unsigned diff_param) {return 0;}
        static ex NoDiff_2P(const ex & x, const ex & y, unsigned diff_param) {return 0;}
        static ex VE_Conjugate(const ex & x, const ex & y) { return VE(x,y).hold(); }
        static ex Diff_ID(const ex & x, unsigned diff_param) {return 1;}
        
    }

    REGISTER_FUNCTION(CV, do_not_evalf_params()) // for use symbol
    REGISTER_FUNCTION(fabs, dummy())
    REGISTER_FUNCTION(PL, do_not_evalf_params())
    REGISTER_FUNCTION(FTX, do_not_evalf_params().derivative_func(NoDiff_2P))
    REGISTER_FUNCTION(CT, do_not_evalf_params().derivative_func(Diff_ID))
    REGISTER_FUNCTION(VE, conjugate_func(VE_Conjugate))
    REGISTER_FUNCTION(VEO, print_func<print_dflt>(print_VEO))
    REGISTER_FUNCTION(epsID, do_not_evalf_params().derivative_func(NoDiff_1P))

}
