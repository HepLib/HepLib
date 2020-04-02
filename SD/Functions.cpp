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
                numeric ex2 = ex_to<numeric>(ex2_in);
                int n2 = floor(ex_to<numeric>(log(ex2)/log10).to_double());
                double d2 = round(ex_to<numeric>(ex2/power(10,n2)).to_double());
                if(d2==10) {d2=1; n2++;}
                c.s << "0(" << d2 << ")";
                if(n2!=0) c.s << "E" << n2;
            } else if(is_zero(ex2_in)) {
                numeric ex1 = ex_to<numeric>(ex1_in);
                bool neg = bool(ex1<0);
                if(neg) ex1 = 0-ex1;
                int n1 = floor(ex_to<numeric>(log(ex1)/log10).to_double());
                int n12 = 20;
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
                numeric ex1 = ex_to<numeric>(ex1_in);
                numeric ex2 = ex_to<numeric>(ex2_in);
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
                    if(n12>20) {
                        n12 = 20;
                        d2 = 0;
                    }
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
                    c.s << sd1 << "(" << d2 << ")";
                    if(n1!=0) c.s << "E" << n1;
                    if(neg) c.s << ")";
                }
            }
            Digits = oDigits;
            return;
        }
        
        static void print_VEO_Old(const ex & ex1_in, const ex & ex2_in, const print_context & c) {
            ex ex1 = ex1_in, ex2 = ex2_in;
            if(abs(ex1) < numeric("1E-30")) ex1 = 0;
            if(abs(ex2) < numeric("1E-30")) ex2 = 0;
            if(ex1==0 || ex2==0) {
                char bf1[128], bf2[128];
                quadmath_snprintf(bf1, sizeof bf1, "%.10QG", CppFormat::ex2q(ex1_in));
                quadmath_snprintf(bf2, sizeof bf2, "%.10QG", CppFormat::ex2q(ex2_in));
                c.s << "(" << bf1 << " +- " << bf2 << ")";
                return;
            }
            int digits = 30;
            if(!ex2.is_zero()) {
                auto ratio = ex_to<numeric>(abs(ex1/ex2));
                digits = floorq(logq(ratio.to_double())/logq(10.Q)) + 2;
                digits = digits > 1 ? digits : 1;
                digits = digits > 30 ? 30 : digits;
            }
            auto oDigits = Digits;
            ostringstream oss;
            try {
                Digits = digits;
                oss << "(" << ex1.evalf();
                Digits = 2;
                oss << " +- " << ex2.evalf() << ")";
                Digits = oDigits;
                c.s << oss.str();
            } catch(...) {
                Digits = oDigits;
                try {
                    char bf1[128], bf2[128];
                    quadmath_snprintf(bf1, sizeof bf1, "%.10QG", CppFormat::ex2q(ex1_in));
                    quadmath_snprintf(bf2, sizeof bf2, "%.10QG", CppFormat::ex2q(ex2_in));
                    c.s << "(" << bf1 << " +- " << bf2 << ")";
                } catch(...) {
                    c.s << Color_Error << "[-NaN-]" << RESET;
                }
            }
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
