/**
 * @file
 * @brief Functions for SD
 */
 
#include "SD.h"
#include <cmath>

namespace HepLib::SD {

    namespace {
        /*-----------------------------------------------------*/
        // Functions used in GiNaC
        /*-----------------------------------------------------*/
        static void print_VEO(const ex & ex1_in, const ex & ex2_in, const print_context & c) {
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
                oss << NN(ex1/power(10,n1));
                string sd1 = oss.str();
                if(sd1.length()<2+n12) {
                    oss.clear();
                    oss.str("");
                    oss << NN(ex1/power(10,n1)+power(10,-4-n12));
                    sd1 = oss.str();
                    if(sd1.length()<2+n12) throw Error("VEO: sd1.length()<2+n12");
                } else {
                    int ci = stoi(sd1.substr(2+n12,1));
                    if(ci>4) {
                        oss.clear();
                        oss.str("");
                        oss << NN(ex1/power(10,n1)+power(10,-n12));
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
                    oss << NN(ex1/power(10,n1));
                    string sd1 = oss.str();
                    if(sd1.length()<3+n12) {
                        oss.clear();
                        oss.str("");
                        oss << NN(ex1/power(10,n1)+power(10,-5-n12));
                        sd1 = oss.str();
                        if(sd1.length()<3+n12) throw Error("VEO2: sd1.length()<3+n12");
                    } else {
                        int ci = stoi(sd1.substr(2+n12,1));
                        if(ci>4) {
                            oss.clear();
                            oss.str("");
                            oss << NN(ex1/power(10,n1)+power(10,-n12));
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
            return;
        }
        
        static void print_VEO2(const ex & ex1_in, const ex & ex2_in, const print_context & c) {
            static ex log10 = log(numeric("10"));
            if(is_zero(ex1_in) && is_zero(ex2_in)) {
                c.s << "00(00)";
            } else if(is_zero(ex1_in)) {
                numeric ex2 = ex_to<numeric>(Rationalize(ex2_in));
                int n2 = floor(ex_to<numeric>(log(ex2)/log10).to_double())-1;
                double d2 = round(ex_to<numeric>(ex2/power(10,n2)).to_double());
                if(d2==100) {d2=10; n2++;}
                c.s << "00(" << d2 << ")";
                if(n2!=0) c.s << "E" << n2;
            } else if(is_zero(ex2_in)) {
                numeric ex1 = ex_to<numeric>(Rationalize(ex1_in));
                bool neg = bool(ex1<0);
                if(neg) ex1 = 0-ex1;
                int n1 = floor(ex_to<numeric>(log(ex1)/log10).to_double());
                int n12 = VEO_Digits;
                ostringstream oss;
                oss << NN(ex1/power(10,n1));
                string sd1 = oss.str();
                if(sd1.length()<2+n12) {
                    oss.clear();
                    oss.str("");
                    oss << NN(ex1/power(10,n1)+power(10,-4-n12));
                    sd1 = oss.str();
                    if(sd1.length()<2+n12) throw Error("VEO2: sd1.length()<2+n12");
                } else {
                    int ci = stoi(sd1.substr(2+n12,1));
                    if(ci>4) {
                        oss.clear();
                        oss.str("");
                        oss << NN(ex1/power(10,n1)+power(10,-n12));
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
                int d2 = round(ex_to<numeric>(ex2/power(10,n2-1)).to_double());
                if(d2==100) {d2=10; n2++;}
                
                if(n2>=n1+2) {
                    if(n2==n1+2 && (ex1>5*power(10,n1-1)))  c.s << "01";
                    else c.s << "00";
                    c.s << "(" << d2 << ")";
                    if(n2!=0) c.s << "E" << n2;
                } else if(n2>n1) {
                    int d1 = round(ex_to<numeric>(ex1/power(10,n2-1)).to_double());
                    if(d1<10)  c.s << "0";
                    c.s << d1;
                    c.s << "(" << d2 << ")";
                    if(n2!=0) c.s << "E" << n2;
                } else {
                    int n12 = n1-n2;
                    if(n12>VEO_Digits) {
                        n12 = VEO_Digits;
                        d2 = 0;
                    }
                    ostringstream oss;
                    oss << NN(ex1/power(10,n1));
                    string sd1 = oss.str();
                    if(sd1.length()<4+n12) {
                        oss.clear();
                        oss.str("");
                        oss << NN(ex1/power(10,n1)+power(10,-6-n12));
                        sd1 = oss.str();
                        if(sd1.length()<4+n12) throw Error("VEO2: sd1.length()<4+n12");
                    } else {
                        int ci = stoi(sd1.substr(3+n12,1));
                        if(ci>4) {
                            oss.clear();
                            oss.str("");
                            oss << NN(ex1/power(10,n1)+power(10,-1-n12));
                            sd1 = oss.str();
                        }
                    }
                    sd1 = sd1.substr(0, 3+n12);
                    if(sd1[sd1.length()-1]=='.') sd1 = sd1.substr(0, sd1.length()-1);
                    if(neg) c.s << "(-";
                    if(sd1.length()==3 && sd1[1]=='.') {
                        c.s << "0." << sd1[0] << sd1[2] << "(" << d2 << ")";
                        if(n1!=-1) c.s << "E" << n1+1;
                    } else {
                        c.s << sd1 << "(" << d2 << ")";
                        if(n1!=0) c.s << "E" << n1;
                    }
                    if(neg) c.s << ")";
                }
            }
            return;
        }
        
    }
    
    namespace {
        static ex conjVE(const ex & x, const ex & y) { return VE(x,y).hold(); }
        static ex zp1D(const ex & x, unsigned diff_param) {return 0;}
        static ex zd1D(const ex & x, const symbol & s) {return 0;} 
        static ex zp2D(const ex & x, const ex & y, unsigned diff_param) {return 0;}
        static ex zd2D(const ex & x, const ex & y, const symbol & s) {return 0;}   
        static ex dCT(const ex & x, const symbol & s) {return x.diff(s);}
        static ex pCT(const ex & x, unsigned diff_param) {return 1;}     
    }
    
    REGISTER_FUNCTION(CV, do_not_evalf_params()) // for use symbol
    REGISTER_FUNCTION(fabs, dummy())
    REGISTER_FUNCTION(PL, do_not_evalf_params().expl_derivative_func(zd1D).derivative_func(zp1D))
    REGISTER_FUNCTION(FTX, do_not_evalf_params().expl_derivative_func(zd2D).derivative_func(zp2D))
    REGISTER_FUNCTION(CT, do_not_evalf_params().expl_derivative_func(dCT).derivative_func(pCT))
    REGISTER_FUNCTION(WRA, do_not_evalf_params().expl_derivative_func(zd1D).derivative_func(zp1D))
    REGISTER_FUNCTION(VE, conjugate_func(conjVE))
    REGISTER_FUNCTION(VEO, print_func<print_dflt>(print_VEO))
    REGISTER_FUNCTION(VEO2, print_func<print_dflt>(print_VEO2))
}
