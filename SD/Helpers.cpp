#include "SD.h"

namespace HepLib {

/*-----------------------------------------------------*/
// IntegratorBase::inDQMP
/*-----------------------------------------------------*/
int IntegratorBase::inDQMP(unsigned xdim, qREAL const *x) {
    if(xdim<=MPXDim) return 3;
    qREAL xmin = 100;
    for(int i=0; i<xdim; i++) {
        if(x[i] < xmin) xmin = x[i];
    }
    if(xmin <= MPXLimit) return 3;
    if(xdim <= QXDim|| xmin < QXLimit) return 2;
    return 1;
}

/*-----------------------------------------------------*/
// Global Functions
/*-----------------------------------------------------*/

vector<ex> get_xy_from(ex pol) {
    exset xyset;
    bool ok = pol.find(x(wild()), xyset);
    if(!ok) {
        ok = pol.find(y(wild()), xyset);
        if(!ok) {
            vector<ex> xys(0);
            return xys;
        }
    }
    vector<ex> xys(xyset.size());
    copy(xyset.begin(), xyset.end(), xys.begin());
    sort(xys.begin(), xys.end(), [&](const auto &a, const auto &b){
        return ex_to<numeric>(normal((b-a)).subs(lst{ x(wild())==wild(), y(wild())==wild() })).is_positive();
    });
    return xys;
}

vector<ex> get_x_from(ex pol) {
    exset xset;
    bool ok = pol.find(x(wild()), xset);
    if(!ok) {
        vector<ex> xs(0);
        return xs;
    }
    vector<ex> xs(xset.size());
    copy(xset.begin(), xset.end(), xs.begin());
    sort(xs.begin(), xs.end(), [&](const auto &a, const auto &b){
        return ex_to<numeric>(normal((b-a)).subs(lst{ x(wild())==wild() })).is_positive();
    });
    return xs;
}

vector<ex> get_y_from(ex pol) {
    exset yset;
    bool ok = pol.find(y(wild()), yset);
    if(!ok) {
        vector<ex> ys(0);
        return ys;
    }
    vector<ex> ys(yset.size());
    copy(yset.begin(), yset.end(), ys.begin());
    sort(ys.begin(), ys.end(), [&](const auto &a, const auto &b){
        return ex_to<numeric>(normal((b-a)).subs(lst{ y(wild())==wild() })).is_positive();
    });
    return ys;
}

vector<ex> get_pl_from(ex pol) {
    exset plset;
    bool ok = pol.find(PL(wild()), plset);
    if(!ok) {
        vector<ex> pls(0);
        return pls;
    }
    vector<ex> pls(plset.size());
    copy(plset.begin(), plset.end(), pls.begin());
    sort(pls.begin(), pls.end(), [&](const auto &a, const auto &b){
        return ex_to<numeric>(normal((b-a)).subs(lst{ PL(wild())==wild() })).is_positive();
    });
    return pls;
}

/*-----------------------------------------------------*/
// VE
/*-----------------------------------------------------*/
ex VESimplify(ex expr, int epN, int epsN) {
    Digits = 50;
    auto ep = SD::ep;
    auto eps = SD::eps;
    auto expr1 = expr.evalf();
    if(expr1.has(eps)) expr1 = mma_series(expr1, eps, epsN);
    expr1 = mma_series(expr1, ep, epN);
    expr1 = expr1.expand();
    ex ret = 0;
    for(int si=expr1.ldegree(eps); si<=epsN; si++) {
    for(int i=expr1.ldegree(ep); i<=epN; i++) {
    
        auto ccRes = expr1.coeff(eps,si).coeff(ep,i).expand();
        lst ccResList;
        if(is_a<add>(ccRes)) {
            for(auto item : ccRes) ccResList.append(item);
        } else {
            ccResList.append(ccRes);
        }
        
        ccRes = 0;
        for(auto item : ccResList) {
            if(is_a<mul>(item)) {
                ex cc = 1, cf = 1;
                for(auto ii : item) {
                    if(is_a<numeric>(ii) || ii.match(VE(wild(1), wild(2)))) {
                        cc *= ii;
                    } else {
                        cf *= ii;
                    }
                }
                ccRes += cc * VF(cf);
            } else if(is_a<numeric>(item) || item.match(VE(wild(1), wild(2)))) {
                ccRes += item*VF(1);
            } else {
                ccRes += VF(item);
            }
        }
    
        
        exset vfs;
        find(ccRes, VF(wild()),vfs);
        
        for(auto vf : vfs) {
            auto tmpIR = ccRes.coeff(vf).expand();
            exset ves;
            tmpIR.find(VE(wild(1), wild(2)), ves);
            auto ntmp = tmpIR.subs(lst{VE(wild(1), wild(2))==0});
            assert(is_a<numeric>(ntmp));
            if(abs(ntmp.evalf())>numeric("1E300")) return SD::NaN;
            ex vIR = ntmp;
            ex eI2 = 0, eR2 = 0;
            for(auto ve : ves) {
                auto co = tmpIR.coeff(ve);
                vIR += co * ve.op(0);
                assert(is_a<numeric>(co));
                numeric nco = ex_to<numeric>(co);
                ex ee = ve.op(1) * ve.op(1);
                eR2 += nco.real_part() * nco.real_part() * ee;
                eI2 += nco.imag_part() * nco.imag_part() * ee;
            }
            ret += VE(ex_to<numeric>(vIR).real_part(), sqrt(eR2)) * pow(eps,si) * pow(ep,i) * vf;
            ret += VE(ex_to<numeric>(vIR).imag_part(), sqrt(eI2)) * pow(eps,si) * pow(ep,i) * vf * I;
        }
    }}
    
    ret = ret.subs(lst{VF(wild())==wild()});
    ret = ret.subs(lst{VE(0,0)==0});
    return ret.collect(lst{eps,ep}, true);
}

ex VEResult(ex expr) {
    return expr.subs(VE(wild(1),wild(2))==VEO(wild(1),wild(2)));
}

/*-----------------------------------------------------*/
// Functions used in GiNaC
/*-----------------------------------------------------*/
static ex NoDiff_1P(const ex & x, unsigned diff_param) {return 0;}
static ex NoDiff_2P(const ex & x, const ex & y, unsigned diff_param) {return 0;}
static ex VE_Conjugate(const ex & x, const ex & y) { return VE(x,y).hold(); }
static ex Diff_ID(const ex & x, unsigned diff_param) {return 1;}

static void print_VEO(const ex & ex1_in, const ex & ex2_in, const print_context & c) {
    ex ex1 = ex1_in, ex2 = ex2_in;
    if(abs(ex1) < numeric("1E-30")) ex1 = 0;
    if(abs(ex2) < numeric("1E-30")) ex2 = 0;
    if(ex1==0 || ex2==0) {
        char bf1[128], bf2[128];
        quadmath_snprintf(bf1, sizeof bf1, "%.10Qg", CppFormat::ex2q(ex1_in));
        quadmath_snprintf(bf2, sizeof bf2, "%.10Qg", CppFormat::ex2q(ex2_in));
        c.s << bf1 << " +- " << bf2;
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
        c.s << RED << "[-NaN-]" << RESET;
    }
}

REGISTER_FUNCTION(fabs, dummy())
REGISTER_FUNCTION(x, dummy())
REGISTER_FUNCTION(y, dummy())
REGISTER_FUNCTION(z, dummy())
REGISTER_FUNCTION(PL, dummy())
REGISTER_FUNCTION(FTX, derivative_func(NoDiff_2P))
REGISTER_FUNCTION(CT, derivative_func(Diff_ID))
REGISTER_FUNCTION(VE, conjugate_func(VE_Conjugate))
REGISTER_FUNCTION(VEO, print_func<print_dflt>(print_VEO))

}
