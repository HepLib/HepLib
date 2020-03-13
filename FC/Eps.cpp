#include "FC.h"

namespace HepLib::FC {

    //-----------------------------------------------------------
    // Pair Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Eps, basic,
        print_func<print_dflt>(&Eps::print).
        print_func<FormFormat>(&Eps::form_print).
        print_func<FCFormat>(&Eps::fc_print)
    )

    Eps::Eps(const Vector &x1, const Vector &x2, const Vector &x3, const Vector &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Vector &x1, const Vector &x2, const Vector &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Vector &x1, const Vector &x2, const Index &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Vector &x1, const Index &x2, const Index &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Index &x1, const Index &x2, const Index &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(vector<Vector> vs, vector<Index> is) {
        int i=0;
        for(auto vi : vs) pis[i++] = vi;
        for(auto ii : is) pis[i++] = ii;
    }

    int Eps::compare_same_type(const basic &other) const {
        const Eps &o = static_cast<const Eps &>(other);
        for(int i=0; i<4; i++) {
            auto c = pis[i].compare(o.pis[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    void Eps::print(const print_dflt &c, unsigned level) const {
        c.s << "𝜀" << "(";
        for(int i=0; i<3; i++) c.s << pis[i] << ",";
        c.s << pis[3] << ")";
    }
    
    void Eps::form_print(const FormFormat &c, unsigned level) const {
        c << "e_(";
        for(int i=0; i<3; i++) c << pis[i] << ",";
        c << pis[3] << ")";
    }
    
    void Eps::fc_print(const FCFormat &c, unsigned level) const {
        c << "LCD[";
        bool first = true;
        for(int i=3; i>=0; i--) {
            if(is_a<Vector>(pis[i]) && first) {
                c << "][";
                first = false;
            }
            c << pis[i];
            if(i!=0) c << ",";
            else c << "]";
        }
    }
    
    size_t Eps::nops() const { return 4; }
    ex Eps::op(size_t i) const {
        return pis[i];
    }
    
    //-----------------------------------------------------------
    // LC function - Levi-Civita
    //-----------------------------------------------------------
    ex LC(ex pi1, ex pi2, ex pi3, ex pi4) {
        bool isEps = true;
        lst pis = lst {pi1, pi2, pi3, pi4};
        for(auto pi : pis) {
            if(!is_a<Vector>(pi) && !is_a<Index>(pi)) {
                isEps = false;
                break;
            }
        }
        
        if(isEps) {
            vector<Vector> vs;
            vector<Index> is;
            ex sign = 1;
            for(int i=0; i<4; i++) {
                if(is_a<Vector>(pis.op(i))) {
                    vs.push_back(ex_to<Vector>(pis.op(i)));
                    sign *= pow(-1, is.size());
                } else is.push_back(ex_to<Index>(pis.op(i)));
            }
            return Eps(vs, is);
        }
        
        for(int i=0; i<4; i++) {
            auto pi = pis.op(i).expand();
            lst alst;
            if(is_a<add>(pi)) {
                for(auto item : pi) alst.append(item);
            } else alst.append(pi);
            
            for(int i=0; i<alst.nops(); i++) {
                if(!is_a<mul>(alst.op(i))) alst.let_op(i) = lst{alst.op(i)};
                ex c=1;
                ex v=1;
                for(auto ii : alst.op(i)) {
                    if(is_a<Vector>(ii) || is_a<Index>(ii)) {
                        if(!is_zero(v-1)) throw Error("Error Found in LC");
                        v = ii;
                    } else c *= ii;
                }
                if(is_zero(v-1)) throw Error("Error Found in LC");
                alst.let_op(i) = lst{c,v};
            }
            pis.let_op(i) = alst;
        }
        
        ex res = 0;
        for(int i0=0; i0<pis.op(0).nops(); i0++)
        for(int i1=0; i1<pis.op(1).nops(); i1++)
        for(int i2=0; i2<pis.op(2).nops(); i2++)
        for(int i3=0; i3<pis.op(3).nops(); i3++) {
            res + LC(get_op(pis,0,i0), get_op(pis,1,i1), get_op(pis,2,i2), get_op(pis,3,i3));
        }
        return res;
    }

}

