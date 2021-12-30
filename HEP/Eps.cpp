/**
 * @file
 * @brief Functions for Levi-Civita Tensor
 */
 
#include "HEP.h"

namespace HepLib {

    //-----------------------------------------------------------
    // Eps Class
    //-----------------------------------------------------------
    GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Eps, basic,
        print_func<print_dflt>(&Eps::print).
        print_func<FormFormat>(&Eps::form_print).
        print_func<FCFormat>(&Eps::fc_print)
    )
    
    DEFAULT_CTOR(Eps)
    GINAC_BIND_UNARCHIVER(Eps);
    IMPLEMENT_HAS(Eps)
    IMPLEMENT_ALL(Eps)

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
    
    bool Eps::is_equal_same_type(const basic & other) const {
        const Eps &o = static_cast<const Eps &>(other);
        for(int i=0; i<4; i++) {
            if(!pis[i].is_equal(o.pis[i])) return false;
        }
        return true;
    }
    
    ex Eps::eval() const {
        if(flags & status_flags::evaluated) return *this;
        bool ok = true;
        int ii = 4;
        for(int i=0; i<4; i++) {
            if(!is_a<Vector>(pis[i]) && !is_a<Index>(pis[i])) ok = false;
            else if(is_a<Vector>(pis[i]) && ii!=4) ok = false; // Vector after Index
            else if(is_a<Index>(pis[i]) && ii==4) ii = i;
            if(!ok) break;
        }
        if(!ok) return LC(pis[0],pis[1],pis[2],pis[3]);
        
        if(isSorted(ii,pis) && isSorted(4-ii,pis+ii)) return this->hold();
        else {
            ex pis2[4];
            for(int i=0; i<4; i++) pis2[i] = pis[i];
            int ac1 = ACSort(ii,pis2);
            int ac2 = ACSort(4-ii,pis2+ii);
            if(ac1 * ac2==0) return 0;
            return ac1 * ac2 * LC(pis2[0],pis2[1],pis2[2],pis2[3]);
        }
    }
    
    void Eps::print(const print_dflt &c, unsigned level) const {
        c.s << "\u03B5" << "(";
        for(int i=0; i<3; i++) c.s << pis[i] << ",";
        c.s << pis[3] << ")";
    }
    
    /**
     * @brief LC in FORM format
     * to make Tr(g5, g1, g2, g3, g4) is the same in both HepLib & FORM, require that
     * LC(a,b,c,d) = i_ * e_(a,b,c,d) ( we use the convention as in FeynCalc, Tr[5,1,2,3,4]=(- i) 4 LC[1,2,3,4])
     * LC is real in HepLib, while e_ is imaginary in FORM.
     * @param c the FormFormat
     * @param level the level
     */
    void Eps::form_print(const FormFormat &c, unsigned level) const {
        c << "(i_*e_(";
        for(int i=0; i<3; i++) c << pis[i] << ",";
        c << pis[3] << "))";
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
            
            if(i==0) c << "]";
            else if(!first || !is_a<Vector>(pis[i-1])) c << ",";
        }
        if(first) c << "[]";
    }
    
    size_t Eps::nops() const { return 4; }
    ex Eps::op(size_t i) const {
        return pis[i];
    }
    ex & Eps::let_op(size_t i) {
        ensure_if_modifiable();
        return pis[i];
    }
    
    void Eps::archive(archive_node & n) const {
        inherited::archive(n);
        for(int i=0; i<4; i++) n.add_ex("pis"+to_string(i), pis[i]);
    }
    
    void Eps::read_archive(const archive_node& n, lst& sym_lst) {
        inherited::read_archive(n, sym_lst);
        for(int i=0; i<4; i++) {
            n.find_ex("pis"+to_string(i), pis[i], sym_lst);
        }
    }
    
    ex Eps::derivative(const symbol & s) const {
        return 0;
    }
    
    /**
     * @brief function similar to LCD in FeynCalc
     * @param pi1 vector/index in the 1st position
     * @param pi2 vector/index in the 2nd position
     * @param pi3 vector/index in the 3rd position
     * @param pi4 vector/index in the 4th position
     * @return expanded/translated to Eps objects
     */
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
            return sign*Eps(vs, is);
        }
        
        for(int i=0; i<4; i++) {
            auto pi = collect_lst(pis.op(i), [](const ex & e)->bool{return Index::has(e) || Vector::has(e);});
            for(auto item : pi) { // check 
                if(!is_a<Vector>(item.op(1)) && !is_a<Index>(item.op(1))) {
                    cout << "pi = " << pi << endl;
                    throw Error("LC Error: there is no Index or Vector.");
                }
            }
            pis.let_op(i) = pi;
        }
        
        ex res = 0;
        for(auto i0 : pis.op(0))
        for(auto i1 : pis.op(1))
        for(auto i2 : pis.op(2))
        for(auto i3 : pis.op(3)) {
            res += i0.op(0)*i1.op(0)*i2.op(0)*i3.op(0) * LC(i0.op(1), i1.op(1), i2.op(1), i3.op(1));
        }
        return res;
    }

}

