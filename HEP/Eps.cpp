/**
 * @file
 * @brief Functions for Levi-Civita Tensor
 */
 
#include "HEP.h"

namespace HepLib {

    //-----------------------------------------------------------
    // Eps Class
    //-----------------------------------------------------------
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Eps, basic,print_func<print_dflt>(&Eps::print).print_func<FormFormat>(&Eps::form_print).print_func<FCFormat>(&Eps::fc_print))
    GiNaC::registered_class_info & Eps::get_class_info_static() { return reg_info; }
    Eps::visitor::~visitor() { }
    Eps * Eps::duplicate() const { Eps * bp = new Eps(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void Eps::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &Eps::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &Eps::get_class_info() { return get_class_info_static(); }
    const char *Eps::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    DEFAULT_CTOR(Eps)
    IMPLEMENT_HAS(Eps)
    IMPLEMENT_ALL(Eps)

    Eps::Eps(const Vector &x1, const Vector &x2, const Vector &x3, const Vector &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Vector &x1, const Vector &x2, const Vector &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Vector &x1, const Vector &x2, const Index &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Vector &x1, const Index &x2, const Index &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(const Index &x1, const Index &x2, const Index &x3, const Index &x4) : pis{x1,x2,x3,x4} { }
    Eps::Eps(vector<Vector> vs, vector<Index> is) {
        pis.resize(vs.size()+is.size());
        int i=0;
        for(auto vi : vs) pis[i++] = vi;
        for(auto ii : is) pis[i++] = ii;
    }
    Eps::Eps(const exvector & pis0) {
        int n = pis0.size();
        pis.resize(n);
        for(int i=0; i<n; i++) pis[i] = pis0[i];
    }

    int Eps::compare_same_type(const basic &other) const {
        if(!is_a<Eps>(other)) throw Error("Eps::compare_same_type");
        const Eps &o = static_cast<const Eps &>(other);
        int n = pis.size();
        int no = o.pis.size();
        if(n!=no) return n>no ? 1 : -1;
        for(int i=0; i<n; i++) {
            auto c = pis[i].compare(o.pis[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    bool Eps::is_equal_same_type(const basic & other) const {
        if(!is_a<Eps>(other)) throw Error("Eps::is_equal_same_type");
        const Eps &o = static_cast<const Eps &>(other);
        int n = pis.size();
        int no = o.pis.size();
        if(n!=no) return false;
        for(int i=0; i<n; i++) {
            if(!pis[i].is_equal(o.pis[i])) return false;
        }
        return true;
    }
    
    ex Eps::eval() const {
        if(flags & status_flags::evaluated) return *this;
        int n = pis.size();
        if(n==4) {
            bool ok = true;
            int ii = 4;
            for(int i=0; i<4; i++) {
                if(!is_a<Vector>(pis[i]) && !is_a<Index>(pis[i])) ok = false;
                else if(is_a<Vector>(pis[i]) && ii!=4) ok = false; // Vector after Index
                else if(is_a<Index>(pis[i]) && ii==4) ii = i;
                if(!ok) break;
            }
            if(!ok) return LC(pis[0],pis[1],pis[2],pis[3]);
            
            const ex* pis_a = pis.data();
            if(isSorted(ii,pis_a) && isSorted(4-ii,pis_a+ii)) return this->hold();
            else {
                ex pis2[4];
                for(int i=0; i<4; i++) pis2[i] = pis[i];
                int ac1 = ACSort(ii,pis2);
                int ac2 = ACSort(4-ii,pis2+ii);
                if(ac1 * ac2==0) return 0;
                return ac1 * ac2 * LC(pis2[0],pis2[1],pis2[2],pis2[3]);
            }
        } else return this->hold();
    }
    
    void Eps::print(const print_dflt &c, unsigned level) const {
        int n = pis.size();
        c.s << "\u03B5" << "(";
        for(int i=0; i<n-1; i++) c.s << pis[i] << ",";
        c.s << pis[n-1] << ")";
    }
    
    /**
     * @brief Eps in FORM format
     * https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527630097.app3
     * to make Tr(g5, g1, g2, g3, g4) is the same in both HepLib & FORM, require that
     * Eps(a,b,c,d) = i_ * e_(a,b,c,d) (we use the convention as in FeynCalc, Tr[5,a,b,c,d]=(- i) 4 Eps(a,b,c,d)=4 eps_(a,b,c,d)), and Eps^{0123}=+1, and g5=i g^{0123}=(-i) Eps(a,b,c,d) gamma(a,b,c,d)/4!
     * Eps is real in HepLib, while e_ is imaginary in FORM.
     * @param c the FormFormat
     * @param level the level
     */
    void Eps::form_print(const FormFormat &c, unsigned level) const {
        int n = pis.size();
        if(n==4) c << "(i_*"; // only multiple i for gamma5
        c << "e_(";
        for(int i=0; i<n-1; i++) c << pis[i] << ",";
        c << pis[n-1] << ")";
        if(n==4) c << ")";
    }
    
    void Eps::fc_print(const FCFormat &c, unsigned level) const {
        int n = pis.size();
        if(n==4) {
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
            //if(first) c << "[]";
        }
    }
    
    size_t Eps::nops() const { return pis.size(); }
    ex Eps::op(size_t i) const {
        return pis[i];
    }
    ex & Eps::let_op(size_t i) {
        ensure_if_modifiable();
        return pis[i];
    }
    
    void Eps::archive(archive_node & n) const {
        inherited::archive(n);
        int nn = pis.size();
        n.add_ex("size", nn);
        for(int i=0; i<nn; i++) n.add_ex("pis"+to_string(i), pis[i]);
    }
    
    void Eps::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        ex nex;
        n.find_ex("size", nex);
        int nn = ex_to<numeric>(nex).to_int();
        if(pis.size()!=nn) pis.resize(nn);
        for(int i=0; i<nn; i++) {
            n.find_ex("pis"+to_string(i), pis[i]);
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

