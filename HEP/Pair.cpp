/**
 * @file
 * @brief Pair class
 */
 
#include "HEP.h"

namespace HepLib {

    //-----------------------------------------------------------
    // Pair Class
    //-----------------------------------------------------------
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Pair, basic,print_func<print_dflt>(&Pair::print).print_func<FormFormat>(&Pair::form_print).print_func<FCFormat>(&Pair::fc_print))
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Eps, basic,print_func<print_dflt>(&Eps::print).print_func<FormFormat>(&Eps::form_print).print_func<FCFormat>(&Eps::fc_print))
    GiNaC::registered_class_info & Pair::get_class_info_static() { return reg_info; }
    Pair::visitor::~visitor() { }
    Pair * Pair::duplicate() const { Pair * bp = new Pair(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void Pair::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &Pair::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &Pair::get_class_info() { return get_class_info_static(); }
    const char *Pair::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    DEFAULT_CTOR(Pair)
    IMPLEMENT_HAS(Pair)
    IMPLEMENT_ALL(Pair)

    /**
     * @brief Pair constructor, ordered internally
     * @param p1 Vector in p1.p2
     * @param p2 Vector  in p1.p2
     * @return p1.p2
     */
    Pair::Pair(const Vector &p1, const Vector &p2) { 
        if(p1.compare(p2)>0) {
            lr[0]=p1;lr[1]=p2; 
        } else {
            lr[0]=p2;lr[1]=p1; 
        }
    }
    
    /**
     * @brief Pair constructor, ordered internally
     * @param i1 Index in i1.i2
     * @param i2 Index  in i1.i2
     * @return i1.i2 i.e., g^{i1,i2}
     */
    Pair::Pair(const Index &i1, const Index &i2) { 
        if(i1.compare(i2)>0) {
            lr[0]=i1;lr[1]=i2; 
        } else {
            lr[0]=i2;lr[1]=i1; 
        }
    }
    
    /**
     * @brief Pair constructor
     * @param p Vector in p^i
     * @param i Index  in p^i
     * @return p^i
     */
    Pair::Pair(const Vector &p, const Index &i)  : lr{p, i} { }
    
    /**
     * @brief Pair constructor
     * @param i Index  in p^i
     * @param p Vector in p^i
     * @return p^i
     */
    Pair::Pair(const Index &i, const Vector &p) : lr{p, i} { }

    int Pair::compare_same_type(const basic &other) const {
        if(!is_a<Pair>(other)) throw Error("Pair::compare_same_type");
        const Pair &o = static_cast<const Pair &>(other);
        int c1 = lr[0].compare(o.lr[0]);
        if(c1!=0) return c1;
        int c2 = lr[1].compare(o.lr[1]);
        return c2;
    }
    
    bool Pair::is_equal_same_type(const basic & other) const {
        if(!is_a<Pair>(other)) throw Error("Pair::is_equal_same_type");
        const Pair &o = static_cast<const Pair &>(other);
        if(!lr[0].is_equal(o.lr[0])) return false;
        return lr[1].is_equal(o.lr[1]);
    }
    
    /**
     * @brief default print function
     * @param c default print format
     * @param level level in print function
     */
    void Pair::print(const print_dflt &c, unsigned level) const {
        c.s << "(" << lr[0] << "." << lr[1] << ")";
    }
    
    /**
     * @brief FormFormat print function
     * @param c FormFormat print format
     * @param level level in print function
     */
    void Pair::form_print(const FormFormat &c, unsigned level) const {
        if(is_a<Vector>(lr[0]) && is_a<Vector>(lr[1])) c << lr[0] << "." << lr[1];
        else if(is_a<Vector>(lr[0]) && is_a<Index>(lr[1])) c << lr[0] << "(" << lr[1] << ")";
        else if(is_a<Index>(lr[0]) && is_a<Index>(lr[1])) c << "d_(" << lr[0] << "," << lr[1] << ")";
    }
    
    /**
     * @brief FCFormat print function
     * @param c FCFormat print format
     * @param level level in print function
     */
    void Pair::fc_print(const FCFormat &c, unsigned level) const {
        if(is_a<Vector>(lr[0]) && is_a<Vector>(lr[1])) c << "SPD[" << lr[0] << "," << lr[1] << "]";
        else if(is_a<Vector>(lr[0]) && is_a<Index>(lr[1])) c << "FVD[" << lr[0] << "," << lr[1] << "]";
        else if(is_a<Index>(lr[0]) && is_a<Index>(lr[1])) {
            auto ii = ex_to<Index>(lr[0]);
            if(ii.type == Index::Type::VD) c << "MTD[" << lr[0] << "," << lr[1] << "]";
            else if(ii.type == Index::Type::CF) c << "SUNFDelta[" << lr[0] << "," << lr[1] << "]";
            else if(ii.type == Index::Type::CA) c << "SUNDelta[" << lr[0] << "," << lr[1] << "]";
            else throw Error("Pair::fc_print unexpected.");
        }
    }
    
    size_t Pair::nops() const { return 2; }
    ex Pair::op(size_t i) const {
        return lr[i];
    }
    ex & Pair::let_op(size_t i) {
        ensure_if_modifiable();
        return lr[i];
    }
    
    ex Pair::eval() const {
        if(flags & status_flags::evaluated) return *this;
        else if(!is_a<Vector>(lr[0]) && !is_a<Index>(lr[0])) return SP(lr[0],lr[1]);
        else if(!is_a<Vector>(lr[1]) && !is_a<Index>(lr[1])) return SP(lr[0],lr[1]);
        else if((is_a<Vector>(lr[0]) && is_a<Vector>(lr[1]))) 
            return Pair(ex_to<Vector>(lr[0]),ex_to<Vector>(lr[1])).hold();
        else if((is_a<Index>(lr[0]) && is_a<Index>(lr[1]))) 
            return Pair(ex_to<Index>(lr[0]),ex_to<Index>(lr[1])).hold();
        else if((is_a<Index>(lr[0]) && is_a<Vector>(lr[1])))
            return Pair(ex_to<Vector>(lr[1]),ex_to<Index>(lr[0])).hold();
        else return this->hold();
    }
    
    void Pair::archive(archive_node & n) const {
        inherited::archive(n);
        for(int i=0; i<2; i++) n.add_ex("lr"+to_string(i), lr[i]);
    }
    
    void Pair::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        for(int i=0; i<2; i++) {
            n.find_ex("lr"+to_string(i), lr[i]);
        }
    }
    
    ex Pair::derivative(const symbol & s) const {
        return 0;
    }
    
    //-----------------------------------------------------------
    // SP function - ScalarProduct
    //-----------------------------------------------------------
    ex SP(const ex & a, bool use_map) { return SP(a,a,use_map); }
    
    /**
     * @brief Function similar to SPD/FVD in FeynCalc
     * @param a the 1st vector/index
     * @param b the 2nd vector/index
     * @param use_map true for apply the SP_map in the function call
     * @return expanded/translated Pair objects
     */
    ex SP(const ex & a, const ex & b, bool use_map) {
        ex ret_sp;
        if(is_a<Vector>(a) && is_a<Vector>(b)) ret_sp = Pair(ex_to<Vector>(a), ex_to<Vector>(b));
        else if(is_a<Vector>(a) && is_a<Index>(b)) ret_sp = Pair(ex_to<Vector>(a), ex_to<Index>(b));
        else if(is_a<Index>(a) && is_a<Vector>(b)) ret_sp = Pair(ex_to<Vector>(b), ex_to<Index>(a));
        else if(is_a<Index>(a) && is_a<Index>(b)) ret_sp = Pair(ex_to<Index>(a), ex_to<Index>(b));
        if(is_a<Pair>(ret_sp)) {
            if(use_map) return ret_sp.subs(SP_map);
            else return ret_sp;
        }

        lst alst, blst;
        auto aex = expand(a);
        auto bex = expand(b);
        if(is_zero(aex) || is_zero(bex)) return 0;
        if(is_a<add>(aex)) {
            for(auto item : aex) alst.append(item);
        } else alst.append(aex);
        for(int i=0; i<alst.nops(); i++) {
            if(!is_a<mul>(alst.op(i))) alst.let_op(i) = lst{alst.op(i)};
            ex c=1;
            ex v=1;
            for(auto ii : alst.op(i)) {
                if(is_a<Vector>(ii) || is_a<Index>(ii)) {
                    if(!is_zero(v-1)) throw Error("Error Found in SP @1 : "+ex2str(alst));
                    v = ii;
                } else c *= ii;
            }
            if(is_zero(v-1)) throw Error("Error Found in SP @2 : "+ex2str(alst));
            alst.let_op(i) = lst{c,v};
        }
        
        if(is_a<add>(bex)) {
            for(auto item : bex) blst.append(item);
        } else blst.append(bex);
        for(int i=0; i<blst.nops(); i++) {
            if(!is_a<mul>(blst.op(i))) blst.let_op(i) = lst{blst.op(i)};
            ex c=1;
            ex v=1;
            for(auto ii : blst.op(i)) {
                if(is_a<Vector>(ii) || is_a<Index>(ii)) {
                    if(!is_zero(v-1)) throw Error("Error Found in SP @3");
                    v = ii;
                } else c *= ii;
            }
            if(is_zero(v-1)) throw Error("Error Found in SP @4");
            blst.let_op(i) = lst{c,v};
        }
        
        ex res = 0;
        for(auto ai : alst) {
            for(auto bi : blst) res += ai.op(0) * bi.op(0) * SP(ai.op(1), bi.op(1), use_map);
        }
        return res;
    }
    
    /**
     * @brief translated the vector dot a.b to a*b, useful in SecDec
     * @param a the 1st vector argument or a symbol
     * @param b the 2nd vector or a symbol
     * @return expression with vector name left only
     */
    ex sp(const ex & a, const ex & b) {
        if(is_a<Vector>(a) && is_a<Vector>(b)) return ex_to<Vector>(a).name * ex_to<Vector>(b).name;
        else throw Error("Error in sp(2).");
    }
    ex sp(const ex & a) {
        if(is_a<Vector>(a)) return ex_to<Vector>(a).name * ex_to<Vector>(a).name;
        else throw Error("Error in sp(1).");
    }
    
    /**
     * @brief return the reference of p1.p2
     * @param p1 the 1st vector/index
     * @param p2 the 2nd vector/index
     * @return the refernce of p1.p2 in SP_map, can be used in assignment
     */
    ex & letSP(const ex &p1, const ex &p2) {
        if(!(is_a<Vector>(p1) || is_a<Index>(p1)) || !(is_a<Vector>(p2) || is_a<Index>(p2)))
            throw Error("Invalide arguments for letSP (2).");
        return SP_map[SP(p1,p2)];
    }
    
    /**
     * @brief return the reference of p.p
     * @param p the vector/index
     * @return the refernce of p.p in SP_map, can be used in assignment
     */
    ex & letSP(const ex &p) {
        if(!is_a<Vector>(p))
            throw Error("Invalide arguments for letSP (1).");
        return SP_map[SP(p)];
    }
    
    /**
     * @brief delete the assignment of p1.p2
     * @param p1 the 1st vector/index
     * @param p2 the 2nd vector/index
     */
    void clearSP(const ex &p1, const ex &p2) {
        if(!(is_a<Vector>(p1) || is_a<Index>(p1)) || !(is_a<Vector>(p2) || is_a<Index>(p2)))
            throw Error("Invalide arguments for resetSP (2).");
        SP_map.erase(SP(p1,p2));
    }
    
    /**
     * @brief delete the assignment of p.p
     * @param p vector/index
     */
    void clearSP(const ex &p) {
        if(!is_a<Vector>(p))
            throw Error("Invalide arguments for resetSP (1).");
        SP_map.erase(SP(p));
    }
    
    /**
     * @brief delete all assignment in SP_map
     */
    void clearSP() {
        SP_map.clear();
    }
    
    /**
     * @brief convert SP(a,b) to sp(a,b)
     * @param exin the input expression
     * @return expression with SP(a,b) replaced to sp(a,b)
     */
    ex SP2sp(const ex & exin) {
        auto ret = exin.subs(SP_map);
        lst vecs = Vector::all(ret);
        exmap spmap;
        for(auto vp1 : vecs) {
            for(auto vp2 : vecs) {
                spmap[SP(vp1,vp2)]=sp(vp1,vp2);
            }
        }
        return ret.subs(spmap);
    }
    
    /**
     * @brief the SP_map with SP(a,b) replaced to sp(a,b)
     * @return the SP_map with SP(a,b) replaced to sp(a,b)
     */
    exmap sp_map() {
        exmap ret;
        for(auto kv : SP_map) {
            auto key = kv.first;
            if(key.nops()!=2 || !is_a<Vector>(key.op(0)) || !is_a<Vector>(key.op(1))) continue;
            auto p1 = ex_to<Vector>(key.op(0));
            auto p2 = ex_to<Vector>(key.op(1));
            ret[p1.name * p2.name] = kv.second.subs(SP_map);
        }
        return ret;
    }

}

