/**
 * @file
 * @brief Basic Functions for FC
 */
 
#include "HEP.h"
#include "flint/ulong_extras.h"

namespace HepLib {

    DEFAULT_CTOR(Index)
    IMPLEMENT_HAS(Index)
    IMPLEMENT_ALL(Index)
    
    DEFAULT_CTOR(Vector)
    IMPLEMENT_HAS(Vector)
    IMPLEMENT_ALL(Vector)
    
    DEFAULT_CTOR(SUNT)
    IMPLEMENT_HAS(SUNT)
    IMPLEMENT_ALL(SUNT)
    
    DEFAULT_CTOR(SUNF)
    IMPLEMENT_HAS(SUNF)
    IMPLEMENT_ALL(SUNF)
    
    DEFAULT_CTOR(SUNF4)
    IMPLEMENT_HAS(SUNF4)
    IMPLEMENT_ALL(SUNF4)

    //-----------------------------------------------------------
    // FormFormat Output
    //-----------------------------------------------------------
    FormFormat::FormFormat(ostream &os, unsigned opt) : print_dflt(os, opt) {}
    FormFormat::FormFormat() : print_dflt(std::cout) {}
    GINAC_IMPLEMENT_PRINT_CONTEXT(FormFormat, print_dflt)
    
    const FormFormat & FormFormat::operator << (const basic & v) const {
        v.print(*this);
        return *this;
    }
    const FormFormat & FormFormat::operator << (const ex & v) const {
        v.print(*this);
        return *this;
    }
    const FormFormat & FormFormat::operator << (const lst & v) const {
        v.print(*this);
        return *this;
    }
    const FormFormat & FormFormat::operator<<(std::ostream& (*v)(std::ostream&)) const {
        s << v;
        return *this;
    }
    
    void FormFormat::power_print(const power & p, const FormFormat & c, unsigned level) {
        if(p.op(1)==2 && !DGamma::has(p)) {
            c << "((" << p.op(0) << ")*(" << p.op(0) << "))";
        } else {
            c << "((" << p.op(0) << ")^(" << p.op(1) << "))";
        }
    }
    
    //-----------------------------------------------------------
    // FCFormat Output
    //-----------------------------------------------------------
    FCFormat::FCFormat(ostream &os, unsigned opt) : print_dflt(os, opt) {}
    FCFormat::FCFormat() : print_dflt(std::cout) {}
    GINAC_IMPLEMENT_PRINT_CONTEXT(FCFormat, print_dflt)
    
    const FCFormat & FCFormat::operator << (const basic & v) const {
        v.print(*this);
        return *this;
    }
    const FCFormat & FCFormat::operator << (const ex & v) const {
        v.print(*this);
        return *this;
    }
    const FCFormat & FCFormat::operator << (const lst & v) const {
        v.print(*this);
        return *this;
    }
    const FCFormat & FCFormat::operator<<(std::ostream& (*v)(std::ostream&)) const {
        s << v;
        return *this;
    }
    
    const FCFormat & FCFormat::operator << (const matrix & mat) const {
        s << "{";
        int nr = mat.rows();
        int nc = mat.cols();
        for(int r=0; r<nr; r++) {
            s << "{";
            for(int c=0; c<nc; c++) {
                mat(r,c).print(*this);
                if(c+1!=nc) s << ",";
            }
            s << "}";
            if(r+1!=nr) s << ",";
        }
        s << "}";
        return *this;
    }
    
    const FCFormat & FCFormat::operator << (const exvector & e) const {
        auto i = e.begin();
        auto vend = e.end();
        if (i==vend) { s << "{}"; return *this; }
        s << "{";
        while (true) {
            i->print(*this);
            ++i;
            if(i==vend) break;
            s << ",";
        }
        s << "}";
        return *this;
    }

    const FCFormat & FCFormat::operator << (const exset & e) const {
        auto i = e.begin();
        auto send = e.end();
        if (i==send) { s << "{}"; return *this; }
        s << "{";
        while (true) {
            i->print(*this);
            ++i;
            if(i==send) break;
            s << ",";
        }
        s << "}";
        return *this;
    }

    const FCFormat & FCFormat::operator << (const exmap & e) const {
        auto i = e.begin();
        auto mend = e.end();
        if (i==mend) { s << "{}"; return *this; }
        s << "{";
        while (true) {
            i->first.print(*this);
            s << "->";
            i->second.print(*this);
            ++i;
            if(i==mend) break;
            s << ",";
        }
        s << "}";
        return *this;
    }
    
    namespace {
        class ncmul_hack : public ncmul { // due to printseq is protected
        public:
            ncmul_hack(ncmul _nm) : ncmul(_nm){ }
            void print(const FCFormat & c, unsigned level) {
                printseq(c, '(', '.', ')', precedence(), level);
            }
        };
        ex mat_conj(const ex & e1, const ex & e2, const ex & e3) {
            return GMat(e1.conjugate(), e3, e2);
        }
    }
    void FCFormat::ncmul_print(const ncmul & nm, const FCFormat & c, unsigned level) {
        ncmul_hack(nm).print(c, level);
    }

    //-----------------------------------------------------------
    // Index Class
    //-----------------------------------------------------------
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Index, basic,print_func<print_context>(&Index::print))
    GiNaC::registered_class_info & Index::get_class_info_static() { return reg_info; }
    Index::visitor::~visitor() { }
    Index * Index::duplicate() const { Index * bp = new Index(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void Index::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &Index::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &Index::get_class_info() { return get_class_info_static(); }
    const char *Index::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    Index::Index(const string &s, const Type t) : name(s), type(t) { }
    int Index::compare_same_type(const basic &other) const {
        if(!is_a<Index>(other)) throw Error("Index::compare_same_type");
        const Index &o = static_cast<const Index &>(other);
        auto ret = name.get_name().compare(o.name.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    bool Index::is_equal_same_type(const basic & other) const {
        if(!is_a<Index>(other)) throw Error("Index::is_equal_same_type");
        const Index &o = static_cast<const Index &>(other);
        return (name.get_name() == o.name.get_name());
    }
    
    void Index::print(const print_context &c, unsigned level) const {
        c.s << name;
    }
    
    Pair Index::operator() (const Index &i) {
        return Pair(*this, i);
    }
    
    Pair Index::operator() (const Vector &p) {
        return Pair(p, *this);
    }
    
    void Index::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_string("name", name.get_name());
        n.add_unsigned("type", type);
    }
    
    void Index::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        string nstr;
        unsigned t;
        n.find_string("name", nstr);
        name = Symbol(nstr);
        n.find_unsigned("type", t);
        type = (Type)t;
    }
    
    ex Index::derivative(const symbol & s) const {
        return 0;
    }
    
    bool Index::hasc(const ex & e) {
        for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) 
            if(is_a<Index>(*i) && ex_to<Index>(*i).type!=Index::Type::VD) return true; 
        return false; 
    }
    
    bool Index::hasv(const ex & e) {
        for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i)
            if(is_a<Index>(*i) && ex_to<Index>(*i).type==Index::Type::VD) return true;
        return false;
    }

    //-----------------------------------------------------------
    // Vector Class
    //-----------------------------------------------------------
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(Vector, basic,print_func<print_context>(&Vector::print))
    GiNaC::registered_class_info & Vector::get_class_info_static() { return reg_info; }
    Vector::visitor::~visitor() { }
    Vector * Vector::duplicate() const { Vector * bp = new Vector(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void Vector::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &Vector::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &Vector::get_class_info() { return get_class_info_static(); }
    const char *Vector::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    Vector::Vector(const string &s) : name(s) { }
    int Vector::compare_same_type(const basic &other) const {
        if(!is_a<Vector>(other)) throw Error("Vector::compare_same_type");
        const Vector &o = static_cast<const Vector &>(other);
        auto ret = name.get_name().compare(o.name.get_name());
        if(ret==0) return 0;
        else if(ret<0) return -1;
        else return 1;
    }
    
    bool Vector::is_equal_same_type(const basic & other) const {
        if(!is_a<Vector>(other)) throw Error("Vector::is_equal_same_type");
        const Vector &o = static_cast<const Vector &>(other);
        return (name.get_name() == o.name.get_name());
    }
    
    void Vector::print(const print_context &c, unsigned level) const {
        c.s << name;
    }
    
    Pair Vector::operator() (const Vector &p) {
        return Pair(*this, p);
    }
    
    Pair Vector::operator() (const Index &i) {
        return Pair(*this, i);
    }
    
    void Vector::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_string("name", name.get_name());
    }
    
    void Vector::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        string nstr;
        unsigned t;
        n.find_string("name", nstr);
        name = Symbol(nstr);
    }
    
    ex Vector::derivative(const symbol & s) const {
        return 0;
    }
    
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(SUNT, basic,print_func<print_dflt>(&SUNT::print).print_func<FormFormat>(&SUNT::form_print).print_func<FCFormat>(&SUNT::fc_print))
    GiNaC::registered_class_info & SUNT::get_class_info_static() { return reg_info; }
    SUNT::visitor::~visitor() { }
    SUNT * SUNT::duplicate() const { SUNT * bp = new SUNT(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void SUNT::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &SUNT::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &SUNT::get_class_info() { return get_class_info_static(); }
    const char *SUNT::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    SUNT::SUNT(ex a, ex i, ex j) : aij{a,i,j} { }
    int SUNT::compare_same_type(const basic &other) const {
        if(!is_a<SUNT>(other)) throw Error("SUNT::compare_same_type");
        const SUNT &o = static_cast<const SUNT &>(other);
        for(int i=0; i<3; i++) {
            auto c = aij[i].compare(o.aij[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    bool SUNT::is_equal_same_type(const basic & other) const {
        if(!is_a<SUNT>(other)) throw Error("SUNT::is_equal_same_type");
        const SUNT &o = static_cast<const SUNT &>(other);
        for(int i=0; i<3; i++) {
            if(!aij[i].is_equal(o.aij[i])) return false;
        }
        return true;
    }
    
    void SUNT::form_print(const FormFormat &c, unsigned level) const {
        if(is_a<lst>(aij[0])) {
            bool first = true;
            for(auto item : aij[0]) {
                if(first) { first=false; c << "T(" << item; }
                else c << "," << item;
            }
        } else c << "T(" << aij[0];
        c << "," << aij[1] << "," << aij[2] << ")";
    }
    
    void SUNT::fc_print(const FCFormat &c, unsigned level) const {
        c << "SUNTF[" << aij[0] << "," << aij[1] << "," << aij[2] << "]";
    }
    
    void SUNT::print(const print_dflt &c, unsigned level) const {
        c.s << "T(" << aij[0] << "," << aij[1] << "," << aij[2] << ")";
    }
    
    size_t SUNT::nops() const { return 3; }
    ex SUNT::op(size_t i) const {
        return aij[i];
    }
    ex & SUNT::let_op(size_t i) {
        ensure_if_modifiable();
        return aij[i];
    }
    
    void SUNT::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("a", aij[0]);
        n.add_ex("i", aij[1]);
        n.add_ex("j", aij[2]);
    }
    
    void SUNT::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        ex o;
        n.find_ex("a", o);
        aij[0] = ex_to<Index>(o);
        n.find_ex("i", o);
        aij[1] = ex_to<Index>(o);
        n.find_ex("j", o);
        aij[2] = ex_to<Index>(o);
    }
    
    ex SUNT::derivative(const symbol & s) const {
        return 0;
    }
    
    ex SUNT::conjugate() const {
        if(!is_a<lst>(aij[0])) return SUNT(aij[0], aij[2], aij[1]);
        lst argv = ex_to<lst>(aij[0]);
        lst as;
        for(auto it=argv.rbegin(); it!=argv.rend(); ++it) as.append(*it);
        return SUNT(as, aij[2], aij[1]);
    }
    
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(SUNF, basic,print_func<print_dflt>(&SUNF::print).print_func<FormFormat>(&SUNF::form_print).print_func<FCFormat>(&SUNF::fc_print))
    GiNaC::registered_class_info & SUNF::get_class_info_static() { return reg_info; }
    SUNF::visitor::~visitor() { }
    SUNF * SUNF::duplicate() const { SUNF * bp = new SUNF(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void SUNF::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &SUNF::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &SUNF::get_class_info() { return get_class_info_static(); }
    const char *SUNF::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    SUNF::SUNF(ex i, ex j, ex k) : ijk{i,j,k} { }
    int SUNF::compare_same_type(const basic &other) const {
        if(!is_a<SUNF>(other)) throw Error("SUNF::compare_same_type");
        const SUNF &o = static_cast<const SUNF &>(other);
        for(int i=0; i<3; i++) {
            auto c = ijk[i].compare(o.ijk[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    bool SUNF::is_equal_same_type(const basic & other) const {
        if(!is_a<SUNF>(other)) throw Error("SUNF::is_equal_same_type");
        const SUNF &o = static_cast<const SUNF &>(other);
        for(int i=0; i<3; i++) {
            if(!ijk[i].is_equal(o.ijk[i])) return false;
        }
        return true;
    }
    
    ex SUNF::eval() const { 
        if(flags & status_flags::evaluated) return *this;
        if(ijk[0].is_equal(ijk[1]) || ijk[1].is_equal(ijk[2]) || ijk[0].is_equal(ijk[2])) return 0;
        bool c01 = ex_less(ijk[0],ijk[1]);
        bool c12 = ex_less(ijk[1],ijk[2]);
        if(c01 && c12) return this->hold();
        bool c02 = ex_less(ijk[0],ijk[2]);
        if(!c01 && c02) return -SUNF(ijk[1],ijk[0],ijk[2]);
        else if(c02 && !c12) return -SUNF(ijk[0],ijk[2],ijk[1]);
        else if(!c02 && c01) return SUNF(ijk[2],ijk[0],ijk[1]);
        else if(c12 && !c02) return SUNF(ijk[1],ijk[2],ijk[0]);
        else if(!c12 && !c01) return -SUNF(ijk[2],ijk[1],ijk[0]);
        else return this->hold();
    }
    
    void SUNF::print(const print_dflt &c, unsigned) const {
        c.s << "f(" << ijk[0] << "," << ijk[1] << "," << ijk[2] << ")";
    }
    
    void SUNF::form_print(const FormFormat &c, unsigned) const {
        c << "f(" << ijk[0] << "," << ijk[1] << "," << ijk[2] << ")";
    }
    
    void SUNF::fc_print(const FCFormat &c, unsigned) const {
        c << "SUNF[" << ijk[0] << "," << ijk[1] << "," << ijk[2] << "]";
    }
    
    size_t SUNF::nops() const { return 3; }
    ex SUNF::op(size_t i) const {
        return ijk[i];
    }
    ex & SUNF::let_op(size_t i) {
        ensure_if_modifiable();
        return ijk[i];
    }
    
    void SUNF::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("i", ijk[0]);
        n.add_ex("j", ijk[1]);
        n.add_ex("k", ijk[2]);
    }
    
    void SUNF::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        ex o;
        n.find_ex("i", o);
        ijk[0] = ex_to<Index>(o);
        n.find_ex("j", o);
        ijk[1] = ex_to<Index>(o);
        n.find_ex("k", o);
        ijk[2] = ex_to<Index>(o);
    }
    
    /**
     * @brief set derivative of SUNF to 0
     * @param s the symbol which the derivative to
     * @return always 0
     */
    ex SUNF::derivative(const symbol & s) const {
        return 0;
    }
    
    //GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(SUNF4, basic,print_func<print_dflt>(&SUNF4::print).print_func<FormFormat>(&SUNF4::form_print).print_func<FCFormat>(&SUNF4::fc_print))
    GiNaC::registered_class_info & SUNF4::get_class_info_static() { return reg_info; }
    SUNF4::visitor::~visitor() { }
    SUNF4 * SUNF4::duplicate() const { SUNF4 * bp = new SUNF4(*this); bp->setflag(GiNaC::status_flags::dynallocated); return bp; }
    void SUNF4::accept(GiNaC::visitor & v) const { if (visitor *p = dynamic_cast<visitor *>(&v)) p->visit(*this); else inherited::accept(v); }
    const GiNaC::registered_class_info &SUNF4::get_class_info() const { return get_class_info_static(); }
    GiNaC::registered_class_info &SUNF4::get_class_info() { return get_class_info_static(); }
    const char *SUNF4::class_name() const { return get_class_info_static().options.get_name(); }
    //GINAC_IMPLEMENT_REGISTERED_CLASS END
    
    SUNF4::SUNF4(ex i, ex j, ex k, ex l) : ijkl{i,j,k,l} { }
    int SUNF4::compare_same_type(const basic &other) const {
        if(!is_a<SUNF4>(other)) throw Error("SUNF4::compare_same_type");
        const SUNF4 &o = static_cast<const SUNF4 &>(other);
        for(int i=0; i<4; i++) {
            auto c = ijkl[i].compare(o.ijkl[i]);
            if(c!=0) return c;
        }
        return 0;
    }
    
    bool SUNF4::is_equal_same_type(const basic & other) const {
        if(!is_a<SUNF4>(other)) throw Error("SUNF4::is_equal_same_type");
        const SUNF4 &o = static_cast<const SUNF4 &>(other);
        for(int i=0; i<4; i++) {
            if(!ijkl[i].is_equal(o.ijkl[i])) return false;
        }
        return true;
    }
    
    /**
     * @brief automatical evaluation of SUNF4
     * @return expression with SUNF4 expanded
     */
    ex SUNF4::eval() const {
        if(flags & status_flags::evaluated) return *this;
        if(ijkl[0].is_equal(ijkl[1]) || ijkl[2].is_equal(ijkl[3])) return 0;
        bool c01 = ex_less(ijkl[0],ijkl[1]);
        bool c23 = ex_less(ijkl[2],ijkl[3]);
        if(c01 && c23) return this->hold(); // 01-23
        else if(!c01 && c23) return -SUNF4(ijkl[1],ijkl[0],ijkl[2],ijkl[3]);
        else if(!c01 && !c23) return SUNF4(ijkl[1],ijkl[0],ijkl[3],ijkl[2]);
        else if(c01 && !c23) return -SUNF4(ijkl[0],ijkl[1],ijkl[3],ijkl[2]);
        else return this->hold();
    }
    
    /**
     * @brief normal priint
     * @param c print_dflat object
     * @param o unsigned option
     */
    void SUNF4::print(const print_dflt &c, unsigned o) const {
        c.s << "f(" << ijkl[0] << "," << ijkl[1] << "," << ijkl[2] << "," << ijkl[3] << ")";
    }
    
    void SUNF4::fc_print(const FCFormat &c, unsigned o) const {
        c << "SUNF[" << ijkl[0] << "," << ijkl[1] << "," << ijkl[2] << "," << ijkl[3] << "]";
    }
    
    /**
     * @brief print the Form Format
     * @param c FormFormat object
     * @param o unsigned option, not used
     */
    void SUNF4::form_print(const FormFormat &c, unsigned o) const {
        c << "f4(" << ijkl[0] << "," << ijkl[1] << "," << ijkl[2] << "," << ijkl[3] << ")";
    }
    
    size_t SUNF4::nops() const { return 4; }
    ex SUNF4::op(size_t i) const {
        return ijkl[i];
    }
    ex & SUNF4::let_op(size_t i) {
        ensure_if_modifiable();
        return ijkl[i];
    }
    
    /**
     * @brief save to archvie
     * @param n the archive_node
     */
    void SUNF4::archive(archive_node & n) const {
        inherited::archive(n);
        n.add_ex("i", ijkl[0]);
        n.add_ex("j", ijkl[1]);
        n.add_ex("k", ijkl[2]);
        n.add_ex("l", ijkl[3]);
    }
    
    /**
     * @brief read from archive
     * @param n the archive_node
     */
    void SUNF4::read_archive(const archive_node& n) {
        inherited::read_archive(n);
        ex o;
        n.find_ex("i", o);
        ijkl[0] = ex_to<Index>(o);
        n.find_ex("j", o);
        ijkl[1] = ex_to<Index>(o);
        n.find_ex("k", o);
        ijkl[2] = ex_to<Index>(o);
        n.find_ex("l", o);
        ijkl[3] = ex_to<Index>(o);
    }
    
    /**
     * @brief set derivative of SUNF4 to 0
     * @param s the variable that derivative to
     * @return 0
     */
    ex SUNF4::derivative(const symbol & s) const {
        return 0;
    }
    
    ex ncmul_expand(const ex & expr) {
        MapFunction inner_expand([](const ex & e, MapFunction & self)->ex{
            if(is_a<add>(e)) {
                ex res = 0;
                for(auto ei : e) res += ncmul_expand(ei);
                return res;
            } else if(is_a<mul>(e) || is_a<ncmul>(e)) {
                lst res = lst{ 1 };
                for(auto ei : e) {
                    ex rei = ncmul_expand(ei);
                    if(!is_a<add>(rei)) rei = lst{ rei };
                    lst ores = res;
                    res = lst{ };
                    for(auto oi : ores) for(auto ri : rei) res.append(oi * ri);
                }
                ex ret = 0;
                for(auto ri : res) ret += ri;
                return ret;
            } else return e.map(self);
        });
        return inner_expand(expr);
    }
    
    /**
     * @brief make contract on matrix, i.e., GMat(a,i1,i2)*GMat(b,i2,i3) -> GMat(a*b,i1,i3)
     * @param expr_in expression contains GMat
     * @return contracted expression
     */
    ex GMatContract(const ex & expr_in) {
        if(!expr_in.has(GMat(w1,w2,w3))) return expr_in;
        
        auto expr = expr_in.subs(pow(GMat(w1,w2,w3),2)==GMat(w1,w2,w3)*GMat(w1,w3,w2));
        auto cv_lst = collect_lst(expr, GMat(w1, w2, w3));
        expr = 0;
        for(auto cv : cv_lst) {
            auto e = cv.op(1);
            if(is_zero(e-1) || e.match(GMat(w1, w2, w3))) {
                if(e.match(GMat(w1, w2, w2))) expr += cv.op(0) * TR(e.op(0));
                else expr += cv.op(0) * e;
                continue;
            } else if(!is_a<mul>(e)) throw Error("GMatContract: collect error: " + ex2str(e));
            
            lst mats;
            for(auto item : e) mats.append(item);
            
            std::map<ex,int,ex_is_less> to_map, from_map;
            std::set<int> todo;
            lst mats_idx;
            for(int i=0; i<mats.nops(); i++) {
                auto item = mats.op(i);
                if(item.op(0).return_type()==return_types::commutative || item.op(0).is_equal(GAS(1))) {
                    mats_idx.append(lst{item,i});
                } else {
                    if(to_map[item.op(1)]!=0 || from_map[item.op(2)]!=0) throw Error("GMatContract: index conflict for "+ex2str(item));
                    to_map[item.op(1)] = i+10; // avoid 0 in map
                    from_map[item.op(2)] = i+10; // avoid 0 in map
                }
                todo.insert(i);
            }
            
            // update to_map/from_map w.r.t mats_idx
            bool checked = false;
            while(true) {
                lst mats_idx2;
                bool ok = true; // double check
                for(int i=0; i<mats_idx.nops(); i++) {
                    auto item = mats_idx.op(i).op(0);
                    int ii = ex_to<numeric>(mats_idx.op(i).op(1)).to_int();
                    if(!checked && 
                        to_map[item.op(1)]==0 && from_map[item.op(2)]==0 && 
                        to_map[item.op(2)]==0 && from_map[item.op(1)]==0) {
                        mats_idx2.append(mats_idx.op(i));
                        continue;
                    }
                    ok = false;
                    checked = false;
                    if(to_map[item.op(1)]==0 && from_map[item.op(2)]==0) {
                        to_map[item.op(1)] = ii+10; // avoid 0 in map
                        from_map[item.op(2)] = ii+10; // avoid 0 in map
                    } else if(to_map[item.op(2)]==0 && from_map[item.op(1)]==0) {
                        to_map[item.op(2)] = ii+10; // avoid 0 in map
                        from_map[item.op(1)] = ii+10; // avoid 0 in map
                        // need to swap the 2nd and 3rd index
                        auto li = get_op(mats, ii, 1); 
                        auto ri = get_op(mats, ii, 2);
                        let_op(mats, ii, 1, ri);
                        let_op(mats, ii, 2, li);
                    } else {
                        throw Error("GMatContract: index conflict (2).");
                    }
                }
                if(mats_idx2.nops()<1) break;
                mats_idx = mats_idx2;
                if(ok) checked=true;
            }
            
            ex retMat = 1;
            while(todo.size()>0) {
                int c = *(todo.begin());
                todo.erase(c);
                ex curMat = mats.op(c).op(0);
                auto li=mats.op(c).op(1);
                auto ri=mats.op(c).op(2);
                while(true) {
                    if(li.is_equal(ri)) {
                        retMat *= TR(curMat);
                        break;
                    }
                    int ti = to_map[ri];
                    int fi = from_map[li];
                    if(ti==0 && fi==0) {
                        retMat *= GMat(curMat, li, ri);
                        break;
                    }
                    if(ti!=0) {
                        auto mat = mats.op(ti-10).op(0);
                        if(curMat.is_equal(GAS(1))) curMat = mat;
                        else if(!mat.is_equal(GAS(1))) curMat = curMat * mat;
                        ri = mats.op(ti-10).op(2);
                        todo.erase(ti-10);
                        continue;
                    }
                    if(fi!=0) {
                        auto mat = mats.op(fi-10).op(0);
                        if(curMat.is_equal(GAS(1))) curMat = mat;
                        else if(!mat.is_equal(GAS(1))) curMat = mat * curMat;
                        li = mats.op(fi-10).op(1);
                        todo.erase(fi-10);
                        continue;
                    }
                }
            }
            expr += cv.op(0) * retMat;
        }
                
        return expr;
    }
    
    ex Contract(const ex & ei) {
        lst idx_lst;
        MapFunction get_idx([&idx_lst](const ex & e, MapFunction & self)->ex{
            if(!Index::has(e) || !Pair::has(e) || e.match(GMat(w1,w2,w3))) return 1; // skip GMat object
            else if(is_a<Pair>(e)) {
                if(is_a<Index>(e.op(0)) || is_a<Index>(e.op(1))) idx_lst.append(e);
                return 1;
            } else return e.map(self);
        });
        get_idx(ei);
        idx_lst.sort();
        idx_lst.unique();
        if(idx_lst.nops()==0) return ei;
        auto cvs = collect_lst(ei, idx_lst);
        ex res = 0;
        for(auto cv : cvs) {
            auto c = cv.op(0);
            auto v = form(cv.op(1)); // contract on itself
            if(!Index::has(v)) res += c * v;
            else {
                if(!is_a<mul>(v)) v = lst{ v };
                exmap repl;
                ex r = 1; // uncontracted remained index
                for(auto vi : v) {
                    if(!is_a<Pair>(vi)) r *= vi; // contract may result in a non-Pair object
                    else if(is_a<Index>(vi.op(1)) && c.has(vi.op(1))) repl[vi.op(1)] = vi.op(0);
                    else if(is_a<Index>(vi.op(0)) && c.has(vi.op(0))) repl[vi.op(0)] = vi.op(1);
                    else r *= vi;
                }
                res += r * c.subs(repl);
            }
        }
        return res.subs(SP_map);
    }
    
    ex GMatExpand(const ex & expr_in) {
        MapFunction inner_expand([&](const ex & e, MapFunction & self)->ex {
            if(!e.has(GMat(w1,w2,w3))) return e;
            else if(e.match(GMat(w1,w2,w3))) {
                auto e0 = e.op(0);
                if(is_a<add>(e0)) {
                    ex res = 0;
                    for(auto item : e0) res += GMatExpand(GMat(item, e.op(1), e.op(2)));
                    return res;
                } else if(is_a<mul>(e0)) {
                    ex c = 1, v = 1;
                    for(auto item : e0) {
                        if(item.return_type()==return_types::commutative) c *= item;
                        else {
                            if(!v.is_equal(1)) {
                                cout << "c=" << c << ", " << "v=" << v << endl;
                                throw Error("GMatExpand: v != 1"); // make sure only one non-commutative object
                            }
                            v = item;
                        }
                    }
                    if(v.is_equal(1)) v = GAS(1); 
                    return c * GMatExpand(GMat(v, e.op(1), e.op(2)));
                } else if(is_a<ncmul>(e0)) { // expand ncmul
                    ex res;
                    bool first = true;
                    for(auto item : e0) {
                        if(first) {
                            res = item;
                            first = false;
                            continue;
                        }
                        ex ncL = res; // previous result
                        if(!is_a<add>(ncL)) ncL = lst{ ncL };
                        ex ncR = item;
                        if(!is_a<add>(ncR)) ncR = lst{ ncR };
                        res = 0; // current result
                        for(auto iL : ncL) for(auto iR : ncR) res += iL * iR;
                    }
                    ex rs = res;
                    res = 0;
                    if(!is_a<add>(rs)) rs = lst{ rs };
                    for(auto item : rs) { // pull out commutative coefficient
                        ex c = 1, v = 1;
                        if(is_a<mul>(item)) {
                            if(item.nops()==1) throw Error("GMatExpand: item.nops == 1"); // make sure
                            for(auto it : item) {
                                if(it.return_type()==return_types::commutative) c *= it;
                                else {
                                    if(!v.is_equal(1)) throw Error("GMatExpand: v != 1"); // make sure only one non-commutative object
                                    v = it;
                                }
                            }
                        } else v = item;
                        
                        while(true) { // recursive replace ɣ.P * ɣ.P -> P^2 and ɣ.mu * ɣ.mu -> d @ v
                            bool to_exit = true;
                            if(is_a<ncmul>(v)) {
                                bool first = true;
                                ex last = 1, vv = 1;
                                for(auto vi : v) {
                                    if(first) {
                                        first = false;
                                        last = vi;
                                    } else {
                                        if(last==vi && is_a<DGamma>(vi)) {
                                            first = true;
                                            last = 1;
                                            if(is_a<Vector>(vi.op(0))) c *= SP(vi.op(0));
                                            else if(is_a<Index>(vi.op(0))) c *= d;
                                            else if(vi.op(0).is_equal(1) || vi.op(0).is_equal(5)) c *= 1; // GAS(1)*GAS(1) = GAS(5)*GAS(5) = 1
                                            else throw Error("GMatExpand: only GAS(i/p/1/5) supported.");
                                            to_exit = false; // need to cycle again
                                        } else {
                                            if(last!=GAS(1) && !last.is_equal(1)) vv = vv * last;
                                            last = vi;
                                        }
                                    }
                                }
                                if(!last.is_equal(1) && last!=GAS(1)) v = vv * last; // check last item
                                else v = vv;
                            }
                            if(to_exit) break;
                        }
                        if(v.is_equal(1)) v = GAS(1); // identity matrix
                        res += c * GMat(v, e.op(1), e.op(2));
                    }
                    return res;
                } else return e;
            } else return e.map(self);
        });
        return inner_expand(expr_in);
    }
    
    ex GMatShift(const ex & expr, const ex & g, bool to_right) {
        if(!expr.has(g)) return expr;
        MapFunction inner_shift([g,to_right](const ex & e, MapFunction & self)->ex{
            if(!e.has(g) || !e.has(GMat(w1,w2,w3))) return e;
            else if(e.match(GMat(w1,w2,w3))) {
                ex eg = e.op(0);
                if(!is_a<ncmul>(eg)) eg = lst{ eg };
                int gi = -1;
                if(to_right) {
                    for(int i=0; i<eg.nops()-1; i++) if(eg.op(i)==g) { gi = i; break; }
                    if(gi==-1 || gi==eg.nops()-1) return e;
                } else {
                    for(int i=eg.nops()-1; i>=0; i--) if(eg.op(i)==g) { gi = i; break; }
                    if(gi==-1 || gi==0) return e;
                }
                int gj = gi + ( to_right ? 1 : -1 );
                ex rem = 1, rem2 = 1;
                for(int i=0; i<eg.nops(); i++) {
                    if(i!=gi && i!=gj) {
                        rem *= eg.op(i);
                        rem2 *= eg.op(i);
                    }
                    if(i==gi) {
                        if(to_right) rem2 *= eg.op(gj)*eg.op(gi);
                        else rem2 *= eg.op(gi)*eg.op(gj);
                    }
                }
                if(eg.op(gi).is_equal(eg.op(gj))) {
                    ex ip = eg.op(gi).op(0);
                    if(rem.is_equal(1)) rem = GAS(1);
                    ex res = GMat(rem, e.op(1), e.op(2));
                    res = GMatShift(res, g, to_right);
                    ex c;
                    if(is_a<Vector>(ip)) c = SP(ip);
                    else if(is_a<Index>(ip)) c = d;
                    else throw Error("GMatShift: only GAS(i/p) supproted.");
                    return c * res;
                }
                if(rem.is_equal(1)) rem = GAS(1);
                if(rem2.is_equal(1)) rem2 = GAS(1);
                ex res = 2*SP(eg.op(gi).op(0), eg.op(gj).op(0)) * GMat(rem, e.op(1), e.op(2));
                res = res - GMat(rem2, e.op(1), e.op(2));
                return GMatShift(res, g, to_right);
            } else return e.map(self);
        });
        ex res = GMatExpand(Contract(expr)); // add Contract & GMatExpand here
        res = collect_ex(res, GMat(w1,w2,w3));
        res = inner_shift(res);
        return res;
    }
    
    ex GMatSimplify(const ex & expr) {
        ex res = GMatContract(expr);
        res = GMatShift(res);
        res = Contract(res);
        return res;
    }
    
    namespace {
        lst shift_12_right(const ex & e) { // return a list of {coeff, gammas}
            if(!is_a<ncmul>(e)) throw Error("input is not a ncmul.");
            if(e.nops()==2) {
                ex e0 = e.op(0);
                if(e.op(0)!=e.op(1)) throw Error("2 items are not equal!");
                else if(is_a<Index>(e0.op(0))) return lst{ lst{ d, 1 }};
                else if(is_a<Vector>(e0.op(0))) return lst{ lst{ SP(e.op(0).op(0)), 1 }};
                else {
                    cout << endl << e << endl;
                    throw Error("shift_12_right: only GAS(i/p) supproted.");
                }
            }
            ex rem = 1;
            int n = e.nops();
            for(int i=2; i<n; i++) rem *= e.op(i);
            lst res = shift_12_right(e.op(0)*rem);
            n = res.nops();
            for(int i=0; i<n; i++) {
                res.let_op(i).let_op(0) = -res.op(i).op(0);
                res.let_op(i).let_op(1) = e.op(1) * res.op(i).op(1);
            }
            if(!is_a<Index>(e.op(0).op(0)) && !is_a<Vector>(e.op(0).op(0))) {
                cout << e << endl;
                throw Error("shift_12_right: not a Vector or Index");
            }
            if(!is_a<Index>(e.op(1).op(0)) && !is_a<Vector>(e.op(1).op(0))) {
                cout << e << endl;
                throw Error("shift_12_right: not a Vector or Index");
            }
            res.append(lst{ 2*SP(e.op(0).op(0), e.op(1).op(0)), rem });
            return res;
        }
    }
    ex GMatShift(const ex & expr) {
        MapFunction inner_shift([](const ex & e, MapFunction & self)->ex{
            if(!e.has(GMat(w1,w2,w3))) return e;
            else if(e.match(GMat(w1,w2,w3))) {
                ex eg = e.op(0);
                if(!is_a<ncmul>(eg)) eg = lst{ eg };
                
                int gi = -1, gj = -1;
                for(int i=0; i<eg.nops(); i++) for(int j=i+1; j<eg.nops(); j++) {
                    if(eg.op(i).is_equal(eg.op(j))) {
                        gi = i;
                        gj = j;
                        goto done;
                    }
                }
                return e;
                done: ;
                
                ex exL = 1, exM=1, exR = 1;
                for(int i=0; i<eg.nops(); i++) {
                    if(i<gi) exL *= eg.op(i);
                    else if(i>gj) exR *= eg.op(i);
                    else exM *= eg.op(i);
                }
                lst cvs = shift_12_right(exM);
                
                ex res = 0;
                for(auto cv : cvs) {
                    ex item = exL*cv.op(1)*exR;
                    if(item.is_equal(1)) item = GAS(1);
                    res += cv.op(0)*GMatShift(GMat(item, e.op(1), e.op(2)));
                }
                
                return res;
            } else return e.map(self);
        });
        ex res = GMatExpand(Contract(expr)); // add Contract & GMatExpand here
        res = collect_ex(res, GMat(w1,w2,w3));
        res = inner_shift(res);
        return res;
    }
    
    namespace {
        void GMat_fc_print(const ex &arg1, const ex &arg2, const ex &arg3, const print_context &c0) {
            auto c = static_cast<const FCFormat &>(c0);
            c << "GMat[" << arg1 << "," << arg2 << "," << arg3 << "]";
        }
    }
    
    REGISTER_FUNCTION(GMat, do_not_evalf_params().print_func<FCFormat>(&GMat_fc_print).conjugate_func(mat_conj).set_return_type(return_types::commutative))
    
    bool IsZero(const ex & e) {
        try {
            exset vs;
            for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
                if(is_a<symbol>(*i) || is_a<Pair>(*i)) vs.insert(*i);
            }
            
            int n = 13;
            for(int i=0; i<5; i++) {
                exmap nsubs;
                for(auto item : vs) {
                    nsubs[item] = ex(1)/n_nth_prime(n);
                    n++;
                }
                ex ret = e.subs(nsubs);
                if(!normal(e).is_zero()) return false;
            }
            
            auto ret = exnormal(e);
            return ret.is_zero();
        } catch(...) { }
        return is_zero(exnormal(e));
    }
    
    ex ToCF(const ex & e) {
        ex res = e;
        bool todo = true;
        while(todo) {
            todo = false;
            auto cvs = collect_lst(res, lst{NF,TF});
            res = 0;
            for(auto cv : cvs) {
                int degTF = cv.op(1).degree(TF);
                int degNF = cv.op(1).ldegree(NF);
                if(degTF>0 && degNF<0) {
                    todo = true;
                    int n = degTF;
                    if(degTF+degNF>0) n = -degNF;
                    res += cv.op(0) * cv.op(1) * pow(TF/NF,-n) * pow(TF*NF-CF,n);
                } else if(degTF>0 && degNF>1) {
                    todo = true;
                    int n = degTF;
                    if(degTF>degNF/2) n = degNF/2;
                    res += cv.op(0) * cv.op(1) * pow(TF*NF*NF,-n) * pow(CF*NF+TF,n);
                } else res += cv.op(0) * cv.op(1);
            }
        }
        return res;
    }
    
    ex ca_neg_pow_sub(const ex & expr) {
        static MapFunction ca_map([](const ex & e, MapFunction & self)->ex {
            if(!e.has(CA)) return e;
            else if(e.match(pow(CA,w)) && e.op(1).info(info_flags::negint)) return pow(CA-2*CF,-e.op(1));
            else return e.map(self);
        });
        return ca_map(expr);
    }
    
    ex ToCACF(const ex & e) { // from FeynCalc
        ex res = e.subs(lst{NA==(NF*NF-1),CF==(NF*NF-1)/(2*NF),TF==ex(1)/2});
        res = exfactor(res);
        // SUNN -> CA
        res = res.subs(NF==CA);
        // (-1+CA^2)->(-2 CA CF)
        res = res.subs(lst{ w*(1-CA)*(1+CA)==-w*2*CA*CF, w*(-1+CA)*(1+CA)==w*2*CA*CF });
        res = res.subs(lst{ w1*pow(1-CA,w2)*pow(1+CA,w2)==w1*pow(-2*CA*CF,w2), w1*pow(-1+CA,w2)*pow(1+CA,w2)==w1*pow(2*CA*CF,w2) });
        // (((2 - CA^2) CF )/CA ) ->(CF (CA - 4 CF))
        res = res.subs(lst{ w*(2-CA*CA)*CF/CA==w*CF*(CA-4*CF), w*(-2+CA*CA)*CF/CA==w*CF*(-CA+4*CF) });
        // (1-CA^2) -> (-2 CA CF)
        res = res.subs(lst{ w*(1-CA)*(1+CA)==-w*2*CA*CF, w*(-1+CA)*(1+CA)==w*2*CA*CF });
        res = res.subs(lst{ w1*pow(1-CA,w2)*pow(1+CA,w2)==-w1*pow(2*CA*CF,w2), w1*pow(-1+CA,w2)*pow(1+CA,w2)==w1*pow(2*CA*CF,w2) });
        // (1/CA)^n -> (CA - 2 CF)^n
		//res = res.subs(lst{ w/CA==w*(CA-2*CF) });
        res = ca_neg_pow_sub(res);
        // ((1 - CA^2)*(CA - 2*CF)) -> (-2*CF)
		res = res.subs(lst{	w*(1-CA)*(1+CA)*(CA-2*CF)==-2*w*CF, w*(-1+CA)*(1+CA)*(CA-2*CF)==2*w*CF });
        // (CA (CA-2 CF)) -> 1
		res = res.subs(lst{	w*CA*(CA-2*CF)==w, w*CA*(-CA+2*CF)==-w });
        res = res.subs(lst{	w1*pow(CA,w2)*pow(CA-2*CF,w2)==w1, w1*pow(CA,w2)*pow(-CA+2*CF,w2)==w*pow(-1,w2) });
        // (CA^2+c)(CA-2CF) -> CA+c(CA-2CF)
        res = res.subs(lst{	w0*(CA*CA+w1)*(CA-2*CF)==w0*(CA+w1*(CA-2*CF)), w0*(CA*CA+w1)*(-CA+2*CF)==-w0*(CA+w1*(CA-2*CF)) });
        res = res.subs(lst{	w0*(CA*CA+w1)*pow(CA-2*CF,w2)==w0*(CA+w1*(CA-2*CF))*pow(CA-2*CF,w2-1), w0*(CA*CA+w1)*pow(-CA+2*CF,w2)==-w0*(CA+w1*(CA-2*CF))*pow(-CA+2*CF,w2-1) });
        return res;
    }
    
    ex HomCACF(const ex & e) {
        ex res = e.subs(lst{NA==(NF*NF-1),CA==NF,CF==(NF*NF-1)/(2*NF),TF==ex(1)/2});
        res = exfactor(res);
        if(!is_a<mul>(res)) res = lst{ res };
        ex c=1, v=1;
        for(auto item : res) {
            if(item.has(NF)) c *= item;
            else v *= item;
        }
        c = collect_ex(c, NF);
        int deg = c.degree(NF);
        int ldeg = -c.ldegree(NF);
        if(ldeg>deg) deg = ldeg;
        if(deg>0) {
            lst vars;
            ex eqn = c;
            for(int i=0; i<=deg; i++) {
                symbol xi;
                eqn -= xi * pow(NF,i) * pow((NF*NF-1)/(2*NF), deg-i);
                vars.append(xi);
            }
            eqn = collect_ex(eqn, NF);
            int nH = eqn.degree(NF);
            int nL = eqn.ldegree(NF);
            lst eqns;
            for(int i=nL; i<=nH; i++) {
                ex cc = eqn.coeff(NF, i);
                if(cc.is_zero()) continue;
                eqns.append(cc==0);
            }
            auto sol = lsolve(eqns, vars);
            if(sol.nops()!=deg+1) {
                cout << "c=" << c << endl;
                cout << "sol=" << sol << endl;
                throw Error("HomCACF: no solution found!");
            }
            c = 0;
            for(int i=0; i<=deg; i++) c += vars.op(i).subs(sol) * pow(CA,i) * pow(CF, deg-i);
        } 
        return c * v;
    }
    
    ex DoColor(const ex & expr, const ex & pref, int method) {
        auto cvs = collect_lst(expr, [](const ex &e)->bool{ return Index::hasc(e); });
        ex res = 0;
        for(auto const & cv : cvs) {
            auto cc = cv.op(0);
            auto vv = cv.op(1);
            if(method==0) vv = HomCACF(form(vv)/pref)*pref;
            else vv = ToCACF(form(vv)/pref)*pref;
            res += cc * vv;
        }
        return res;
    }
}

