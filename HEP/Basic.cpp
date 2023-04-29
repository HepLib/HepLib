/**
 * @file
 * @brief Basic Functions for FC
 */
 
#include "HEP.h"
#include "cln/cln.h"

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
            return Matrix(e1.conjugate(), e3, e2);
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
    
    /**
     * @brief make contract on matrix, i.e., Matrix(a,i1,i2)*Matrix(b,i2,i3) -> Matrix(a*b,i1,i3)
     * @param expr_in expression contains Matrix
     * @return contracted expression
     */
    ex MatrixContract(const ex & expr_in) {
        if(!expr_in.has(Matrix(w1,w2,w3))) return expr_in;
        
        auto expr = expr_in.subs(pow(Matrix(w1,w2,w3),2)==Matrix(w1,w2,w3)*Matrix(w1,w3,w2));
        auto cv_lst = collect_lst(expr, Matrix(w1, w2, w3));
        expr = 0;
        for(auto cv : cv_lst) {
            auto e = cv.op(1);
            if(is_zero(e-1) || e.match(Matrix(w1, w2, w3))) {
                expr += cv.op(0) * e;
                continue;
            } else if(!is_a<mul>(e)) throw Error("MatrixContract: collect error: " + ex2str(e));
            
            lst mats;
            for(auto item : e) mats.append(item);
            
            std::map<ex,int,ex_is_less> to_map, from_map;
            std::set<int> todo;
            lst mats_idx;
            for(int i=0; i<mats.nops(); i++) {
                auto item = mats.op(i);
                if(item.op(0).return_type()==return_types::commutative || item.op(0).is_equal(GAS(1))  || item.op(0).is_equal(color_ONE())) {
                    mats_idx.append(lst{item,i});
                } else {
                    if(to_map[item.op(1)]!=0 || from_map[item.op(2)]!=0) throw Error("MatrixContract: index conflict for "+ex2str(item));
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
                        throw Error("MatrixContract: index conflict (2).");
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
                    if(is_zero(li-ri)) {
                        retMat *= TR(curMat);
                        break;
                    }
                    int ti = to_map[ri];
                    int fi = from_map[li];
                    if(ti==0 && fi==0) {
                        retMat *= Matrix(curMat, li, ri);
                        break;
                    }
                    if(ti!=0) {
                        auto mat = mats.op(ti-10).op(0);
                        if(is_zero(curMat-GAS(1))) curMat = mat;
                        else if(!is_zero(mat-GAS(1))) curMat = curMat * mat;
                        ri = mats.op(ti-10).op(2);
                        todo.erase(ti-10);
                        continue;
                    }
                    if(fi!=0) {
                        auto mat = mats.op(fi-10).op(0);
                        if(is_zero(curMat-GAS(1))) curMat = mat;
                        else if(!is_zero(mat-GAS(1))) curMat = mat * curMat;
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
    
    namespace {
        void Matrix_fc_print(const ex &arg1, const ex &arg2, const ex &arg3, const print_context &c0) {
            auto c = static_cast<const FCFormat &>(c0);
            c << "Matrix[" << arg1 << "," << arg2 << "," << arg3 << "]";
        }
    }
    
    REGISTER_FUNCTION(Matrix, do_not_evalf_params().print_func<FCFormat>(&Matrix_fc_print).conjugate_func(mat_conj).set_return_type(return_types::commutative))
    
    bool IsZero(const ex & e) {
        exset vs;
        for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
            if(is_a<symbol>(*i) || is_a<Pair>(*i)) vs.insert(*i);
        }
        
        int n = 3;
        for(int i=0; i<5; i++) {
            exmap nsubs;
            for(auto item : vs) {
                nsubs[item] = ex(1)/numeric(cln::nextprobprime(n));
                n++;
            }
            ex ret = e.subs(nsubs);
            if(!normal(e).is_zero()) return false;
        }
        
        auto ret = normal_fermat(e);
        return ret.is_zero();
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
        return res;
    }
    
    ex DoColor(const ex & expr, const ex & pref) {
        auto cvs = collect_lst(expr, [](const ex &e)->bool{ return Index::hasc(e); });
        ex res = 0;
        for(auto const & cv : cvs) {
            auto cc = cv.op(0);
            auto vv = cv.op(1);
            vv = ToCACF(form(vv)/pref)*pref;
            res += cc * vv;
        }
        return res;
    }
}

