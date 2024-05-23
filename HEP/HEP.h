/**
 * @file 
 * @brief HEP header file
 */
 
#pragma once

#include "BASIC.h"
#include "IBP.h"

namespace HepLib {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    extern exmap SP_map;
    
    extern int form_trace_mode;
    extern const int form_trace_auto;
    extern const int form_trace_all;
    extern const int form_trace_each_all;
    extern const int form_trace_each_each;
    
    extern int form_expand_mode;
    extern const int form_expand_none;
    extern const int form_expand_tr;
    extern const int form_expand_ci;
    extern const int form_expand_li;
    extern const int form_expand_all;
    
    extern bool Apart_using_fermat;
    extern bool form_using_su3;
    extern bool form_using_dim4;
    extern bool form_using_gamma5;
    
    class Index;
    class Vector;
    class Pair;
        
    /**
     * @brief class for FormFormat Output
     */
    class FormFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(FormFormat, print_dflt)
    public:
        FormFormat(ostream &os, unsigned opt=0);
        static void power_print(const power & p, const FormFormat & c, unsigned level=0);
        
        template<class T> const FormFormat & operator << (const T & v) const {
            s << v;
            return *this;
        };
        const FormFormat & operator << (const basic & v) const;
        const FormFormat & operator << (const ex & v) const;
        const FormFormat & operator << (const lst & v) const;
        const FormFormat & operator<<(std::ostream& (*v)(std::ostream&)) const;
        
        #ifndef DOXYGEN_SKIP
        class _init {
            public: _init();
        };
    private:
        static _init FormFormat_init;
        #endif
    };
        
    /**
     * @brief class for FCFormat Output
     */
    class FCFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(FCFormat, print_dflt)
    public:
        FCFormat(ostream &os, unsigned opt=0);
        static void ncmul_print(const ncmul & p, const FCFormat & c, unsigned level=0);
        
        template<class T> const FCFormat & operator << (const T & v) const {
            s << v;
            return *this;
        };
        const FCFormat & operator << (const basic & v) const;
        const FCFormat & operator << (const ex & v) const;
        const FCFormat & operator << (const lst & v) const;
        const FCFormat & operator<<(std::ostream& (*v)(std::ostream&)) const;
        
        const FCFormat & operator << (const matrix & v) const;
        const FCFormat & operator << (const exvector & v) const;
        const FCFormat & operator << (const exmap & v) const;
        const FCFormat & operator << (const exset & v) const;
        
        #ifndef DOXYGEN_SKIP
        class _init {
            public: _init();
        };
    private:
        static _init FCFormat_init;
        #endif
    };
    extern FCFormat fcout;
    
    /**
     * @brief class for index object
     */
    class Index : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(Index, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::exmap Dimension;
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const Index &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        Index(); // classname
        Index * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        enum Type {VD, CF, CA};
        Index(const string &s, const Type type=Type::VD);
        Pair operator() (const Index & i);
        Pair operator() (const Vector & p);
        Symbol name;
        Type type;
        void print(const print_context &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static bool hasc(const ex &e);
        static bool hasv(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;        
    };
    
    /**
     * @brief class for vector object
     */
    class Vector : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(Vector, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const Vector &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        Vector(); // classname
        Vector * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        Vector(const string &s);
        Pair operator() (const Vector & p);
        Pair operator() (const Index & mu);
        Symbol name;
        void print(const print_context &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;        
    };
    
    /**
     * @brief class for SUNT object
     */
    class SUNT : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(SUNT, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const SUNT &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        SUNT(); // classname
        SUNT * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        SUNT(ex a, ex i, ex j);
        ex aij[3]; // Index
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void print(const print_dflt &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        ex conjugate() const override;
        bool is_equal_same_type(const basic & other) const override;        
    };
    
    /**
     * @brief class for SUNF object
     */
    class SUNF : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(SUNF, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const SUNF &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        SUNF(); // classname
        SUNF * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        SUNF(ex i, ex j, ex k);
        ex ijk[3]; // Index
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        ex eval() const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;        
    };
    
    
    /**
     * @brief class for SUNF4 object
     */
    class SUNF4 : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(SUNF4, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const SUNF4 &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        SUNF4(); // classname
        SUNF4 * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        SUNF4(ex i, ex j, ex k, ex l);
        ex ijkl[4]; // Index
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        ex eval() const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;        
    };
    
    /**
     * @brief class for Pair object
     */
    class Pair : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(Pair, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const Pair &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        Pair(); // classname
        Pair * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        Pair(const Vector &p1, const Vector &p2);
        Pair(const Index &i1, const Index &i2);
        Pair(const Vector &p, const Index &i);
        Pair(const Index &i, const Vector &p);
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        ex eval() const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;
    private:
        ex lr[2];
    };
    
    ex SP(const ex &a, bool use_map=false);
    ex SP(const ex &a, const ex &b, bool use_map=false);
    ex sp(const ex & a, const ex & b);
    ex sp(const ex & a);
    ex& letSP(const ex &p1, const ex &p2);
    ex& letSP(const ex &p);
    void clearSP(const ex &p1, const ex &p2);
    void clearSP(const ex &p);
    void clearSP();
    ex SP2sp(const ex & exin);
    exmap sp_map();
    
    /**
     * @brief class for Levi-Civita object
     * https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527630097.app3
     * to make Tr(g5, g1, g2, g3, g4) is the same in both HepLib & FORM, require that
     * Eps(a,b,c,d) = i_ * e_(a,b,c,d) (we use the convention as in FeynCalc, Tr[5,a,b,c,d]=(- i) 4 Eps(a,b,c,d)=4 eps_(a,b,c,d)), and Eps^{0123}=+1, and g5=i g^{0123}=(-i) Eps(a,b,c,d) gamma(a,b,c,d)/4!
     * Eps is real in HepLib, while e_ is imaginary in FORM.
     */
    class Eps : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(Eps, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const Eps &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        Eps(); // classname
        Eps * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        exvector pis;
        Eps(const Vector &p1, const Vector &p2, const Vector &p3, const Vector &p4);
        Eps(const Vector &p1, const Vector &p2, const Vector &p3, const Index &i1);
        Eps(const Vector &p1, const Vector &p2, const Index &i1, const Index &i2);
        Eps(const Vector &p1, const Index &i1, const Index &i2, const Index &i3);
        Eps(const Index &i1, const Index &i2, const Index &i3, const Index &i4);
        Eps(vector<Vector> vs, vector<Index> is);
        Eps(const exvector & pis0);
        size_t nops() const override;
        ex op(size_t i) const override;
        ex & let_op(size_t i) override;
        ex eval() const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;
    };
    ex LC(ex pi1, ex pi2, ex pi3, ex pi4);
    
    /**
     * @brief class for Dirac Gamma object
     */
    class DGamma : public basic {
    //GINAC_DECLARE_REGISTERED_CLASS(DGamma, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const DGamma &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        DGamma(); // classname
        DGamma * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        ex pi;
        unsigned rl;
        DGamma(const Vector &p, unsigned rl=0);
        DGamma(const Index &i, unsigned rl=0);
        DGamma(int int_1567, unsigned _rl=0);
        DGamma(const DGamma &g, unsigned _rl);
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        return_type_t return_type_tinfo() const override;
        unsigned return_type() const override { return return_types::noncommutative; }
        bool match_same_type(const basic & other) const override;
        unsigned get_rl();
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        ex eval() const override;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        ex conjugate() const override;
        bool is_equal_same_type(const basic & other) const override;
    };
    
    //-----------------------------------------------------------
    // TR/GAS functions
    //-----------------------------------------------------------
    DECLARE_FUNCTION_3P(Matrix)
    DECLARE_FUNCTION_1P(TR)
    DECLARE_FUNCTION_1P(TTR)
    DECLARE_FUNCTION_1P(HF)
    
    inline ex GAS(const Vector &p, unsigned rl=0) { return DGamma(p,rl); }
    inline ex GAS(const Index &i, unsigned rl=0) { return DGamma(i,rl); }
    ex GAS(const ex &expr, unsigned rl=0);
    
    // Form, TIR, Apart
    ex charge_conjugate(const ex &);
    ex form(const ex &expr, int verb=0);
    ex UnContract(const ex expr, const lst &loop_ps, const lst &ext_ps=lst{}); // Eps/DGamma always uncontract
    ex TIR(const ex &expr_in, const lst &loop_ps, const lst &ext_ps);
    ex MatrixContract(const ex & expr_in);
    ex Apart(const matrix & mat);
    ex Apart(const ex &expr_in, const lst &vars, exmap sgnmap={});
    ex Apart(const ex &expr_in, const lst &loops, const lst & extmoms, exmap sgnmap={});
    ex ApartIR2ex(const ex & expr_in);
    ex ApartIR2F(const ex & expr_in);
    ex F2ex(const ex & expr_in);
    ex ApartIRC(const ex & expr_in);
    void ApartIBP(exvector &io_vec, int IBPmethod, const lst & loops, const lst & exts,
        const lst & cut_props=lst{}, std::function<lst(const IBP &, const ex &)> uf=LoopUF);
    inline void ApartIBP(int IBPmethod, exvector &io_vec, const lst & loops, const lst & exts,
        const lst & cut_props=lst{}, std::function<lst(const IBP &, const ex &)> uf=LoopUF) {
        return ApartIBP(io_vec, IBPmethod, loops, exts, cut_props, uf);
    }
    exmap ApartRules(const exvector &airs, bool irc=true);
        
    struct AIOption {
        int pn_sector = 0; // treat each sector as a spearate problem when total denominators > pn_sector, set 0 to disable this feature
        bool ap_rules = true; // minimize the total number of ibp problems
        int IBPmethod = 1; // 0
        lst Internal; // Internal for Apart/IBP
        lst External; // External for Apart/IBP
        lst DSP; // DSP for IBP
        exmap smap; // Sign Map for Apart
        lst Cut; // Cut Propagator. optional
        lst CSP; // SP in Cut, to be cleared. optional
        lst ISP; // SP for IBP. optional
        bool CutFirst = true;
        bool keep0F = false; // keep 0 exponent in F
        int NIBP = 0;
        ex apart1 = 0; // set Apart(1,{x,y,...}) to apart1
        lst pat = { F(w1,w2), gs, nL, nH };
        std::function<ex(const ex &, const ex &)> cv = nullptr;
        string SaveDir = ""; // save temporary result, and restart from it
        std::function<lst(const IBP &, const ex &)> UF = LoopUF;
        void init_smap() { for(auto li : Internal) smap[SP(li)] = 1; }
    };
    void ApartIBP(exvector &io_vec, AIOption aip);
    
    bool IsZero(const ex & e);
    
    #ifndef DOXYGEN_SKIP
    
    class ApartIR1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function ApartIR(const T1 & p1) {
        return GiNaC::function(ApartIR1_SERIAL::serial, ex(p1));
    }
    
    class ApartIR2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function ApartIR(const T1 & p1, const T2 & p2) {
        return GiNaC::function(ApartIR2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    #endif
    
    ex ToCF(const ex & e);
    ex ToCACF(const ex & e);
    ex HomCACF(const ex & e);
    ex DoColor(const ex & e, const ex & pref=1, int method=0);
    ex A0(const ex m2, int n=1, const ex d=4-2*ep);
    
        
}

