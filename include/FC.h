#pragma once

#define DEFAULT_CTOR(classname) \
classname::classname() { setflag(status_flags::evaluated | status_flags::expanded); }

#define IMPLEMENT_HAS(classname) \
bool classname::has(const ex &e) { \
    for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) if(is_a<classname>(*i)) return true; \
    return false; \
}


#include "Basic.h"
#include "Process.h"

namespace HepLib::FC {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    extern realsymbol D;
    extern realsymbol CA;
    extern realsymbol CF;
    extern realsymbol NA;
    extern realsymbol NF;
    extern realsymbol gs;
    extern exmap sp_map;
    
    class Index;
    class Vector;
    class Pair;
    
    //-----------------------------------------------------------
    // Filed/Propagator/Vertex Function
    //-----------------------------------------------------------
    DECLARE_FUNCTION_3P(Propagator)
    DECLARE_FUNCTION_3P(InField)
    DECLARE_FUNCTION_3P(OutField)
    DECLARE_FUNCTION_3P(Matrix)
    
    // Field function 2-3
    class Field2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function Field(const T1 & p1, const T2 & p2) {
        return GiNaC::function(Field2_SERIAL::serial, ex(p1), ex(p2));
    }
    class Field3_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3>
    inline GiNaC::function Field(const T1 & p1, const T2 & p2, const T3 & p3) {
        return GiNaC::function(Field3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
    }
    
    // Vertex function 2-6
    class Vertex2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function Vertex(const T1 & p1, const T2 & p2) {
        return GiNaC::function(Vertex2_SERIAL::serial, ex(p1), ex(p2));
    }
    class Vertex3_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3>
    inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3) {
        return GiNaC::function(Vertex3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
    }
    class Vertex4_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4>
    inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4) {
        return GiNaC::function(Vertex4_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4));
    }
    class Vertex5_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5) {
        return GiNaC::function(Vertex5_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5));
    }
    class Vertex6_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6) {
        return GiNaC::function(Vertex6_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5), ex(p6));
    }
    
    
    //-----------------------------------------------------------
    // Qgraf Class
    //-----------------------------------------------------------
    class Qgraf {
    public:
        string Output;
        string Model;
        string In;
        string Out;
        int Loops;
        string Options;
        vector<string> Others;
        ex Amplitudes(symtab st);
        
        static Index LI(ex fn);
        static Index DI(ex fn);
        static Index TI(ex fn);
        static Index FI(ex fn);
        static Index CI(ex fn);
        static Index AI(ex fn);
        
        static ex QuarkPropagator(ex e, ex m=0);
        static ex GluonPropagator(ex e);
        static ex GhostPropagator(ex e);
        static ex q2gVertex(ex e);
        static ex g3Vertex(ex e);
        static ex g4Vertex(ex e);
        static ex gh2gVertex(ex e);
        
    };
    
    //-----------------------------------------------------------
    // FormFormat Output
    //-----------------------------------------------------------
    class FormFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(FormFormat, print_dflt)
    public:
        FormFormat(ostream &os, unsigned opt=0);
        static void power_print(const power & p, const FormFormat & c, unsigned level=0);
        OUT_FORMAT_DECLARE(FormFormat)
        // static initialization
        class _init {
            public: _init();
        };
    private:
        static _init FormFormat_init;
    };
        
    //-----------------------------------------------------------
    // FCFormat Output
    //-----------------------------------------------------------
    class FCFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(FCFormat, print_dflt)
    public:
        FCFormat(ostream &os, unsigned opt=0);
        static void ncmul_print(const ncmul & p, const FCFormat & c, unsigned level=0);
        OUT_FORMAT_DECLARE(FCFormat)
        // static initialization
        class _init {
            public: _init();
        };
    private:
        static _init FCFormat_init;
    };
    
    class Error : public exception {
    public:
        string msg;
        const char * what() const throw ();
        Error(const char * _msg);
    };
    
    //-----------------------------------------------------------
    // Index Class
    //-----------------------------------------------------------
    class Index : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Index, basic)
    public:
        enum Type {VD, CF, CA};
        Index(const string &s, const Type type=Type::VD);
        Pair operator() (const Index & i);
        Pair operator() (const Vector & p);
        symbol name;
        Type type;
        void print(const print_context &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(Index);
    
    //-----------------------------------------------------------
    // Vector Class
    //-----------------------------------------------------------
    class Vector : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Vector, basic)
    public:
        Vector(const string &s);
        Pair operator() (const Vector & p);
        Pair operator() (const Index & mu);
        symbol name;
        void print(const print_context &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(Vector);
    
    //-----------------------------------------------------------
    // SUNT/SUNF Class
    //-----------------------------------------------------------
    class SUNT : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(SUNT, basic)
    public:
        SUNT(Index i, Index j, Index a);
        Index ija[3];
        size_t nops() const override;
        ex op(size_t i) const override;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void print(const print_dflt &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        ex derivative(const symbol & s) const override;
        ex conjugate() const override;
    };
    GINAC_DECLARE_UNARCHIVER(SUNT);
    
    class SUNF : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(SUNF, basic)
    public:
        SUNF(Index i, Index j, Index k);
        Index ijk[3];
        size_t nops() const override;
        ex op(size_t i) const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(SUNF);
    
    ex SUNF4(Index i, Index j, Index k, Index l);
    ex SUNSimplify(const ex & inexpr);
    
    //-----------------------------------------------------------
    // Pair Class
    //-----------------------------------------------------------
    class Pair : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Pair, basic)
    public:
        Pair(const Vector &p1, const Vector &p2);
        Pair(const Index &i1, const Index &i2);
        Pair(const Vector &p, const Index &i);
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        ex derivative(const symbol & s) const override;
        void operator = (const ex & e);
    private:
        ex lr[2];
    };
    GINAC_DECLARE_UNARCHIVER(Pair);
    ex SP(ex a);
    ex SP(ex a, ex b);
    void letSP(ex a, ex b, ex ab);
    void letSP(ex a, ex a2);
    void resetSP();
    
    //-----------------------------------------------------------
    // Eps Class
    //-----------------------------------------------------------
    class Eps : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Eps, basic)
    public:
        Eps(const Vector &p1, const Vector &p2, const Vector &p3, const Vector &p4);
        Eps(const Vector &p1, const Vector &p2, const Vector &p3, const Index &i1);
        Eps(const Vector &p1, const Vector &p2, const Index &i1, const Index &i2);
        Eps(const Vector &p1, const Index &i1, const Index &i2, const Index &i3);
        Eps(const Index &i1, const Index &i2, const Index &i3, const Index &i4);
        Eps(vector<Vector> vs, vector<Index> is);
        size_t nops() const override;
        ex op(size_t i) const override;
        ex & let_op(size_t i) override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        ex derivative(const symbol & s) const override;
    private:
        ex pis[4];
    };
    GINAC_DECLARE_UNARCHIVER(Eps);
    ex LC(ex pi1, ex pi2, ex pi3, ex pi4);
    
    //-----------------------------------------------------------
    // DiracGamma Class
    //-----------------------------------------------------------
    class DiracGamma : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(DiracGamma, basic)
    public:
        DiracGamma(const Vector &p, unsigned rl=0);
        DiracGamma(const Index &i, unsigned rl=0);
        DiracGamma(int int_1567, unsigned _rl=0);
        DiracGamma(const DiracGamma &g, unsigned _rl);
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        return_type_t return_type_tinfo() const override;
        unsigned return_type() const override { return return_types::noncommutative; }
        bool match_same_type(const basic & other) const override;
        unsigned get_rl();
        size_t nops() const override;
        ex op(size_t i) const override;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        ex derivative(const symbol & s) const override;
        ex conjugate() const override;
    private:
        ex pi;
        unsigned rl;
    };
    GINAC_DECLARE_UNARCHIVER(DiracGamma);
    
    //-----------------------------------------------------------
    // TR/GAS functions
    //-----------------------------------------------------------
    DECLARE_FUNCTION_1P(TR)
    
    inline ex GAS(const Vector &p) { return DiracGamma(p); }
    inline ex GAS(const Index &i) { return DiracGamma(i); }
    ex GAS(ex expr);
    
    // Form, TIR
    ex form(const ex &expr, bool verb=false, bool all=true);
    ex TIR(const lst &vis, const lst &eps);
    
    
    //-----------------------------------------------------------
    // Quarkonium
    //-----------------------------------------------------------
    namespace Quarkonium {
        enum IO {In, Out};
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu);
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu);
        ex ColorProj(Index i, Index j, Index a);
        ex ColorProj();
    }
    
    
}
