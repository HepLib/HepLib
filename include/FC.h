/**
 * @file 
 * @brief FC header file
 */
 
#pragma once

#include "Basic.h"
#include "Process.h"
#include "IBP.h"

/**
 * @brief namespace for Index, Vector, DiracGamma, etc.
 */
namespace HepLib::FC {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    using namespace HepLib::IBP;
    
    extern const Symbol NA;
    extern const Symbol NF;
    extern const Symbol gs;
    extern const Symbol as;
    extern const Symbol mu;
    extern const Symbol nL;
    extern const Symbol nH;
    extern exmap SP_map;
    
    extern int trace_method;
    extern const int trace_all;
    extern const int trace_each_all;
    extern const int trace_each_each;
    
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
        OUT_FORMAT_DECLARE(FormFormat)
        /**
         * @brief inner class for some static initializations
         */
        class _init {
            public: _init();
        };
    private:
        static _init FormFormat_init;
    };
        
    /**
     * @brief class for FCFormat Output
     */
    class FCFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(FCFormat, print_dflt)
    public:
        FCFormat(ostream &os, unsigned opt=0);
        static void ncmul_print(const ncmul & p, const FCFormat & c, unsigned level=0);
        OUT_FORMAT_DECLARE(FCFormat)
        /**
         * @brief inner class for some static initializations
         */
        class _init {
            public: _init();
        };
    private:
        static _init FCFormat_init;
    };
    extern FCFormat FCout;
    
    /**
     * @brief class for index object
     */
    class Index : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Index, basic)
    public:
        enum Type {VD, CF, CA};
        Index(const string &s, const Type type=Type::VD);
        Pair operator() (const Index & i);
        Pair operator() (const Vector & p);
        Symbol name;
        Type type;
        void print(const print_context &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static bool hasc(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(Index);
    
    /**
     * @brief class for vector object
     */
    class Vector : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Vector, basic)
    public:
        Vector(const string &s);
        Pair operator() (const Vector & p);
        Pair operator() (const Index & mu);
        Symbol name;
        void print(const print_context &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(Vector);
    
    /**
     * @brief class for SUNT object
     */
    class SUNT : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(SUNT, basic)
    public:
        SUNT(ex i, ex j, ex a);
        ex ija[3]; // Index
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void print(const print_dflt &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        ex conjugate() const override;
    };
    GINAC_DECLARE_UNARCHIVER(SUNT);
    
    /**
     * @brief class for SUNF object
     */
    class SUNF : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(SUNF, basic)
    public:
        SUNF(ex i, ex j, ex k);
        ex ijk[3]; // Index
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        ex eval() const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(SUNF);
    
    /**
     * @brief class for SUNF4 object
     */
    class SUNF4 : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(SUNF4, basic)
    public:
        SUNF4(ex i, ex j, ex k, ex l);
        ex ijkl[4]; // Index
        size_t nops() const override;
        ex op(size_t i) const override;
        ex& let_op(size_t i) override;
        ex eval() const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(SUNF4);
    
    /**
     * @brief class for Pair object
     */
    class Pair : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Pair, basic)
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
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
    private:
        ex lr[2];
    };
    GINAC_DECLARE_UNARCHIVER(Pair);
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
     */
    class Eps : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(Eps, basic)
    public:
        ex pis[4];
        Eps(const Vector &p1, const Vector &p2, const Vector &p3, const Vector &p4);
        Eps(const Vector &p1, const Vector &p2, const Vector &p3, const Index &i1);
        Eps(const Vector &p1, const Vector &p2, const Index &i1, const Index &i2);
        Eps(const Vector &p1, const Index &i1, const Index &i2, const Index &i3);
        Eps(const Index &i1, const Index &i2, const Index &i3, const Index &i4);
        Eps(vector<Vector> vs, vector<Index> is);
        size_t nops() const override;
        ex op(size_t i) const override;
        ex & let_op(size_t i) override;
        ex eval() const override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void form_print(const FormFormat &c, unsigned level = 0) const;
        void fc_print(const FCFormat &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
    };
    GINAC_DECLARE_UNARCHIVER(Eps);
    ex LC(ex pi1, ex pi2, ex pi3, ex pi4);
    
    /**
     * @brief class for DiracGamma object
     */
    class DiracGamma : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(DiracGamma, basic)
    public:
        ex pi;
        unsigned rl;
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
        ex& let_op(size_t i) override;
        ex eval() const override;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        ex derivative(const symbol & s) const override;
        ex conjugate() const override;
    };
    GINAC_DECLARE_UNARCHIVER(DiracGamma);
    
    //-----------------------------------------------------------
    // TR/GAS functions
    //-----------------------------------------------------------
    DECLARE_FUNCTION_1P(TR)
    DECLARE_FUNCTION_1P(HF)
    
    inline ex GAS(const Vector &p, unsigned rl=0) { return DiracGamma(p,rl); }
    inline ex GAS(const Index &i, unsigned rl=0) { return DiracGamma(i,rl); }
    ex GAS(ex expr, unsigned rl=0);
    
    // Form, TIR, Apart
    ex charge_conjugate(const ex &);
    ex form(const ex &expr, bool all_in_one=false, int verb=0);
    ex TIR(const ex &expr_in, const lst &loop_ps, const lst &ext_ps);
    ex MatrixContract(const ex & expr_in);
    ex Apart(const ex &expr_in, const lst &vars, exmap sign_map=exmap());
    ex Apart(const ex &expr_in, const lst &loops, const lst & extmoms);
    ex ApartIR2ex(const ex & expr_in);
    ex ApartIR2F(const ex & expr_in);
    ex F2ex(const ex & expr_in);
    ex ApartIRC(const ex & expr_in, const ex & cut_props=lst{});
    void ApartIBP(int IBPmethod, exvector &air_vec, const lst & loops_exts=lst{}, const lst & cut_props=lst{}, std::function<lst(const Base &, const ex &)> uf=IBP::LoopUF);
    
    /**
     * @brief ApartIR function with 1 argument
     */
    class ApartIR1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function ApartIR(const T1 & p1) {
        return GiNaC::function(ApartIR1_SERIAL::serial, ex(p1));
    }
    
    /**
     * @brief ApartIR function with 2 arguments
     */
    class ApartIR2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function ApartIR(const T1 & p1, const T2 & p2) {
        return GiNaC::function(ApartIR2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    /**
     * @brief interface to qgraf program
     */
    namespace Qgraf {
        
        //-----------------------------------------------------------
        // Filed/Propagator/Vertex Function
        //-----------------------------------------------------------
        DECLARE_FUNCTION_3P(Propagator)
        DECLARE_FUNCTION_3P(InField)
        DECLARE_FUNCTION_3P(OutField)
        DECLARE_FUNCTION_3P(Matrix)
        
        /**
         * @brief Field function with 2 arguments
         */
        class Field2_SERIAL { public: static unsigned serial; };
        template<typename T1, typename T2>
        inline GiNaC::function Field(const T1 & p1, const T2 & p2) {
            return GiNaC::function(Field2_SERIAL::serial, ex(p1), ex(p2));
        }
        
        /**
         * @brief Field function with 3 arguments
         */
        class Field3_SERIAL { public: static unsigned serial; };
        template<typename T1, typename T2, typename T3>
        inline GiNaC::function Field(const T1 & p1, const T2 & p2, const T3 & p3) {
            return GiNaC::function(Field3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
        }
        
        /**
         * @brief Vertex function with 2 arguments
         */
        class Vertex2_SERIAL { public: static unsigned serial; };
        template<typename T1, typename T2>
        inline GiNaC::function Vertex(const T1 & p1, const T2 & p2) {
            return GiNaC::function(Vertex2_SERIAL::serial, ex(p1), ex(p2));
        }
        
        /**
         * @brief Vertex function with 3 arguments
         */
        class Vertex3_SERIAL { public: static unsigned serial; };
        template<typename T1, typename T2, typename T3>
        inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3) {
            return GiNaC::function(Vertex3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
        }
        
        /**
         * @brief Vertex function with 4 arguments
         */
        class Vertex4_SERIAL { public: static unsigned serial; };
        template<typename T1, typename T2, typename T3, typename T4>
        inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4) {
            return GiNaC::function(Vertex4_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4));
        }
        
        /**
         * @brief Vertex function with 5 arguments
         */
        class Vertex5_SERIAL { public: static unsigned serial; };
        template<typename T1, typename T2, typename T3, typename T4, typename T5>
        inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5) {
            return GiNaC::function(Vertex5_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5));
        }
        
        /**
         * @brief Vertex function with 6 arguments
         */
        class Vertex6_SERIAL { public: static unsigned serial; };
        template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
        inline GiNaC::function Vertex(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5, const T6 & p6) {
            return GiNaC::function(Vertex6_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5), ex(p6));
        }
        
        /**
         * @brief main interface to qgraf program
         */
        class Process {
        public:
            string Model;
            string In;
            string Out;
            string LoopPrefix = "q";
            int Loops;
            string Options;
            vector<string> Others;
            lst Amplitudes(symtab st, bool debug=false);
            /**
             * @brief inner class for some initialization
             */
            class _init {
                public: _init();
            };
        private:
            static _init Process_init;
        };
        
        extern const string Style;
        lst TopoLines(const ex & amp);
        void DrawPDF(const lst & amps, string fn, bool debug=false);
        vector<lst> ShrinkCut(ex amp, lst prop, int n=1);
        bool HasLoop(ex amp, lst prop);
        extern map<ex,string,ex_is_less> LineTeX;
        extern map<ex,string,ex_is_less> VerTeX;
        extern map<ex,string,ex_is_less> InOutTeX;
        
        Index LI(ex fn);
        Index DI(ex fn);
        Index TI(ex fn);
        Index FI(ex fn);
        Index CI(ex fn);
        Index AI(ex fn);
        Index RLI(ex fn);
        Index RDI(ex fn);
        Index RTI(ex fn);
        Index RFI(ex fn);
        Index RCI(ex fn);
        Index RAI(ex fn);
        
        ex QuarkPropagator(ex e, ex m=0);
        ex GluonPropagator(ex e);
        ex GhostPropagator(ex e);
        ex q2gVertex(ex e);
        ex g3Vertex(ex e);
        ex g4Vertex(ex e);
        ex gh2gVertex(ex e);
        
        ex IndexCC(ex e, bool all=true);
        ex GluonFFV(ex e, ex n);
        ex QuarkFFV(ex e, ex n);
        
        ex eikonalPropagator(ex e, ex n, int mode); // 0 for gluon, others for quark/anti-quark
        ex eikonalPropagatorR(ex e, ex n, int mode); // right side from cut
        ex eikonalVertex(ex e, ex n, int mode); // 0 for gluon, 1 for quark, 2 for anti-quark, in<0 & out>0
        ex eikonalVertexR(ex e, ex n, int mode); // right side from cut
            
        ex GluonSum(int qi);
        ex QuarkSum(int qi, ex p, ex m);
        ex AntiQuarkSum(int qi, ex p, ex m);
        ex GhostSum(int qi);
        ex AntiGhostSum(int qi);
        ex J1Sum(int qi, ex p);
    };
    
    /**
     * @brief namespace for functions helpful in Heavy Quarkonium
     */
    namespace Quarkonium {
        enum IO {In, Out};
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu);
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu);
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu, int i, int j);
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu, int i, int j);
        ex ColorProj(int i, int j, Index a);
        ex ColorProj(int i, int j);
        
        ex S1L1Proj(ex si, ex qi, ex p);
        ex S1L1Proj(ex si, ex qi, ex mu, ex p);
        ex S1L1Proj(ex si, ex qi, ex mu1, ex mu2, ex p);
        ex S1L2Proj(ex si, ex qi1, ex qi2, ex mu, ex p);
        ex S1L2Proj(ex si, ex qi1, ex qi2, ex mu1, ex mu2, ex p);
        ex S1L1Sum(ex si, ex siR, ex qi, ex qiR, ex p, int J);
        
        ex LProj(const ex &expr_in, const lst &pqi, string prefix="lpj");
        
        ex Gamma5(const string pre, int start=1);
        
        ex DoPS(lst moms, ex amp, int si=-1, ex q2=1);
        ex nPS(int n, ex q2=1);
    }
    
    /**
     * @brief namespace for Renormalization Constant
     */
    namespace RC {
        ex Z2(string name, ex m, int loop=2);
        ex Z2(Symbol n, ex m, int loop=2);
        ex Zm(ex m, int loop=2);
        ex asBare(int loop=2);
        ex asLO();
    }
        
}

