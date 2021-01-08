/**
 * @file 
 * @brief QGRAF header file
 */
 
#pragma once

#include "HEP.h"

/**
 * @brief namespace for Index, Vector, DGamma, etc.
 */
namespace HepLib::Qgraf {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    //-----------------------------------------------------------
    // Filed/Propagator/Vertex Function
    //-----------------------------------------------------------
    DECLARE_FUNCTION_3P(Propagator)
    DECLARE_FUNCTION_3P(InField)
    DECLARE_FUNCTION_3P(OutField)
    
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
        static string Style;
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
    
    ex QuarkPropagator(ex e, ex m=0, bool color=true);
    ex GluonPropagator(ex e, bool color=true);
    ex GhostPropagator(ex e, bool color=true);
    ex q2gVertex(ex e, bool color=true);
    ex g3Vertex(ex e);
    ex g4Vertex(ex e);
    ex gh2gVertex(ex e, bool color=true);
    
    ex IndexL2R(ex e, bool all=true);
    ex IndexCC(ex e, bool all=true);
    ex GluonFFV(ex e, ex n);
    ex QuarkFFV(ex e, ex n);
    
    ex eikonalPropagator(ex e, ex n, int mode); // 0 for gluon, others for quark/anti-quark
    ex eikonalPropagatorR(ex e, ex n, int mode); // right side from cut
    ex eikonalVertex(ex e, ex n, int mode); // 0 for gluon, 1 for quark, 2 for anti-quark, in<0 & out>0
    ex eikonalVertexR(ex e, ex n, int mode); // right side from cut
        
    ex GluonSumL(int qi, bool color=true);
    ex QuarkSumL(int qi, ex p, ex m, bool color=true);
    ex AntiQuarkSumL(int qi, ex p, ex m, bool color=true);
    ex GhostSumL(int qi);
    ex AntiGhostSumL(int qi);
    ex J1SumL(int qi, ex p);
    
    ex GluonSumR(int qi, bool color=true);
    ex QuarkSumR(int qi, ex p, ex m, bool color=true);
    ex AntiQuarkSumR(int qi, ex p, ex m, bool color=true);
    ex GhostSumR(int qi);
    ex AntiGhostSumR(int qi);
    ex J1SumL(int qi, ex p);
    
    inline ex GluonSum(int qi, bool color=true) { return GluonSumL(qi,color); };
    inline ex QuarkSum(int qi, ex p, ex m, bool color=true) { return QuarkSumL(qi,p,m,color); };
    inline ex AntiQuarkSum(int qi, ex p, ex m, bool color=true) { return AntiQuarkSumL(qi,p,m,color); };
    inline ex GhostSum(int qi) { return GhostSumL(qi); };
    inline ex AntiGhostSum(int qi) { return AntiGhostSumL(qi); };
    inline ex J1Sum(int qi, ex p) { return J1SumL(qi,p); };
        
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
        ex Zas(int loop);
    }
        
}

