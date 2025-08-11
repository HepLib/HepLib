/**
 * @file 
 * @brief QGRAF header file
 */
 
#pragma once

#include "HEP.h"

/**
 * @brief namespace for generating Feynman diagrams or amplitudes.
 */
namespace HepLib::QGRAF {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    //-----------------------------------------------------------
    // Filed/Propagator/Vertex Function
    //-----------------------------------------------------------
    DECLARE_FUNCTION_3P(Propagator)
    DECLARE_FUNCTION_3P(InField)
    DECLARE_FUNCTION_3P(OutField)
    
    #ifndef DOXYGEN_SKIP
    
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
    
    #endif
    
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
        lst Amplitudes(symtab st);
        
        #ifndef DOXYGEN_SKIP
        class _init {
            public: _init();
        };
    private:
        static _init Process_init;
        #endif
    };
    
    lst TopoLines(const ex & amp);
    void DrawPDF(const lst & amps, string fn, int nr=4, bool single_page=false);
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
    
    // Propagator @ SM https://arxiv.org/abs/1209.6213v2
    ex QuarkPropagator(const ex & e, const ex & m=0);
    ex LeptonPropagator(const ex & e, const ex & m=0);
    ex GluonPropagator(const ex & e, const ex & xi=1);
    ex GluonGhostPropagator(const ex & e, const ex & xi=1, const ex & eta_G=1);
    ex AZWPropagator(const ex & e, const ex & m=0, const ex & xi=1);
    ex AZWGhostPropagator(const ex & e, const ex & m=0, const ex & xi=1, const ex & eta_G=1);
    ex ScalarPropagator(const ex & e, const ex & m, const ex & xi=1);
    
    // Vertex @ SM https://arxiv.org/abs/1209.6213v2
    ex QuarkGluonVertex(const ex & e, const ex & eta_s=-1);
    ex Gluon3Vertex(const ex & e, const ex & eta_s=-1);
    ex Gluon4Vertex(const ex & e, const ex & eta_s=-1);
    ex GhostGluonVertex(const ex & e, const ex & eta_s=-1, const ex & eta_G=1);
    ex QuarkAVertex(const ex & e, const ex & eq, const ex & eta_e=-1);
    ex LeptonAVertex(const ex & e, const ex & eq, const ex & eta_e=-1);
    
    ex IndexL2R(ex e, bool all=true);
    ex IndexCC(ex e, bool all=true);
    
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
    
    class Models {
    public:
        static string GluonFF;
        static string QuarkFF;
        static string SM;
        static string QCD;
        
        static ex FeynRulesSM(const ex & amp, const ex & xi=1);
        static lst FeynRulesSM(const lst & amps, const ex & xi=1);
        static ex FeynRulesQCD(const ex & amp, const ex & xi=1);
        static lst FeynRulesQCD(const lst & amps, const ex & xi=1);
    };
    
}

