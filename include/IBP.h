/**
 * @file 
 * @brief IBP header file
 */
 
#pragma once

#include "Basic.h"

namespace HepLib::IBP {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    extern const Symbol d;
    
    DECLARE_FUNCTION_1P(a)

    class FIRE {
    public:
        lst Internal;
        lst External;
        lst Replacements;
        lst Propagators;
        lst Integrals;
        lst Cuts;
        int pos_pref = 1;
        lst mi_pref;
        
        string WorkingDir;
        int ProblemNumber;
        
        int ProblemDimension;
        lst MasterIntegrals;
        lst Rules;
        lst Pairs;
        
        void Export();
        void Run();
        void Import();
        void Reduce();
        
        static lst SortPermutation(const ex & in_expr, const lst & xs);
        static lst LoopUF(const FIRE & fire, const ex & corner);
        static pair<exmap,lst> FindRules(vector<FIRE> &fs, bool mi=true, std::function<lst(const FIRE &, const ex &)> uf=LoopUF);
        static pair<exmap,lst> FindRules(vector<FIRE*> &fs, bool mi=true, std::function<lst(const FIRE &, const ex &)> uf=LoopUF); 
        
        static int Version;
        static ex VectorDimension;
        
        static lst UF(const ex & ps, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs); 
    private:
        vector<exmap> IBPs;
        
    };
            
    class Laporta {
    private:
        map<ex,lst,ex_is_less> ccF;
        map<ex,lst,ex_is_less> ccI;
        map<lst,ex,ex_is_less> ccFinv;
        map<ex,lst,ex_is_less> ccImax;
    public:
        
        lst F2lst(ex f);
        lst I2lst(ex ibp);
        
        lst Internal;
        lst External;
        lst Replacements;
        lst Propagators;
        lst Integrals;
        lst Cuts;
        
        int ProblemDimension;
        lst MasterIntegrals;
        lst Rules;
        lst Pairs;
                
        void Prepare(lst loop, lst ext, lst prop, lst repl);
        void Generate(vector<lst> seeds);
        void Generate2(lst seed);
        void Reduce();
        
        lst preIBPs;
        lst IBPs;
        static ex collectF(ex);
    };


}
