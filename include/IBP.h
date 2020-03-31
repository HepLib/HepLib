#pragma once

#include "Basic.h"

namespace HepLib::IBP {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    extern const Symbol d;
    
    DECLARE_FUNCTION_1P(a)
    DECLARE_FUNCTION_2P(F)

    class FIRE {
    public:
        lst Internal;
        lst External;
        lst Variables;
        lst Replacements;
        lst Propagators;
        exvector Integrals;
        vector<exmap> IBPs;
        ex UF(ex corner);
        
        string WorkingDir;
        int ProblemNumber;
        
        int Dimension;
        exvector MasterIntegrals;
        exmap Rules;
        
        void Reduce();
        
        static exmap FindRules(vector<FIRE> fs, bool mi=true); 
        
    };


}
