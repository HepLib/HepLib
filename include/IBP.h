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
        lst Replacements;
        lst Propagators;
        lst Integrals;
        ex UF(ex corner);
        
        string WorkingDir;
        int ProblemNumber;
        int Version = 6;
        
        int Dimension;
        lst MasterIntegrals;
        lst Rules;
        
        void Reduce();
        
        static lst FindRules(vector<FIRE> &fs, bool mi=true);
        static lst FindRules(vector<FIRE*> &fs, bool mi=true); 
    private:
        vector<exmap> IBPs;
        
    };


}
