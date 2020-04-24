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
    DECLARE_FUNCTION_2P(F)

    class FIRE {
    public:
        lst Internal;
        lst External;
        lst Replacements;
        lst Propagators;
        lst Integrals;
        lst Cuts;
        ex UF(const ex & corner) const;        
        ex VectorDimension = d;
        
        string WorkingDir;
        int ProblemNumber;
        
        int ProblemDimension;
        lst MasterIntegrals;
        lst Rules;
        
        void Reduce();
        
        static pair<exmap,lst> FindRules(vector<FIRE> &fs, bool mi=true);
        static pair<exmap,lst> FindRules(vector<FIRE*> &fs, bool mi=true); 
        static int Version;
    private:
        vector<exmap> IBPs;
        
    };


}
