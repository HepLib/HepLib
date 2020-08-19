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
    
    class Base {
    public:
        lst Internal;
        lst External;
        lst Replacements;
        lst Propagators;
        lst Integrals;
        lst Cuts;
        
        string WorkingDir;
        int ProblemNumber = 0;
        ex VectorDimension = d;
        
        int ProblemDimension;
        lst MasterIntegrals;
        lst Rules;
        lst Pairs;
        
        virtual void Export() = 0;
        virtual void Run() = 0;
        virtual void Import() = 0;
        void Reduce();
        
    };

    class FIRE : public Base {
    public:
    
        void Export() override;
        void Run() override;
        void Import() override;
    
        int pos_pref = 1;
        lst mi_pref;
        
        static int Version;
        static int Threads;
        
    };
            
    class KIRA : public Base {
    public:
        
        int SeedRuns = 2;
        
        void Export() override;
        void Run() override;
        void Import() override;
        
    };
    
    
    lst SortPermutation(const ex & in_expr, const lst & xs);
    lst LoopUF(const Base & fire, const ex & corner);
    pair<exmap,lst> FindRules(vector<Base*> &fs, bool mi=true, std::function<lst(const Base &, const ex &)> uf=LoopUF); 
    lst UF(const ex & ps, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs);


}
