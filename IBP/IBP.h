/**
 * @file 
 * @brief IBP header file
 */
 
#pragma once

#include "BASIC.h"

/**
 * @brief namespace for IBP reduction
 */
namespace HepLib::IBP {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    extern const Symbol d;
    DECLARE_FUNCTION_1P(a)
    
    /**
     * @brief Base class for IBP reduction
     */
    class Base {
    public:
        lst Internal;
        lst External;
        lst Replacements;
        lst Propagators;
        lst Integrals; // lst of index lst
        lst Cuts; // index start from 1
        lst DSP; // { {q1,q1}, {q1,p}, ... } Diff SP
        lst ISP; // { q1*q1, q1*p } Independent SP
        
        string WorkingDir;
        int ProblemNumber = 0;
        lst PIntegrals;
        lst MIntegrals;
        lst Rules;
        bool reCut = true;
        
        virtual void Export() { throw Error("Export() not implemented!"); };
        virtual void Run() { throw Error("Run() not implemented!"); };
        virtual void Import() { throw Error("Import() not implemented!"); };
        void Reduce();
        
    };

    /**
     * @brief IBP reduction using FIRE program
     */
    class FIRE : public Base {
    public:
    
        void Export() override;
        void Run() override;
        void Import() override;
        int pos_pref = 1;
        static int Version;
        static int Threads;
    };
    
    /**
     * @brief IBP reduction using KIRA program
     */
    class KIRA : public Base {
    public:
    
        static string KArgs;
        
        int ra = 2;
        int sa = 3;
        int dmax = -1;
        
        void Export() override;
        void Run() override;
        void Import() override;
        
    private:
        string Fout(const ex & expr);
        ex Fin(const string & expr);
        map<ex,unsigned long long,ex_is_less> i2w;
        map<unsigned long long,ex> w2i;
    };

    
    /**
     * @brief IBP reduction using KIRA program with user-defined equations
     */
    class UKIRA : public Base {
    public:
    
        static int Rounds;
        static string KArgs; // check kira --help
        
        bool using_uw = true;
        map<int,ex> Shift;
        
        void Export() override;
        void Run() override;
        void Import() override;
        
        int ra = 2;
        int sa = 2;
        int rap = 1;
        int sap = 1;
        int sort_option = 0;
        int seed_option = 0;
        
    private:
        lst ibps;
        
        string Fout(const ex & expr);
        ex Fin(const string & expr);
        map<ex,unsigned long long,ex_is_less> i2w;
        map<unsigned long long,ex> w2i;
    };
    
    class Laporta : public Base {
    
    public:

        static int Rounds;
        
        bool using_uw = true;
        map<int,ex> Shift;
        
        void Export() override;
        void Run() override;
        void Import() override;
        
        int ra = 1;
        int sa = 1;
        int rap = 1;
        int sap = 1;
        int sort_option = 0;
        int seed_option = 0;
        
    private:
        lst ibps;
        int Round = 0;
        lst _Integrals;
        lst _Rules;
        lst RIntegrals;
        
        string Fout(const ex & expr);
        ex Fin(const string & expr);
        map<ex,long long,ex_is_less> i2w;
        map<long long,ex> w2i;
        
     };
    
    lst SortPermutation(const ex & in_expr, const lst & xs);
    lst LoopUF(const Base & fire, const ex & corner);
    lst UF(const ex & ps, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs);
    pair<exmap,lst> FindRules(vector<Base*> fs, bool mi=true, std::function<lst(const Base &, const ex &)> uf=LoopUF);
    inline pair<exmap,lst> FindRules(Base& ibp, bool mi=true, std::function<lst(const Base &, const ex &)> uf=LoopUF) {
        vector<Base*> fs;
        fs.push_back(&ibp);
        return FindRules(fs, mi, uf);
    }


}
