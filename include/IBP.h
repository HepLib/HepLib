/**
 * @file 
 * @brief IBP header file
 */
 
#pragma once

#include "Basic.h"

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
        
        string WorkingDir;
        int ProblemNumber = 0;
        
        lst MasterIntegrals;
        lst Rules;
        lst Pairs;
        
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
    
        int ProblemDimension;
        int pos_pref = 1;
        lst mi_pref;
        
        static int Version;
        static int Threads;
    };
    
    /**
     * @brief IBP reduction using KIRA program
     */
    class KIRA : public Base {
    public:
    
        string cmd_args = "";
        lst mi_pref;
        
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
        
        bool use_weight = true;
        string cmd_args = "";
        map<int,ex> Shift;
        lst mi_pref;
        
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
        int Round = 0;
        lst _Integrals;
        lst _Rules;
        lst RIntegrals;
        
        string Fout(const ex & expr);
        ex Fin(const string & expr);
        map<ex,unsigned long long,ex_is_less> i2w;
        map<unsigned long long,ex> w2i;
    };
    
    class Laporta : public Base {
    
    private:
         map<ex,lst,ex_is_less> ccF;
         map<ex,lst,ex_is_less> ccI;
         map<lst,ex,ex_is_less> ccFinv;
         map<ex,lst,ex_is_less> ccImax;

    public:

        lst F2lst(ex f);
        lst I2lst(ex ibp);

        void Prepare(lst loop, lst ext, lst prop, lst repl);
        void Generate(vector<lst> seeds);
        void Generate2(lst seed);
        void Reduce();

        lst preIBPs;
        lst IBPs;
        static ex collectF(ex);
     };
    
    lst SortPermutation(const ex & in_expr, const lst & xs);
    lst LoopUF(const Base & fire, const ex & corner);
    pair<exmap,lst> FindRules(vector<Base*> &fs, bool mi=true, std::function<lst(const Base &, const ex &)> uf=LoopUF); 
    lst UF(const ex & ps, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs);


}
