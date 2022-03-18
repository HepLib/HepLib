/**
 * @file 
 * @brief IBP header file
 */
 
#pragma once

#include "BASIC.h"

/**
 * @brief namespace for IBP reduction
 */
namespace HepLib {

    using namespace std;
    using namespace GiNaC;
    using namespace HepLib;
    
    DECLARE_FUNCTION_1P(a)
    
    /**
     * @brief IBP base class for IBP reduction
     */
    class IBP {
    public:
        lst Internal;
        lst External;
        lst Replacements;
        lst Propagators;
        lst Integrals; // lst of index lst
        lst Cuts; // index start from 1
        lst DSP; // { {q1,q1}, {q1,p}, ... } Diff SP
        lst ISP; // { q1*q1, q1*p } Independent SP
        map<int,ex> Shift; // index start from 1
        bool reCut = false;
        string WorkingDir;
        int ProblemNumber = 0;
        lst PIntegrals;
        
        lst MIntegrals;
        lst Rules;
        bool IsAlwaysZero = false;
        
        virtual void Export() { throw Error("Export() not implemented!"); };
        virtual void Run() { throw Error("Run() not implemented!"); };
        virtual void Import() { throw Error("Import() not implemented!"); };
        pair<exmap,lst> FindRules(bool mi=true);
        bool IsZero(ex sector);
        void Reduce();
        void Export(string garfn); // Export to .gar
        void Import(string garfn); // Import from .gar
        ex TO(); // to single list for output
        void FROM(ex s); // from a single expression 
        exmap SP2Pn();
        exmap Dinv(const lst & ns);
        ex D(const ex & x, const lst &ns);
        void RM(bool keep_start_config=false);
        
        static void ReShare(const vector<IBP*> & fs);
    };

    /**
     * @brief IBP reduction using FIRE program
     */
    class FIRE : public IBP {
    public:
        void Export() override;
        void Run() override;
        void Import() override;
        int pos_pref = 2;
        static int Version;
        static int Threads;
        static int fThreads;
        static int lThreads;
        static int sThreads;
        static void RRTables(const string & filename, int pnum);
        static void ThieleTables(const string & filename, int si, int ei);
    };
    
    /**
     * @brief IBP reduction using KIRA program
     */
    class KIRA : public IBP {
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
    class UKIRA : public IBP {
    public:
    
        static string KArgs; // check kira --help
        
        bool using_uw = true;
        
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
        
        string Fout(const ex & expr);
        ex Fin(const string & expr);
        map<ex,unsigned long long,ex_is_less> i2w;
        map<unsigned long long,ex> w2i;
    };
    
    /**
     * @brief IBP reduction using Direc method
     */
    class Direct : public IBP {
    public:
        class Condition {
        public:
            vector<pair<int,int>> cs; // (1st==-1 for <=, 1st==1 for >=, 1st==0 for ==) 2nd
            inline bool IsOK(ex ns) {
                if(ns.nops() != cs.size()) return false;
                for(int i=0; cs.size()-i>0; i++) {
                    if(cs[i].first==-1 && ns.op(i)>cs[i].second) return false;
                    else if(cs[i].first==1 && ns.op(i)<cs[i].second) return false;
                    else if(cs[i].first==0 && !ns.op(i).is_equal(cs[i].second)) return false;
                }
                return true;
            }
            inline ex cs2ex() {
                lst ret;
                for(auto item : cs) ret.append(lst{item.first, item.second});
                return ret;
            }
            inline void ex2cs(ex e) {
                cs.clear();
                for(auto item : e) {
                    cs.push_back(make_pair(ex2int(item.op(0)), ex2int(item.op(1))));
                }
            }
        };
        
        void Export() override;
        void Run() override;
        void Import() override;
        
    private:
        lst ibps;
        vector<pair<Condition,ex>> ConSolVec; // the conditional solution vector
    };
    
    class Laporta : public IBP {
    
    public:
        
        bool using_uw = true;
        
        void Export() override;
        void Run() override;
        void Import() override;
        
        int ra = 1;
        int sa = 1;
        int rap = 1;
        int sap = 1;
        int sort_option = 0;
        int seed_option = 1;
        
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
    
    extern exmap MapPreSP;
    exmap SortPermutation(const ex & in_expr, const lst & xs);
    lst LoopUF(const IBP & fire, const ex & corner);
    lst UF(const ex & ps, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs);
    pair<exmap,lst> FindRules(vector<IBP*> fs, bool mi=true, std::function<lst(const IBP &, const ex &)> uf=LoopUF);
    inline pair<exmap,lst> FindRules(IBP& ibp, bool mi=true, std::function<lst(const IBP &, const ex &)> uf=LoopUF) {
        vector<IBP*> fs;
        fs.push_back(&ibp);
        return FindRules(fs, mi, uf);
    }
    
    ex GPolynomial(const IBP & IBP);
    void GPermutation(const ex & uf, const lst & xs);
}
