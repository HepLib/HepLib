/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once

#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <chrono>
#include <fstream>
#include <tuple>
#include <list>
#include <omp.h>
#include <sys/wait.h>

#define DEFAULT_CTOR(classname) \
classname::classname() { setflag(status_flags::evaluated | status_flags::expanded); }

#define IMPLEMENT_HAS(classname) \
bool classname::has(const ex &e) { \
    for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) if(is_a<classname>(*i)) return true; \
    return false; \
}

#define IMPLEMENT_ALL(classname) \
lst classname::all(const ex &e) { \
    lst ret;\
    for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) if(is_a<classname>(*i)) ret.append(*i); \
    ret.sort(); \
    ret.unique(); \
    return ret; \
}

/*-----------------------------------------------------*/
// operator << Macro for GiNaC Output Format
/*-----------------------------------------------------*/
#define OUT_FORMAT_DECLARE(classname) \
    template<class T> const classname & operator << (const T & v) const { \
        s << v; \
        return *this; \
    }; \
    const classname & operator << (const basic & v) const; \
    const classname & operator << (const ex & v) const; \
    const classname & operator<<(std::ostream& (*v)(std::ostream&)) const;

#define OUT_FORMAT_IMPLEMENT(classname) \
    const classname & classname::operator << (const basic & v) const { \
        v.print(*this); \
        return *this; \
    } \
    const classname & classname::operator << (const ex & v) const { \
        v.print(*this); \
        return *this; \
    } \
    const classname & classname::operator<<(std::ostream& (*v)(std::ostream&)) const { \
        s << v; \
        return *this; \
    }

namespace HepLib {

    using namespace std;
    using namespace GiNaC;
    
    /*-----------------------------------------------------*/
    // Terminal Color
    /*-----------------------------------------------------*/
    #define RESET   "\033[0m"
    #define BLACK   "\033[30m"      /* Black */
    #define RED     "\033[31m"      /* Red */
    #define GREEN   "\033[32m"      /* Green */
    #define YELLOW  "\033[33m"      /* Yellow */
    #define BLUE    "\033[34m"      /* Blue */
    #define MAGENTA "\033[35m"      /* Magenta */
    #define CYAN    "\033[36m"      /* Cyan */
    #define WHITE   "\033[37m"      /* White */
    #define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
    #define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
    #define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
    #define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
    #define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
    #define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
    #define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
    #define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
    
    //-----------------------------------------------------------
    // Symbol Class
    //-----------------------------------------------------------
    class Symbol : public symbol {
    GINAC_DECLARE_REGISTERED_CLASS(Symbol, symbol)
    public:
        // is_real=false for pure imaginary
        // check=true to if the name exist then error
        Symbol(const string &s, bool is_real=true, bool check=false);
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        
        ex eval() const override; // for performance reasons
        ex evalf() const override; // for performance reasons
        ex conjugate() const override;
        ex real_part() const override;
        ex imag_part() const override;
        bool isReal = true;
        
        static bool has(const ex &e);
        static lst all(const ex &e);
        static std::map<std::string, ex> Table;
        static exmap AssignMap;
        static void Assign(const Symbol & s, const ex & v);
        static void Assign(const string & str, const ex & v);
        static void clearAssign(const Symbol &s);
        static void clearAssign(const string &str);
        static void clearAssign();
    };
    GINAC_DECLARE_UNARCHIVER(Symbol);
    
    
    /*-----------------------------------------------------*/
    // Eroor Class
    /*-----------------------------------------------------*/
    class Error : public exception {
    public:
        string msg;
        const char * what() const throw ();
        Error(string _msg);
    };

    /*-----------------------------------------------------*/
    // Global Symbol
    /*-----------------------------------------------------*/
    const symbol & get_symbol(const string & s, bool check=false);

    /*-----------------------------------------------------*/
    // split
    /*-----------------------------------------------------*/
    vector<std::string> split(const string& s, char delimiter);

    /*-----------------------------------------------------*/
    // Helper Classes
    /*-----------------------------------------------------*/
    class MatHelper {
    public:
        static bool has_zero_row(const matrix &mat);
        static bool is_zero_row(const matrix &mat, int r);
        static vector<int> zero_row_index(const matrix &mat);
        static matrix remove_zero_rows(const matrix &mat);
        static matrix sub(matrix mat, int r, int nr, int c, int nc);
        template <typename F> static void map_inplace(matrix &m, F f);
        template <typename F> static matrix map(const matrix &m, F f);
    };

    template <typename F>
    void MatHelper::map_inplace(matrix &m, F f) {
        for (unsigned i = 0; i < m.nops(); i++) {
            m.let_op(i) = f(m.op(i));
        }
    }

    template <typename F>
    matrix MatHelper::map(const matrix &m, F f) {
        matrix r(m.rows(), m.cols());
        for (unsigned i = 0; i < m.nops(); i++) {
            r.let_op(i) = f(m.op(i));
        }
        return r;
    }

    class lstHelper {
    public:
        template <typename F> static void map_inplace(lst &m, F f);
        template <typename F> static lst map(const lst &m, F f);
        static ex sum(lst m);
    };

    template <typename F>
    void lstHelper::map_inplace(lst &m, F f) {
        for (unsigned i = 0; i < m.nops(); i++) {
            m.let_op(i) = f(m.op(i));
        }
    }

    template <typename F>
    lst lstHelper::map(const lst &m, F f) {
        lst r = m;
        for (unsigned i = 0; i < m.nops(); i++) {
            r.let_op(i) = f(m.op(i));
        }
        return r;
    }

    /*-----------------------------------------------------*/
    // Global Functions
    /*-----------------------------------------------------*/
    string now(bool use_date = true);
    lst gather_symbols(const ex & e);
    lst gather_symbols(const vector<ex> & ve);

    inline bool file_exists(const char* fn) {
        return (access(fn,F_OK)!=-1);
    }
    
    inline ex subs_naive(const ex & expr, const ex & repls) {
        return subs(expr, repls, subs_options::no_pattern);
    };
    inline ex subs_naive(const ex & expr, const exmap & repls) {
        return subs(expr, repls, subs_options::no_pattern);
    };
    
    inline ex subs_all(const ex & expr, const ex & repls) {
        return subs(expr, repls, subs_options::algebraic);
    };
    inline ex subs_all(const ex & expr, const exmap & repls) {
        return subs(expr, repls, subs_options::algebraic);
    };
    
    /*-----------------------------------------------------*/
    // vector : GiNaC_Parallel
    /*-----------------------------------------------------*/
    extern lst GiNaC_archive_Symbols;
    void GiNaC_archive_Symbols_from(ex);
    void GiNaC_archive_Symbols_from(vector<ex>);
    vector<ex> GiNaC_Parallel(
        int ntotal, int nbatch,
        std::function<ex(int)> f,
        const string & key = "",
        bool rm = true,
        int prtlvl = 0
    );
    inline vector<ex> GiNaC_Parallel(
        int ntotal,
        std::function<ex(int)> f,
        const string & key = "",
        bool rm = true,
        int prtlvl = 0
    ) { return  GiNaC_Parallel(ntotal, 0, f, key, rm, prtlvl); }
    
    /*-----------------------------------------------------*/
    // Helpers
    /*-----------------------------------------------------*/
    string RunOS(const string & cmd);
    void garRead(const string &garfn, map<string, ex> &resMap);
    ex garRead(const string &garfn, const char* key);
    ex garRead(const string &garfn);
    void garWrite(const string &garfn, const map<string, ex> &resMap);
    void garWrite(const string &garfn, const ex & res);
    ex str2ex(const string &expr, symtab stab);
    ex str2ex(const string &expr);
    lst str2lst(const string &expr, symtab stab);
    lst str2lst(const string &expr);
    string file2str(string filename);
    ex file2ex(string filename);
    string ex2str(const ex &expr);
    lst exvec2lst(const exvector & exvec);
    exvector lst2exvec(const lst & alst);
    lst xlst(int ei);
    lst xlst(int bi, int ei);
    
    void let_op_append(ex & ex_in, const ex item);
    void let_op_prepend(ex & ex_in, const ex item);
    void let_op_remove_last(ex & ex_in);
    void let_op_remove_first(ex & ex_in);

    void let_op_append(ex & ex_in, int index, const ex item);
    void let_op_prepend(ex & ex_in, int index, const ex item);
    void let_op_remove_last(ex & ex_in, int index);
    void let_op_remove_first(ex & ex_in, int index);
    void let_op_append(lst & ex_in, int index, const ex item);
    void let_op_prepend(lst & ex_in, int index, const ex item);
    void let_op_remove_last(lst & ex_in, int index);
    void let_op_remove_first(lst & ex_in, int index);

    void let_op_append(ex & ex_in, int index1, int index2, const ex item);
    void let_op_prepend(ex & ex_in, int index1, int index2, const ex item);
    void let_op_remove_last(ex & ex_in, int index1, int index2);
    void let_op_remove_first(ex & ex_in, int index1, int index2);
    void let_op_append(lst & ex_in, int index1, int index2, const ex item);
    void let_op_prepend(lst & ex_in, int index1, int index2, const ex item);
    void let_op_remove_last(lst & ex_in, int index1, int index2);
    void let_op_remove_first(lst & ex_in, int index1, int index2);

    void let_op(ex &ex_in, int index1, int index2, const ex item);
    void let_op(lst &ex_in, int index1, int index2, const ex item);
    void let_op(ex &ex_in, int index1, int index2, int index3, const ex item);
    void let_op(lst &ex_in, int index1, int index2, int index3, const ex item);

    ex get_op(const ex ex_in, int index1, int index2);
    ex get_op(const lst ex_in, int index1, int index2);
    ex get_op(const ex ex_in, int index1, int index2, int index3);
    ex get_op(const lst ex_in, int index1, int index2, int index3);

    /*-----------------------------------------------------*/
    // Series at s=0 similar to Mathematica
    /*-----------------------------------------------------*/
    ex mma_series(ex const & expr, const symbol &s, int sn);
    
    ex mma_expand(const ex &expr, std::function<bool(const ex &)>, int depth=0);
    ex mma_expand(ex const &expr, lst const &pats, int depth=0);
    ex mma_expand(ex const &expr, ex const &pat, int depth=0);
    
    ex mma_collect(const ex &expr, std::function<bool(const ex &)>, bool ccf=false, bool cvf=false);
    ex mma_collect(const ex &expr, lst const &pats, bool ccf=false, bool cvf=false);
    ex mma_collect(const ex &expr, ex const &pat, bool ccf=false, bool cvf=false);
    
    ex mma_diff(ex const expr, ex const xp, unsigned nth=1, bool expand=false);
    

    /*-----------------------------------------------------*/
    // Evalf
    /*-----------------------------------------------------*/
    ex Evalf(ex);

    /*-----------------------------------------------------*/
    // xPositive
    /*-----------------------------------------------------*/
    bool xPositive(ex const expr);
    int xSign(ex const expr);

    /*-----------------------------------------------------*/
    // Global object wildcard
    /*-----------------------------------------------------*/
    extern ex w, w0, w1, w2, w3, w4, w5;
    extern string InstallPrefix;
    extern const Symbol iEpsilon;
    extern const Symbol ep;
    extern const Symbol D;
    extern int Verbose;
    extern int ParallelProcess;
    extern pid_t PID;
    
    /*-----------------------------------------------------*/
    // Global Colors
    /*-----------------------------------------------------*/
    extern const char* Color_Error;
    extern const char* Color_Warn;
    extern const char* Color_HighLight;

    /*-----------------------------------------------------*/
    // Customized GiNaC Function
    /*-----------------------------------------------------*/
    DECLARE_FUNCTION_1P(coCF)
    DECLARE_FUNCTION_1P(coVF)
    
    DECLARE_FUNCTION_1P(x)
    DECLARE_FUNCTION_1P(y)
    DECLARE_FUNCTION_1P(z)
    
    // F wrapper function upto 2 arguments
    class F1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function F(const T1 & p1) {
        return GiNaC::function(F1_SERIAL::serial, ex(p1));
    }
    
    class F2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function F(const T1 & p1, const T2 & p2) {
        return GiNaC::function(F2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    // WF wrapper function upto 5 arguments
    class WF1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function WF(const T1 & p1) {
        return GiNaC::function(WF1_SERIAL::serial, ex(p1));
    }
    
    class WF2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2) {
        return GiNaC::function(WF2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    class WF3_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2, const T3 & p3) {
        return GiNaC::function(WF3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
    }
    
    class WF4_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4) {
        return GiNaC::function(WF4_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4));
    }
    
    class WF5_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5) {
        return GiNaC::function(WF5_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5));
    }
    
    // iWF internal wrapper function upto 5 arguments
    class iWF1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function iWF(const T1 & p1) {
        return GiNaC::function(iWF1_SERIAL::serial, ex(p1));
    }
    
    class iWF2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2) {
        return GiNaC::function(iWF2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    class iWF3_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2, const T3 & p3) {
        return GiNaC::function(iWF3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
    }
    
    class iWF4_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4) {
        return GiNaC::function(iWF4_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4));
    }
    
    class iWF5_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5) {
        return GiNaC::function(iWF5_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5));
    }
    
    /*-----------------------------------------------------*/
    // MapFunction Class
    /*-----------------------------------------------------*/
    class MapFunction : public map_function {
    public:
        ex operator()(const ex &e);
        MapFunction(std::function<ex(const ex &, MapFunction &)>);
    private:
        std::function<ex(const ex &, MapFunction &)> Function;
    };
    
    
    /*-----------------------------------------------------*/
    // Parser Class
    /*-----------------------------------------------------*/
    class Parser {
    public:
        prototype_table FTable;
        symtab STable;
        ex Read(string instr);
        ex ReadFile(string filename);
        Parser();
        Parser(symtab st);
    };
    
    /*-----------------------------------------------------*/
    // isFunction
    /*-----------------------------------------------------*/
    inline bool isFunction(const ex &e, string func_name) {
        return is_a<GiNaC::function>(e) && ex_to<GiNaC::function>(e).get_name()==func_name;
    }
    inline bool isFunction(const ex &e, string func_name, int nargs) {
        return is_a<GiNaC::function>(e) && ex_to<GiNaC::function>(e).get_name()==func_name && e.nops()==nargs;
    }
    
    /*-----------------------------------------------------*/
    // string Functions
    /*-----------------------------------------------------*/
    void string_replace_all(string &str, const string &from, const string &to);
    void string_trim(string &str);
    
    void Combinations(int n, int m, std::function<void(const int*)> f);
    void CombinationsR(int n, int m, std::function<void(const int*)> f);
    void Permutations(int n, std::function<void(const int*)> f);
    void Permutations(int n, int m, std::function<void(const int*)> f);
    void PermutationsR(int n, int m, std::function<void(const int*)> f);
    
    /*-----------------------------------------------------*/
    // Rationalize
    /*-----------------------------------------------------*/
    extern MapFunction Rationalize;
    
    /*-----------------------------------------------------*/
    // LeafCount & sort
    /*-----------------------------------------------------*/
    long long int LeafCount(const ex & e);
    bool ex_less(const ex &a, const ex &b);
    void sort_lst(lst & ilst, bool less=true);
    void sort_lst_by(lst & ilst, int n, bool less=true);
    void sort_vec(exvector & ivec, bool less=true);
    void sort_vec_by(exvector & ivec, int n, bool less=true);
    
    /*-----------------------------------------------------*/
    // Other Functions
    /*-----------------------------------------------------*/
    inline void append_to(const exvector & exv, lst & alst) { for(auto item : exv) alst.append(item); }
    inline void append_to(const lst & alst, exvector & exv) { for(auto item : alst) exv.push_back(item); }
    
    
}
