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
#include <sys/stat.h>

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
    
/**
 * @brief Extension to GiNaC
 */
namespace GiNaC {
    ex ginac_factor(const ex& poly, unsigned options=0);
}

/**
 * @brief HepLib namespace
 */
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
    
    /**
     * @brief class extended to GiNaC symbol class, represent a positive symbol
     */
    class Symbol : public symbol {
    GINAC_DECLARE_REGISTERED_CLASS(Symbol, symbol)
    public:
        Symbol(const string &s);
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        
        ex eval() const override; // for performance reasons
        ex evalf() const override; // for performance reasons
        ex conjugate() const override;
        ex real_part() const override;
        ex imag_part() const override;
        
        unsigned get_domain() const override { return domain::positive; }
        
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
    
    /**
     * @brief class extended to GiNaC symbol class, pure imaginary symbol 
     */
    class iSymbol : public symbol {
    GINAC_DECLARE_REGISTERED_CLASS(iSymbol, symbol)
    public:
        iSymbol(const string &s);
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        
        unsigned get_domain() const override { return domain::complex; }
        
        ex eval() const override; // for performance reasons
        ex evalf() const override; // for performance reasons
        ex conjugate() const override;
        ex real_part() const override;
        ex imag_part() const override;
        
        static bool has(const ex &e);
        static lst all(const ex &e);
        static std::map<std::string, ex> Table;
    };
    GINAC_DECLARE_UNARCHIVER(iSymbol);
    
    
    /**
     * @brief class used to wrap error message 
     */
    class Error : public exception {
    public:
        string msg;
        const char * what() const throw ();
        Error(string _msg);
    };

    /*-----------------------------------------------------*/
    // Global Symbol
    /*-----------------------------------------------------*/
    const symbol & get_symbol(const string & s);

    /*-----------------------------------------------------*/
    // split
    /*-----------------------------------------------------*/
    vector<std::string> split(const string& s, char delimiter);

    /**
     * @brief class as Matrix Helper 
     */
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

    /**
     * @brief class as lst Helper 
     */
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

    inline bool file_exists(string fn) {
        return (access(fn.c_str(),F_OK)!=-1);
    }
    
    inline bool key_exists(const exmap &map, const ex & key) {
        return (map.find(key)!=map.end());
    }
    
    inline bool dir_exists(string dir) {
        struct stat buffer;
        return (stat(dir.c_str(), &buffer)==0);
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
        const string & pre = "  "
    );
    inline vector<ex> GiNaC_Parallel(
        int ntotal,
        std::function<ex(int)> f,
        const string & key = "",
        bool rm = true,
        const string & pre = "  "
    ) { return  GiNaC_Parallel(ntotal, 0, f, key, rm, pre); }
    
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
    vector<string> file2strvec(string filename, bool skip_empty=false);
    ex file2ex(string filename);
    ex file2ex(string filename, symtab st);
    int ex2int(ex);
    void ex2file(const ex &, string filename);
    void ex2file(string filename, const ex &);
    string ex2str(const ex &expr);
    ex q2ex(__float128);
    __float128 ex2q(ex);
    lst exvec2lst(const exvector & exvec);
    exvector lst2exvec(const lst & alst);
    lst xlst(int ei);
    lst xlst(int bi, int ei);
    int CpuCores();
    
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
    
    pair<ex,epvector> mma_expand(const ex &expr, std::function<bool(const ex &)>, int depth);
    ex mma_expand(const ex &expr, std::function<bool(const ex &)>);
    ex mma_expand(ex const &expr, lst const &pats);
    ex mma_expand(ex const &expr, ex const &pat);
    
    ex mma_collect(const ex &expr, std::function<bool(const ex &)>, bool ccf=false, bool cvf=false);
    ex mma_collect(const ex &expr, lst const &pats, bool ccf=false, bool cvf=false);
    ex mma_collect(const ex &expr, ex const &pat, bool ccf=false, bool cvf=false);
    
    lst mma_collect_lst(const ex &expr, std::function<bool(const ex &)>);
    lst mma_collect_lst(const ex &expr, lst const &pats);
    lst mma_collect_lst(const ex &expr, ex const &pat);
    
    ex mma_diff(ex const expr, ex const xp, unsigned nth=1, bool expand=false);
    
    extern bool fermat_use_array;
    ex fermat_numer_denom(const ex & expr, bool factor=false);
    ex fermat_normal(const ex & expr, bool factor=false);
    
    ex form_factor(const ex & expr);
    
    enum FactorMethod { GiNaC, FORM };
    ex exfactor(const ex & expr, FactorMethod fm = FORM);
    
    ex collect_factors(const ex & expr);
    
    /*-----------------------------------------------------*/
    // EvalF/D/Q/MP
    /*-----------------------------------------------------*/
    ex EvalF(ex);
    ex EvalL(ex);
    ex EvalQ(ex);
    ex EvalMP(ex);
    extern int NNDigits;
    ex NN(ex expr,int digits=NNDigits);

    /*-----------------------------------------------------*/
    // xPositive
    /*-----------------------------------------------------*/
    bool xPositive(ex const expr);
    int xSign(ex const expr);

    /*-----------------------------------------------------*/
    // Global object wildcard
    /*-----------------------------------------------------*/
    extern ex w, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9;
    extern string InstallPrefix;
    extern string INC_FLAGS;
    extern string LIB_FLAGS;
    extern const iSymbol iEpsilon;
    extern const ex iEpsilonN;
    extern const Symbol ep;
    extern const Symbol D;
    extern int Verbose;
    extern int ParallelProcess;

    /*-----------------------------------------------------*/
    // Global Colors
    /*-----------------------------------------------------*/
    extern const char* ErrColor;
    extern const char* WarnColor;
    extern const char* Color_HighLight;

    /*-----------------------------------------------------*/
    // Customized GiNaC Function
    /*-----------------------------------------------------*/
    DECLARE_FUNCTION_1P(coCF)
    DECLARE_FUNCTION_1P(coVF)    
    DECLARE_FUNCTION_1P(x)
    DECLARE_FUNCTION_1P(y)
    DECLARE_FUNCTION_1P(z)
    
    /**
     * @brief F function with 1 argument  
     */
    class F1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function F(const T1 & p1) {
        return GiNaC::function(F1_SERIAL::serial, ex(p1));
    }
    
    /**
     * @brief F function with 2 arguments
     */
    class F2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function F(const T1 & p1, const T2 & p2) {
        return GiNaC::function(F2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    /**
     * @brief WF function with 1 argument  
     */
    class WF1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function WF(const T1 & p1) {
        return GiNaC::function(WF1_SERIAL::serial, ex(p1));
    }
    
    /**
     * @brief WF function with 2 arguments
     */
    class WF2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2) {
        return GiNaC::function(WF2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    /**
     * @brief WF function with 3 arguments
     */
    class WF3_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2, const T3 & p3) {
        return GiNaC::function(WF3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
    }
    
    /**
     * @brief WF function with 4 arguments
     */
    class WF4_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4) {
        return GiNaC::function(WF4_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4));
    }
    
    /**
     * @brief WF function with 5 arguments
     */
    class WF5_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    inline GiNaC::function WF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5) {
        return GiNaC::function(WF5_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5));
    }
    
    /**
     * @brief inner iWF function with 1 argument
     */
    class iWF1_SERIAL { public: static unsigned serial; };
    template<typename T1>
    inline GiNaC::function iWF(const T1 & p1) {
        return GiNaC::function(iWF1_SERIAL::serial, ex(p1));
    }
    
    /**
     * @brief inner iWF function with 2 arguments
     */
    class iWF2_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2) {
        return GiNaC::function(iWF2_SERIAL::serial, ex(p1), ex(p2));
    }
    
    /**
     * @brief inner iWF function with 3 arguments
     */
    class iWF3_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2, const T3 & p3) {
        return GiNaC::function(iWF3_SERIAL::serial, ex(p1), ex(p2), ex(p3));
    }
    
    /**
     * @brief inner iWF function with 4 arguments
     */
    class iWF4_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4) {
        return GiNaC::function(iWF4_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4));
    }
    
    /**
     * @brief inner iWF function with 5 arguments
     */
    class iWF5_SERIAL { public: static unsigned serial; };
    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    inline GiNaC::function iWF(const T1 & p1, const T2 & p2, const T3 & p3, const T4 & p4, const T5 & p5) {
        return GiNaC::function(iWF5_SERIAL::serial, ex(p1), ex(p2), ex(p3), ex(p4), ex(p5));
    }
    
    /**
     * @brief class to wrap map_function of GiNaC
     */
    class MapFunction : public map_function {
    public:
        ex operator()(const ex &e);
        MapFunction(std::function<ex(const ex &, MapFunction &)>);
        static ex subs(const ex & expr, const ex & pat, std::function<ex(const ex &)> f);
    private:
        std::function<ex(const ex &, MapFunction &)> Function;
    };
    
    
    /**
     * @brief class to parse for string or file, helpful with Symbol class
     */
    class Parser {
    public:
        prototype_table FTable;
        symtab STable;
        ex Read(string instr,bool s2S=true);
        ex ReadFile(string filename,bool s2S=true);
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
        return is_a<GiNaC::function>(e) && ex_to<GiNaC::function>(e).get_name()==func_name && (e.nops()-nargs)==0;
    }
    
    /*-----------------------------------------------------*/
    // hasFunction
    /*-----------------------------------------------------*/
    inline bool has_function(const ex & expr) {
        for(const_preorder_iterator i = expr.preorder_begin(); i != expr.preorder_end(); ++i) {
            if(is_a<GiNaC::function>(*i)) return true;
        }
        return false;
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
    bool isSorted(const lst & exs);
    bool isSorted(int n, const ex exs[]);
    int ACSort(lst & exs);
    int ACSort(int n, ex exs[]);
    
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
    inline lst CoPat(const ex & e, std::function<bool(const ex &)> f) {  
        if(is_a<mul>(e)) {
            ex cc=1, vv=1;
            for(auto item : e) {
                if(f(item)) vv *= item;
                else cc *= item;
            }
            return lst{cc,vv};
        } else if(f(e)) return lst{1,e};
        else return lst{e,1};
    }
    
    /**
     * @brief XIntegral Class, preface to SecDec
     */
    class XIntegral : public basic {
    GINAC_DECLARE_REGISTERED_CLASS(XIntegral, basic)
    public:
        ex Functions=lst{};
        ex Exponents=lst{};
        ex Deltas=lst{};
        size_t nops() const override;
        ex op(size_t i) const override;
        ex & let_op(size_t i) override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n, lst& sym_lst) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        XIntegral(ex fed);
        XIntegral(ex loops, ex ps, ex ns);
    };
    GINAC_DECLARE_UNARCHIVER(XIntegral);
    
    // global init class
    class _global_init {
    public:
        class _init {
            public: _init();
        };
    private:
        static _init init_object;
    };
    
    /*-----------------------------------------------------*/
    // Interface to FORM & Fermat
    /*-----------------------------------------------------*/
    
    /**
     * @brief interface to communicate with Fermat program
     */
    class Fermat {
    public:
        static int buffer_size;
        string Sentinel = "---EOF---";
        void Init(string fer_path="fer64");
        string Execute(string);
        void Exit();
        ~Fermat();
                
    private:
        bool inited = false;
        bool exited = false;
        int P2C[2];
        int C2P[2];
        pid_t fpid = 0;
        pid_t pid = 0;
    };
    
    /**
     * @brief interface to communicate with Form program
     */
    class Form {
    public:
        static int buffer_size;
        string Sentinel = "---EOF---";
        string Prompt = "***EOF***";
        void Init(string form_path="form");
        string Execute(string script, const string & out_var="[o]");
        void Exit();
        ~Form();
    
    private:
        bool inited = false;
        bool exited = false;
        int io[2][2];
        int stdo[2];
        pid_t fpid = 0;
        pid_t pid = 0;
    };
    
}
