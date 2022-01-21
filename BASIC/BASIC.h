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
classname::classname() { setflag(status_flags::evaluated | status_flags::expanded | status_flags::hf_expanded); }

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
    
    extern string Version;
    
    inline bool has_any(ex expr, lst ps) {
        for(auto pi : ps) if(expr.has(pi)) return true;
        return false;
    }
    inline bool has_all(ex expr, lst ps) {
        for(auto pi : ps) if(!expr.has(pi)) return false;
        return true;
    }
    
    inline bool is_equal_any(ex expr, lst ps) {
        for(auto pi : ps) if(expr.is_equal(pi)) return true;
        return false;
    }
    
    inline bool match_any(ex expr, lst ps) {
        for(auto pi : ps) if(expr.match(pi)) return true;
        return false;
    }
    inline bool match_all(ex expr, lst ps) {
        for(auto pi : ps) if(!expr.match(pi)) return false;
        return true;
    }
    
    // Terminal Color
    #define RESET "\033[0m"
    #define BLACK "\033[30m"
    #define RED "\033[31m"
    #define GREEN "\033[32m"
    #define YELLOW "\033[33m"
    #define BLUE "\033[34m"
    #define MAGENTA "\033[35m"
    #define CYAN "\033[36m"
    #define WHITE "\033[37m"
    #define BOLDBLACK "\033[1m\033[30m"
    #define BOLDRED "\033[1m\033[31m"
    #define BOLDGREEN "\033[1m\033[32m"
    #define BOLDYELLOW "\033[1m\033[33m"
    #define BOLDBLUE "\033[1m\033[34m"
    #define BOLDMAGENTA "\033[1m\033[35m"
    #define BOLDCYAN "\033[1m\033[36m"
    #define BOLDWHITE "\033[1m\033[37m"
    
    /**
     * @brief class extended to GiNaC symbol class, represent a positive symbol
     */
    class Symbol : public symbol {
    //GINAC_DECLARE_REGISTERED_CLASS(Symbol, symbol)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const Symbol &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef symbol inherited; // supername
        Symbol(); // classname
        Symbol * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        explicit Symbol(const string &s);
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        
        //bool info(unsigned inf) const override; 
        ex eval() const override; // for performance reasons
        ex evalf() const override; // for performance reasons
        ex series(const relational & s, int order, unsigned options = 0) const override;
        ex subs(const exmap & m, unsigned options = 0) const override { return subs_one_level(m, options); } // overwrites basic::subs() for performance reasons
        //ex normal(exmap & repl, exmap & rev_lookup, lst & modifier) const override;
        //ex to_rational(exmap & repl) const override;
        //ex to_polynomial(exmap & repl) const override;
        //bool is_polynomial(const ex & var) const override;
    
        ex conjugate() const override;
        ex real_part() const override;
        ex imag_part() const override;
        void set_name(string n);
        
        unsigned get_domain() const override { return domain::positive; }
        
        static bool has(const ex &e);
        static lst all(const ex &e);
        static std::map<std::string, ex> Table;
        
        void set(const ex & v) const;
        void unset() const;
        
        unsigned get_serial() { return serial; }
        
        static exmap vmap;
        static void set(const Symbol & s, const ex & v);
        static void set(const string & str, const ex & v);
        static void unset(const Symbol &s);
        static void unset(const string &str);
        static void unset_all();
        static ex set_all(const ex & expr);   
    protected:
        unsigned calchash() const override;
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;
    };
    
    /**
     * @brief class extended to GiNaC symbol class, pure imaginary symbol 
     */
    class iSymbol : public symbol {
    //GINAC_DECLARE_REGISTERED_CLASS(iSymbol, symbol)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const iSymbol &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef symbol inherited; // supername
        iSymbol(); // classname
        iSymbol * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        iSymbol(const string &s);
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        
        unsigned get_domain() const override { return domain::complex; }
        
        //bool info(unsigned inf) const override; 
        ex eval() const override; // for performance reasons
        ex evalf() const override; // for performance reasons
        ex series(const relational & s, int order, unsigned options = 0) const override;
        ex subs(const exmap & m, unsigned options = 0) const override { return subs_one_level(m, options); } // overwrites basic::subs() for performance reasons
        //ex normal(exmap & repl, exmap & rev_lookup, lst & modifier) const override;
        //ex to_rational(exmap & repl) const override;
        //ex to_polynomial(exmap & repl) const override;
        //bool is_polynomial(const ex & var) const override;
        
        ex conjugate() const override;
        ex real_part() const override;
        ex imag_part() const override;
        void set_name(string n);
        
        static bool has(const ex &e);
        static lst all(const ex &e);
        static std::map<std::string, ex> Table; 
        
    protected:
        unsigned calchash() const override;
        ex derivative(const symbol & s) const override;
        bool is_equal_same_type(const basic & other) const override;       
    };
    
    
    /**
     * @brief class used to wrap error message 
     */
    class Error : public exception {
    public:
        string msg;
        const char * what() const throw ();
        Error(string _msg);
    };

    vector<std::string> split(const string& s, char delimiter);


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
    lst gather_symbols(const exvector & ve);

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
        
    inline ex subs_all(const ex & expr, const ex & repls) {
        return subs(expr, repls, subs_options::algebraic);
    };
    inline ex subs_all(const ex & expr, const exmap & repls) {
        return subs(expr, repls, subs_options::algebraic);
    };
    
    /*-----------------------------------------------------*/
    // vector : GiNaC_Parallel
    /*-----------------------------------------------------*/
    extern int GiNaC_Parallel_Process;
    extern map<string, int> GiNaC_Parallel_NP;
    extern int GiNaC_Parallel_Batch;
    extern map<string,int> GiNaC_Parallel_NB;
    extern map<string,bool> GiNaC_Parallel_RM;
    extern map<string,string> GiNaC_Parallel_PRE;
    exvector GiNaC_Parallel(int ntotal, std::function<ex(int)> f, const string & key = "");
    
    /*-----------------------------------------------------*/
    // Helpers
    /*-----------------------------------------------------*/
    string RunOS(const string & cmd);
    void garRead(const string &garfn, map<string, ex> &resMap);
    ex garRead(const string &garfn, const char* key);
    ex garRead(const string &garfn);
    void garWrite(const string &garfn, const map<string, ex> &resMap);
    inline void garWrite(const map<string, ex> &resMap, const string &garfn) { garWrite(garfn,resMap); }
    void garWrite(const string &garfn, const ex & res);
    inline void garWrite(const ex & res, const string &garfn) { garWrite(garfn,res); }
    
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
    string ex2str(const exvector &expr);
    string ex2str(const exmap &expr);
    string ex2str(const exset &expr);
    inline string in2str(int i) { return to_string(i); }
    #ifdef _USE_FLOAT128
    ex q2ex(__float128);
    __float128 ex2q(ex);
    #else
    ex q2ex(long double);
    long double ex2q(ex);
    #endif
    lst vec2lst(const exvector & ev);
    exvector lst2vec(const lst & alst);
    lst add2lst(const ex & expr);
    lst mul2lst(const ex & expr);
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
    ex series_ex(ex const & expr, const symbol &s, int sn);
        
    ex expand_ex(const ex &expr, std::function<bool(const ex &)>);
    
    inline ex expand_ex(const ex &expr) {
        return expand_ex(expr, [](const ex & e)->bool{ return true; });
    }
    
    inline ex expand_ex(ex const &expr, lst const &pats) {
        return expand_ex(expr, [pats](const ex & e)->bool {
            for(auto pat : pats) { if(e.has(pat)) return true; }
            return false;
        });
    }
    
    inline ex expand_ex(ex const &expr, ex const &pat) {
        return expand_ex(expr, [pat](const ex & e)->bool { return e.has(pat); });
    }
    
    ex collect_ex(const ex &expr, std::function<bool(const ex &)>, bool cf=false, bool vf=false, int opt=0);
    
    inline ex collect_ex(const ex &expr, lst const &pats, bool cf=false, bool vf=false, int opt=0) {
        return collect_ex(expr, [pats](const ex & e)->bool {
            for(auto pat : pats) { if(e.has(pat)) return true; }
            return false;
        }, cf, vf, opt);
    }
    
    inline ex collect_ex(const ex &expr, ex const &pat, bool cf=false, bool vf=false, int opt=0) {
        return collect_ex(expr, [pat](const ex & e)->bool { return e.has(pat); }, cf, vf, opt);
    }
    
    inline ex collect_o(const ex &expr, std::function<bool(const ex &)> f, int opt=0) {
        return collect_ex(expr, f, false, false, opt);
    }
    inline ex collect_o(const ex &expr, lst const &pats, int opt=0) {
        return collect_ex(expr, pats, false, false, opt);
    }
    inline ex collect_o(const ex &expr, ex const &pat, int opt=0) {
        return collect_ex(expr, pat, false, false, opt);
    }
    
    inline ex collect_c(const ex &expr, std::function<bool(const ex &)> f, int opt=0) {
        return collect_ex(expr, f, true, false, opt);
    }
    inline ex collect_c(const ex &expr, lst const &pats, int opt=0) {
        return collect_ex(expr, pats, true, false, opt);
    }
    inline ex collect_c(const ex &expr, ex const &pat, int opt=0) {
        return collect_ex(expr, pat, true, false, opt);
    }

    inline ex collect_v(const ex &expr, std::function<bool(const ex &)> f, int opt=0) {
        return collect_ex(expr, f, false, true, opt);
    }
    inline ex collect_v(const ex &expr, lst const &pats, int opt=0) {
        return collect_ex(expr, pats, false, true, opt);
    }
    inline ex collect_v(const ex &expr, ex const &pat, int opt=0) {
        return collect_ex(expr, pat, false, true, opt);
    }
    
    inline ex collect_cv(const ex &expr, std::function<bool(const ex &)> f, int opt=0) {
        return collect_ex(expr, f, true, true, opt);
    }
    inline ex collect_cv(const ex &expr, lst const &pats, int opt=0) {
        return collect_ex(expr, pats, true, true, opt);
    }
    inline ex collect_cv(const ex &expr, ex const &pat, int opt=0) {
        return collect_ex(expr, pat, true, true, opt);
    }

    lst collect_lst(const ex &expr, std::function<bool(const ex &)>, int opt=0);
    
    inline lst collect_lst(const ex &expr, lst const &pats, int opt=0) {
        return collect_lst(expr, [pats](const ex & e)->bool {
            for(auto pat : pats) { if(e.has(pat)) return true; }
            return false;
        }, opt);
    }
    
    inline lst collect_lst(const ex &expr, ex const &pat, int opt=0) {
        return collect_lst(expr, [pat](const ex & e)->bool { return e.has(pat); }, opt);
    } 
    
    ex diff_ex(ex const expr, ex const xp, unsigned nth=1, bool expand=false);
    
    extern bool using_cache;
    extern long long cache_limit;
    extern bool fermat_using_array;
    ex fermat_eval(const ex & expr);
    ex numer_denom_fermat(const ex & expr, bool dfactor=false);
    ex numer_fermat(const ex & expr);
    inline ex fermat_numer_denom(const ex & expr, bool dfactor=false) { return numer_denom_fermat(expr,dfactor); }
    
    extern map<ex,long long,ex_is_less> fermat_weight;
    ex normal_fermat(const ex & expr, bool dfactor=false);
    inline ex fermat_normal(const ex & expr, bool dfactor=false) { return normal_fermat(expr,dfactor); }
    
    ex form_eval(const ex & expr);
    ex factor_form(const ex & expr, bool nd=true);
    inline ex form_factor(const ex & expr, bool nd=true) { return factor_form(expr,nd); }
    
    ex exfactor(const ex & expr, int opt = 1);
    ex exnormal(const ex & expr, int opt = -5);
    ex exnd(const ex & expr, int opt = 1);
    
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
    extern bool Debug;

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
    
    #ifndef DOXYGEN_SKIP
    
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
    
    #endif
    
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
    // sort
    /*-----------------------------------------------------*/
    long long node_number(const ex & expr, int level=0);
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
    //GINAC_DECLARE_REGISTERED_CLASS(XIntegral, basic)
    private:
        static GiNaC::registered_class_info reg_info;
    public:
        static GiNaC::registered_class_info &get_class_info_static();
        class visitor {
        public:
            virtual void visit(const XIntegral &) = 0; // classname
            virtual ~visitor();
        };
        template<class B, typename... Args> friend B & dynallocate(Args &&... args);
        typedef basic inherited; // supername
        XIntegral(); // classname
        XIntegral * duplicate() const override; // classname
        void accept(GiNaC::visitor & v) const override;
        const GiNaC::registered_class_info &get_class_info() const override;
        GiNaC::registered_class_info &get_class_info() override;
        const char *class_name() const override;
    protected:
        int compare_same_type(const GiNaC::basic & other) const override;
    // GINAC_DECLARE_REGISTERED_CLASS END
    
    public:
        ex Functions=lst{};
        ex Exponents=lst{};
        ex Deltas=lst{};
        size_t nops() const override;
        ex op(size_t i) const override;
        ex & let_op(size_t i) override;
        void print(const print_dflt &c, unsigned level = 0) const;
        void archive(archive_node & n) const override;
        void read_archive(const archive_node& n) override;
        static bool has(const ex &e);
        static lst all(const ex &e);
        XIntegral(ex fed);
        XIntegral(ex loops, ex ps, ex ns);        
    };
    
    #ifndef DOXYGEN_SKIP
    class _global_init {
    public:
        class _init {
            public: _init();
        };
    private:
        static _init init_object;
    };
    #endif
    
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
    
        
    /**
     * @brief class for Mathematica Format Output
     */
    class MMAFormat : public print_dflt {
        GINAC_DECLARE_PRINT_CONTEXT(MMAFormat, print_dflt)
    public:
        MMAFormat(ostream &os, unsigned opt=0);
        
        template<class T> const MMAFormat & operator << (const T & v) const {
            s << v;
            return *this;
        };
        const MMAFormat & operator << (const basic & v) const;
        const MMAFormat & operator << (const ex & v) const;
        const MMAFormat & operator << (const lst & v) const;
        const MMAFormat & operator<<(std::ostream& (*v)(std::ostream&)) const;
        const MMAFormat & operator << (const matrix & v) const;
        const MMAFormat & operator << (const exvector & v) const;
        const MMAFormat & operator << (const exmap & v) const;
        const MMAFormat & operator << (const exset & v) const;
    };
    extern MMAFormat mout;
    
    void garWrite(const exvector &exv, string garfn);
    inline void garWrite(string garfn, const exvector &exv) { garWrite(exv,garfn); }
    void garRead(exvector &exv, string garfn);
    inline void garRead(string garfn, exvector &exv) { garRead(exv, garfn); }
    ex add_collect_normal(const exvector &exv, lst const &pats);
    
    class Server {
    public:
        int Round = 3;
        int Port = 0;
        unsigned Total;
        void Start();
        int Verbose = 1;
        string Skip = "[ID].run";
        string DL;
        string FUNC;
        
        static string Next(string sip, string sport);
    };
    
    void ReShare(const ex & e);
    void ReShare(const lst & e);
    void ReShare(const ex & e1, const ex & e2);
    void ReShare(const ex & e1, const ex & e2, const ex & e3);
    void ReShare(const exvector & ev);
    void ReShare(const exvector & ev1, const exvector & ev2);
}

typedef void (*RUN)(std::string dir_id);
