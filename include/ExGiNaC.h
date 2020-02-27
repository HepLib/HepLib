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

namespace HepLib {

    using namespace GiNaC;
    using namespace std;

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

    /*-----------------------------------------------------*/
    // Global Symbol
    /*-----------------------------------------------------*/
    const symbol & get_symbol(const string & s);

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

    /*-----------------------------------------------------*/
    // vector : GiNaC_Parallel
    /*-----------------------------------------------------*/
    extern lst GiNaC_archive_Symbols;
    void GiNaC_archive_Symbols_from(ex);
    void GiNaC_archive_Symbols_from(vector<ex>);
    vector<ex> GiNaC_Parallel(
        int nproc,
        vector<ex> const &invec,
        std::function<ex(ex const &, int)> f,
        const char* key = NULL,
        int verb = 0,
        bool rm = true,
        int prtlvl = 0
    );
    
    /*-----------------------------------------------------*/
    // Helpers
    /*-----------------------------------------------------*/
    string RunOS(const char * cmd);
    void garRead(const char *garfn, map<string, ex> &resMap);
    ex garRead(const char *garfn, const char* key);
    ex garResult(const char *garfn);
    ex str2ex(const char *expr, symtab stab);
    lst str2lst(const char *expr, symtab stab);
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
    ex mma_series(ex const expr, symbol const s, int sn);
    bool has_pats(ex const &item, lst const &pats);
    ex mma_expand(ex const &expr, lst const &pats, int depth=0);
    ex mma_expand(ex const &expr, ex const &pat, int depth=0);
    ex mma_collect(ex const expr, lst const pat, bool ccf=false, bool cvf=false);
    ex mma_collect(ex const expr, ex const pat, bool ccf=false, bool cvf=false);
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
    
    /*-----------------------------------------------------*/
    // Global Colors
    /*-----------------------------------------------------*/
    extern const char* Color_Error;
    extern const char* Color_Warn;
    extern const char* Color_HighLight;

    /*-----------------------------------------------------*/
    // Customized GiNaC Function
    /*-----------------------------------------------------*/
    DECLARE_FUNCTION_1P(CCF)
    DECLARE_FUNCTION_1P(CVF)

    DECLARE_FUNCTION_1P(VF)
    DECLARE_FUNCTION_1P(VF1)
    DECLARE_FUNCTION_2P(VF2)
    DECLARE_FUNCTION_3P(VF3)

    DECLARE_FUNCTION_1P(FF) // not used internally, for user use only
    DECLARE_FUNCTION_2P(CV) // not used internally, for user use only

    DECLARE_FUNCTION_1P(x)
    DECLARE_FUNCTION_1P(y)
    DECLARE_FUNCTION_1P(z)
}
