#pragma once

#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <assert.h>
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

/*********************************************************/
// Global Symbol
/*********************************************************/
const symbol & get_symbol(const string & s);

/*********************************************************/
// split
/*********************************************************/
vector<std::string> split(const string& s, char delimiter);

/*********************************************************/
// Helper Classes
/*********************************************************/
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
    static lst sub(lst m, int s, int n);
    static ex sum(lst m);
    static lst subs(lst m, ex r);
    static lst subs(lst m, exmap r);
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



/*********************************************************/
// Global Functions
/*********************************************************/
string now(bool use_date = true);
lst gather_symbols(const ex & e);
lst gather_symbols(const vector<ex> & ve);
lst gather_symbols(const vector<pair<lst,lst>> & ve);
lst gather_symbols(const vector<pair<ex,ex>> & ve);

inline bool file_exists(const char* fn) {
    return (access(fn,F_OK)!=-1);
}

/*********************************************************/
// vector : GiNaC_Parallel
/*********************************************************/
template <typename F, typename T>
vector<ex> GiNaC_Parallel(int nproc, lst syms, vector<T> const &invec, F f, const char* key = NULL, int verb = 0, bool rm = true, int prtlvl = 0) {
    auto ppid = getpid();
    int para_run = 0;
    int para_max_run = nproc<0 ? omp_get_max_threads() : nproc;
    ostringstream cmd;
    cmd << "mkdir -p " << ppid;
    system(cmd.str().c_str());
    
    int total = invec.size();
    for(auto item : invec) {
        para_run++;
        if(verb > 1) {
            cout << "\r  ";
            for(int pi=0;pi<prtlvl;pi++) cout << "   ";
            cout << "\\--Evaluating ["<<para_run<<"/"<<total<<"] ... "<< flush;
        }

        if(para_max_run>0) {
            auto pid = fork();
            if (pid < 0) perror("fork() error");
            else if (pid != 0) {
                if(para_run >= para_max_run) wait(NULL);
                continue;
            }
        }
        
        try {
            auto res = f(item, para_run);
            archive a;
            a.archive_ex(res, "res");
            a.archive_ex(19790923, "c");
            ostringstream garfn;
            if(key == NULL) garfn << ppid << "/" << para_run << ".gar";
            else garfn << ppid << "/" << para_run << "." << key << ".gar";
            ofstream out(garfn.str().c_str());
            out << a;
            out.close();
        } catch(...) {
            if(para_max_run>0) exit(0);
            throw runtime_error("error in GiNaC_Parallel.");
        }
        if(para_max_run>0) exit(0);
    }
    
    auto cpid = getpid();
    if(cpid!=ppid) exit(0); // make sure
    if(para_max_run>0) while (wait(NULL) != -1) { }
    if(verb > 1 && para_run > 0) cout << "@" << now(false) << endl;

    vector<ex> ovec;
    for(int i=1; i<=para_run; i++) {
        if(verb > 1) {
            if(key == NULL) {
                cout << "\r  ";
                for(int pi=0; pi<prtlvl; pi++) cout << "   ";
                cout << "\\--Reading *.gar ["<<i<<"/"<<para_run<<"] ... "<< flush;
            } else {
                cout << "\r  ";
                for(int pi=0;pi<prtlvl;pi++) cout << "   ";
                cout << "\\--Reading *."<<key<<".gar ["<<i<<"/"<<para_run<<"] ... "<< flush;
            }
        }
        
        archive a;
        ostringstream garfn;
        if(key == NULL) garfn << ppid << "/" << i << ".gar";
        else garfn << ppid << "/" << i << "." << key << ".gar";
        ifstream in(garfn.str().c_str());
        in >> a;
        in.close();
        remove(garfn.str().c_str());
        auto c = a.unarchive_ex(syms, "c");
        auto res = a.unarchive_ex(syms, "res");
        if(c!=19790923) throw runtime_error("*.gar error!");
        ovec.push_back(res);
    }
    
    cmd.clear();
    cmd.str("");
    cmd << "rm -r " << ppid;
    if(rm) system(cmd.str().c_str());
    if(verb > 1 && para_run > 0) cout << "@" << now(false) << endl;
    return ovec;
}



/*********************************************************/
// lst : GiNaC_Parallel
/*********************************************************/
template <typename F>
lst GiNaC_Parallel(int nproc, lst syms, lst const &ilst, F f, const char* key = NULL, int verb = 0, bool rm = true, int prtlvl = 0) {
    auto ppid = getpid();
    int para_run = 0;
    int para_max_run = nproc<0 ? omp_get_max_threads() : nproc;
    ostringstream cmd;
    cmd << "mkdir -p " << ppid;
    system(cmd.str().c_str());
    
    int total = ilst.nops();
    for(int i=0; i<total; i++) {
        auto item = ilst.op(i);
        para_run++;
        if(verb > 1) {
            cout << "\r  ";
            for(int pi=0;pi<prtlvl;pi++) cout << "   ";
            cout << "\\--Evaluating ["<<para_run<<"/"<<total<<"] ... "<< flush;
        }
        
        if(para_max_run>0) {
            auto pid = fork();
            if (pid < 0) perror("fork() error");
            else if (pid != 0) {
                if(para_run >= para_max_run) wait(NULL);
                continue;
            }
        }
        
        try {
            auto res = f(item, para_run);
            archive a;
            a.archive_ex(res, "res");
            a.archive_ex(19790923, "c");
            ostringstream garfn;
            if(key == NULL) garfn << ppid << "/" << para_run << ".gar";
            else garfn << ppid << "/" << para_run << "." << key << ".gar";
            ofstream out(garfn.str().c_str());
            out << a;
            out.close();
        } catch(...) {
            if(para_max_run>0) exit(0);
            throw runtime_error("error in GiNaC_Parallel.");
        }
        if(para_max_run>0) exit(0);
    }
    
    auto cpid = getpid();
    if(cpid!=ppid) exit(0); // make sure
    if(para_max_run>0) while (wait(NULL) != -1) { }
    if(verb > 1 && para_run > 0) cout << "@" << now(false) << endl;

    lst olst;
    for(int i=1; i<=para_run; i++) {
        if(verb > 1) {
            if(key == NULL) {
                cout << "\r  ";
                for(int pi=0; pi<prtlvl; pi++) cout << "   ";
                cout << "\\--Reading *.gar ["<<i<<"/"<<para_run<<"] ... "<< flush;
            } else {
                cout << "\r  ";
                for(int pi=0;pi<prtlvl;pi++) cout << "   ";
                cout << "\\--Reading *."<<key<<".gar ["<<i<<"/"<<para_run<<"] ... "<< flush;
            }
        }
        
        archive a;
        ostringstream garfn;
        if(key == NULL) garfn << ppid << "/" << i << ".gar";
        else garfn << ppid << "/" << i << "." << key << ".gar";
        ifstream in(garfn.str().c_str());
        in >> a;
        in.close();
        remove(garfn.str().c_str());
        auto c = a.unarchive_ex(syms, "c");
        auto res = a.unarchive_ex(syms, "res");
        if(c!=19790923) throw runtime_error("*.gar error!");
        olst.append(res);
    }
    
    cmd.clear();
    cmd.str("");
    cmd << "rm -r " << ppid;
    if(rm) system(cmd.str().c_str());
    if(verb > 1 && para_run > 0) cout << "@" << now(false) << endl;
    return olst;
}

template <typename FUN>
ex GiNaC_Replace(ex expr, const ex &pat, FUN f) {
    exset ps;
    expr.find(pat, ps);
    lst repl;
    for(auto pi : ps) {
        auto pif = f(pi);
        repl.append(pi==pif);
    }
    return expr.subs(repl);
}

string RunOS(const char * cmd);

ex garResult(const char *garfn, lst syms);

/*********************************************************/
// Series at s=0 similar to Mathematica
/*********************************************************/
ex mma_series(ex expr, symbol s, int sn);

/*********************************************************/
// Customized GiNaC Function
/*********************************************************/
DECLARE_FUNCTION_1P(VF)
DECLARE_FUNCTION_1P(VF1)
DECLARE_FUNCTION_2P(VF2)
DECLARE_FUNCTION_3P(VF3)

/*********************************************************/
// Terminal Color
/*********************************************************/
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

}
