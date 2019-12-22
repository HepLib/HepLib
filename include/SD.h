#pragma once

#include "ExGiNaC.h"

#include <dlfcn.h>

#include <string>
#include <fcntl.h>
#include <signal.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <sstream>
#include <ios>
#include <regex>

extern "C" {
    #include <quadmath.h>
}

namespace HepLib {

/*-----------------------------------------------------*/
// SD Input
/*-----------------------------------------------------*/
struct FeynmanParameter {
    lst LoopMomenta;
    lst tLoopMomenta;
    lst Propagators;
    lst Exponents;
    exmap lReplacements;
    exmap tReplacements;
    exmap nReplacements;
    ex Prefactor = 1;
    bool isQuasi = false;
    bool isAsy = false;
};

struct XIntegrand {
    lst Functions;
    lst Exponents;
    exmap nReplacements;
    vector<lst> Deltas;
    bool isAsy = false;
};

/*-----------------------------------------------------*/
// Global Functions
/*-----------------------------------------------------*/
vector<ex> get_xy_from(ex pol);
vector<ex> get_x_from(ex pol);
vector<ex> get_y_from(ex pol);
vector<ex> get_pl_from(ex pol);

/*-----------------------------------------------------*/
// SecDec Classes
/*-----------------------------------------------------*/
class SecDecBase {
public:
    virtual vector<exmap> x2y(const ex &xpol) =0;
    vector<exmap> x2y(const lst &xpols);
};

class SecDecG : public SecDecBase {
public:
    vector<exmap> x2y(const ex &xpol) override;
    static vector<vector<int>> RunQHull(const matrix &pts);
private:
    vector<matrix> ZeroFaces(const matrix &pts);
    matrix NormalVectors(const vector<matrix> &zfs);
    matrix DualCone(const matrix &pts);
    vector<vector<int>> QHull(const matrix &dc, int dim);
    vector<matrix> SimplexifyR(const matrix &dc, int dim);
    vector<matrix> Simplexify(const matrix &dc, int dim);
    vector<matrix> SimplexCones(matrix pts);
};

/*-----------------------------------------------------*/
// Integrator Classes
/*-----------------------------------------------------*/
typedef __float128 qREAL;
typedef __complex128 qCOMPLEX;

class IntegratorBase {
public:
    typedef int (*SD_Type) (const unsigned int xn, const qREAL xx[], const unsigned int yn, qREAL y[], const qREAL pl[], const qREAL las[]);
    typedef qREAL (*FT_Type) (const qREAL xx[], const qREAL pl[]);
    virtual ex Integrate() =0;
    virtual int inDQMP(qREAL const *x);
    
    FT_Type FT = NULL;
    SD_Type Integrand = NULL;
    SD_Type IntegrandQ = NULL;
    SD_Type IntegrandMP = NULL;
    const qREAL* Lambda;
    const qREAL* Parameter;
    int XDim;
    
    long long RunMAX = 100;
    long long RunPTS = 10000;
    qREAL EpsAbs = 1E-5;
    qREAL EpsRel = 0;
    int ReIm = 3; // 1-Re, 2-Im, 3-ReIm
    int Verbose = 0;
    int NANMax = 250;
    int nNAN = 0;
    
    int DQMP = 0;
    int QXDim = 2;
    int MPXDim = 1;
    qREAL QXLimit = 1.Q-4;
    qREAL MPXLimit = 1.Q-8;
    qREAL QFLimit = 1.Q-3;
    qREAL MPFLimit = 1.Q-5;
    
    bool UseCpp = true;
    long long NEval = 0;
    int MPDigits = 50;
};

class HCubature : public IntegratorBase {
public:
    static int Wrapper(unsigned int xdim, long long npts, const qREAL *x, void *fdata, unsigned int ydim, qREAL *y);
    
    typedef void (*PrintHookerType) (qREAL*, qREAL*, long long int*, void *);
    
    virtual ex Integrate() override;
    static void DefaultPrintHooker(qREAL*, qREAL*, long long int*, void*);
    PrintHookerType PrintHooker = DefaultPrintHooker;
    long long MaxPTS;
    bool use_last = false;
    
private:
    qREAL LastResult[2];
    qREAL LastAbsErr[2];
    long long lastNRUN = 0;
    int lastnNAN = 0;
    int LastState = 0;
};

class CUBA : public IntegratorBase {
public:
    
    enum METHOD { VEGAS, CUHRE };
    METHOD Method = CUHRE;
    int VERBOSE = 0;
    
    // CUHRE Parameters
    int CUHRE_KEY = 0;
    
    // VEGAS Parameters
    int VEGAS_SEED = 0;
    long long VEGAS_NSTART = 1000;
    long long VEGAS_NINCREASE = 1000;
    long long VEGAS_NBATCH = 1000;
    
    static int Wrapper(const int *xdim, const qREAL *x, const int *ydim, qREAL *y, void *fdata);
    long long MaxPTS;
    virtual ex Integrate() override;
    
private:
    qREAL LastResult[2];
    qREAL LastAbsErr[2];
};

/*-----------------------------------------------------*/
// Minimize Classes
/*-----------------------------------------------------*/
typedef long double dREAL;
class MinimizeBase {
public:
    typedef dREAL (*FunctionType)(int nvars, dREAL* x, dREAL* pl, dREAL *las);
    virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false, int TryPTS=0, int SavePTS=0) =0;
    dREAL ZeroValue = -1E-20;
    virtual void Minimize(int nvars, FunctionType func, dREAL *ip)=0;
    virtual void ForceStop()=0;
};

class HookeJeeves : public MinimizeBase {
public:
    virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false, int TryPTS=0, int SavePTS=0) override;
    bool Exit = false;
    virtual void Minimize(int nvars, FunctionType func, dREAL *ip) override;
    virtual void ForceStop() override;
    
private:
    dREAL best_nearby(dREAL* delta, dREAL* point, dREAL prevbest, int nvars);
    int hooke(int nvars, dREAL* startpt, dREAL* endpt, dREAL rho, dREAL epsilon, int itermax);
    dREAL ObjectWrapper(int nvars, dREAL* x);
    FunctionType ObjectFunction;
    dREAL UpperBound[50];
    dREAL LowerBound[50];
    dREAL *PL;
    dREAL *LAS;
};

class MinUit : public MinimizeBase {
public:    
    virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false, int TryPTS=0, int SavePTS=0) override;
    virtual void Minimize(int nvars, FunctionType func, dREAL *IP) override;
    virtual void ForceStop() override;
};

/*-----------------------------------------------------*/
// CppFormat Class
/*-----------------------------------------------------*/
class CppFormat : public print_csrc_cl_N {
    GINAC_DECLARE_PRINT_CONTEXT(CppFormat, print_csrc_cl_N)
public:
    CppFormat(ostream &os, const char* s = "L", unsigned opt = 0);
    static ex q2ex(qREAL);
    static qREAL ex2q(ex);
    const char* suffix;
    static void QPrint(qREAL qr);
    const char* MQuote = "\"";
private:
    static void print_integer(const CppFormat & c, const cln::cl_I & x);
    static void print_real(const CppFormat & c, const cln::cl_R & x);
    static void print_numeric(const numeric & p, const CppFormat & c, unsigned level);
};

/*-----------------------------------------------------*/
// VE
/*-----------------------------------------------------*/
ex VESimplify(ex expr, int epN = 0, int epsN = 0);
ex VEResult(ex expr);

/*-----------------------------------------------------*/
// ErrMin with HookeJeeves
/*-----------------------------------------------------*/
class ErrMin {
public:
    static int Verbose;
    static IntegratorBase *Integrator;
    static qREAL *paras;
    static dREAL *lambda;
    static dREAL err_max;
    static dREAL err_min;
    static long long MaxRND;
    static long long RunRND;
    static MinimizeBase *miner;
    static dREAL hjRHO;
    static ex lastResErr;
    static dREAL IntError(int nvars, dREAL *las, dREAL *n1, dREAL *n2);
};

/*-----------------------------------------------------*/
// SD Class
/*-----------------------------------------------------*/
class SD {

public:
    static const symbol iEpsilon;
    static const symbol ep;
    static const symbol eps;
    static const symbol vs;
    static const symbol vz;
    static const realsymbol NaN;
    static bool use_dlclose;
    static bool debug;
    
    int ParallelProcess = -1;
    lst ParallelSymbols = lst{ ep, eps, vs, vz, iEpsilon };
    
    int epN = 0;
    int epsN = 0;
    int sN = 0;
    int Verbose = 1;
    int PoleRequested = -100;
    exmap nReplacements;
    vector<pair<lst, lst>> FunExp;
    vector<lst> Deltas;
    vector<ex> Integrands;
    vector<pair<ex, ex>> expResult;
    SecDecBase *SecDec = NULL;
    IntegratorBase *Integrator = NULL;
    MinimizeBase *Minimizer = NULL;
    ex ResultError;
    const char * CFLAGS = "";
    bool IsZero = false;
    bool CheckEnd = true;
    //bool use_CCF = true;
    bool use_ErrMin = false;
    bool use_las = false;
    bool save_las = false;
    bool use_IBF = false;
    bool use_ff = false; // use FindMinimum in F-term
    bool use_exp = false; // use exp in contour deformation
    bool use_MP = true;
    bool use_FT = true;
    int MPDigits = 50; // digits in mpREAL for MP
    lst BisectionPoints = lst { ex(1)/13, ex(1)/19, ex(1)/29, ex(1)/59, ex(1)/41, ex(1)/37, ex(1)/43, ex(1)/53  };
    
    map<int, numeric> Parameter; // used Contours and Integrates, use PL in Prepares part
    
    // used in Contours
    dREAL CTMax = 50;
    int CTTryPTS = 3;
    int CTSavePTS = 5;
    
    long long TryPTS = 500000;
    long long LambdaSplit = 5;
    qREAL LambdaMax = 10;
    int CTry = 2;
    int CTryLeft = 1;
    int CTryRight = 1;
    dREAL CTryRightRatio = 1.5; 
    
    long long RunMAX = 20;
    long long RunPTS = 500000;
    qREAL EpsAbs = 1E-4;
    int ReIm = 3; // 1-Re, 2-Im, 3-ReIm
    
    void Initialize(FeynmanParameter fpi);
    void Initialize(XIntegrand xint);
    void Normalizes();
    void Scalelesses(bool verb=false);
    void SDPrepares();
    void EpsEpExpands();
    void RemoveDeltas();
    void XReOrders();
    void XTogethers();
    void XExpands();
    void KillPowers(bool repeat=true);
    bool IsBad(ex f, vector<exmap> vmap);
    vector<pair<lst, lst>> AutoEnd(pair<lst, lst> po_ex);
    void CIPrepares(const char* key = NULL);
    void Contours(const char * key = NULL, const char *pkey = NULL);
    void Integrates(const char* key = NULL, const char *pkey = NULL, int kid=0);
    void Evaluate(FeynmanParameter fpi, const char *key = NULL);
    void Evaluate(XIntegrand xint, const char *key = NULL);
    
    static ex PrefactorFIESTA(int nLoop);
    ex VEResult();
    void VEPrint(bool endlQ=true);
    double FindMinimum(ex expr, bool compare0 = false);
    static ex Factor(const ex expr);
    static ex PExpand(ex xpol, bool delta=true);
    static int PRank(matrix m);
        
private:
    vector<lst> DS(const pair<lst, lst> po_ex);
    pair<lst, lst> Normalize(const pair<lst, lst> &input);
    static int epRank(ex);
    static int epsRank(ex);
    static int vsRank(ex);
    void DoAsy();

    void CompileMatDet();
    vector<lst> ciResult;
    lst FT_N_XN; // list of { ft, n, xn }
    exmap LambdaMap;
};

/*-----------------------------------------------------*/
// Common SubExpression Parser
/*-----------------------------------------------------*/
class cseParser {
public:
    ex Parse(ex expr, bool reset=true);
    const char* oc = "o";
    int on();
    vector<pair<int, ex>> os();
private:
    map<ex, ex, ex_is_less> ex_var_map;
    int no = 0;
    vector<pair<int, ex>> o_ex_vec;
};

/*-----------------------------------------------------*/
// Customized GiNaC Function
/*-----------------------------------------------------*/
DECLARE_FUNCTION_1P(fabs)
DECLARE_FUNCTION_1P(x)
DECLARE_FUNCTION_1P(y)
DECLARE_FUNCTION_1P(z)
DECLARE_FUNCTION_1P(PL)
DECLARE_FUNCTION_1P(CT)
DECLARE_FUNCTION_2P(FTX)
DECLARE_FUNCTION_2P(VE)
DECLARE_FUNCTION_2P(VEO)

}

