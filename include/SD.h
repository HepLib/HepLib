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

/*********************************************************/
// SD Input
/*********************************************************/
struct FeynmanParameter {
    lst LoopMomenta;
    lst tLoopMomenta;
    lst Propagators;
    lst Exponents;
    exmap lReplacements;
    exmap tReplacements;
    exmap nReplacements;
    ex Prefactor = 1;
};

struct XIntegrand {
    lst Functions;
    lst Exponents;
    exmap nReplacements;
    vector<lst> Deltas;
    ex Prefactor = 1;
};

/*********************************************************/
// Global Functions
/*********************************************************/
vector<ex> get_xy_from(ex pol);
vector<ex> get_x_from(ex pol);
vector<ex> get_y_from(ex pol);
vector<ex> get_pl_from(ex pol);
vector<pair<lst, lst>> diff_wrt(const pair<lst, lst> &input, ex xi);

/*********************************************************/
// SecDec Classes
/*********************************************************/
class SecDecBase {
public:
    virtual vector<exmap> x2y(const ex &xpol) =0;
    vector<exmap> x2y(const lst &xpols);
};

class SecDecG : public SecDecBase {
public:
    vector<exmap> x2y(const ex &xpol) override;
private:
    vector<vector<int>> RunQHull(const matrix &pts);
    vector<matrix> ZeroFaces(const matrix &pts);
    matrix NormalVectors(const vector<matrix> &zfs);
    matrix DualCone(const matrix &pts);
    vector<vector<int>> QHull(const matrix &dc, int dim);
    vector<matrix> SimplexifyR(const matrix &dc, int dim);
    vector<matrix> Simplexify(const matrix &dc, int dim);
    vector<matrix> SimplexCones(matrix pts);
};

class SecDecX : public SecDecBase {
public:
    vector<exmap> x2y(const ex &xpol) override;
private:
    bool PointOver(ex xx, ex yy);
    ex OnlyLowPoints(ex xx);
    ex NormShift(ex xx);
    ex MakeOneStep(ex xx);
    ex MakeOneStep(ex xx, ex facet);
    ex MakeOneStep(ex xx, ex facet, ex graphas);
    ex FindSD(ex xxx);
    ex FindSD(ex xxx, bool extra);
};

/*********************************************************/
// Integrator Classes
/*********************************************************/
typedef __float128 qREAL;
typedef __complex128 qCOMPLEX;

class IntegratorBase {
public:
    typedef int (*SD_Type) (const unsigned int xn, const qREAL xx[], const unsigned int yn, qREAL y[], const qREAL pl[], const qREAL las[]);
    virtual ex Integrate(unsigned int xn, SD_Type fp, SD_Type fpQ, const qREAL pl[], const qREAL las[]) =0;
    
    long long RunMAX = 100;
    long long RunPTS = 100000;
    qREAL EpsAbs = 1E-5;
    qREAL EpsRel = 0;
    int ReIm = 3; // 1-Re, 2-Im, 3-ReIm
    int Verbose = 0;
    bool UseQ = false;
};

class HCubature : public IntegratorBase {
public:
    static bool useQ(unsigned xdim, qREAL const *x);
    static int Wrapper(unsigned int xdim, long long npts, const qREAL *x, void *fdata, unsigned int ydim, qREAL *y);
    
    typedef void (*PrintHookerType) (qREAL*, qREAL*, long long int*, void *);
    
    virtual ex Integrate(unsigned int, SD_Type, SD_Type, const qREAL[], const qREAL[]) override;
    static void DefaultPrintHooker(qREAL*, qREAL*, long long int*, void*);
    PrintHookerType PrintHooker = DefaultPrintHooker;
    long long MaxPTS;
    bool use_last = false;
    qREAL XN = 1.Q;
    
private:
    SD_Type Integrand;
    SD_Type IntegrandQ;
    const qREAL* Lambda;
    const qREAL* Parameter;
    qREAL LastResult[2];
    qREAL LastAbsErr[2];
    int LastState = 0;
};

/*********************************************************/
// Minimize Classes
/*********************************************************/
typedef long double dREAL;
class MinimizeBase {
public:
    typedef dREAL (*FunctionType)(int nvars, dREAL* x, dREAL* pl, dREAL *las);
    virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false) =0;
    dREAL ZeroValue = -1E-20;
    virtual void Minimize(int nvars, FunctionType func, dREAL *ip)=0;
    virtual void ForceStop()=0;
};

class HookeJeeves : public MinimizeBase {
public:
    virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false) override;
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
    virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false) override;
    virtual void Minimize(int nvars, FunctionType func, dREAL *IP) override;
    virtual void ForceStop() override;
};

class ExOpt : public MinimizeBase {
public:
    virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false) override;
    virtual void Minimize(int nvars, FunctionType func, dREAL *IP) override;
    virtual void ForceStop() override;
    bool Compare0 = false;
    dREAL UB[50];
    dREAL LB[50];
    FunctionType ObjectFunction;
    dREAL *PL;
    dREAL *LAS;
};

/*********************************************************/
// CppFormat Class
/*********************************************************/
class CppFormat : public print_csrc_cl_N {
    GINAC_DECLARE_PRINT_CONTEXT(CppFormat, print_csrc_cl_N)
public:
    CppFormat(ostream &os, const char* s = "L", unsigned opt = 0);
    static ex q2ex(qREAL);
    static qREAL ex2q(ex);
    const char* suffix;
    static void QPrint(qREAL qr);
private:
    static void print_integer(const CppFormat & c, const cln::cl_I & x);
    static void print_real(const CppFormat & c, const cln::cl_R & x);
    static void print_numeric(const numeric & p, const CppFormat & c, unsigned level);
};

/*********************************************************/
// VE
/*********************************************************/
ex VESimplify(ex expr, int epN = 0, int epsN = 0);
ex VEResult(ex expr);

/*********************************************************/
// Wrapper for use_cpp = false, i.e. use GiNaC
/*********************************************************/
class GWrapper {
public:
    static void InitMinFunction(ex minF, vector<ex> xs, int ri = 1);
    static dREAL MinFunction(int nvars, dREAL* x, dREAL* pl, dREAL *las);
    
    static void InitIntFunction(ex intF, vector<ex> xs);
    static int IntFunction(const unsigned int xn, const qREAL xx[], const unsigned int yn, qREAL y[], const qREAL pl[], const qREAL las[]);
    
    static ex Lambda;
    
private:
    static ex MinF;
    static vector<ex> Xs;
    static ex IntF;
    static int ReIm;
};

/*********************************************************/
// ILWrapper with HookeJeeves
/*********************************************************/
class ILWrapper {
public:
    static int Verbose;
    static unsigned int xsize;
    static IntegratorBase::SD_Type fp;
    static IntegratorBase::SD_Type fpQ;
    static IntegratorBase *Integrator;
    static qREAL *paras;
    static dREAL *lambda;
    static dREAL err_max;
    static dREAL err_min;
    static long long MaxPTS;
    static long long RunPTS;
    static MinimizeBase *miner;
    static dREAL hjRHO;
    static ex lastResErr;
    static dREAL IntError(int nvars, dREAL *las, dREAL *n1, dREAL *n2);
};

/*********************************************************/
// SD Class
/*********************************************************/
class SD {

public:
    static const symbol iEpsilon;
    static const symbol ep;
    static const symbol eps;
    static const realsymbol NaN;
    static bool use_dlclose;
    static bool debug;
    
    int ParallelProcess = -1;
    lst ParallelSymbols = lst{ ep, eps, iEpsilon };
    
    int epN = 0;
    int epsN = 0;
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
    bool CheckF1 = false;
    bool use_CCF = true;
    bool use_cpp = true;
    bool use_ilwrapper = false;
    lst BisectionPoints = lst { ex(1)/11, ex(1)/13, ex(1)/17, ex(1)/19, ex(1)/23  };
    
    map<int, numeric> Parameter;
    map<int, numeric> ParameterUB;
    map<int, numeric> ParameterLB;
    
    long long TryPTS = 500000;
    long long LambdaSplit = 5;
    qREAL LambdaMax = 100;
    int CTry = 1;
    int CTryLeft = 1;
    int CTryRight = 1;
    dREAL CTryRightRatio = 1.5;
    
    long long RunMAX = 20;
    long long RunPTS = 500000;
    qREAL EpsAbs = 1E-5;
    int ReIm = 3; // 1-Re, 2-Im, 3-ReIm
    
    void Initialize(FeynmanParameter fpi);
    void Initialize(XIntegrand xint);
    void Normalizes();
    void Scalelesses();
    void SDPrepares();
    void EpsEpExpands();
    void RemoveDeltas();
    void XReOrders();
    bool IsBadF1(ex f, vector<exmap> vmap);
    vector<pair<lst, lst>> AutoF1(pair<lst, lst> po_ex);
    void CIPrepares(const char* key = NULL);
    void Contours(const char * key = NULL, const char *pkey = NULL);
    void Integrates(const char* key = NULL, const char *pkey = NULL, int kid=0);
    void Evaluate(FeynmanParameter fpi);
    void Evaluate(XIntegrand xint);
    
    static ex PrefactorFIESTA(int nLoop);
    ex VEResult();
    double FindMinimum(ex expr, bool compare0 = false);
        
private:
    vector<pair<exmap, ex>> SDPrepare(const pair<lst, lst> po_ex);
    pair<lst, lst> Normalize(const pair<lst, lst> &input);
    static int epRank(ex);
    static int epsRank(ex);

    void CompileMatDet();
    vector<lst> ciResult;
    lst FT_N_NX;
    exmap LambdaMap;
};

/*********************************************************/
// Customized GiNaC Function
/*********************************************************/
DECLARE_FUNCTION_1P(fabs)
DECLARE_FUNCTION_1P(x)
DECLARE_FUNCTION_1P(y)
DECLARE_FUNCTION_1P(z)
DECLARE_FUNCTION_1P(t)
DECLARE_FUNCTION_1P(PL)
DECLARE_FUNCTION_1P(CT)
DECLARE_FUNCTION_2P(FTX)
DECLARE_FUNCTION_2P(VE)
DECLARE_FUNCTION_2P(VEO)

}

