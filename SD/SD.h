/**
 * @file 
 * @brief SecDec header file
 */
 
#pragma once

#include "BASIC.h"

#include <dlfcn.h>

#include <string>
#include <signal.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <sstream>
#include <ios>
#include <regex>

extern "C" {
    #include <quadmath.h>
}

/**
 * @brief namespace for Numerical integration with Sector Decomposition method
 */
namespace HepLib::SD {
    
    using namespace HepLib;
    
    /*-----------------------------------------------------*/
    // Global Variables
    /*-----------------------------------------------------*/
    
    extern const Symbol eps;
    extern const Symbol vs;
    extern const Symbol vz;
    extern const Symbol epz;
    extern const Symbol NaN;
    
    /*-----------------------------------------------------*/
    // Global Functions
    /*-----------------------------------------------------*/
    vector<ex> get_xy_from(ex pol);
    vector<ex> get_x_from(ex pol);
    vector<ex> get_y_from(ex pol);
    vector<ex> get_z_from(ex pol);
    vector<ex> get_pl_from(ex pol);
    int epRank(ex);
    int epsRank(ex);
    int vsRank(ex);
    int x_free_index(ex expr);
    int y_free_index(ex expr);
    ex Factor(const ex expr);
    ex FactorOutX(const ex expr);
    ex exp_simplify(const ex);
    ex pow_simplify(const ex);
    ex xyz_pow_simplify(const ex expr);
    
    /*-----------------------------------------------------*/
    // Customized GiNaC Function
    /*-----------------------------------------------------*/
    DECLARE_FUNCTION_1P(fabs)
    DECLARE_FUNCTION_1P(PL)
    DECLARE_FUNCTION_1P(CT)
    DECLARE_FUNCTION_2P(FTX)
    DECLARE_FUNCTION_2P(VE)
    DECLARE_FUNCTION_2P(VEO)
    DECLARE_FUNCTION_2P(VEO2)
    DECLARE_FUNCTION_1P(epsID)
    DECLARE_FUNCTION_2P(CV) // not used internally, for user use only
    DECLARE_FUNCTION_1P(WRA) // for wick rotation
    extern int VEO_Digits;
    
    /**
     * @brief wrap parameters for loop integrals
     */
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

    /**
     * @brief wrap parameters for generic parameter integrals
     */
    struct XIntegrand {
        lst Functions;
        lst Exponents;
        exmap nReplacements;
        lst Deltas;
        bool isAsy = false;
    };

    /**
     * @brief base class of SecDec
     */
    class SecDecBase {
    public:
        virtual vector<exmap> x2y(const ex &xpol) =0;
        vector<exmap> x2y(const lst &xpols);
        static bool VerifySD(vector<exmap> vmap, bool quick = true);
        static ex XMonomials(const ex & expr);
    };

    /**
     * @brief SecDec by geometric method
     */
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

    /**
     * @brief base for numerical integrator
     */
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
        
        time_t StartTimer; // used internally
        
        long long RunMAX = 100;
        long long RunPTS = 100000;
        long long MinPTS = 100000;
        long long RunTime = 0;
        qREAL EpsAbs = 1E-5;
        qREAL EpsRel = 0;
        int ReIm = 3; // 1-Re, 2-Im, 3-ReIm
        int NANMax = 250;
        int nNAN = 0;
        
        int DQMP = 0;
        int QXDim = 0;
        int MPXDim = 0;
        qREAL QXLimit = 1E-6Q;
        qREAL MPXLimit = 1E-8Q;
        qREAL QFLimit = -1;
        qREAL MPFLimit = -1;
        
        bool UseCpp = true;
        long long NEval = 0;
        int MPDigits = 80;
    };

    /**
     * @brief numerical integrator using HCubature
     */
    class HCubature : public IntegratorBase {
    public:
        static int Wrapper(unsigned int xdim, long long npts, const qREAL *x, void *fdata, unsigned int ydim, qREAL *y);
        
        typedef void (*PrintHookerType) (qREAL*, qREAL*, long long int*, void *);
        
        virtual ex Integrate() override;
        static void DefaultPrintHooker(qREAL*, qREAL*, long long int*, void*);
        PrintHookerType PrintHooker = DefaultPrintHooker;
        long long MaxPTS;
        bool use_last = false;
        long long lastNRUN = 0;
        qREAL LastResult[2];
        qREAL LastAbsErr[2];
        int lastnNAN = 0;
        int LastState = 0;
    };

    /**
     * @brief numerical integrator using CUBA
     */
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


    typedef long double dREAL;
    /**
     * @brief base for class to minimize a function
     */
    class MinimizeBase {
    public:
        typedef dREAL (*FunctionType)(int nvars, dREAL* x, dREAL* pl, dREAL *las);
        virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false, int TryPTS=0, int SavePTS=0) =0;
        dREAL ZeroValue = -1E-20;
        virtual void Minimize(int nvars, FunctionType func, dREAL *ip)=0;
        virtual void ForceStop()=0;
    };

    /**
     * @brief class to minimize a function using HookeJeeves
     */
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

    /**
     * @brief class to minimize a function using MinUit
     */
    class MinUit : public MinimizeBase {
    public:
        virtual dREAL FindMinimum(int nvars, FunctionType func, dREAL *PL = NULL, dREAL *las = NULL, dREAL *UB = NULL, dREAL *LB = NULL, dREAL *IP = NULL, bool compare0 = false, int TryPTS=0, int SavePTS=0) override;
        virtual void Minimize(int nvars, FunctionType func, dREAL *IP) override;
        virtual void ForceStop() override;
    };

    /**
     * @brief class to export GiNaC expression to cpp format
     */
    class CppFormat : public print_csrc_cl_N {
        GINAC_DECLARE_PRINT_CONTEXT(CppFormat, print_csrc_cl_N)
    public:
        CppFormat(ostream &os, const string & s = "L", unsigned opt = 0);
        string suffix;
        string MQuote = "\"";
        
        template<class T> const CppFormat & operator << (const T & v) const {
            s << v;
            return *this;
        };
        const CppFormat & operator << (const basic & v) const;
        const CppFormat & operator << (const ex & v) const;
        const CppFormat & operator << (const lst & v) const;
        const CppFormat & operator<<(std::ostream& (*v)(std::ostream&)) const;
        
        #ifndef DOXYGEN_SKIP
        class _init {
            public: _init();
        };
        #endif
    private:
        #ifndef DOXYGEN_SKIP
        static _init CppFormat_init;
        #endif
        static void print_integer(const CppFormat & c, const cln::cl_I & x);
        static void print_real(const CppFormat & c, const cln::cl_R & x);
        static void print_numeric(const numeric & p, const CppFormat & c, unsigned level);
    };

    /*-----------------------------------------------------*/
    // VE
    /*-----------------------------------------------------*/
    ex VESimplify(ex expr, int epN = 0, int epsN = 0);
    ex VEResult(ex expr);
    ex VEResult2(ex expr); // keep two digits in error
    ex VEMaxErr(ex expr);

    /**
     * @brief ErrMin with HookeJeeves
     */
    class ErrMin {
    public:
        static IntegratorBase *Integrator;
        static qREAL *paras;
        static dREAL err_max;
        static dREAL err_min;
        static long long MaxRND;
        static long long RunRND;
        static MinimizeBase *miner;
        static dREAL *lambda;
        static dREAL hjRHO;
        static ex lastResErr;
        static dREAL IntError(int nvars, dREAL *las, dREAL *n1, dREAL *n2);
    };
    
    /**
     * @brief class to manipulate with Cheng-Wu theorem
     */
    class ChengWu {
    public:
        static bool isProjective(const ex fe, const ex delta);
        static void Projectivize(ex &fe, const ex delta, const ex xsum=0);
        static void Scalelize(ex &fe, const lst xs, const ex cy);
        static void Scalelize(ex &fe, const ex xi, const ex cy);
        static vector<ex> Binarize(ex const fe, ex const eqn);
        static void Binarize(ex const fe, ex const eqn, vector<ex> & ovec);
        static bool isLinearizable(const ex ft, const ex delta, lst & xcs);
        static void Linearize(const lst xcs, ex & fe, ex & ft);
        static bool isPartilizable(const ex ft, const ex delta, lst &xcs, int mode=0);
        static void Partilize(const lst xcs, const lst delta, const ex fe, exvector & ret_lst);
        
        static exvector Evaluate(const ex & fe);
        static exvector WickRotation(const exvector & fe_vec);
        static exvector Apply(const vector<ex> & fe_vec, const ex & ft=0);
        inline static exvector Apply(const ex & fe, const ex & ft=0) {
            exvector fe_vec;
            fe_vec.push_back(fe);
            return Apply(fe_vec, ft);
        }
    };

    /**
     * @brief SecDec the main class to use Sector Decompostion method
     */
    class SecDec {

    public:
        static bool use_dlclose;
        static bool debug;
        static string cpp;
        
        int epN = 0;
        int epsN = 0;
        int vsN = 0;
        int PoleRequested = -5;
        bool vs_before_ep = false;
        exmap nReplacements;
        vector<ex> FunExp; // each item : { {f1,f2,...}, {n1,n2,...}, { delta_list1, delta_list2 } }
        vector<ex> Integrands;
        vector<ex> expResult;
        SecDecBase *SecDec = NULL;
        IntegratorBase *Integrator = NULL;
        MinimizeBase *Minimizer = NULL;
        ex ResultError;
        bool IsZero = false;
        bool CheckEnd = false;
        //bool use_CCF = false;
        bool use_ErrMin = false;
        bool use_las = false;
        bool save_las = false;
        bool use_IBF = false;
        bool use_XReOrders = false;
        int MPDigits = 80; // digits in mpREAL for MP
        lst BisectionPoints = lst { ex(1)/13, ex(1)/19, ex(1)/29, ex(1)/59, ex(1)/41, ex(1)/37, ex(1)/43, ex(1)/53  };
        
        map<int, numeric> Parameter; // used Contours and Integrates, use PL in Prepares part
        
        // used in Contours
        bool CTMaxF = true;
        dREAL CTLaMax = 10; // CTLaMax<0 for explict REAL mode
        int CTTryPTS = 3;
        int CTSavePTS = 3;
        
        long long TryPTS = 500000;
        long long LambdaSplit = 5;
        qREAL IntLaMax = 50;
        int CTry = 1;
        int CTryLeft = 1;
        int CTryRight = 1;
        dREAL CTryRightRatio = 1.5;
        int soLimit = 10000;
        
        long long RunMAX = 20;
        long long RunPTS = 500000;
        map<int, long long> MinPTS;
        qREAL EpsAbs = 1E-4;
        int ReIm = 3; // 1-Re, 2-Im, 3-ReIm
        
        void Initialize(FeynmanParameter fpi);
        void Initialize(XIntegrand xint);
        void Normalizes();
        void Scalelesses();
        void SDPrepares();
        void EpsEpExpands();
        void RemoveDeltas();
        void XReOrders();
        void XTogethers();
        void XExpands();
        void KillPowers(int bits=1+2);
        bool IsBad(ex f, vector<exmap> vmap);
        vector<ex> AutoEnd(ex po_ex);
        void CIPrepares(const string & key = "");
        void Contours(const string & key = "", const string & pkey = "");
        void Integrates(const string & key="", const string & pkey="", int kid=0);
        void ReIntegrates(const string & key, const string & pkey, qREAL err);
        void Evaluate(FeynmanParameter fpi, const string & key = "");
        void Evaluate(XIntegrand xint, const string & key = "");
        void Evaluate(vector<ex> FunExp, const string & key = "");
        void MB();
        void XEnd();
        void ChengWu(const ex & ft=0);
        
        static bool VerifySD(vector<exmap> vmap, bool quick = true);
        static ex XRefined(ex const & ft);
        static lst XRefined_lst(ex const & ft);
        static ex PrefactorFIESTA(int nLoop);
        ex VEResult();
        void VEPrint(bool endlQ=true);
        static ex PExpand(ex xpol, bool delta=true);
        static int PRank(matrix m);
        static ex ContinuousWRA(ex expr_in, int nc=15);
        
        #ifndef DOXYGEN_SKIP
        class _init {
            public: _init();
        };
        #endif
        ~SecDec();
                
    private:
        vector<ex> DS(const ex po_ex);
        lst Normalize(const ex &input);
        void DoAsy();
        bool KillPowerD(ex fe, int kpi);
        bool KillPower(ex fe, int kpi, int bits);

        void CompileMatDet();
        vector<lst> ciResult;
        lst FT_N_XN; // list of { ft, n, xn }
        exmap LambdaMap;
        
        #ifndef DOXYGEN_SKIP
        static _init SD_init;
        #endif
    };

    /**
     * @brief class for  Common SubExpression Parser
     */
    class cseParser {
    public:
        ex Parse(ex expr, bool reset=true);
        string oc = "o";
        int on();
        vector<pair<int, ex>> os();
    private:
        map<ex, ex, ex_is_less> ex_var_map;
        int no = 0;
        vector<pair<int, ex>> o_ex_vec;
    };
    
}


