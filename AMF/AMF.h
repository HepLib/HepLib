/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once

#include "BASIC.h"
#include "IBP.h"
#include "DE.h"
#include "exFlint.h"

namespace HepLib {

    namespace {
        ex matrix_den_lcm(const matrix & m) {
            auto verb = Verbose;
            Verbose = 0;
            auto den_vec = GiNaC_Parallel(m.nops(), [&m](int idx) {
                auto den = exnormal(m.op(idx)).denom();
                return den;
            });
            Verbose = verb;
            
            symbol s("s"); // since n*(a+b) will expand in GiNaC
            exmap pn_map;
            for(int i=0; i<den_vec.size(); i++) {
                auto den = den_vec[i];
                den = exfactor(s*den);
                if(!is_a<mul>(den)) den = lst{den};
                for(auto item : den) {
                    if(item.info(info_flags::integer)) continue;
                    ex p = item;
                    ex n = 1;
                    if(item.match(pow(w1,w2)) && item.op(1).info(info_flags::integer)) {
                        p = item.op(0);
                        n = item.op(1);
                    }
                    auto itr = pn_map.find(p);
                    if(itr==pn_map.end() || itr->second<n) pn_map[p] = n;
                }
            }
            ex res = 1;
            for(auto kv : pn_map) res *= pow(kv.first, kv.second);
            res = res.subs(s==1);
            
            return res;
        }
        
        bool isComplete(const lst & prop, const lst & lps, const lst & eps) {
            int nl = lps.nops();
            int ne = eps.nops();
            int nle = nl*(nl+1)/2 + nl*ne;
            int n = prop.nops();
            if(nle!=n) return false;
            matrix mat(n, n);
            for(int c=0; c<n; c++) {
                int r = 0;
                auto cc = prop.op(c).expand();
                for(int i=0; i<nl; i++) {
                    mat(r, c) = cc.coeff(lps.op(i), 2);
                    r++;
                    auto ci = cc.coeff(lps.op(i));
                    for(int j=i+1; j<nl; j++) {
                        mat(r, c) = ci.coeff(lps.op(j));
                        r++;
                    }
                    for(int j=0; j<ne; j++) {
                        mat(r, c) = ci.coeff(eps.op(j));
                        r++;
                    }
                    cc = cc.coeff(lps.op(i), 0);
                }
            }
            return mat.rank()==n;
        }
        
        lst express_by(const ex & e, const lst & prop, const lst & lps, const lst & eps) {
            int nl = lps.nops();
            int ne = eps.nops();
            int nle = nl*(nl+1)/2 + nl*ne;
            int n = prop.nops();
            matrix mat(n,n), cv(n,1), rv(n+1,1);
            for(int c=0; c<=n; c++) {
                int r = 0;
                ex cc;
                if(c==n) cc = e.expand();
                else cc = prop.op(c).expand();
                for(int i=0; i<nl; i++) {
                    if(c==n) cv(r,0) = cc.coeff(lps.op(i), 2);
                    else mat(r,c) = cc.coeff(lps.op(i), 2);
                    r++;
                    auto ci = cc.coeff(lps.op(i));
                    for(int j=i+1; j<nl; j++) {
                        if(c==n) cv(r,0) = ci.coeff(lps.op(j));
                        else mat(r,c) = ci.coeff(lps.op(j));
                        r++;
                    }
                    for(int j=0; j<ne; j++) {
                        if(c==n) cv(r,0) = ci.coeff(eps.op(j));
                        else mat(r,c) = ci.coeff(eps.op(j));
                        r++;
                    }
                    cc = cc.coeff(lps.op(i), 0);
                }
                rv(c,0) = cc;
            }
            cv = mat.inverse(solve_algo::gauss).mul(cv);
            lst ret;
            ex rem = 0;
            for(int i=0; i<n; i++) {
                ret.append(cv(i,0));
                rem += rv(i,0) * cv(i,0);
            }
            rem = rv(n,0) - rem;
            ret.append(rem);
            
            return ret;
        }
    }
    
    typedef vector<vector<vector<vector<fmpz_poly_q_struct*>>>> abrc_fmpz_poly_q_t; // M[a][b][r][c]
    typedef vector<vector<vector<vector<fmpz_poly_struct*>>>> abrc_fmpz_poly_t; // M[a][b][r][c]
    typedef vector<vector<map<int,vector<vector<fmpq_mat_struct*>>>>> abikn_fmpq_mat_t; // U[a][b][ila][k][n]
    typedef vector<map<int,vector<vector<fmpq_mat_struct*>>>> aikn_fmpq_mat_t; // I[a][ila][k][n]
        
    typedef vector<vector<vector<vector<gr_poly_struct*>>>> abrc_gr_poly_t; // M[a][b][r][c]
    typedef vector<vector<map<int,vector<vector<gr_mat_struct*>>>>> abikn_gr_mat_t; // U[a][b][ila][k][n]
    typedef vector<map<int,vector<vector<gr_mat_struct*>>>> aikn_gr_mat_t; // I[a][ila][k][n]
    
    class DEX {
    public:
        static int Threads;
        
        vector<ex> las; // ila to la in ex
        
        DEX(const symbol & x, const string & pre = "  ");
        ~DEX();
        void init(const matrix & m);
        matrix T();
        matrix M();
        void fuchsify();
        void clear();
        
        /* U-Series expansion - Rational */
        umat_t series(int xn, const lst & slas={}); // U[la][k][n]
        //matrix series(xn, x0, slas={}); // non-exist due to log(x0) is not rational
                
        /* U-Series expansion - GR */
        //umat_t series(int xn, gr_ctx_ct ctx, const lst & slas={}); // U[la][k][n]
        matrix series(int xn, gr_ptr z0, gr_ctx_t ctx, const lst & slas={});
        
        /* I-Series expansion - Rational */
        imat_t series(int xn, const matrix & m, const lst & slas={}); // imat -> umat.m
        //matrix series(xn, cb, x0, slas={}); // NOT-exist due to log(x0) is not rational
            
        /* I-Series expansion - GR */
        //imat_t series(int xn, const matrix & m, const lst & slas={});
        matrix series(int xn, const matrix & m, gr_ptr z0, gr_ctx_t ctx, const lst & slas={}); // imat -> umat.m
        
        /* Taylor expansion - Rational */
        vector<matrix> taylor(int xn, const matrix I0, const ex & x0); // I[a]
        //matrix taylor(xn, I0, x0, dx); // NOT-exist due to log(x0) is not rational
        
        /* Taylor expansion - GR */
        void taylor(int xn, gr_mat_t imat, gr_ptr z0, gr_ptr dz, gr_ctx_t ctx, const string & es=""); // imat as in and out
        
        /* Note: for Taylor expansion, one should note the initial I0 in which
                 the T matrix has to be taken into considerrations,
                 since there may be a permutation even without calling fuchsify */
        
    private:
        string pre;
        bool fuchsified = false;
        bool taylor_inited = false;
        bool gr_taylor_inited = false;
        bool taylor_clearq = true; // only save QMat, rational Mat will be cleared
        int N; // DE dimension
        const symbol & x;
        vector<matrix> Ts; // sequence of T matrix
        vector<pair<int,int>> bs; // each block start_index,size
        
        abrc_fmpz_poly_q_t Mat; // DE block matrix
        map<ex,int,ex_is_less> la2i; // convert la to ila
        vector<fmpq*> qlas; // ila to la in fmpq_t, the inner vector's length is 1
        vector<vector<map<int,int>>> UK; // UK[a][b][ila]: K+1 for each block U matrix
        vector<map<int,int>> IK; // IK[a][ila]: K+1 for each block I matrix
        abikn_fmpq_mat_t U0; // U0[a][b][ila][k][0]
        int nlas = -1;
        int kmmax = -1;
        
        abrc_fmpz_poly_t QMat; // Mat for taylor
        vector<fmpz_poly_struct*> QD; // QD[br] : denominator
        
        gr_ctx_t _ctx_;
        abrc_gr_poly_t GMat; // Mat for ntaylor
        vector<gr_poly_struct*> GD; // GD[br] : denominator
        
        // U: Rational versions
        void series(abikn_fmpq_mat_t & U, int xn, const vector<fmpq*> & qslas);
        void all_series(abikn_fmpq_mat_t & U, int xn, const vector<fmpq*> & qslas);
        void ab_series(abikn_fmpq_mat_t & U, int xn, const vector<fmpq*> & qslas);
            
        // U: GR version
        void ab_series(abikn_gr_mat_t & U, int xn, gr_ctx_t ctx, const vector<fmpq*> & qslas);
        
        // I: Rational version
        void series(aikn_fmpq_mat_t & I, int xn, aikn_fmpq_mat_t & In0, int nc, const vector<fmpq*> & qslas);
        
        // I: GR versions
        void all_series(aikn_gr_mat_t & I, int xn, aikn_gr_mat_t & In0, int nc, gr_ctx_t ctx, const vector<fmpq*> & qslas);
        void a_series(aikn_gr_mat_t & I, int xn, aikn_gr_mat_t & In0, int nc, gr_ctx_t ctx, const vector<fmpq*> & qslas);
        
        // T: Rational version
        void taylor(vector<vector<fmpq_mat_struct*>> & I, int xn, const matrix I0, const ex & x0);
        
        // T: GR versions
        void an_taylor(vector<vector<gr_mat_struct*>> & I, int xn, gr_mat_t imat, gr_ptr x0, gr_ctx_t ctx, const string & es="");
        void a_taylor(vector<vector<gr_mat_struct*>> & I, int xn, gr_mat_t imat, gr_ptr x0, gr_ctx_t ctx, const string & es="");
    
    };
    
    class AMF {
    public:
        virtual ~AMF() { }
        AMF(const ex & nd) : d0(nd), x(iet) { }
        AMF(const AMF & amf) : x(iet) {
            using_FR = amf.using_FR;
            Propagator = amf.Propagator;
            iPropagator = amf.iPropagator;
            Replacement = amf.Replacement;
            Internal = amf.Internal;
            External = amf.External;
            xn = amf.xn;
            xxn = amf.xxn;
            dp = amf.dp;
            d0 = amf.d0;
            rr = amf.rr;
            T1 = amf.T1;
            T2 = amf.T2;
            LT1 = amf.LT1;
            LT2 = amf.LT2;
            TP = amf.TP;
            LEN = amf.LEN;
            PosPref = amf.PosPref;
            pre = amf.pre;
        }
        AMF &operator=(const AMF &) = delete;
        
        bool using_FR = true;
        lst Integral;
        lst iIntegral; // Integral for iPropagator
        lst NIntegral;
        lst iPropagator; // internal Propagator, maybe different from original Propagator
        lst Propagator;
        lst Replacement;
        lst Internal;
        lst External;
        lst Rules;
        lst TopTopo;
        lst MIntegral;
        lst NMIntegral;
        lst xMIntegral;
        lst NxMIntegral;
        int xn = 200;
        int xxn = 20; // ExterXOrder
        int dp = 200;
        int rr = 2;
        int T1 = CpuCores()/2;
        int T2 = CpuCores()/2;
        int LT1 = 1;
        int LT2 = 1;
        int TP = CpuCores()/2;
        int LEN = 50;
        int PosPref = 1;
        string pre = "  ";
        
        virtual void operator()() { throw Error("Not implemented!"); }
    protected:
        ex d0;
        const symbol & x;
        lst MIxMI;
    };

    class Single : public AMF {
    public:
        Single(const ex & nd) : AMF(nd) { }
        Single(const AMF & amf) : AMF(amf) { }
        Single(const Single & amf) : AMF(amf) { }
        Single &operator=(const Single &) = delete;
                
        void operator()() override;
        static bool isVacuum(int nl, int np);
        static ex Vacuum(int nl, int np);
        
        void UsePosition(int i) { prop_idx = i; }
        
    private:
        bool init_de = false;
        int prop_idx = -1;
        matrix Mat; // original DE matrix
        lst pts;
        
        void InitDE();
        void Poles();
        matrix oo(); // return I at last pts
        void oo2o(matrix & iC); // input I at last pts, update I at first pts
        matrix o(); // return U at first pts
        
        lst apart(const ex & lps, const ex & eps, const lst & expr);
        ex region(const lst & ls, const lst & es, const lst & ps, const lst & ns);
    };
    
    class Loop : public AMF {
    public:
        Loop(const ex & nd) : AMF(nd) { }
        Loop(const AMF & amf) : AMF(amf) { }
        Loop(const Loop & amf) : AMF(amf) { }
        Loop &operator=(const Loop &) = delete;
                
        void operator()() override;
        static bool isVacuum(int nl, int np);
        static ex Vacuum(int nl, int np);
        
    private:
        bool init_de = false;
        ex eLoop;
        matrix Mat; // original DE matrix
        lst pts;
        
        void InitDE();
        void Poles();
        matrix oo(); // return I at last pts
        void oo2o(matrix & iC); // input I at last pts, update I at first pts
        matrix o(); // return U at first pts
        
        lst apart(const ex & lps, const ex & eps, const lst & expr);
        ex region(const lst & ls, const lst & es, const lst & ps, const lst & ns);
    };
    
    class All : public AMF {
    public:
        All(const ex & nd) : AMF(nd) { }
        
        lst ooMIntegral;
        bool using_fBC = false;
        void operator()() override;
        static bool isVacuum(int nl, int np);
        static ex Vacuum(int nl, int np);
        ex Vacuum(int nl, int np, const ex & nd);
        
    private:
        bool init_de = false;
        matrix Mat; // original DE matrix
        lst pts;
        matrix fBC;
        
        void InitDE();
        void Poles();
        matrix oo(); // return I at last pts
        void oo2o(matrix & iC); // input I at last pts, update I at first pts
        matrix o(); // return U at first pts
        
    };
    
    class Fit {
    public:
        enum class Modes { Single, All, Loop };
        lst Integral;
        lst NIntegral;
        lst Propagator;
        lst Replacement;
        lst Internal;
        lst External;
        Modes Mode = Modes::Single;
        
        void operator()(int goal, int order, int rr=2);
        void Parallel(int goal, int order, int rr=2);
            
        lst operator()(int dp, int xn, int rr, lst eps);
        lst Parallel(int dp, int xn, int rr, lst eps);
        
        ex GenerateNumericalConfig(int goal, int order);
        static ex PolyFit(const lst & xs, const lst & ys, int lp);
    };
    

}
