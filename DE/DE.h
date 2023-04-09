/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once

#include "BASIC.h"
#include "IBP.h"
#include "FlintArb.h"

namespace HepLib {

    typedef map<ex, unsigned, ex_is_less> ev_am_t;
    ev_am_t ev_am(const matrix & mat);
    typedef vector<pair<ex,int>> jcf_t;
    pair<matrix,jcf_t> jordan(const matrix &m);
    pair<matrix,jcf_t> jordan(const matrix &m, ev_am_t & ev2am);
    matrix proj_mat(const pair<matrix,matrix> & a01);
    matrix normal(const matrix & mat);
    matrix exnormal(const matrix & mat);
    int prank(const matrix & mat, const symbol &x);
    ex xpow(const ex & e, const ex & x);
    void xpow(matrix & mat, const ex & x);
    void subs(matrix & mat, const ex & s, unsigned opt);
    matrix a0_mat(const matrix & m, const symbol &x, int pr=17790923);
    pair<matrix,matrix> a01_mat(const matrix & m, const symbol &x, int pr=17790923);
    matrix matrix_diff(const matrix &mat, const symbol &x);
    matrix transform(const matrix &m, const matrix &t, const symbol &x);
    matrix x2y(const matrix & mat, const ex & y, const symbol & x);
    
    typedef vector<vector<map<ex,vector<vector<matrix>>,ex_is_less>>> block_umat_t; // U[a][b][la][k][n]
    typedef map<ex,vector<vector<matrix>>,ex_is_less> umat_t; // U[la][k][n]  
    typedef vector<vector<map<ex,vector<matrix>,ex_is_less>>> block_umat_n0_t; // U[a][b][la][k] (n=0)
    typedef vector<map<ex,vector<vector<matrix>>,ex_is_less>> block_imat_t; // I[a][la][k][n]
    typedef map<ex,vector<vector<matrix>>,ex_is_less> imat_t; // I[la][k][n]
    typedef vector<map<ex,vector<matrix>,ex_is_less>> block_imat_n0_t; // I[a][la][k] (n=0)
    typedef vector<map<ex,matrix,ex_is_less>> block_imat_kn0_t; // // I[a][la] (k=n=0)
    
    class DE {
    public:        
        DE(const symbol & x);
        void init(const matrix & m);
        virtual vector<vector<matrix>> m2b(const matrix &);
        virtual matrix b2m(vector<vector<matrix>> &);
        virtual vector<matrix> c2b(const matrix &);
        virtual matrix b2c(vector<matrix> &);
        map<ex,vector<vector<matrix>>,ex_is_less> b2m(const block_umat_t & bu);
        virtual matrix b2m(block_umat_t & bu, const ex & x);
        virtual matrix b2m(block_imat_t & bi, const ex & x);
        virtual matrix T();
        virtual matrix M();
        virtual void fuchsify();
        virtual block_umat_t series(int xn); // U[a][b][la][k][n]
        virtual block_imat_t series(int xn, block_imat_n0_t & In0); // I[a][la][k][n] & In0[a][la][k]
        virtual block_imat_t series(int xn, const vector<matrix> & bc);
    //protected:
        int N; // DE dimension
        const symbol & x;
        vector<matrix> Ts; // sequence of T matrix
        matrix Mat; // internal DE matrix
        vector<pair<int,int>> bs; // each block start_index,size
        vector<vector<map<ex,int,ex_is_less>>> UK; // UK[a][b][la]: K+1 for each block U matrix
        vector<map<ex,int,ex_is_less>> IK; // IK[a][la]: K+1 for each block I matrix
        block_umat_n0_t U0; // U0[a][b][la][k]: U0[la][k] for each block
    };
        
    typedef vector<vector<vector<vector<fmpz_poly_q_t>>>> block_mat_fmpz_poly_q_t; // M[a][b][r][c]
    typedef vector<vector<vector<vector<fmpz_poly_t>>>> block_mat_fmpz_poly_t; // M[a][b][r][c]
    typedef vector<vector<vector<vector<acb_poly_t>>>> block_mat_acb_poly_t; // M[a][b][r][c]
    typedef vector<vector<map<int,vector<vector<fmpq_mat_t>>>>> block_umat_fmpq_mat_t; // U[a][b][ila][k][n]
    typedef vector<vector<map<int,vector<vector<acb_mat_t>>>>> block_umat_acb_mat_t; // U[a][b][ila][k][n]
    typedef vector<map<int,vector<vector<fmpq_mat_t>>>> block_imat_fmpq_mat_t; // I[a][ila][k][n]
    typedef vector<map<int,vector<vector<acb_mat_t>>>> block_imat_acb_mat_t; // I[a][ila][k][n]
    
    class DEX {
    public:        
        DEX(const symbol & x);
        ~DEX();
        void init(const matrix & m);
        matrix T();
        matrix M();
        void fuchsify();
        void clear();
        
        /* U-Series expansion - rational */
        umat_t series(int xn, const lst & slas={}); // U[la][k][n]
        //matrix series(xn, x0, slas={}); // non-exist due to log(x0)
        
        /* U-Series expansion - acb */
        umat_t series(int xn, slong dp, const lst & slas={}); // U[la][k][n]
        matrix series(int xn, slong dp, const ex & x0, const lst & slas={});
        
        /* I-Series expansion - rational */
        imat_t series(int xn, const matrix & m, const lst & slas={});
        //matrix series(xn, cb, x0, slas={}); // non-exist due to log(x0)
    
        /* I-Series expansion - acb */
        imat_t series(int xn, const matrix & m, slong dp, const lst & slas={});
        matrix series(int xn, const matrix & m, slong dp, const ex & x0, const lst & slas={});
        
        /* Taylor expansion */
        vector<matrix> taylor(int xn, const matrix I0, const ex & x0); // I[n]
        //matrix taylor(xn, I0, x0, dx); // non-exist due to log(x0)
        vector<matrix> taylor(int xn, const matrix I0, const ex & x0, slong dp); // I[n]
        matrix taylor(int xn, const matrix I0, const ex & x0, slong dp, const ex & dx);
        
        /* for Taylor expansion, one should note the inital I0 in which the T matrix has to be taken into considerrations, since there may be a permutation even without calling fuchsify */
        
    private:
        void series(block_umat_fmpq_mat_t & U, int xn, const vector<fmpq_t> & qslas); 
        void series(block_umat_acb_mat_t & U, int xn, slong dp, const vector<fmpq_t> & qslas);
        void series(block_imat_fmpq_mat_t & I, int xn, block_imat_fmpq_mat_t & In0, int nc, const vector<fmpq_t> & qslas); 
        void series(block_imat_acb_mat_t & I, int xn, block_imat_acb_mat_t & In0, int nc, slong dp, const vector<fmpq_t> & qslas); 
        void taylor(vector<vector<fmpq_mat_t>> & I, int xn, const matrix I0, const ex & x0); 
        void taylor(vector<vector<acb_mat_t>> & I, int xn, const matrix I0, const ex & x0, slong dp); 
    public:
        
    //protected:
        bool fuchsified = false;
        bool taylor_inited = false;
        bool ntaylor_inited = false;
        bool taylor_clearq = true; // only save Q/MatTMat, rational Mat will be cleared
        int N; // DE dimension
        slong rel_fp = 100;
        slong abs_fp = 333;
        const symbol & x;
        vector<matrix> Ts; // sequence of T matrix
        vector<pair<int,int>> bs; // each block start_index,size
        
    //private:
        block_mat_fmpz_poly_q_t Mat; // DE block matrix
        map<ex,int,ex_is_less> la2i; // convert la to ila
        vector<vector<fmpq_t>> qlas; // ila to la in fmpq_t, the inner vector's length is 1
        vector<ex> las; // ila to la in ex
        vector<vector<map<int,int>>> UK; // UK[a][b][ila]: K+1 for each block U matrix
        vector<map<int,int>> IK; // IK[a][ila]: K+1 for each block I matrix
        block_umat_fmpq_mat_t U0; // U0[a][b][ila][k][0]
        int nlas = -1;
        int kmmax = -1;
        
        block_mat_fmpz_poly_t QMat; // Mat for taylor
        vector<vector<fmpz_poly_t>> QD; // QD[br][0] : denominator
        block_mat_acb_poly_t TMat; // Mat for ntaylor
        vector<vector<acb_poly_t>> TD; // TD[br][0] : denominator
    };
    matrix u2m(const umat_t & umat, const ex & x);
    matrix i2m(const imat_t & imat, const ex & x);
    matrix i2m(const vector<matrix> & tmat, const ex & x);

    
    class AMF {
    public:
        AMF(IBP & ibp);
        void InitDE();
        void Poles(const ex & rr);
        void ExportDE(const string fn);
        void ImportDE(const string fn);
        
        const symbol & x;
        IBP & ibp;
        lst Rules;
        lst MIntegrals;
        matrix Mat; // original DE matrix
        lst _MIntegrals; // MI @ BC
        lst pts;
        ex d0;
        
        lst Evaluate(int xn, int dp=500);
//        lst FitEps(const lst & eps, int xn, int dp=500, int lp=-1, int nproc=1);
//        lst FitEps(int epn, int xn, int dp, int nproc=1);
        
        static ex Vacuum(int nl, int np);
    };
    matrix PolynomialFit(const exvector & xs, const exvector & ys, unsigned int k, int k0);


}
