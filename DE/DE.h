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
    int prank(const matrix & m);
    matrix a0_mat(const matrix & m, const symbol &x, int pr=17790923);
    pair<matrix,matrix> a01_mat(const matrix & m, const symbol &x, int pr=17790923);
    matrix matrix_diff(const matrix &mat, const symbol &x);
    matrix transform(const matrix &m, const matrix &t, const symbol &x);
    matrix x2y(const matrix & mat, const ex & y, const symbol & x);
    
    typedef vector<vector<map<ex,vector<vector<matrix>>,ex_is_less>>> block_umat_t; // U[a][b][la][k][n]  
    typedef vector<vector<map<ex,vector<matrix>,ex_is_less>>> block_umat_n0_t; // U[a][b][la][k] (n=0)
    typedef vector<map<ex,vector<vector<matrix>>,ex_is_less>> block_imat_t; // I[a][la][k][n]
    typedef vector<map<ex,vector<matrix>,ex_is_less>> block_imat_n0_t; // I[a][la][k] (n=0)
    typedef vector<map<ex,matrix,ex_is_less>> block_imat_kn0_t; // // I[a][la] (k=n=0)
    
    class DE {
    public:        
        DE(const symbol & x);
        void init(const matrix & m);
        virtual vector<vector<matrix>> m2b(const matrix &);
        virtual matrix b2m(vector<vector<matrix>> &);
        virtual matrix d2m(vector<vector<matrix>> &);
        virtual vector<matrix> c2b(const matrix &);
        virtual matrix b2c(vector<matrix> &);
        virtual matrix u2mat(block_umat_t & bu, const ex & x);
        virtual matrix i2mat(block_imat_t & bi, const ex & x);
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
    typedef vector<vector<map<int,vector<vector<fmpq_mat_t>>>>> block_umat_fmpq_mat_t; // U[a][b][ila][k][n]
    typedef vector<map<int,vector<vector<fmpq_mat_t>>>> block_imat_fmpq_mat_t; // I[a][ila][k][n]
    typedef vector<map<int,vector<vector<acb_mat_t>>>> block_imat_acb_mat_t; // I[a][ila][k][n]
    
    class DEX {
    public:        
        DEX(const symbol & x);
        ~DEX();
        void init(const matrix & m);
        virtual vector<vector<matrix>> m2b(const matrix &);
        virtual matrix b2m(vector<vector<matrix>> &);
        virtual matrix d2m(vector<vector<matrix>> &);
        virtual vector<matrix> c2b(const matrix &);
        virtual matrix b2c(vector<matrix> &);
        virtual matrix u2mat(block_umat_t & bu, const ex & x);
        virtual matrix i2mat(block_imat_t & bi, const ex & x);
        virtual matrix T();
        virtual matrix M();
        virtual void fuchsify();
        virtual block_umat_t series(int xn); // U[a][b][la][k][n]
        virtual block_imat_t series(int xn, const matrix & m, slong dp=-1);
        virtual block_imat_t series(int xn, const vector<matrix> & bc, slong dp=-1);
        virtual block_imat_t series(int xn, block_imat_fmpq_mat_t & In0, int nc); // I[a][la][k][n] & In0[a][la][k][0]
        virtual block_imat_t series(int xn, block_imat_acb_mat_t & In0, int nc, slong dp); // I[a][la][k][n] & In0[a][la][k][0]
        void clear();
    //protected:
        bool fuchsified = false;
        int N; // DE dimension
        const symbol & x;
        vector<matrix> Ts; // sequence of T matrix
        vector<pair<int,int>> bs; // each block start_index,size
        
    private:
        block_mat_fmpz_poly_q_t Mat; // DE block matrix
        map<ex,int,ex_is_less> la2i; // convert la to ila
        vector<vector<fmpq_t>> qlas; // ila to la in fmpq_t, the inner vector's length is 1
        vector<ex> las; // ila to la in ex
        vector<vector<map<int,int>>> UK; // UK[a][b][ila]: K+1 for each block U matrix
        vector<map<int,int>> IK; // IK[a][ila]: K+1 for each block I matrix
        block_umat_fmpq_mat_t U0; // U0[a][b][ila][k][0]
    };

    
    class AMF {
    public:
        AMF(IBP & ibp);
        void InitDE();
        
        const symbol & x;
        IBP & ibp;
        lst Rules;
        lst MIntegrals;
        matrix Mat; // original DE matrix
        lst _MIntegrals; // MI @ BC
    };


}
