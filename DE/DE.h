/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once

#include "BASIC.h"
#include "IBP.h"
#include "exFlint.h"

namespace HepLib {

    typedef map<ex, unsigned, ex_is_less> ev_am_t;
    ev_am_t ev_am(const matrix & mat);
    typedef vector<pair<ex,int>> jcf_t;
    pair<matrix,jcf_t> jordan(const matrix &m);
    pair<matrix,jcf_t> jordan(const matrix &m, ev_am_t & ev2am);
    matrix proj_mat(const pair<matrix,matrix> & a01);
    matrix normal(const matrix & mat);
    matrix exnormal(const matrix & mat);
    int prank(const matrix & mat, const symbol &x, bool para=false);
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
        virtual ~DE() { }
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
    
    matrix u2m(const umat_t & umat, const ex & x);
    matrix i2m(const imat_t & imat, const ex & x);
    matrix i2m(const vector<matrix> & tmat, const ex & x);

}
