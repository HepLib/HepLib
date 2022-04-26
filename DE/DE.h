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
        
    class DE {
    public:
        DE(const symbol & x);
        void init(matrix m);
        virtual matrix T();
        virtual void fuchsify();
    protected:
        const symbol & x;
        vector<matrix> Ts; // sequence of T matrix
        matrix Mat; // internal DE matrix
        vector<pair<int,int>> bs; // each block start_index,size
        int N; // DE dimension
    };


}
