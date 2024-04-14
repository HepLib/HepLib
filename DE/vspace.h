/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once
#include "BASIC.h"
#include "IBP.h"
#include "exFlint.h"

namespace HepLib::ED {

    // https://github.com/magv/fuchsia.cpp 
    class matrix_hack : public matrix {
    public:
        void append_rows(const matrix &src);
        exvector &mvec();
        void resize(unsigned nrows);
    };
    
    struct vspace {
        matrix basis;
        vspace(unsigned n);
        vspace(const matrix &b);
        unsigned dim() const;
        matrix basis_col(unsigned i) const;
        matrix basis_row(unsigned i) const;
        const matrix basis_cols() const;
        const matrix &basis_rows() const;
        bool contains(const matrix &v) const;
        void add_rows(const matrix &v);
        void normalize();
    };
    
    vspace nullspace(const matrix &m);
    matrix exnormal(const matrix &m);
    vector<matrix> dual_basis(const matrix &u);
    vector<matrix> dual_basis_left(const matrix &v);
    vector<matrix> eigenvectors(const matrix &m, const ex &eval);
    vector<matrix> eigenvectors_left(const matrix &m, const ex &eval);
    void rescale_submatrix(matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc);
    void echelon_form_gauss(matrix &m);

    matrix matrix_solve_left(const matrix &m, const matrix &vars, const matrix &rhs);
    matrix matrix_cut(const matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc);
    
    //=*********************************************************************=
    
    typedef map<ex, unsigned, ex_is_less> ev_am_t;
    typedef vector<pair<ex,int>> jcf_t;
    ev_am_t ev_am(const matrix & mat);
    pair<matrix,jcf_t> jordan(const matrix &m);
    pair<matrix,jcf_t> jordan(const matrix &m, ev_am_t & ev2am);
    matrix proj_mat(const matrix & a0, const matrix & a1);

}
