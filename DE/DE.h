/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once

#include "BASIC.h"
#include "IBP.h"

namespace HepLib::DE {

    extern Symbol iet;

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
    matrix normal(const matrix &m);
    vector<matrix> dual_basis(const matrix &u);
    vector<matrix> dual_basis_left(const matrix &v);
    vector<matrix> eigenvectors(const matrix &m, const ex &eval);
    vector<matrix> eigenvectors_left(const matrix &m, const ex &eval);
    void rescale_submatrix(matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc);
    void echelon_form_gauss(matrix &m);

    matrix imatrix(unsigned n);
    map<ex, unsigned, ex_is_less> eigenvalues(const matrix &m);
    
    int prank(const matrix & mat, const symbol &x);
    pair<matrix,matrix> a01_matrix(const matrix &mat, const symbol &x, int pr=19790923);
    matrix a0_matrix(const matrix &mat, const symbol &x, int pr=19790923);
    matrix with_balance_t(const matrix &m, const matrix &P, const symbol &x);
    matrix with_balance_ti(const matrix &m, const matrix &P, const symbol &x);
    
    matrix matrix_solve_left(const matrix &m, const matrix &vars, const matrix &rhs);
    void matrix_map_inplace(matrix &m, std::function<ex(const ex &)> f);
    matrix matrix_map(const matrix &m, std::function<ex(const ex &)> f);
    matrix matrix_cut(const matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc);
    
    // Dx J = M.J & Dx J' = M'.J' with J = T.J' & M' = Ti.M.T - Ti.Dx T
    matrix matrix_diff(const matrix &mat, const symbol &x, const int n=1);
    matrix transform(const matrix &m, const matrix &t, const symbol &x);
    matrix transform(const matrix &m, const matrix &tinverse, const matrix &t, const symbol &x);
    
    bool is_jordan_form(const matrix & mat);
    pair<matrix, vector<pair<ex,int>>> jordan(const matrix &m);
    
    extern int NDigits;
    
    class BJF {
    public:
        BJF(matrix _A, ex _b, int K);
        matrix operator()(int i, int j);
    private:
        matrix A;
        matrix B;
        ex b;
    };

    class BJFinv {
    public:
        BJFinv(matrix _A, ex _b, int K);
        matrix operator()(int i, int j);
    private:
        matrix A;
        matrix U;
        map<int,matrix> Ain;
        map<pair<int,int>,matrix> cache;
        ex b;
    };
        
    class Base {
    public:
        Base(const Base & b);
        Base(const symbol & x);
        Base(const matrix & m, const symbol & x);
        matrix Mat;
        const symbol & x;
        vector<matrix> Ts;
        void Apply(const matrix & t);
        void Apply(const lst & diag); // t is diagnal matrix
        void x2y(const ex & y); // final expression still in x
        void Fuchsify();
        void Shear();
        matrix Series(const ex & x0, const int xn=0);
        void Normalize();
        void info();
        void xpow();
        matrix MatT();
    private:
        
    };
    
    class AMFlow {
    public:
        Base o;
        Base oo;
        AMFlow(IBP::Base & ibp);
        void InitDE();
        matrix FSS(const int xn=0);
        void Scale();
        
    //private:
        IBP::Base & ibp;
        lst Rules;
        lst MIntegrals;
    };

}
