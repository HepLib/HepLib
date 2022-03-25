/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once

#include "BASIC.h"
#include "IBP.h"

namespace HepLib {

namespace EoD {

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
    
}

    ex matrix_norm(const matrix & mat, unsigned opt=0);
    ex matrix_den_lcm(const matrix & mat);
    ex xpow(const ex & e, const ex & x);
    void xpow(matrix & mat, const ex & x);
    void subs(matrix & mat, const ex & s, unsigned opt);
    
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
        
    class DE {
    public:
        int WDigits = -1;
        const symbol & x;
        
        DE(const DE & b);
        DE(const symbol & x);
        DE(const matrix & m, const symbol & x);
        DE(const symbol & x, const matrix & m);
        
        void Apply(const matrix & t, bool st=true); // st=true to save t to Ts
        void Apply(const lst & diag, bool st=true); // T is diagnal matrix
        void x2y(const ex & y); // final expression still in x
        void Fuchsify();
        void Shear();
        matrix Series(const ex & x0=0, const unsigned int xn=0, const lst & las={}); // C & C0 matrix
        matrix Taylor(const ex & x0, const ex & dx, const unsigned int xn=0);
        void Normalize();
        void info();
        void xpow();
        void subs(const ex & sub, unsigned opt=0);
        void subs(const exmap & sub, unsigned opt=0);
        void subs(const lst & l, const lst & r, unsigned opt=0);
        matrix MatT();
        void Reset();
        lst EigenValues();
        
    private:
        symbol scn; // cn
        symbol a; // alpha
        // used in Taylor function
        bool taylor_inited = false;
        vector<matrix> tT;
        int ts;
        // used in Series function
        bool series_inited = false;
        vector<vector<vector<vector<matrix>>>> sT; // sT[la][cm][i][j]
        int ss;
        lst las;
        map<ex, unsigned, ex_is_less> laK;
        map<ex, vector<matrix>, ex_is_less> laC0;
    protected:
        vector<matrix> Ts;
        matrix Mat;
    };
    
    class DESS : public DE {
    public:
        DESS(IBP & ibp, const symbol & x);
        void InitDE();
        
    //private:
        const symbol & x;
        IBP & ibp;
        lst Rules;
        lst MIntegrals;
    };
    
    class AMF { // DE @ origin
    public:
        unsigned int xN = 150;
        ex d0 = d;
        int WDigits = -1;
        
        AMF(IBP & ibp);
        void InitDE();
        lst Evaluate();
        lst FitEps(const lst & eps, int lp=0, bool parallel=true);
        lst FitEps(int goal, int order, bool parallel=true);
        
    //private:
        const symbol & x;
        IBP & ibp;
        lst Rules;
        lst MIntegrals;
        lst pts;
        matrix Mat; // original DE matrix
        
        //get iet1 by expansion around regular point iet2
        matrix RU(const ex & iet1, const ex & iet2); 
    };
    
    matrix PolynomialFit(const exvector & xs, const exvector & ys, unsigned int k, int k0=0);

}
