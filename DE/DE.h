/**
 * @file 
 * @brief Basic header file
 */
 
#pragma once

#include "BASIC.h"
#include "IBP.h"
#include "FLINT.h"

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
    
    bool is_sheared_form(const matrix & mat);
    pair<matrix, vector<pair<ex,int>>> jordan(const matrix &m);
    pair<matrix, vector<pair<ex,int>>> jordan(const matrix &m, const map<ex, unsigned, ex_is_less> & eval2almul);
}

    class BJF {
    public:
        BJF(matrix _A, ex _b, int K);
        matrix operator()(int i, int j);
    private:
        matrix A;
        matrix B;
        ex b;
    };

    class invBJF {
    public:
        invBJF(matrix _A, ex _b, int K);
        matrix operator()(int i, int j);
    private:
        matrix A;
        matrix U;
        map<int,matrix> Ain;
        map<pair<int,int>,matrix> cache;
        ex b;
    };
        
    class SeriesT { // Series T-Matrix
    public:
        vector<vector<vector<matrix>>> T; // T[cm][i][j] with a = lambda+cn
        int s; // s in DESS
        int kmax;
        exvector las; // lambda list
        map<ex, unsigned, ex_is_less> K; // K_lambda: K[la] 
        map<ex, vector<matrix>, ex_is_less> C0; // C0[la] C_0 matrix for la
        bool inited = false;
        void Resize();
        void Reset();
    };
    
    class TaylorT { // Taylor T-Matrix
    public:
        bool inited = false;
        vector<matrix> T; // T[cm] the single T in Taylor expsion
        int s; // s in DESS
    };
    
    typedef map<ex,vector<vector<matrix>>,ex_is_less> CMatrix; // C[la][k][n] : coefficient of x^la*log(x)^k/k!
    matrix C2Mat(const CMatrix & cmat, const ex & x0);
        
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
        
        exvector las;
        CMatrix Series(const int xN=0); 
        matrix Series(const ex & x0=0, const int xN=0, const lst & las={}); // C matrix
        matrix Taylor(const ex & x0, const ex & dx, const int xN=0);
        void Normalize();
        void info();
        void xpow();
        void subs(const ex & sub, unsigned opt=0);
        void subs(const exmap & sub, unsigned opt=0);
        void subs(const lst & l, const lst & r, unsigned opt=0);
        matrix MatT();
        void Reset();
        
    private:
        TaylorT TT;
        SeriesT ST;
        void STInit();
    protected:
        vector<matrix> Ts;
        matrix Mat;
        symbol a; // alpha
        matrix a0;
        bool fuchsified = false;
        bool sheared = false;
    };
    
    class NDE {
    public:
    
        int WDigits = -1;
        const symbol & x;
        
        NDE(const NDE & b);
        NDE(const symbol & x);
        NDE(const matrix & m, const symbol & x);
        NDE(const symbol & x, const matrix & m);
        
        void Apply(const matrix & t, bool st=true); // st=true to save t to Ts
        void Apply(const lst & diag, bool st=true); // T is diagnal matrix
        void x2y(const ex & y); // final expression still in x
        void Fuchsify();
        void Shear();
        
        exvector las;
        CMatrix Series(const int xN=0, const lst & las={}); 
        matrix Series(const ex & x0=0, const int xN=0, const lst & las={}); // C matrix
        matrix Taylor(const ex & x0, const ex & dx, const int xN=0);
        void info();
        void xpow();
        void subs(const ex & sub, unsigned opt=0);
        void subs(const exmap & sub, unsigned opt=0);
        void subs(const lst & l, const lst & r, unsigned opt=0);
        matrix MatT();
        void Reset();
        
    private:
        TaylorT TT;
        SeriesT ST;
    protected:
        vector<matrix> Ts;
        matrix Mat;
        symbol a; // alpha
        matrix a0;
        MX mx;
        bool fuchsified = false;
        bool sheared = false;
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
        lst NEvaluate();
        lst FitEps(const lst & eps, int lp=0, bool parallel=true);
        lst FitEps(int goal, int order, bool parallel=true);
        static ex Vacuum(int nl, int np);
        
        void ExportDE(const string fn);
        void ImportDE(const string fn);
        
    //private:
        const symbol & x;
        IBP & ibp;
        lst Rules;
        lst MIntegrals;
        lst pts;
        matrix Mat; // original DE matrix
        lst _MIntegrals; // MI @ BC
        
        //get iet1 by expansion around regular point iet2
        matrix RU(const ex & iet1, const ex & iet2, NDE & de); 
        matrix RU(const ex & iet1, const ex & iet2); 
    };
    
    matrix PolynomialFit(const exvector & xs, const exvector & ys, unsigned int k, int k0=0);
    ex matrix_norm(const matrix & mat, unsigned opt=0);
    ex matrix_den_lcm(const matrix & mat);
    ex xpow(const ex & e, const ex & x);
    void xpow(matrix & mat, const ex & x);
    void subs(matrix & mat, const ex & s, unsigned opt);

}
