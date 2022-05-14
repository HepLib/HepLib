#pragma once

#include "flint/fmpz_poly_q.h"
#include "flint/fmpz_poly_factor.h"
#include "flint/fmpz_poly_mat.h"
#include "flint/fmpz_mpoly_factor.h"
#include "flint/fmpq_mpoly_factor.h"
#include "acb_mat.h"
#include "bool_mat.h"

#include "BASIC.h"

namespace HepLib {
    
    slong fp2dp(slong fp);
    slong dp2fp(slong dp);
    //=*********************************************************************=
    void _to_(fmpz_poly_mat_t m, const matrix & mat);
    void _to_(fmpz_poly_q_t f, const ex & e);
    void _to_(fmpz_poly_t f, const ex & e);
    void _to_(fmpz_t f, const ex & e);
    void _to_(fmpz_mat_t m, const matrix & mat);
    void _to_(fmpq_t f, const ex & e);
    void _to_(fmpq_mat_t m, const matrix & mat);
    void _to_(arb_t r, const ex & expr, slong fp);
    void _to_(arb_mat_t m, const matrix & mat, slong fp);
    void _to_(acb_t z, const ex & expr, slong fp);
    void _to_(acb_mat_t m, const matrix & mat, slong fp);
    void _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx, const ex & e);
    void _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx, const ex & e);
    //=*********************************************************************=
    matrix _to_(const ex & x, fmpz_poly_mat_t m);
    ex _to_(const ex & x, fmpz_poly_q_t f);
    ex _to_(const ex & x, fmpz_poly_t f);
    ex _to_(fmpz_t f);
    matrix _to_(fmpz_mat_t m);
    ex _to_(fmpq_t q);
    matrix _to_(fmpq_mat_t m);
    ex _to_(arb_t r, slong fp);
    matrix _to_(arb_mat_t m, slong fp);
    ex _to_(acb_t z, slong fp);
    matrix _to_(acb_mat_t m, slong fp);
    ex _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx);
    ex _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx);
    //=*********************************************************************=
    class MX { 
    public:
        MX();
        MX(const matrix &m);
        MX(const MX & mx);
        MX(const vector<vector<fmpz_poly_q_t>> & M);
        ~MX();
        void init(const matrix & m);
        void init(const MX & mx);
        void init(const vector<vector<fmpz_poly_q_t>> & M);
        void clear();
        MX & add(const MX & mx);
        MX & sub(const MX & mx);
        MX & mul(const MX & mx);
        MX & mul_left(const MX & mx);
        MX & add(const matrix & m);
        MX & sub(const matrix & m);
        MX & mul(const matrix & m);
        MX & mul_left(const matrix & m);
        MX & dx();
        MX & scale(const ex & s);
        MX & scale(fmpz_poly_q_t f);
        MX & scale(fmpz_poly_t f);
        MX & balance(const matrix & P);
        MX & transform(const matrix & t, const matrix & ti);
        MX & shift(const ex & x0);
        matrix operator()(const ex & x);
        void operator()(vector<vector<fmpz_poly_q_t>> & M);
        void operator()(vector<vector<acb_poly_t>> & M, slong fp);
        int prank();
        int degree();
        void series(int xn);
        int denlcm(fmpz_poly_t dl); // return degree of dl and M got updated
        matrix coeff(int i);
        void coeff(fmpq_mat_t m, int i);
        void coeff(acb_mat_t m, int i, slong fp);
        matrix a0(); // if pr<0, we use pr=0
        pair<matrix,matrix> a01(); // if pr<0, we use pr=0
        
        int nr=-1; // rows
        int nc=-1; // cols
        vector<vector<fmpz_poly_q_t>> M;
    };
    //=*********************************************************************=
    ex _factor_(const ex & x,fmpz_poly_t f);
    ex _factor_(const lst & xs,fmpz_mpoly_t f,fmpz_mpoly_ctx_t ctx);
    ex _factor_(const lst & xs,fmpq_mpoly_t f,fmpq_mpoly_ctx_t ctx);
    ex factor_flint(const ex & e);
    ex factor_fpq(const ex & e);
    ex normal_flint(const ex & expr);
    matrix normal_flint(const matrix & mat);
    
    ex den_lcm(const ex & e);
    
}
