#pragma once

#include "BASIC.h"

#include "flint/fmpz_poly_q.h"
#include "flint/fmpz_mpoly_q.h"
#include "flint/fmpq_poly.h"
#include "flint/fmpz_poly_factor.h"
#include "flint/fmpz_poly_mat.h"
#include "flint/fmpz_poly_q.h"
#include "flint/fmpz_mpoly_factor.h"
#include "flint/fmpq_mpoly_factor.h"
#include "flint/fmpz_mat.h"
#include "flint/fmpq_mat.h"
#include "flint/acb_mat.h"
#include "flint/bool_mat.h"
#include "flint/arb_fmpz_poly.h"
#include "flint/acb_poly.h"
#include "flint/acf.h"
#include "flint/gr_mat.h"
#include "flint/gr_special.h"



namespace HepLib {
    
    slong fp2dp(slong fp);
    slong dp2fp(slong dp);
    //=*********************************************************************=
    void _to_(acf_t z, acb_t zb);
    void _to_(acb_t zb, acf_t z);
    void _to_(gr_ptr z, const ex & e, gr_ctx_t ctx); // assume gr is acf
    inline void _to_(gr_ptr z, gr_ctx_t ctx, const ex & e) { _to_(z, e, ctx); }
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
    void _to_(gr_mat_t m, const matrix & mat, gr_ctx_t ctx); // gr is acf
    void _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx, const ex & e);
    void _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx, const ex & e);
    void _to_(const lst & xs, fmpz_mpoly_q_t f, fmpz_mpoly_ctx_t ctx, const ex & e);
    //=*********************************************************************=
    ex _to_(gr_ptr z, gr_ctx_t ctx); // assume gr is acf
    matrix _to_(const ex & x, fmpz_poly_mat_t m);
    ex _to_(const ex & x, fmpz_poly_q_t f);
    ex _to_(const ex & x, fmpz_poly_t f);
    ex _to_(const ex & x, fmpq_poly_t f);
    ex _to_(fmpz_t f);
    matrix _to_(fmpz_mat_t m);
    ex _to_(fmpq_t q);
    matrix _to_(fmpq_mat_t m);
    ex _to_(arb_t r, slong fp);
    matrix _to_(arb_mat_t m, slong fp);
    ex _to_(acb_t z, slong fp);
    matrix _to_(acb_mat_t m, slong fp);
    matrix _to_(gr_mat_t m, gr_ctx_t ctx); // gr is acf
    ex _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx);
    ex _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx);
    ex _to_(const lst & xs, fmpz_mpoly_q_t f, fmpz_mpoly_ctx_t ctx);
    //=*********************************************************************=
    class MX { 
    public:
        MX();
        MX(const matrix &m);
        MX(const MX & mx);
        MX(const vector<vector<fmpz_poly_q_struct*>> & M);
        ~MX();
        void init(const matrix & m);
        void init(const MX & mx);
        void init(const vector<vector<fmpz_poly_q_struct*>> & M);
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
        void operator()(vector<vector<fmpz_poly_struct*>> & M);
        void operator()(vector<vector<fmpz_poly_q_struct*>> & M);
        void operator()(vector<vector<acb_poly_struct*>> & M, slong fp);
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
        vector<vector<fmpz_poly_q_struct*>> M;
    };
    //=*********************************************************************=
    class MQ { 
    public:
        MQ();
        ~MQ();
        void init(const vector<vector<fmpz_poly_q_struct*>> & M);
        void clear();
        void scale(fmpz_poly_q_t f);
        void scale(fmpz_poly_t f);
        int degree();
        int denlcm(fmpz_poly_t dl); // return degree of dl and M got updated
        void coeff(fmpq_mat_t m, int i);
        void coeff(acb_mat_t m, int i, slong fp);
        void coeff(gr_mat_t m, int i, gr_ctx_t ctx);
        void operator()(vector<vector<fmpz_poly_struct*>> & M);
        void operator()(vector<vector<fmpz_poly_q_struct*>> & M);
        void operator()(vector<vector<acb_poly_struct*>> & M, slong fp);
        void operator()(vector<vector<gr_poly_struct*>> & M, gr_ctx_t ctx);
    private:
        int nr=-1; // rows
        int nc=-1; // cols
        vector<vector<fmpz_poly_q_struct*>> M;
    };
    //=*********************************************************************=
    ex _factor_(const ex & x,fmpz_poly_t f);
    ex _factor_(const lst & xs,fmpz_mpoly_t f,fmpz_mpoly_ctx_t ctx);
    ex _factor_(const lst & xs,fmpq_mpoly_t f,fmpq_mpoly_ctx_t ctx);
    //ex factor_flint(const ex & e); // moved to BASIC.h    
    //ex normal_flint(const ex & expr, int opt=o_flint); // moved to BASIC.h
    
    matrix normal_flint(const matrix & mat);
    lst poly_roots(const ex & poly, slong fp);
    
    ex den_lcm(const ex & e);
    
}
