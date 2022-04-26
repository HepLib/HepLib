#pragma once

#include "flint/fmpz_poly_q.h"
#include "flint/fmpz_poly_factor.h"
#include "flint/fmpz_poly_mat.h"
#include "flint/acb_mat.h"
#include "flint/bool_mat.h"
#include "flint/fmpz_mpoly_factor.h"
#include "flint/fmpq_mpoly_factor.h"

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
    ex _factor_(const ex & x,fmpz_poly_t f);
    ex _factor_(const lst & xs,fmpz_mpoly_t f,fmpz_mpoly_ctx_t ctx);
    ex _factor_(const lst & xs,fmpq_mpoly_t f,fmpq_mpoly_ctx_t ctx);
    ex factor_flint(const ex & e);
    ex factor_fpq(const ex & e);
    ex normal_flint(const ex & expr);
    
    ex den_lcm(const ex & e);
    
}
