#pragma once

#include "flint/fmpz_poly_q.h"
#include "flint/fmpz_poly_factor.h"
#include "flint/fmpz_poly_mat.h"
#include "flint/acb_mat.h"

#include "BASIC.h"

namespace HepLib {
    
    slong cln_ceiling(const ex & e);
    slong fp2dp(slong fp);
    slong dp2fp(slong dp);
    void ex_to_flint(fmpz_poly_q_t fr, const ex & expr);
    ex flint_to_ex(fmpz_poly_q_t fr, const ex & x);
    ex flint_to_ex(fmpz_poly_t fr, const ex & x);
    ex flint_to_ex(fmpz_t fr, const ex & x);
    ex factor_flint(const ex & expr);
    ex den_lcm_flint(const ex & expr);
    ex normal_flint(const ex & expr);
    
    void ex_to_fmpq(fmpq_t z, const ex & expr);
    void mat_to_fmpq(fmpq_mat_t m, const matrix & mat);
    ex fmpq_to_ex(fmpq_t z);
    matrix fmpq_to_mat(fmpq_mat_t m);
    
    void ex_to_arb(arb_t z, const ex & expr, slong dp);
    void mat_to_arb(arb_mat_t m, const matrix & mat, slong dp);
    ex arb_to_ex(arb_t z, slong dp);
    matrix arb_to_mat(arb_mat_t m, slong dp);
    
    void ex_to_acb(acb_t z, const ex & expr, slong dp);
    void mat_to_acb(acb_mat_t m, const matrix & mat, slong dp);
    ex acb_to_ex(acb_t z, slong dp);
    matrix acb_to_mat(acb_mat_t m, slong dp);
    
    slong get_rel_err(arb_t r);
    slong get_rel_err(acb_t z);
    slong get_rel_err(arb_mat_t m);
    slong get_rel_err(acb_mat_t m);
    
    extern slong error_check_dp;
    bool error_check(arb_t r);
    bool error_check(acb_t r);
    
    class MX { 
        friend class NDE;
    public:
        MX();
        MX(const matrix &m);
        MX(const MX & mx);
        ~MX();
        void init(const matrix & m);
        void init(const MX & mx);
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
        MX & balance(const matrix & P);
        MX & transform(const matrix & t, const matrix & ti);
        MX & shift(const ex & x0);
        matrix operator()(const ex & x);
        int prank();
        int degree();
        matrix coeff(int i, bool evalf=false);
        void lcm(bool init_mat=true); // now M is x*M, Q is its den's lcm
        matrix a0();
        pair<matrix,matrix> a01();
        
        int n=-1; // rows & cols
        int s=-1; // s in DESS
        void Q(arb_t Qm, slong m, slong dp);
        void B(arb_mat_t Bm, slong m, arb_t a, slong dp);
        void Q(acb_t Qm, slong m, slong dp);
        void B(acb_mat_t Bm, slong m, acb_t a, slong dp);
        void Q(fmpq_t Qm, slong m);
        void B(fmpq_mat_t Bm, slong m, fmpq_t a);
    private:
        vector<vector<fmpz_poly_q_t>> M;
        fmpz_poly_t Qx;
        vector<fmpz_mat_t> QxM; // coefficient matrix of Q.x.M
        vector<fmpz_t> Qs; // coefficient of Q
        int _s_ = -1; // _s_>0 need to clear QxM & Qs
    };
        
}
