
#include "exFlint.h"
#include <cmath>

namespace HepLib {
    
    namespace {
        inline lst syms(const ex & e) {
            exset ss;
            for(const_preorder_iterator i=e.preorder_begin(); i!=e.preorder_end(); ++i) 
                if(is_a<symbol>(*i)) ss.insert(*i); 
            lst ls;
            for(auto item : ss) ls.append(item);
            return ls;
        }
        
    }
    
    //=*********************************************************************=
    
    slong fp2dp(slong fp) {
        static long double r = 0.30102999566398119521373889472449302676818988146211L;
        return ceil(r*fp);
    }
        
    slong dp2fp(slong dp) {
        static long double r = 3.3219280948873623478703194294893901758648313930246L;
        return ceil(r*dp);
    }
    
    void _to_(acf_t z, acb_t zb) {
        arf_set(acf_realref(z), arb_midref(acb_realref(zb)));
        arf_set(acf_imagref(z), arb_midref(acb_imagref(zb)));
    }

    void _to_(acb_t zb, acf_t z) {
        arb_t rr, ii;
        arb_init(rr);
        arb_init(ii);
        arb_set_arf(rr, acf_realref(z));
        arb_set_arf(ii, acf_imagref(z));
        acb_set_arb_arb(zb, rr, ii);
        arb_clear(rr);
        arb_clear(ii);
    }
    
    void _to_(gr_ptr z, const ex & expr, gr_ctx_t ctx) { // assume z is acf_t
        int status = GR_SUCCESS;
        if(is_a<numeric>(expr)) {
            auto n = ex_to<numeric>(expr);
            if(n.is_integer()) {
                fmpz_t q;
                fmpz_init(q);
                _to_(q, n);
                status |= gr_set_fmpz(z, q, ctx);
                fmpz_clear(q);
            } else if(n.is_rational()) {
                fmpq_t q;
                fmpq_init(q);
                _to_(q, n);
                status |= gr_set_fmpq(z, q, ctx);
                fmpq_clear(q);
            } else if(n.is_cinteger()) {
                fmpz_t q1, q2;
                fmpz_init(q1);
                fmpz_init(q2);
                _to_(q1, real(n));
                _to_(q2, imag(n));
                acf_struct* zz = (acf_struct*)z;
                arf_set_fmpz(acf_realref(zz), q1);
                arf_set_fmpz(acf_imagref(zz), q2);
                fmpz_clear(q1);
                fmpz_clear(q2);
            } else { // back to acb version
                slong fp;
                status |= gr_ctx_get_real_prec(&fp, ctx);
                acb_t zb;
                acb_init(zb);
                _to_(zb, expr, fp);
                _to_((acf_struct*)z, zb);
                acb_clear(zb);
            }
        } else if(expr.match(sin(w))) { // sin
            _to_(z, expr.op(0), ctx);
            status |= gr_sin(z, z, ctx);
        } else if(expr.match(cos(w))) { // cos
            _to_(z, expr.op(0), ctx);
            status |= gr_cos(z, z, ctx);
        } else if(expr.match(tan(w))) { // tan
            _to_(z, expr.op(0), ctx);
            status |= gr_tan(z, z, ctx);
        } else if(expr.match(asin(w))) { // asin
            _to_(z, expr.op(0), ctx);
            status |= gr_asin(z, z, ctx);
        } else if(expr.match(acos(w))) { // acos
            _to_(z, expr.op(0), ctx);
            status |= gr_acos(z, z, ctx);
        } else if(expr.match(sin(w))) { // atan
            _to_(z, expr.op(0), ctx);
            status |= gr_atan(z, z, ctx);
        } else if(expr.match(exp(w))) { // exp
            _to_(z, expr.op(0), ctx);
            status |= gr_exp(z, z, ctx);
        } else if(expr.match(log(w))) { // log
            _to_(z, expr.op(0), ctx);
            status |= gr_log(z, z, ctx);
        } else if(expr.match(sqrt(w))) { // sqrt
            _to_(z, expr.op(0), ctx);
            status |= gr_sqrt(z, z, ctx);
        } else if(expr.match(tgamma(w))) { // gamma
            // it seems gr NOT supported yet
            //_to_(z, expr.op(0), ctx);
            //status |= gr_gamma(z, z, ctx);
            slong fp;
            status |= gr_ctx_get_real_prec(&fp, ctx);
            acb_t zb;
            acb_init(zb);
            _to_(zb, expr.op(0), 2*fp);
            acb_gamma(zb, zb, 2*fp);
            _to_((acf_struct*)z, zb);
            acb_clear(zb);
        } else if(is_a<power>(expr) && expr.op(1).info(info_flags::integer)) {
            slong n = ex_to<numeric>(expr.op(1)).to_int();
            _to_(z, expr.op(0), ctx);
            status |= gr_pow_si(z, z, n, ctx);
        } else if(is_a<add>(expr) || is_a<mul>(expr)) {
            bool isa = is_a<add>(expr);
            if(isa) status |= gr_zero(z, ctx);
            else status |= gr_one(z, ctx);
            gr_ptr t = gr_heap_init(ctx);
            for(auto const & item : expr) {
                _to_(t, item, ctx);
                if(isa) status |= gr_add(z, z, t, ctx);
                else status |= gr_mul(z, z, t, ctx);
            }
            gr_heap_clear(t, ctx);
        } else if(expr==Pi) {
            status |= gr_pi(z, ctx);
        } else if(expr==Euler) {
            status |= gr_euler(z, ctx);
        } else {
            cout << "input: " << expr << endl;
            throw Error("_to_gr Not Supported Yet!");
        }
        if(status != GR_SUCCESS) {
            cout << "input: " << expr << endl;
            throw Error("_to_: status != GR_SUCCESS!");
        }
    }
    
    ex _to_(gr_ptr z, gr_ctx_t ctx) { // assume gr is a acf, use acb version
        int status = GR_SUCCESS;
        slong fp;
        status |= gr_ctx_get_real_prec(&fp, ctx);
        acb_t zb;
        acb_init(zb);
        _to_(zb, (acf_struct*)z);
        ex res = _to_(zb, fp);
        acb_clear(zb);
        if(status != GR_SUCCESS) throw Error("_to_: status != GR_SUCCESS!");
        return res;
    }
    
    void _to_(fmpz_poly_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(fmpz_poly_mat_entry(m,r,c), mat(r,c));
        }
    }
    
    void _to_(fmpz_poly_q_t f, const ex & e) {
        if(syms(e).nops()>1) {
            cout << endl << e << endl;
            throw Error("_to_(fmpz_poly_q_t, ex): >=2 variables found.");
        }
        if(is_a<symbol>(e)) {
            if(fmpz_poly_q_set_str(f, "2  0 1/1  1")) throw Error("_to_, fmpz_poly_q_set_str failed!");
            return;
        } else if(e.info(info_flags::rational)) {
            auto nd = e.numer_denom();
            ostringstream oss;
            oss << "1  " << nd.op(0) << "/1  " << nd.op(1);
            string str = oss.str();
            if(fmpz_poly_q_set_str(f, str.c_str())) {
                cout << "str: " << oss.str() << endl;
                throw Error("_to_(fmpz_poly_q_t, ex) failed.");
            }
            return;
        } else if(is_a<add>(e)) {
            fmpz_poly_q_zero(f);
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            for(auto item : e) {
                _to_(fi, item);
                fmpz_poly_q_add(f, f, fi);
            }
            fmpz_poly_q_clear(fi);
            return;
        } else if(is_a<mul>(e)) {
            fmpz_poly_q_one(f);
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            for(auto item : e) {
                _to_(fi, item);
                fmpz_poly_q_mul(f, f, fi);
            }
            fmpz_poly_q_clear(fi);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::posint)) {
            ulong n = ex_to<numeric>(e.op(1)).to_int();
            _to_(f, e.op(0));
            fmpz_poly_q_pow(f, f, n);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::negint)) {
            ulong n = -ex_to<numeric>(e.op(1)).to_int();
            _to_(f, e.op(0));
            fmpz_poly_q_pow(f, f, n);
            fmpz_poly_q_inv(f,f);
            return;
        } else {
            cout << "expr = " << e << endl;
            throw Error("_to_(fmpz_poly_q_t, ex) NOT supported region.");
        }
    }
    
    void _to_(fmpz_poly_t f, const ex & e) {
        if(syms(e).nops()>1) throw Error("_to_(fmpz_poly_t, ex): >=2 variables found.");
        if(is_a<symbol>(e)) {
            if(fmpz_poly_set_str(f, "2  0 1")) throw Error("_to_: fmpz_poly_set_str failed.");
            return;
        } else if(e.info(info_flags::integer)) {
            ostringstream oss;
            oss << "1  " << e;
            string str = oss.str();
            if(fmpz_poly_set_str(f, str.c_str())) {
                cout << "str: " << oss.str() << endl;
                throw Error("fmpz_poly_set_str error.");
            }
            return;
        } else if(is_a<add>(e)) {
            fmpz_poly_zero(f);
            fmpz_poly_t fi;
            fmpz_poly_init(fi);
            for(auto item : e) {
                _to_(fi, item);
                fmpz_poly_add(f, f, fi);
            }
            fmpz_poly_clear(fi);
            return;
        } else if(is_a<mul>(e)) {
            fmpz_poly_one(f);
            fmpz_poly_t fi;
            fmpz_poly_init(fi);
            for(auto item : e) {
                _to_(fi, item);
                fmpz_poly_mul(f, f, fi);
            }
            fmpz_poly_clear(fi);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::posint)) {
            ulong n = ex_to<numeric>(e.op(1)).to_int();
            _to_(f, e.op(0));
            fmpz_poly_pow(f, f, n);
            return;
        } else {
            cout << "expr = " << e << endl;
            throw Error("_to_(fmpz_poly_t, ex) NOT supported region");
        }
    }
    
    void _to_(fmpz_t q, const ex & expr) {
        auto str = ex2str(expr);
        if(fmpz_set_str(q, str.c_str(), 10)) {
            cout << "str: " << str << endl;
            throw Error("fmpz_set_str error.");
        }
    }
    
    void _to_(fmpq_t q, const ex & expr) {
        auto str = ex2str(expr);
        if(fmpq_set_str(q, str.c_str(), 10)) {
            cout << "str: " << str << endl;
            throw Error("fmpq_set_str error.");
        }
    }
    
    void _to_(fmpz_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(fmpz_mat_entry(m,r,c), mat(r,c));
        }
    }
    
    void _to_(fmpq_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(fmpq_mat_entry(m,r,c), mat(r,c));
        }
    }
    
    void _to_v0_(arb_t r, const ex & expr, slong fp) {
        set_precision(fp2dp(fp)+10);
        auto ne = expr.evalf();
        reset_precision();
        if(!is_a<numeric>(ne)) {
            cout << endl << "ne = " << ne << endl;
            throw Error("_to_(arb_t, ex, slong): NOT a number");
        }
        auto str = ex2str(ne);
        if(arb_set_str(r, str.c_str(), fp)) {
            cout << "str: " << str << endl;
            throw Error("arb_set_str error.");
        }
    }
    
    void _to_(arb_t z, const ex & expr, slong fp) { // only work for rational expression
        if(is_a<numeric>(expr)) {
            auto n = ex_to<numeric>(expr);
            if(n.is_integer()) {
                fmpz_t q;
                fmpz_init(q);
                _to_(q, n);
                arb_set_fmpz(z, q);
                fmpz_clear(q);
            } else if(n.is_rational()) {
                fmpq_t q;
                fmpq_init(q);
                _to_(q, n);
                arb_set_fmpq(z, q, fp);
                fmpq_clear(q);
            } else {
                string str = ex2str(n);
                if(arb_set_str(z, str.c_str(), fp)) {
                    cout << "str: " << str << endl;
                    throw Error("arb_set_str error.");
                }
            }
        } else if(expr.match(sin(w))) { // sin
            _to_(z, expr.op(0), fp);
            arb_sin(z, z, fp);
        } else if(expr.match(cos(w))) { // cos
            _to_(z, expr.op(0), fp);
            arb_cos(z, z, fp);
        } else if(expr.match(tan(w))) { // tan
            _to_(z, expr.op(0), fp);
            arb_tan(z, z, fp);
        } else if(expr.match(asin(w))) { // asin
            _to_(z, expr.op(0), fp);
            arb_asin(z, z, fp);
        } else if(expr.match(acos(w))) { // acos
            _to_(z, expr.op(0), fp);
            arb_acos(z, z, fp);
        } else if(expr.match(sin(w))) { // atan
            _to_(z, expr.op(0), fp);
            arb_atan(z, z, fp);
        } else if(expr.match(exp(w))) { // exp
            _to_(z, expr.op(0), fp);
            arb_exp(z, z, fp);
        } else if(expr.match(log(w))) { // log
            _to_(z, expr.op(0), fp);
            arb_log(z, z, fp);
        } else if(expr.match(sqrt(w))) { // sqrt
            _to_(z, expr.op(0), fp);
            arb_sqrt(z, z, fp);
        } else if(expr.match(tgamma(w))) { // gamma
            _to_(z, expr.op(0), fp);
            arb_gamma(z, z, fp);
        } else if(is_a<add>(expr) || is_a<mul>(expr)) {
            bool isa = is_a<add>(expr);
            if(isa) arb_zero(z);
            else arb_one(z);
            arb_t t;
            arb_init(t);
            for(auto const & item : expr) {
                _to_(t, item, fp);
                if(isa) arb_add(z, z, t, fp);
                else arb_mul(z, z, t, fp);
            }
            arb_clear(t);
        } else {
            cout << "input: " << expr << endl;
            throw Error("_to_arb Not Supported Yet!");
        }
    }
    
    void _to_(arb_mat_t m, const matrix & mat, slong fp) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(arb_mat_entry(m,r,c), mat(r,c), fp);
        }
    }
    
    void _to_v0_(acb_t z, const ex & expr, slong fp) { // old version
        set_precision(fp2dp(fp)+10);
        auto ne = expr.evalf();
        reset_precision();
        if(!is_a<numeric>(ne)) {
            cout << endl << "ne = " << ne << endl;
            throw Error("_to_(acb_t, ex, slong): NOT a number");
        }
        
        arb_t rb, ib;
        arb_init(rb);
        arb_init(ib);
        auto nre = ne.real_part();
        auto nie = ne.imag_part();
        auto rstr = ex2str(nre);
        if(arb_set_str(rb, rstr.c_str(), fp)) {
            cout << "str: " << rstr << endl;
            throw Error("arb_set_str error.");
        }
        auto istr = ex2str(nie);
        if(arb_set_str(ib, istr.c_str(), fp)) {
            cout << "str: " << istr << endl;
            throw Error("arb_set_str error.");
        }
        acb_set_arb_arb(z, rb, ib);
        arb_clear(rb);
        arb_clear(ib);
    }
    
    void _to_(acb_t z, const ex & expr, slong fp) { // only work for rational expression
        if(is_a<numeric>(expr)) {
            auto n = ex_to<numeric>(expr);
            if(n.is_integer()) {
                fmpz_t q;
                fmpz_init(q);
                _to_(q, n);
                acb_set_fmpz(z, q);
                fmpz_clear(q);
            } else if(n.is_rational()) {
                fmpq_t q;
                fmpq_init(q);
                _to_(q, n);
                acb_set_fmpq(z, q, fp);
                fmpq_clear(q);
            } else if(n.is_cinteger()) {
                fmpz_t q1, q2;
                fmpz_init(q1);
                fmpz_init(q2);
                _to_(q1, real(n));
                _to_(q2, imag(n));
                acb_set_fmpz_fmpz(z, q1, q2);
                fmpz_clear(q1);
                fmpz_clear(q2);
            } else if(n.is_crational()) {
                fmpq_t q1, q2;
                fmpq_init(q1);
                fmpq_init(q2);
                _to_(q1, real(n));
                _to_(q2, imag(n));
                arb_set_fmpq(acb_realref(z), q1, fp);
                arb_set_fmpq(acb_imagref(z), q2, fp);
                fmpq_clear(q1);
                fmpq_clear(q2);
            } else {
                arb_t rb, ib;
                arb_init(rb);
                arb_init(ib);
                _to_(rb, real(n), fp);
                _to_(ib, imag(n), fp);
                acb_set_arb_arb(z, rb, ib);
                arb_clear(rb);
                arb_clear(ib);
            }
        } else if(expr.match(sin(w))) { // sin
            _to_(z, expr.op(0), fp);
            acb_sin(z, z, fp);
        } else if(expr.match(cos(w))) { // cos
            _to_(z, expr.op(0), fp);
            acb_cos(z, z, fp);
        } else if(expr.match(tan(w))) { // tan
            _to_(z, expr.op(0), fp);
            acb_tan(z, z, fp);
        } else if(expr.match(asin(w))) { // asin
            _to_(z, expr.op(0), fp);
            acb_asin(z, z, fp);
        } else if(expr.match(acos(w))) { // acos
            _to_(z, expr.op(0), fp);
            acb_acos(z, z, fp);
        } else if(expr.match(sin(w))) { // atan
            _to_(z, expr.op(0), fp);
            acb_atan(z, z, fp);
        } else if(expr.match(exp(w))) { // exp
            _to_(z, expr.op(0), fp);
            acb_exp(z, z, fp);
        } else if(expr.match(log(w))) { // log
            _to_(z, expr.op(0), fp);
            acb_log(z, z, fp);
        } else if(expr.match(sqrt(w))) { // sqrt
            _to_(z, expr.op(0), fp);
            acb_sqrt(z, z, fp);
        } else if(expr.match(tgamma(w))) { // gamma
            _to_(z, expr.op(0), fp);
            acb_gamma(z, z, fp);
        } else if(is_a<add>(expr) || is_a<mul>(expr)) {
            bool isa = is_a<add>(expr);
            if(isa) acb_zero(z);
            else acb_one(z);
            acb_t t;
            acb_init(t);
            for(auto const & item : expr) {
                _to_(t, item, fp);
                if(isa) acb_add(z, z, t, fp);
                else acb_mul(z, z, t, fp);
            }
            acb_clear(t);
        } else {
            cout << "input: " << expr << endl;
            throw Error("_to_acb Not Supported Yet!");
        }
    }
    
    void _to_(acb_mat_t m, const matrix & mat, slong fp) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(acb_mat_entry(m,r,c), mat(r,c), fp);
        }
    }
    
    void _to_(gr_mat_t m, const matrix & mat, gr_ctx_t ctx) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(gr_mat_entry_ptr(m,r,c,ctx), mat(r,c), ctx);
        }
    }
    
    void _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx, const ex & e) {
        string vars[xs.nops()];
        const char* cvars[xs.nops()];
        for(int i=0; i<xs.nops(); i++) {
            vars[i] = ex2str(xs.op(i));
            cvars[i] = vars[i].c_str();
        }
        fmpz_mpoly_ctx_init(ctx, xs.nops(), ORD_LEX);
        fmpz_mpoly_init(f, ctx);
        string str = ex2str(e);
        if(fmpz_mpoly_set_str_pretty(f, str.c_str(), cvars, ctx)) {
            cout << "str: " << e << endl;
            throw Error("fmpz_mpoly_set_str_pretty error.");
        }
        // call fmpz_mpoly_ctx_clear() and fmpz_mpoly_clear() in user side
    }
    
    void _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx, const ex & e) {
        string vars[xs.nops()];
        const char* cvars[xs.nops()];
        for(int i=0; i<xs.nops(); i++) {
            vars[i] = ex2str(xs.op(i));
            cvars[i] = vars[i].c_str();
        }
        fmpq_mpoly_ctx_init(ctx, xs.nops(), ORD_LEX);
        fmpq_mpoly_init(f, ctx);
        string str = ex2str(e);
        if(fmpq_mpoly_set_str_pretty(f, str.c_str(), cvars, ctx)) {
            cout << "str: " << e << endl;
            throw Error("fmpq_mpoly_set_str_pretty error.");
        }
        // call fmpz_mpoly_ctx_clear() and fmpz_mpoly_clear() in user side
    }
        
    //=*********************************************************************=
    
    matrix _to_(const ex & x, fmpz_poly_mat_t m) {
        auto nr = fmpz_poly_mat_nrows(m);
        auto nc = fmpz_poly_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(x, fmpz_poly_mat_entry(m,r,c));
        }
        return mat;
    }
    
    ex _to_(const ex & x, fmpz_poly_q_t f) {
        ex res_den = 0;
        auto den = fmpz_poly_q_denref(f);
        auto xdn = fmpz_poly_length(den);
        fmpz_t z;
        fmpz_init(z);
        for(slong i=0; i<xdn; i++) {
            fmpz_poly_get_coeff_fmpz(z,den,i);
            auto str = fmpz_get_str(NULL,10,z);
            res_den += numeric(str) * pow(x,i);
            flint_free(str);
        }
        ex res_num = 0;
        auto num = fmpz_poly_q_numref(f);
        auto xnn = fmpz_poly_length(num);
        for(slong i=0; i<xnn; i++) {
            fmpz_poly_get_coeff_fmpz(z,num,i);
            auto str = fmpz_get_str(NULL,10,z);
            if(xdn==1) res_num += numeric(str)/res_den * pow(x,i);
            else res_num += numeric(str) * pow(x,i);
            flint_free(str);
        }
        fmpz_clear(z);
        
        if(xdn==1) return res_num;
        return res_num/res_den;
    }

    ex _to_(const ex & x, fmpz_poly_t f) {
        ex res = 0;
        auto xn = fmpz_poly_length(f);
        fmpz_t z;
        fmpz_init(z);
        for(slong i=0; i<xn; i++) {
            fmpz_poly_get_coeff_fmpz(z,f,i);
            auto str = fmpz_get_str(NULL,10,z);
            res += numeric(str) * pow(x,i);
            flint_free(str);
        }
        fmpz_clear(z);
        return res;
    }
    
    ex _to_(const ex & x, fmpq_poly_t f) {
        ex res = 0;
        auto xn = fmpq_poly_length(f);
        fmpq_t z;
        fmpq_init(z);
        for(slong i=0; i<xn; i++) {
            fmpq_poly_get_coeff_fmpq(z,f,i);
            auto str = fmpq_get_str(NULL,10,z);
            res += numeric(str) * pow(x,i);
            flint_free(str);
        }
        fmpq_clear(z);
        return res;
    }
    
    ex _to_(fmpz_t f) {
        auto str = fmpz_get_str(NULL,10,f);
        auto res = numeric(str);
        flint_free(str);
        return res;
    }
    
    matrix _to_(fmpz_mat_t m) {
        auto nr = fmpz_mat_nrows(m);
        auto nc = fmpz_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(fmpz_mat_entry(m,r,c));
        }
        return mat;
    }
    
    ex _to_(fmpq_t q) {        
        ex n = _to_(fmpq_numref(q));
        ex d = _to_(fmpq_denref(q));
        ex res = n/d;
        return res;
    }
    
    matrix _to_(fmpq_mat_t m) {
        auto nr = fmpq_mat_nrows(m);
        auto nc = fmpq_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(fmpq_mat_entry(m,r,c));
        }
        return mat;
    }
    
    ex _to_(arb_t r, slong fp) {
        auto dp = fp2dp(fp);
        auto cstr = arb_get_str(r,fp,ARB_STR_NO_RADIUS);
        string str = cstr;
        flint_free(cstr);
        if(string_contain(str,".") || string_contain(str,"e")) str += "_"+to_string(dp);
        ex n = numeric(cln::cl_R(str.c_str()));
        return n;
    }
    
    matrix _to_(arb_mat_t m, slong dp) {
        auto nr = arb_mat_nrows(m);
        auto nc = arb_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(arb_mat_entry(m,r,c), dp);
        }
        return mat;
    }
    
    ex _to_(acb_t z, slong fp) {
        auto dp = fp2dp(fp);
        arb_t ri;
        arb_init(ri);
        acb_get_real(ri,z);
        auto str = arb_get_str(ri,fp,ARB_STR_NO_RADIUS);
        string rstr = str;
        flint_free(str);
        if(string_contain(rstr,".") || string_contain(rstr,"e")) rstr += "_"+to_string(dp);
        acb_get_imag(ri,z);
        str = arb_get_str(ri,fp,ARB_STR_NO_RADIUS);
        string istr = str;
        flint_free(str);
        if(string_contain(istr,".") || string_contain(istr,"e")) istr += "_"+to_string(dp);
        arb_clear(ri);
        ex n = numeric(cln::complex(cln::cl_R(rstr.c_str()),cln::cl_R(istr.c_str())));
        return n;
    }
    
    matrix _to_(acb_mat_t m, slong fp) {
        auto nr = acb_mat_nrows(m);
        auto nc = acb_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(acb_mat_entry(m,r,c), fp);
        }
        return mat;
    }
    
    matrix _to_(gr_mat_t m, gr_ctx_t ctx) {
        auto nr = acb_mat_nrows(m);
        auto nc = acb_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(gr_mat_entry_ptr(m,r,c,ctx), ctx);
        }
        return mat;
    }
    
    ex _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx) {
        symtab st;
        string vars[xs.nops()];
        const char* cvars[xs.nops()];
        for(int i=0; i<xs.nops(); i++) {
            vars[i] = ex2str(xs.op(i));
            cvars[i] = vars[i].c_str();
            st[cvars[i]] = xs.op(i);
        }
        auto cstr = fmpz_mpoly_get_str_pretty(f, cvars, ctx);
        string str(cstr);
        flint_free(cstr);
        return str2ex(str,st);
    }
    
    ex _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx) {
        symtab st;
        string vars[xs.nops()];
        const char* cvars[xs.nops()];
        for(int i=0; i<xs.nops(); i++) {
            vars[i] = ex2str(xs.op(i));
            cvars[i] = vars[i].c_str();
            st[cvars[i]] = xs.op(i);
        }
        auto cstr = fmpq_mpoly_get_str_pretty(f, cvars, ctx);
        string str(cstr);
        flint_free(cstr);
        return str2ex(str,st);
    }
    
    ex factor_flint(const ex & expr_in, bool nd) {
        if(nd) {
            exmap map_rat;
            auto expr = expr_in.to_rational(map_rat);
            auto nd = numer_denom(expr);
            expr = factor_flint(nd.op(0), false) / factor_flint(nd.op(1), false);
            return expr.subs(map_rat,nopat);
        } else {
            exmap map_rat;
            ex res;
            auto expr = expr_in.to_polynomial(map_rat);
            auto xs = syms(expr);
            if(xs.nops()<1) return expr;
            if(xs.nops()==1) {
                auto x = xs.op(0);
                fmpz_poly_t f;
                fmpz_poly_init(f);
                _to_(f,expr);
                res = _factor_(x,f);
                fmpz_poly_clear(f);
            } else {
                fmpz_mpoly_t f;
                fmpz_mpoly_ctx_t ctx;
                _to_(xs,f,ctx,expr); // call fmpz_mpoly_ctx_init() & fmpz_mpoly_init()
                res = _factor_(xs,f,ctx);
                fmpz_mpoly_clear(f,ctx);
                fmpz_mpoly_ctx_clear(ctx);
            }
            res = res.subs(map_rat,nopat);
            return res;
        }
    }
    
    ex _factor_(const ex & x, fmpz_poly_t fp) {
        fmpz_poly_factor_t fs;
        fmpz_poly_factor_init(fs);
        fmpz_poly_factor_zassenhaus(fs,fp);
        
        fmpz_poly_t f;
        fmpz_poly_init(f);
        fmpz_poly_set_fmpz(f,&fs->c);
        ex res = _to_(x,f);
        
        for (int i=0; i<fs->num; i++) {
            fmpz_poly_set(f,fs->p+i);
            ex fx = _to_(x,f);
            if(fs->exp[i]==1) res *= fx;
            else res *= pow(fx,fs->exp[i]);
        }
        
        fmpz_poly_clear(f);
        fmpz_poly_factor_clear(fs);
        return res;
    }
    
    ex _factor_(const lst & xs, fmpz_mpoly_t fp, fmpz_mpoly_ctx_t ctx) {
        fmpz_mpoly_factor_t fs;
        fmpz_mpoly_factor_init(fs, ctx);
        if(!fmpz_mpoly_factor(fs, fp, ctx)) {
            flint_abort();
            throw Error("fmpz_mpoly_factor failed.");
        }
        
        fmpz_t c;
        fmpz_init(c);
        fmpz_mpoly_factor_get_constant_fmpz(c,fs,ctx);
        ex res = _to_(c);
        fmpz_clear(c);
        
        fmpz_mpoly_t f;
        fmpz_mpoly_init(f,ctx);
        for (int i=0; i<fs->num; i++) {
            fmpz_mpoly_set(f,fs->poly+i,ctx);
            ex fx = _to_(xs,f,ctx);
            if(fs->exp[i]==1) res *= fx;
            else res *= pow(fx,fs->exp[i]);
        }
        
        fmpz_mpoly_clear(f,ctx);
        fmpz_mpoly_factor_clear(fs,ctx);
        return res;
    }
    
    ex _factor_(const lst & xs, fmpq_mpoly_t fp, fmpq_mpoly_ctx_t ctx) {
        fmpq_mpoly_factor_t fs;
        fmpq_mpoly_factor_init(fs, ctx);
        if(!fmpq_mpoly_factor(fs, fp, ctx)) {
            flint_abort();
            throw Error("fmpz_mpoly_factor failed.");
        }
        
        fmpq_t c;
        fmpq_init(c);
        fmpq_mpoly_factor_get_constant_fmpq(c,fs,ctx);
        ex res = _to_(c);
        fmpq_clear(c);
        
        fmpq_mpoly_t f;
        fmpq_mpoly_init(f,ctx);
        for (int i=0; i<fs->num; i++) {
            fmpq_mpoly_set(f,fs->poly+i,ctx);
            ex fx = _to_(xs,f,ctx);
            if(fs->exp[i]==1) res *= fx;
            else res *= pow(fx,fs->exp[i]);
        }
        
        fmpq_mpoly_clear(f,ctx);
        fmpq_mpoly_factor_clear(fs,ctx);
        return res;
    }
    
    inline void _to_(const lst & xs, fmpz_mpoly_q_t f, fmpz_mpoly_ctx_t ctx, const ex & e) {
        if(e.info(info_flags::rational)) {
            fmpq_t fq;
            fmpq_init(fq);
            _to_(fq,e);
            fmpz_mpoly_q_set_fmpq(f,fq,ctx);
            fmpq_clear(fq);
            return;
        } else if(is_a<add>(e)) {
            fmpz_mpoly_q_zero(f,ctx);
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            for(auto item : e) {
                _to_(xs,fi,ctx,item);
                fmpz_mpoly_q_add(f, f, fi, ctx);
            }
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(is_a<mul>(e)) {
            fmpz_mpoly_q_one(f,ctx);
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            for(auto item : e) {
                _to_(xs,fi,ctx,item);
                fmpz_mpoly_q_mul(f, f, fi, ctx);
            }
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::posint)) {
            ulong n = ex_to<numeric>(e.op(1)).to_int();
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            _to_(xs, fi, ctx, e.op(0));
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(f), fmpz_mpoly_q_numref(fi), n, ctx);
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(f), fmpz_mpoly_q_denref(fi), n, ctx);
            fmpz_mpoly_q_canonicalise(f,ctx);
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::negint)) {
            ulong n = -ex_to<numeric>(e.op(1)).to_int();
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            _to_(xs, fi, ctx, e.op(0));
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(f), fmpz_mpoly_q_numref(fi), n, ctx);
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(f), fmpz_mpoly_q_denref(fi), n, ctx);
            fmpz_mpoly_q_canonicalise(f,ctx);
            fmpz_mpoly_q_inv(f,f,ctx);
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(e.is_polynomial(xs)) {
            string vars[xs.nops()];
            const char* cvars[xs.nops()];
            for(int i=0; i<xs.nops(); i++) {
                vars[i] = ex2str(xs.op(i));
                cvars[i] = vars[i].c_str();
            }
            fmpz_mpoly_q_one(f,ctx);
            string es = ex2str(e);
            if(fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_numref(f), es.c_str(), cvars, ctx)) {
                cout << e << endl;
                cout << xs << endl;
                throw Error("fmpz_mpoly_set_str_pretty error.");
            }
            return;
        } else {
            cout << "expr = " << e << endl;
            throw Error("_to_(lst, fmpz_mpoly_q_t, fmpz_mpoly_ctx_t, ex): Not supported region");
        }
    }
    
    ex normal_flint(const ex & expr_in, int opt) {
        exmap map_rat;
        ex res;
        auto expr = expr_in.to_rational(map_rat);
        auto xs = syms(expr);
        if(xs.nops()<1) return expr;
        if(xs.nops()==1) {
            auto x = xs.op(0);
            auto sx = ex2str(x);
            symtab st;
            st[sx] = x;
            
            if(opt==o_flint) {
                fmpz_poly_q_t f;
                fmpz_poly_q_init(f);
                _to_(f, expr);
                auto cstr = fmpz_poly_get_str_pretty(fmpz_poly_q_numref(f), sx.c_str());
                string nstr(cstr);
                flint_free(cstr);
                cstr = fmpz_poly_get_str_pretty(fmpz_poly_q_denref(f), sx.c_str());
                string dstr(cstr);
                flint_free(cstr);
                fmpz_poly_q_clear(f);
                res = str2ex(nstr,st)/str2ex(dstr,st);
            } else if(opt==o_flintf) {
                fmpz_poly_q_t f;
                fmpz_poly_q_init(f);
                _to_(f, expr);
                auto num = _factor_(x, fmpz_poly_q_numref(f));
                auto den = _factor_(x, fmpz_poly_q_denref(f));
                fmpz_poly_q_clear(f);
                res = num/den;
            } else if(opt==o_flintfD) {
                fmpz_poly_q_t f;
                fmpz_poly_q_init(f);
                _to_(f, expr);
                auto cstr = fmpz_poly_get_str_pretty(fmpz_poly_q_numref(f), sx.c_str());
                string nstr(cstr);
                flint_free(cstr);
                auto num = str2ex(nstr,st);
                auto den = _factor_(x, fmpz_poly_q_denref(f));
                fmpz_poly_q_clear(f);
                res = num/den;
            } else throw Error("normal_flint: unsupported option.");
        } else {
            ex e = expr;
            symtab st;
            string vars[xs.nops()];
            const char* cvars[xs.nops()];
            for(int i=0; i<xs.nops(); i++) {
                vars[i] = ex2str(xs.op(i));
                cvars[i] = vars[i].c_str();
                st[cvars[i]] = xs.op(i);
            }
            if(opt==o_flint) {
                fmpz_mpoly_q_t f;
                fmpz_mpoly_ctx_t ctx;
                fmpz_mpoly_ctx_init(ctx, xs.nops(), ORD_LEX);
                fmpz_mpoly_q_init(f, ctx);
                _to_(xs,f,ctx,e);
                auto cstr = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), cvars, ctx);
                string nstr(cstr);
                flint_free(cstr);
                cstr = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(f), cvars, ctx);
                string dstr(cstr);
                flint_free(cstr);
                fmpz_mpoly_q_clear(f, ctx);
                fmpz_mpoly_ctx_clear(ctx);
                res = str2ex(nstr,st)/str2ex(dstr,st);
            } else if(opt==o_flintf) {
                fmpz_mpoly_q_t f;
                fmpz_mpoly_ctx_t ctx;
                fmpz_mpoly_ctx_init(ctx, xs.nops(), ORD_LEX);
                fmpz_mpoly_q_init(f, ctx);
                _to_(xs,f,ctx,e);
                auto num = _factor_(xs, fmpz_mpoly_q_numref(f), ctx);
                auto den = _factor_(xs, fmpz_mpoly_q_denref(f), ctx);
                fmpz_mpoly_q_clear(f, ctx);
                fmpz_mpoly_ctx_clear(ctx);
                res = num/den;
            } else if(opt==o_flintfD) {
                fmpz_mpoly_q_t f;
                fmpz_mpoly_ctx_t ctx;
                fmpz_mpoly_ctx_init(ctx, xs.nops(), ORD_LEX);
                fmpz_mpoly_q_init(f, ctx);
                _to_(xs,f,ctx,e);
                auto cstr = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), cvars, ctx);
                string nstr(cstr);
                flint_free(cstr);
                auto num = str2ex(nstr,st);
                auto den = _factor_(xs, fmpz_mpoly_q_denref(f), ctx);
                fmpz_mpoly_q_clear(f, ctx);
                fmpz_mpoly_ctx_clear(ctx);
                res = num/den;
            } else throw Error("normal_flint: unsupported option.");
        }
        res = res.subs(map_rat,nopat);
        return res;
    }
    
    matrix normal_flint(const matrix & mat) {
        matrix m = mat;
        for(int i=0; i<m.nops(); i++) m.let_op(i) = normal_flint(m.op(i));
        return m;
    }   
    
    ex den_lcm(const ex & expr) {
        auto xs = syms(expr);
        if(xs.nops()!=1) throw Error("0 or >1 symbols found.");
        auto x = xs.op(0);
        auto str = ex2str(x);
        auto sx = str.c_str();
        symtab st;
        st[sx] = x;
        
        fmpz_poly_t lcm;
        fmpz_poly_init(lcm);
        if(fmpz_poly_set_str(lcm, "1  1")) throw Error("denlcm, fmpz_poly_set_str failed!");
        
        fmpz_poly_q_t f;
        fmpz_poly_q_init(f);
        
        for(int i=0; i<expr.nops(); i++) {
            _to_(f, expr.op(i));
            fmpz_poly_lcm(lcm, lcm, fmpz_poly_q_denref(f));
        }
        
        ex res = _to_(x,lcm);
        fmpz_poly_q_clear(f);
        fmpz_poly_clear(lcm);
        return res;
    }
    
    lst poly_roots(const ex & pex, slong fp) {
        lst root_lst;
        fmpz_poly_t poly;
        fmpz_poly_init(poly);
        _to_(poly,pex);
        fmpz_poly_factor_t fac;
        fmpz_poly_factor_init(fac);
        fmpz_poly_factor_squarefree(fac, poly);
        for(int i=0; i<fac->num; i++) {
            auto deg = fmpz_poly_degree(fac->p + i);
            auto roots = _acb_vec_init(deg);
            arb_fmpz_poly_complex_roots(roots, fac->p + i, 0, fp);
            for(int j = 0; j < deg; j++) root_lst.append(_to_(roots+j,fp));
            _acb_vec_clear(roots, deg);
        }
        fmpz_poly_factor_clear(fac);
        fmpz_poly_clear(poly);
        return root_lst;
    }
    
}
