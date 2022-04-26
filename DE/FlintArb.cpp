
#include "FlintArb.h"
#include "cln/cln.h"

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
        
        inline int ldegree(fmpz_poly_t f) {
            fmpz_t z;
            fmpz_init(z);
            auto n = fmpz_poly_length(f);
            int res = n;
            for(int i=0; i<n; i++) {
                fmpz_poly_get_coeff_fmpz(z,f,i);
                if(!fmpz_is_zero(z)) { res = i; break; }
            }
            fmpz_clear(z);
            return res;
        }
    }
    
    //=*********************************************************************=
    
    slong fp2dp(slong fp) {
        static cln::cl_R r("0.30102999566398119521373889472449302676818988146211_50");
        return cln::cl_I_to_long(cln::ceiling1(r*fp));
    }
        
    slong dp2fp(slong dp) {
        static cln::cl_R r("3.3219280948873623478703194294893901758648313930246_50");
        return cln::cl_I_to_long(cln::ceiling1(r*dp));
    }
    
    void _to_(fmpz_poly_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.rows();
        fmpz_poly_t p;
        fmpz_poly_init(p);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(p, mat(r,c));
            fmpz_poly_set(fmpz_poly_mat_entry(m,r,c), p);
        }
        fmpz_poly_clear(p);
    }
    
    void _to_(fmpz_poly_q_t f, const ex & e) {
        if(is_a<symbol>(e)) {
            fmpz_poly_q_set_str(f, "2  0 1/1  1");
            return;
        } else if(e.info(info_flags::rational)) {
            auto nd = e.numer_denom();
            ostringstream oss;
            oss << "1  " << nd.op(0) << "/1  " << nd.op(1);
            fmpz_poly_q_set_str(f, oss.str().c_str());
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
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            _to_(fi, e.op(0));
            fmpz_poly_q_pow(f, fi, n);
            fmpz_poly_q_clear(fi);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::negint)) {
            ulong n = -ex_to<numeric>(e.op(1)).to_int();
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            _to_(fi, e.op(0));
            fmpz_poly_q_pow(f, fi, n);
            fmpz_poly_q_inv(f,f);
            fmpz_poly_q_clear(fi);
            return;
        } else {
            cout << "expr = " << e << endl;
            throw Error("ex_to_fmpz_poly_q_t Not supported region");
        }
    }
    
    void _to_(fmpz_poly_t f, const ex & e) {
        if(is_a<symbol>(e)) {
            fmpz_poly_set_str(f, "2  0 1");
            return;
        } else if(e.info(info_flags::integer)) {
            ostringstream oss;
            oss << "1  " << e;
            fmpz_poly_set_str(f, oss.str().c_str());
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
            fmpz_poly_t fi;
            fmpz_poly_init(fi);
            _to_(fi, e.op(0));
            fmpz_poly_pow(f, fi, n);
            fmpz_poly_clear(fi);
            return;
        } else {
            cout << "expr = " << e << endl;
            throw Error("ex_to_fmpz_poly_t Not supported region");
        }
    }
    
    void _to_(fmpz_t q, const ex & expr) {
        auto str = ex2str(expr);
        fmpz_set_str(q,str.c_str(),10);
    }
    
    void _to_(fmpq_t q, const ex & expr) {
        auto str = ex2str(expr);
        fmpq_set_str(q,str.c_str(),10);
    }
    
    void _to_(fmpz_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.rows();
        fmpz_t z;
        fmpz_init(z);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(z, mat(r,c));
            fmpz_set(fmpz_mat_entry(m,r,c), z);
        }
        fmpz_clear(z);
    }
    
    void _to_(fmpq_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.rows();
        fmpq_t q;
        fmpq_init(q);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(q, mat(r,c));
            fmpq_set(fmpq_mat_entry(m,r,c), q);
        }
        fmpq_clear(q);
    }
    
    void _to_(arb_t r, const ex & expr, slong fp) {
        set_precision(fp2dp(fp)+10);
        auto ne = expr.evalf();
        reset_precision();
        if(!is_a<numeric>(ne)) {
            cout << endl << "ne = " << ne << endl;
            throw Error("to_arb_t: NOT a number");
        }
        auto str = ex2str(ne);
        arb_set_str(r, str.c_str(), fp);
    }
    
    void _to_(arb_mat_t m, const matrix & mat, slong fp) {
        auto nr = mat.rows();
        auto nc = mat.rows();
        arb_t z;
        arb_init(z);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(z, mat(r,c), fp);
            arb_set(arb_mat_entry(m,r,c), z);
        }
        arb_clear(z);
    }
    
    void _to_(acb_t z, const ex & expr, slong fp) {
        set_precision(fp2dp(fp)+10);
        auto ne = expr.evalf();
        reset_precision();
        if(!is_a<numeric>(ne)) {
            cout << endl << "ne = " << ne << endl;
            throw Error("to_acb_t: NOT a number");
        }
        
        arb_t rb, ib;
        arb_init(rb);
        arb_init(ib);
        auto nre = ne.real_part();
        auto nie = ne.imag_part();
        auto rstr = ex2str(nre);
        arb_set_str(rb, rstr.c_str(), fp);
        auto istr = ex2str(nie);
        arb_set_str(ib, istr.c_str(), fp);
        acb_set_arb_arb(z, rb, ib);
        arb_clear(rb);
        arb_clear(ib);
    }
    
    void _to_(acb_mat_t m, const matrix & mat, slong fp) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(acb_mat_entry(m,r,c), mat(r,c), fp);
        }
    }
    
    void _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx, const ex & e) {
        exmap x2x;
        for(int i=0; i<xs.nops(); i++) x2x[xs.op(i)] = Symbol("x"+to_string(i+1));
        fmpz_mpoly_ctx_init(ctx, xs.nops(), ORD_LEX);
        fmpz_mpoly_init(f, ctx);
        if(fmpz_mpoly_set_str_pretty(f, ex2str(e.subs(x2x,nopat)).c_str(), NULL, ctx))
            throw Error("ex_to_fmpz_mpoly_t failed.");
    }
    
    void _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx, const ex & e) {
        exmap x2x;
        for(int i=0; i<xs.nops(); i++) x2x[xs.op(i)] = Symbol("x"+to_string(i+1));
        fmpq_mpoly_ctx_init(ctx, xs.nops(), ORD_LEX);
        fmpq_mpoly_init(f, ctx);
        if(fmpq_mpoly_set_str_pretty(f, ex2str(e.subs(x2x,nopat)).c_str(), NULL, ctx))
            throw Error("ex_to_fmpq_mpoly_t failed.");
    }
    
    //=*********************************************************************=
    
    matrix _to_(const ex & x, fmpz_poly_mat_t m) {
        auto nr = fmpz_poly_mat_nrows(m);
        auto nc = fmpz_poly_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(x,fmpz_poly_mat_entry(m,r,c));
        }
        return mat;
    }
    
    ex _to_(const ex & x, fmpz_poly_q_t f) {
        if(!fmpz_poly_q_is_canonical(f)) fmpz_poly_q_canonicalise(f);
        
        ex res_den = 0;
        auto den = fmpz_poly_q_denref(f);
        auto xdn = fmpz_poly_length(den);
        for(slong i=0; i<xdn; i++) {
            fmpz_t z;
            fmpz_init(z);
            fmpz_poly_get_coeff_fmpz(z,den,i);
            auto str = fmpz_get_str(NULL,10,z);
            res_den += numeric(str) * pow(x,i);
            fmpz_clear(z);
            flint_free(str);
        }
        
        ex res_num = 0;
        auto num = fmpz_poly_q_numref(f);
        auto xnn = fmpz_poly_length(num);
        for(slong i=0; i<xnn; i++) {
            fmpz_t z;
            fmpz_init(z);
            fmpz_poly_get_coeff_fmpz(z,num,i);
            auto str = fmpz_get_str(NULL,10,z);
            if(xdn==1) res_num += numeric(str)/res_den * pow(x,i);
            else res_num += numeric(str) * pow(x,i);
            fmpz_clear(z);
            flint_free(str);
        }
        if(xdn==1) return res_num;
        return res_num/res_den;
    }

    ex _to_(const ex & x, fmpz_poly_t f) {
        ex res = 0;
        auto xn = fmpz_poly_length(f);
        for(slong i=0; i<xn; i++) {
            fmpz_t z;
            fmpz_init(z);
            fmpz_poly_get_coeff_fmpz(z,f,i);
            auto str = fmpz_get_str(NULL,10,z);
            res += numeric(str) * pow(x,i);
            fmpz_clear(z);
            flint_free(str);
        }
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
        auto cstr = fmpq_get_str(NULL,10,q);
        ex n = numeric(cstr);
        flint_free(cstr);
        return n;
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
    
    matrix _to_(acb_mat_t m, slong dp) {
        auto nr = acb_mat_nrows(m);
        auto nc = acb_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(acb_mat_entry(m,r,c), dp);
        }
        return mat;
    }
    
    ex _to_(const lst & xs, fmpz_mpoly_t f, fmpz_mpoly_ctx_t ctx) {
        symtab x2x;
        string vars[xs.nops()];
        const char* cvars[xs.nops()];
        for(int i=0; i<xs.nops(); i++) {
            vars[i] = "x"+to_string(i+1);
            cvars[i] = vars[i].c_str();
            x2x[cvars[i]] = xs.op(i);
        }
        auto cstr = fmpz_mpoly_get_str_pretty(f, cvars, ctx);
        string str(cstr);
        flint_free(cstr);
        return str2ex(str,x2x); 
    }
    
    ex _to_(const lst & xs, fmpq_mpoly_t f, fmpq_mpoly_ctx_t ctx) {
        symtab x2x;
        string vars[xs.nops()];
        const char* cvars[xs.nops()];
        for(int i=0; i<xs.nops(); i++) {
            vars[i] = "x"+to_string(i+1);
            cvars[i] = vars[i].c_str();
            x2x[cvars[i]] = xs.op(i);
        }
        auto cstr = fmpq_mpoly_get_str_pretty(f, cvars, ctx);
        string str(cstr);
        flint_free(cstr);
        return str2ex(str,x2x); 
    }
        
    //=*********************************************************************=
    
    ex factor_flint(const ex & expr) {
        auto xs = syms(expr);
        if(xs.nops()<1) return expr;
        ex res;
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
            _to_(xs,f,ctx,expr);
            res = _factor_(xs,f,ctx);
            fmpz_mpoly_clear(f,ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }
        return res;
    }
    
    ex factor_fpq(const ex & expr) {
        auto xs = syms(expr);
        if(xs.nops()<1) return expr;
        fmpq_mpoly_t f;
        fmpq_mpoly_ctx_t ctx;
        _to_(xs,f,ctx,expr);
        ex res = _factor_(xs,f,ctx);
        fmpq_mpoly_clear(f,ctx);
        fmpq_mpoly_ctx_clear(ctx);
        return res;
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
    
    ex normal_flint(const ex & expr) {
        auto xs = syms(expr);
        if(xs.nops()<1) return expr;
        if(xs.nops()>1) throw Error(">=2 variables found.");
        auto x = xs.op(0);
        auto sx = ex2str(x).c_str();
        symtab st;
        st[sx] = x;
        fmpz_poly_q_t f;
        fmpz_poly_q_init(f);
        _to_(f, expr);
        auto res = _to_(x,f);
        fmpz_poly_q_clear(f);
        return res;
    }
    
    ex den_lcm(const ex & expr) {
        auto xs = syms(expr);
        if(xs.nops()!=1) throw Error("0 or >1 symbols found.");
        auto x = xs.op(0);
        auto sx = ex2str(x).c_str();
        symtab st;
        st[sx] = x;
        
        fmpz_poly_t lcm;
        fmpz_poly_init(lcm);
        fmpz_poly_set_str(lcm, "1  1");
        
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
    
    
}
