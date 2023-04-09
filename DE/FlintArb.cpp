
#include "FlintArb.h"
#include "cln/cln.h"
#include "fmpz_mpoly_q.h"

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
        
        inline int degree(fmpz_poly_t f) { return fmpz_poly_degree(f); }
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
    
    slong cln_ceiling(const ex & e) {
        auto n = ex_to<numeric>(e).to_cl_N();
        return cln::cl_I_to_long(cln::ceiling1(cln::realpart(n))); 
    }
    
    void _to_(fmpz_poly_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        fmpz_poly_t p;
        fmpz_poly_init(p);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(p, mat(r,c));
            fmpz_poly_set(fmpz_poly_mat_entry(m,r,c), p);
        }
        fmpz_poly_clear(p);
    }
    
    void _to_(fmpz_poly_q_t f, const ex & e) {
        if(syms(e).nops()>1) throw Error(">=2 variables found.");
        if(is_a<symbol>(e)) {
            fmpz_poly_q_set_str(f, "2  0 1/1  1");
            return;
        } else if(e.info(info_flags::rational)) {
            auto nd = e.numer_denom();
            ostringstream oss;
            oss << "1  " << nd.op(0) << "/1  " << nd.op(1);
            if(fmpz_poly_q_set_str(f, oss.str().c_str())) {
                cout << "str: " << oss.str() << endl;
                throw Error("_to_(fmpz_poly_q_t f, const ex & e)");
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
        if(syms(e).nops()>1) throw Error(">=2 variables found.");
        if(is_a<symbol>(e)) {
            fmpz_poly_set_str(f, "2  0 1");
            return;
        } else if(e.info(info_flags::integer)) {
            ostringstream oss;
            oss << "1  " << e;
            if(fmpz_poly_set_str(f, oss.str().c_str())) {
                cout << "str: " << oss.str() << endl;
                throw Error("void _to_(fmpz_poly_t f, const ex & e)");
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
        if(fmpz_set_str(q,str.c_str(),10)) {
            cout << "str: " << str << endl;
            throw Error("fmpz_set_str error.");
        }
    }
    
    void _to_(fmpq_t q, const ex & expr) {
        auto str = ex2str(expr);
        if(fmpq_set_str(q,str.c_str(),10)) {
            cout << "str: " << str << endl;
            throw Error("fmpq_set_str error.");
        }
    }
    
    void _to_(fmpz_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.cols();
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
        auto nc = mat.cols();
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
        if(arb_set_str(r, str.c_str(), fp)) {
            cout << "str: " << str << endl;
            throw Error("arb_set_str error.");
        }
    }
    
    void _to_(arb_mat_t m, const matrix & mat, slong fp) {
        auto nr = mat.rows();
        auto nc = mat.cols();
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
    
    void _to_(acb_mat_t m, const matrix & mat, slong fp) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            _to_(acb_mat_entry(m,r,c), mat(r,c), fp);
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
        if(fmpz_mpoly_set_str_pretty(f, ex2str(e).c_str(), cvars, ctx)) {
            cout << "str: " << e << endl;
            throw Error("fmpz_mpoly_set_str_pretty error.");
        }
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
        if(fmpq_mpoly_set_str_pretty(f, ex2str(e).c_str(), cvars, ctx)) {
            cout << "str: " << e << endl;
            throw Error("fmpq_mpoly_set_str_pretty error.");
        }
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
    
    matrix _to_(acb_mat_t m, slong fp) {
        auto nr = acb_mat_nrows(m);
        auto nc = acb_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = _to_(acb_mat_entry(m,r,c), fp);
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
    
    //=*********************************************************************=
        
    MX::MX() { }
    MX::MX(const matrix & m) { init(m); }
    MX::MX(const MX & mx) { init(mx); }
    MX::MX(const vector<vector<fmpz_poly_q_t>> & m) { init(m); }
    
    MX::~MX() { clear(); nr=-1; nc=-1; }
        
    void MX::clear() {
        if(nr<1||nc<1) return;
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) fmpz_poly_q_clear(M[r][c]);
        }
        nr=-1; nc=-1;
    }
    
    void MX::init(const matrix & m) {
        clear();
        nr = m.rows();
        nc = m.cols();
        if(nr<1||nc<1) return;
        M.resize(nr);
        for(int r=0; r<nr; r++) {
            M[r] = vector<fmpz_poly_q_t>(nc);
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_init(M[r][c]);
                _to_(M[r][c],m(r,c));
            }
        }
    }
    
    void MX::init(const MX & mx) {
        clear();
        if(mx.nr<1||mx.nc<1) return;
        nr = mx.nr;
        nc = mx.nc;
        M.resize(nr);
        for(int r=0; r<nr; r++) {
            M[r] = vector<fmpz_poly_q_t>(nc);
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_init(M[r][c]);
                fmpz_poly_q_set(M[r][c],mx.M[r][c]);
            }
        }
    }
    
    void MX::init(const vector<vector<fmpz_poly_q_t>> & m) {
        clear();
        nr = m.size();
        if(nr<1) return;
        nc = m[0].size();
        if(nc<1) return;
        M.resize(nr);
        for(int r=0; r<nr; r++) {
            M[r] = vector<fmpz_poly_q_t>(nc);
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_init(M[r][c]);
                fmpz_poly_q_set(M[r][c],m[r][c]);
            }
        }
    }
    
    MX & MX::add(const MX & mx) {
        #pragma omp parallel for schedule(runtime) collapse(2)
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_add(M[r][c],M[r][c],mx.M[r][c]);
                flint_cleanup();
            }
        }
        return *this;
    }
    
    MX & MX::sub(const MX & mx) {
        #pragma omp parallel for schedule(runtime) collapse(2)
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_sub(M[r][c],M[r][c],mx.M[r][c]);
                flint_cleanup();
            }
        }
        return *this;
    }
    
    MX & MX::mul(const MX & mx) { 
        if(this!=&mx && mx.nc==nc) {
            vector<fmpz_poly_q_t> row(nc); // keep a row
            for(int c=0; c<nc; c++) fmpz_poly_q_init(row[c]);
            for(int r=0; r<nr; r++) {
                for(int c=0; c<nc; c++) fmpz_poly_q_set(row[c], M[r][c]);
                #pragma omp parallel for schedule(runtime)
                for(int c=0; c<nc; c++) {
                    fmpz_poly_q_zero(M[r][c]);
                    for(int i=0; i<nc; i++) fmpz_poly_q_addmul(M[r][c],row[i],mx.M[i][c]);
                    flint_cleanup();
                }
            }
            for(int c=0; c<nc; c++) fmpz_poly_q_clear(row[c]);
        } else {
            vector<vector<fmpz_poly_q_t>> mat(nr);
            int nc2 = mx.nc;
            #pragma omp parallel for schedule(runtime)
            for(int r=0; r<nr; r++) {
                mat[r] = vector<fmpz_poly_q_t>(nc2);
                for(int c=0; c<nc2; c++) {
                    fmpz_poly_q_init(mat[r][c]);
                    fmpz_poly_q_zero(mat[r][c]);
                    for(int i=0; i<nc; i++) fmpz_poly_q_addmul(mat[r][c],M[r][i],mx.M[i][c]);
                }
                flint_cleanup();
            }
            clear();
            init(mat);
            for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) fmpz_poly_q_clear(mat[r][c]);
        }
        return *this;
    }
    
    MX & MX::mul_left(const MX & mx) {
        if(this!=&mx && mx.nr==nr) {
            vector<fmpz_poly_q_t> col(nr);
            for(int r=0; r<nr; r++) fmpz_poly_q_init(col[r]);
            for(int c=0; c<nc; c++) {
                for(int r=0; r<nr; r++) fmpz_poly_q_set(col[r], M[r][c]);
                #pragma omp parallel for schedule(runtime)
                for(int r=0; r<nr; r++) {
                    fmpz_poly_q_zero(M[r][c]);
                    for(int i=0; i<nr; i++) fmpz_poly_q_addmul(M[r][c],mx.M[r][i],col[i]);
                    flint_cleanup();
                }
            }
            for(int r=0; r<nr; r++) fmpz_poly_q_clear(col[r]);
        } else {
            int nr2 = mx.nr;
            vector<vector<fmpz_poly_q_t>> mat(nr2);
            #pragma omp parallel for schedule(runtime)
            for(int r=0; r<nr2; r++) {
                mat[r] = vector<fmpz_poly_q_t>(nc);
                for(int c=0; c<nc; c++) {
                    fmpz_poly_q_init(mat[r][c]);
                    fmpz_poly_q_zero(mat[r][c]);
                    for(int i=0; i<nr; i++) fmpz_poly_q_addmul(mat[r][c],mx.M[r][i],M[i][c]);
                }
                flint_cleanup();
            }
            clear();
            init(mat);
            for(int r=0; r<nr2; r++) for(int c=0; c<nc; c++) fmpz_poly_q_clear(mat[r][c]);
        }
        return *this;
    }
    
    MX & MX::add(const matrix & m) { MX mx(m); return add(mx); }
    MX & MX::sub(const matrix & m) { MX mx(m); return sub(mx); }
    MX & MX::mul(const matrix & m) { MX mx(m); return mul(mx); }
    MX & MX::mul_left(const matrix & m) { MX mx(m); return mul_left(mx); }
    
    MX & MX::scale(const ex & s) {
        fmpz_poly_q_t f;
        fmpz_poly_q_init(f);
        _to_(f,s);
        #pragma omp parallel for schedule(runtime) collapse(2)
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_mul(M[r][c],M[r][c],f);
                flint_cleanup();
            }
        }
        fmpz_poly_q_clear(f);
        return *this;
    }
    
    MX & MX::scale(fmpz_poly_q_t f) {
        #pragma omp parallel for schedule(runtime) collapse(2)
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_mul(M[r][c],M[r][c],f);
                flint_cleanup();
            }
        }
        return *this;
    }
    
    MX & MX::scale(fmpz_poly_t f) {
        #pragma omp parallel for schedule(runtime) collapse(2)
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto nr = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_mul(nr,nr,f);
                fmpz_poly_q_canonicalise(M[r][c]);
                flint_cleanup();
            }
        }
        return *this;
    }
    
    MX & MX::dx() {
        #pragma omp parallel for schedule(runtime) collapse(2)
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_derivative(M[r][c],M[r][c]);
                flint_cleanup();
            }
        }
        return *this;
    }
        
    MX & MX::balance(const matrix & P) {
        if(nr!=nc) throw Error("nr!=nc");
        static symbol x("x");
        MX mx1(P);
        MX mx2(mx1);
        mx1.scale(x);
        mx2.scale(1/x);
        auto coP = ex_to<matrix>(unit_matrix(nr)).sub(P);
        MX mx3(coP);
        auto mx4(mx3);
        mx3.sub(mx1);
        mx4.sub(mx2);
        return mul_left(mx3).mul(mx4).add(mx2);
    }

    MX & MX::transform(const matrix & t, const matrix & ti) {
        if(nr!=nc) throw Error("nr!=nc");
        MX tX(t), tiX(ti);
        return mul(tX).sub(tX.dx()).mul_left(tiX);
    }
    
    MX & MX::shift(const ex & x0) {
        if(!x0.info(info_flags::rational)) {
            cout << endl << "x0 = " << x0 << endl;
            throw Error("x0 is not rational.");
        }
        auto nd = normal(-x0).numer_denom();
        fmpz_t qn, qd;
        fmpz_init(qn); 
        fmpz_init(qd);
        _to_(qn,nd.op(0));
        _to_(qd,nd.op(1));
        #pragma omp parallel for schedule(runtime)
        for(int r=0; r<nr; r++) {
            fmpz_t z;
            fmpz_init(z);
            for(int c=0; c<nc; c++) {
                auto num = fmpz_poly_q_numref(M[r][c]);
                auto nn = fmpz_poly_length(num);
                auto den = fmpz_poly_q_denref(M[r][c]);
                auto dn = fmpz_poly_length(den);
                slong nr = nn;
                if(nn>dn) {
                    nr = nn;
                    fmpz_pow_ui(z,qd,nn-dn);
                    fmpz_poly_scalar_mul_fmpz(den,den,z);
                } else if(dn>nn) {
                    nr = dn;
                    fmpz_pow_ui(z,qd,dn-nn);
                    fmpz_poly_scalar_mul_fmpz(num,num,z);
                }
                for(int i=0; i<nn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(num, i);
                    if(ci==NULL) throw Error("ci == NULL");
                    fmpz_pow_ui(z,qd,nn-i);
                    fmpz_mul(ci, ci, z);
                }
                for(int i=0; i<dn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(den, i);
                    fmpz_pow_ui(z,qd,dn-i);
                    fmpz_mul(ci, ci, z);
                }
                fmpz* rs;
                rs = _fmpz_vec_init(nr);
                for(int i=0; i<nr; i++) fmpz_set(rs+i,qn);
                _fmpz_poly_newton_to_monomial(num->coeffs, rs, nn);
                _fmpz_poly_newton_to_monomial(den->coeffs, rs, dn);
                for(int i=0; i<nn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(num, i);
                    fmpz_pow_ui(z,qd,i);
                    fmpz_mul(ci, ci, z);
                }
                for(int i=0; i<dn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(den, i);
                    fmpz_pow_ui(z,qd,i);
                    fmpz_mul(ci, ci, z);
                }
                _fmpz_vec_clear(rs,nr);
                fmpz_poly_q_canonicalise(M[r][c]);
            }
            fmpz_clear(z);
            flint_cleanup();
        }
        fmpz_clear(qn); 
        fmpz_clear(qd);
        return *this;
    }
    
    matrix MX::operator()(const ex & x) {
        if(nr<1||nc<1) return matrix();
        matrix rmat(nr,nc);
        symtab st;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            rmat(r,c) = _to_(x,M[r][c]);
        }
        return rmat;
    }
    
    void MX::operator()(vector<vector<fmpz_poly_t>> & m) {
        int nr_m = m.size();
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            fmpz_poly_set(m[r][c],fmpz_poly_q_numref(M[r][c]));
        }
    }
    
    void MX::operator()(vector<vector<fmpz_poly_q_t>> & m) {
        int nr_m = m.size();
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            fmpz_poly_q_set(m[r][c],M[r][c]);
        }
    }
    
    void MX::operator()(vector<vector<acb_poly_t>> & m, slong fp) {
        int nr_m = m.size();
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            acb_poly_set_fmpz_poly(m[r][c],fmpz_poly_q_numref(M[r][c]),fp);
        }
    }
    
    int MX::prank() {
        int pr;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            fmpz_poly_q_canonicalise(M[r][c]);
            auto ipr = ldegree(fmpz_poly_q_numref(M[r][c]));
            if(ipr==0) ipr = -ldegree(fmpz_poly_q_denref(M[r][c]));
            if(r==0 && c==0) pr = ipr;
            else if(ipr<pr) pr = ipr;
        }
        return -1-pr;
    }
    
    matrix MX::a0() {
        fmpz_t z;
        fmpz_init(z);
        auto pr = prank();
        matrix m(nr,nc);
        if(pr<0) return m;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto num = fmpz_poly_q_numref(M[r][c]);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto ipr_num = ldegree(num);
            auto ipr_den = ldegree(den);
            if(ipr_den-ipr_num<pr+1) m(r,c)=0;
            else if(ipr_num>0) {
                if(ipr_den!=0) throw Error("something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num);
                m(r,c) = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,0);
                m(r,c) /= _to_(z);
            } else {
                if(ipr_num!=0) throw Error("something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,0);
                m(r,c) = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,ipr_den);
                m(r,c) /= _to_(z);
            }
        }
        fmpz_clear(z);
        return m;
    }
    
    pair<matrix,matrix> MX::a01() {
        fmpz_t z;
        fmpz_init(z);
        auto pr = prank();
        if(pr<0) pr = 0;
        matrix m0(nr,nc), m1(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto num = fmpz_poly_q_numref(M[r][c]);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto ipr_num = ldegree(num);
            auto ipr_den = ldegree(den);
            if(ipr_den-ipr_num<pr) { m0(r,c)=0; m1(r,c)=0; continue; }
            ex n0, n1, d0, d1;
            if(ipr_num>0) {
                if(ipr_den!=0) throw Error("1: something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num);
                n0 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num+1);
                n1 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,0);
                d0 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,1);
                d1 = _to_(z);
            } else {
                if(ipr_num!=0) throw Error("2: something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,0);
                n0 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,num,1);
                n1 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,ipr_den);
                d0 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,ipr_den+1);
                d1 = _to_(z);
            }
            if(ipr_den-ipr_num==pr) { 
                m0(r,c)=0; 
                m1(r,c)=n0/d0; 
                continue; 
            } else {
                if(ipr_den-ipr_num!=pr+1) throw Error("3: something is wrong here.");
                m0(r,c)=n0/d0; 
                m1(r,c)=(d0*n1-d1*n0)/(d0*d0); 
            }
        }
        fmpz_clear(z);
        return make_pair(m0,m1);
    }
    
    int MX::degree() {
        int deg;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            if(fmpz_poly_length(den)!=1) {
                cout << endl;
                fmpz_poly_q_print_pretty(M[r][c],"x"); 
                cout << endl;
                throw Error("degree: denominator is NOT constant");
            }
            auto cdeg = fmpz_poly_degree(num);
            if(r==0 && c==0) deg = cdeg;
            else if(deg<cdeg) deg = cdeg;
        }
        return deg;
    }
    
    void MX::series(int xn) {
        #pragma omp parallel for schedule(runtime)
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            fmpq_poly_t den_q, num_q;
            fmpq_poly_init(den_q);
            fmpq_poly_init(num_q);
            fmpq_poly_set_fmpz_poly(den_q,den);
            fmpq_poly_set_fmpz_poly(num_q,num);
            auto dn = ldegree(den);
            if(dn>0) fmpq_poly_shift_right(den_q,den_q,dn);
            fmpq_poly_div_series(num_q,num_q,den_q,xn+dn+1);
            fmpq_poly_get_numerator(num,num_q);
            fmpz_t z;
            fmpz_init(z);
            fmpq_poly_get_denominator(z,num_q);
            fmpz_poly_set_fmpz(den,z);
            if(dn>0) fmpz_poly_shift_left(den,den,dn);
            fmpz_clear(z);
            fmpq_poly_clear(den_q);
            fmpq_poly_clear(num_q);
            fmpz_poly_q_canonicalise(M[r][c]);
            flint_cleanup();
        }
    }
    
    int MX::denlcm(fmpz_poly_t dl) { // M will be updated
        vector<fmpz_poly_t> lcms(nr);
        #pragma omp parallel for schedule(runtime)
        for(int r=0; r<nr; r++) {
            fmpz_poly_init(lcms[r]);
            fmpz_poly_set_str(lcms[r], "1  1");
            for(int c=0; c<nc; c++) {
                fmpz_poly_lcm(lcms[r], lcms[r], fmpz_poly_q_denref(M[r][c]));
            }
            flint_cleanup();
        }
        fmpz_poly_set_str(dl, "1  1");
        for(int r=0; r<nr; r++) {
            fmpz_poly_lcm(dl,dl,lcms[r]);
            fmpz_poly_clear(lcms[r]);
        }
        #pragma omp parallel for schedule(runtime) collapse(2)
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto num = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_mul(num, num, dl);
                fmpz_poly_q_canonicalise(M[r][c]);
                flint_cleanup();
            }
        }
        return fmpz_poly_degree(dl);
    }
    
    matrix MX::coeff(int i) {
        fmpz_t z;
        fmpz_init(z);
        matrix rmat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            slong ldeg = ldegree(den);
            if(fmpz_poly_degree(den)!=ldeg) throw Error("coeff: not a monomial in denominator.");
            fmpz_poly_get_coeff_fmpz(z,den,ldeg);
            ex _den = _to_(z);
            if(i+ldeg>=0) {
                fmpz_poly_get_coeff_fmpz(z,num,(slong)(i+ldeg));
                ex _num = _to_(z);
                rmat(r,c) = _num/_den;
            } else rmat(r,c) = 0;
        }
        fmpz_clear(z);
        return rmat;
    }
    
    void MX::coeff(fmpq_mat_t m, int i) {
        fmpz_t zn,zd;
        fmpz_init(zn);
        fmpz_init(zd);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto item = fmpq_mat_entry(m,r,c);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            slong ldeg = ldegree(den);
            if(fmpz_poly_degree(den)!=ldeg) throw Error("coeff: not a monomial in denominator.");
            fmpz_poly_get_coeff_fmpz(zd,den,ldeg);
            if(fmpz_is_zero(zd)) throw Error("denominator is zero.");
            if(i+ldeg>=0) {
                fmpz_poly_get_coeff_fmpz(zn,num,(slong)(i+ldeg));
                fmpq_set_fmpz_frac(item,zn,zd);
            } else fmpq_zero(item);
        }
        fmpz_clear(zn);
        fmpz_clear(zd);
    }
    
    void MX::coeff(acb_mat_t m, int i, slong fp) {
        fmpz_t zn,zd;
        fmpz_init(zn);
        fmpz_init(zd);
        fmpq_t q;
        fmpq_init(q);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto item = acb_mat_entry(m,r,c);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            slong ldeg = ldegree(den);
            if(fmpz_poly_degree(den)!=ldeg) throw Error("coeff: not a monomial in denominator.");
            fmpz_poly_get_coeff_fmpz(zd,den,ldeg);
            if(fmpz_is_zero(zd)) throw Error("denominator is zero.");
            if(i+ldeg>=0) {
                fmpz_poly_get_coeff_fmpz(zn,num,(slong)(i+ldeg));
                fmpq_set_fmpz_frac(q,zn,zd);
                acb_set_fmpq(item,q,fp);
            } else acb_zero(item);
        }
        fmpz_clear(zn);
        fmpz_clear(zd);
        fmpq_clear(q);
    }
    
    //=*********************************************************************=
    
    MQ::MQ() { }    
    MQ::~MQ() { clear(); }
        
    void MQ::clear() {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) fmpz_poly_q_clear(M[r][c]);
        }
        nr=-1; nc=-1;
    }
        
    void MQ::init(const vector<vector<fmpz_poly_q_t>> & m) {
        clear();
        nr = m.size();
        if(nr<1) throw Error("MQ::init, nr<1");
        nc = m[0].size();
        if(nc<1) throw Error("MQ::init, nc<1");
        M.resize(nr);
        for(int r=0; r<nr; r++) {
            M[r] = vector<fmpz_poly_q_t>(nc);
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_init(M[r][c]);
                fmpz_poly_q_set(M[r][c],m[r][c]);
            }
        }
    }
    
    void MQ::scale(fmpz_poly_q_t f) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_mul(M[r][c],M[r][c],f);
            }
        }
    }
    
    void MQ::scale(fmpz_poly_t f) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto nr = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_mul(nr,nr,f);
                fmpz_poly_q_canonicalise(M[r][c]);
            }
        }
    }
    
    int MQ::degree() {
        int deg;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            if(fmpz_poly_length(den)!=1) {
                cout << endl;
                fmpz_poly_q_print_pretty(M[r][c],"x"); 
                cout << endl;
                throw Error("degree: denominator is NOT constant");
            }
            auto cdeg = fmpz_poly_degree(num);
            if(r==0 && c==0) deg = cdeg;
            else if(deg<cdeg) deg = cdeg;
        }
        return deg;
    }
    
    int MQ::denlcm(fmpz_poly_t dl) { // M will be updated
        vector<fmpz_poly_t> lcms(nr);
        for(int r=0; r<nr; r++) {
            fmpz_poly_init(lcms[r]);
            fmpz_poly_set_str(lcms[r], "1  1");
            for(int c=0; c<nc; c++) {
                fmpz_poly_lcm(lcms[r], lcms[r], fmpz_poly_q_denref(M[r][c]));
            }
        }
        fmpz_poly_set_str(dl, "1  1");
        for(int r=0; r<nr; r++) {
            fmpz_poly_lcm(dl,dl,lcms[r]);
            fmpz_poly_clear(lcms[r]);
        }
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto num = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_mul(num, num, dl);
                fmpz_poly_q_canonicalise(M[r][c]);
            }
        }
        return fmpz_poly_degree(dl);
    }
    
    void MQ::coeff(fmpq_mat_t m, int i) {
        fmpz_t zn,zd;
        fmpz_init(zn);
        fmpz_init(zd);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto item = fmpq_mat_entry(m,r,c);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            slong ldeg = ldegree(den);
            if(fmpz_poly_degree(den)!=ldeg) throw Error("coeff: not a monomial in denominator.");
            fmpz_poly_get_coeff_fmpz(zd,den,ldeg);
            if(fmpz_is_zero(zd)) throw Error("denominator is zero.");
            if(i+ldeg>=0) {
                fmpz_poly_get_coeff_fmpz(zn,num,(slong)(i+ldeg));
                fmpq_set_fmpz_frac(item,zn,zd);
            } else fmpq_zero(item);
        }
        fmpz_clear(zn);
        fmpz_clear(zd);
    }
    
    void MQ::coeff(acb_mat_t m, int i, slong fp) {
        fmpz_t zn,zd;
        fmpz_init(zn);
        fmpz_init(zd);
        fmpq_t q;
        fmpq_init(q);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto item = acb_mat_entry(m,r,c);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            slong ldeg = ldegree(den);
            if(fmpz_poly_degree(den)!=ldeg) throw Error("coeff: not a monomial in denominator.");
            fmpz_poly_get_coeff_fmpz(zd,den,ldeg);
            if(fmpz_is_zero(zd)) throw Error("denominator is zero.");
            if(i+ldeg>=0) {
                fmpz_poly_get_coeff_fmpz(zn,num,(slong)(i+ldeg));
                fmpq_set_fmpz_frac(q,zn,zd);
                acb_set_fmpq(item,q,fp);
            } else acb_zero(item);
        }
        fmpz_clear(zn);
        fmpz_clear(zd);
        fmpq_clear(q);
    }
    
    void MQ::operator()(vector<vector<fmpz_poly_t>> & m) {
        int nr_m = m.size();
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            fmpz_poly_set(m[r][c],fmpz_poly_q_numref(M[r][c]));
        }
    }
    
    void MQ::operator()(vector<vector<fmpz_poly_q_t>> & m) {
        int nr_m = m.size();
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            fmpz_poly_q_set(m[r][c],M[r][c]);
        }
    }
    
    void MQ::operator()(vector<vector<acb_poly_t>> & m, slong fp) {
        int nr_m = m.size();
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            acb_poly_set_fmpz_poly(m[r][c],fmpz_poly_q_numref(M[r][c]),fp);
        }
    }
    
    //=*********************************************************************=
    
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
                _to_(xs,f,ctx,expr);
                res = _factor_(xs,f,ctx);
                fmpz_mpoly_clear(f,ctx);
                fmpz_mpoly_ctx_clear(ctx);
            }
            res = res.subs(map_rat,nopat);
            return res;
        }
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
    
    inline void _to_q_(const lst & xs, fmpz_mpoly_q_t f, fmpz_mpoly_ctx_t ctx, const ex & e) {
        if(is_a<add>(e)) {
            fmpz_mpoly_q_zero(f,ctx);
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            for(auto item : e) {
                _to_q_(xs,fi,ctx,item);
                fmpz_mpoly_q_add(f, f, fi, ctx);
            }
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(is_a<mul>(e)) {
            fmpz_mpoly_q_one(f,ctx);
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            for(auto item : e) {
                _to_q_(xs,fi,ctx,item);
                fmpz_mpoly_q_mul(f, f, fi, ctx);
            }
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::posint)) {
            ulong n = ex_to<numeric>(e.op(1)).to_int();
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            _to_q_(xs, fi, ctx, e.op(0));
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(f), fmpz_mpoly_q_numref(fi), n, ctx);
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(f), fmpz_mpoly_q_denref(fi), n, ctx);
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(is_a<power>(e) && e.op(1).info(info_flags::negint)) {
            ulong n = -ex_to<numeric>(e.op(1)).to_int();
            fmpz_mpoly_q_t fi;
            fmpz_mpoly_q_init(fi,ctx);
            _to_q_(xs, fi, ctx, e.op(0));
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(f), fmpz_mpoly_q_numref(fi), n, ctx);
            fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(f), fmpz_mpoly_q_denref(fi), n, ctx);
            fmpz_mpoly_q_inv(f,f,ctx);
            fmpz_mpoly_q_clear(fi,ctx);
            return;
        } else if(e.info(info_flags::rational)) {
            fmpq_t fq;
            fmpq_init(fq);
            _to_(fq,e);
            fmpz_mpoly_q_set_fmpq(f,fq,ctx);
            fmpq_clear(fq);
            return;
        } else if(e.is_polynomial(xs)) {
            string vars[xs.nops()];
            const char* cvars[xs.nops()];
            for(int i=0; i<xs.nops(); i++) {
                vars[i] = ex2str(xs.op(i));
                cvars[i] = vars[i].c_str();
            }
            fmpz_mpoly_q_one(f,ctx);
            if(fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_numref(f), ex2str(e).c_str(), cvars, ctx)) {
                cout << e << endl;
                cout << xs << endl;
                throw Error("_to_q_ failed.");
            }
            return;
        } else {
            cout << "expr = " << e << endl;
            throw Error("_to_q_ Not supported region");
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
                _to_q_(xs,f,ctx,e);
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
                _to_q_(xs,f,ctx,e);
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
                _to_q_(xs,f,ctx,e);
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
