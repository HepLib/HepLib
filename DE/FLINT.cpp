
#include "FLINT.h"
#include "cln/cln.h"

namespace HepLib {
    
    namespace {
        inline ex syms(const ex & e) {
            exset ret;
            for(const_preorder_iterator i=e.preorder_begin(); i!=e.preorder_end(); ++i) 
                if(is_a<symbol>(*i)) ret.insert(*i); 
            if(ret.size()<1) return 0;
            if(ret.size()>1) throw Error("more symbols found!");
            return *(ret.begin());
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
    
    slong n2exp(const ex & e) {
        if(is_zero(e)) return 0;
        auto cn = ex_to<numeric>(abs(e)).to_cl_N();
        auto rcn = cln::realpart(cn);
        return cln::cl_I_to_long(cln::ceiling1(cln::ln(abs(rcn))/cln::ln(2)));
    }
    
    slong get_rel_err(arb_t r) {
        return fp2dp(arb_rel_error_bits(r));
    }
    slong get_rel_err(acb_t z) {
        arb_t t;
        arb_init(t);
        acb_get_real(t,z);
        auto s1 = get_rel_err(t);
        acb_get_imag(t,z);
        auto s2 = get_rel_err(t);
        if(s1>s2) return s1;
        arb_clear(t);
        return s2;
    }
    slong get_rel_err(arb_mat_t m) {
        slong s = -ARF_PREC_EXACT;
        auto nr = arb_mat_nrows(m);
        auto nc = arb_mat_ncols(m);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto si = get_rel_err(arb_mat_entry(m,r,c));
            if(si>s) s = si;
        }
        return s;
    }
    slong get_rel_err(acb_mat_t m) {
        slong s = -ARF_PREC_EXACT;
        auto nr = arb_mat_nrows(m);
        auto nc = arb_mat_ncols(m);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto si = get_rel_err(acb_mat_entry(m,r,c));
            if(si>s) s = si;
        }
        return s;
    }
    
    void ex_to_flint(fmpz_poly_q_t f, const ex & expr) {
        if(is_a<symbol>(expr)) {
            fmpz_poly_q_set_str(f, "2  0 1/1  1");
            return;
        } else if(is_a<numeric>(expr)) {
            auto nd = expr.numer_denom();
            ostringstream oss;
            oss << "1  " << nd.op(0) << "/1  " << nd.op(1);
            fmpz_poly_q_set_str(f, oss.str().c_str());
            return;
        } else if(is_a<add>(expr)) {
            fmpz_poly_q_zero(f);
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            for(auto item : expr) {
                ex_to_flint(fi, item);
                fmpz_poly_q_add(f, f, fi);
            }
            fmpz_poly_q_clear(fi);
            return;
        } else if(is_a<mul>(expr)) {
            fmpz_poly_q_one(f);
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            for(auto item : expr) {
                ex_to_flint(fi, item);
                fmpz_poly_q_mul(f, f, fi);
            }
            fmpz_poly_q_clear(fi);
            return;
        } else if(is_a<power>(expr) && expr.op(1)>0) {
            ulong n = ex_to<numeric>(expr.op(1)).to_int();
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            ex_to_flint(fi, expr.op(0));
            fmpz_poly_q_pow(f, fi, n);
            fmpz_poly_q_clear(fi);
            return;
        } else if(is_a<power>(expr) && expr.op(1)<0) {
            ulong n = -ex_to<numeric>(expr.op(1)).to_int();
            fmpz_poly_q_t fi;
            fmpz_poly_q_init(fi);
            ex_to_flint(fi, expr.op(0));
            fmpz_poly_q_pow(f, fi, n);
            fmpz_poly_q_inv(f,f);
            fmpz_poly_q_clear(fi);
            return;
        } else throw Error("Not supported region");
    }
        
    void ex_to_fmpq(fmpq_t q, const ex & expr) {
        auto str = ex2str(expr);
        fmpq_set_str(q,str.c_str(),10);
    }
    
    void ex_to_fmpz(fmpz_t q, const ex & expr) {
        auto str = ex2str(expr);
        fmpz_set_str(q,str.c_str(),10);
    }
    
    ex fmpq_to_ex(fmpq_t q) {
        auto cstr = fmpq_get_str(NULL,10,q);
        ex n = numeric(cstr);
        flint_free(cstr);
        return n;
    }
    
    void mat_to_fmpq(fmpq_mat_t m, const matrix & mat) {
        auto nr = mat.rows();
        auto nc = mat.rows();
        fmpq_t q;
        fmpq_init(q);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            ex_to_fmpq(q, mat(r,c));
            fmpq_set(fmpq_mat_entry(m,r,c), q);
        }
        fmpq_clear(q);
    }
    
    matrix fmpq_to_mat(fmpq_mat_t m) {
        auto nr = fmpq_mat_nrows(m);
        auto nc = fmpq_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = fmpq_to_ex(fmpq_mat_entry(m,r,c));
        }
        return mat;
    }
    
    void ex_to_arb(arb_t r, const ex & expr, slong dp) {
        set_precision(dp);
        auto ne = expr.evalf();
        if(!is_a<numeric>(ne)) {
            cout << endl << "ne = " << ne << endl;
            throw Error("ne is NOT numeric");
        }
        auto fp = dp2fp(dp);
        auto str = ex2str(ne);
        arb_set_str(r, str.c_str(), fp);
        reset_precision();
    }
    
    ex arb_to_ex(arb_t r, slong dp) {
        set_precision(dp);
        if(!error_check(r)) throw Error("too low precision.");
        auto fp = dp2fp(dp);
        auto cstr = arb_get_str(r,fp,ARB_STR_NO_RADIUS);
        string str = cstr;
        flint_free(cstr);
        if(string_contain(str,".") || string_contain(str,"e")) str += "_"+to_string(dp);
        ex n = numeric(cln::cl_R(str.c_str()));
        reset_precision();
        return n;
    }
    
    void mat_to_arb(arb_mat_t m, const matrix & mat, slong dp) {
        auto nr = mat.rows();
        auto nc = mat.rows();
        arb_t z;
        arb_init(z);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            ex_to_arb(z, mat(r,c), dp);
            arb_set(arb_mat_entry(m,r,c), z);
        }
        arb_clear(z);
    }
    
    matrix arb_to_mat(arb_mat_t m, slong dp) {
        auto nr = arb_mat_nrows(m);
        auto nc = arb_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = arb_to_ex(arb_mat_entry(m,r,c), dp);
        }
        return mat;
    }
    
    bool error_check(arb_t r) {
        auto dp = error_check_dp;
        auto fp = dp2fp(dp);
        if(get_rel_err(r)>-dp) {
            mag_t mag;
            arb_get_mag(mag,r);
            if(mag_cmp_2exp_si(mag,-fp)<0) return true;
            cout << endl; arb_printd(r,30); cout << endl;
            mag_clear(mag);
            return false;
        }
        return true;
    }
    
    void ex_to_acb(acb_t z, const ex & expr, slong dp) {
        set_precision(dp+10);
        auto ne = expr.evalf();
        if(!is_a<numeric>(ne)) {
            cout << endl << "ne = " << ne << endl;
            throw Error("ne is NOT numeric");
        }
        arb_t rb, ib;
        arb_init(rb);
        arb_init(ib);
        auto nre = ne.real_part();
        auto nie = ne.imag_part();
        auto fp = dp2fp(dp);
        auto rstr = ex2str(nre);
        arb_set_str(rb, rstr.c_str(), fp);
        auto istr = ex2str(nie);
        arb_set_str(ib, istr.c_str(), fp);
        acb_set_arb_arb(z, rb, ib);
        arb_clear(rb);
        arb_clear(ib);
        reset_precision();
    }
    
    ex acb_to_ex(acb_t z, slong dp) {
        set_precision(dp);
        auto fp = dp2fp(dp);
        arb_t ri;
        arb_init(ri);
        acb_get_real(ri,z);
        if(!error_check(ri)) throw Error("too low precision.");
        auto str = arb_get_str(ri,fp,ARB_STR_NO_RADIUS);
        string rstr = str;
        flint_free(str);
        if(string_contain(rstr,".") || string_contain(rstr,"e")) rstr += "_"+to_string(dp);
        acb_get_imag(ri,z);
        if(!error_check(ri)) throw Error("too low precision.");
        str = arb_get_str(ri,fp,ARB_STR_NO_RADIUS);
        string istr = str;
        flint_free(str);
        if(string_contain(istr,".") || string_contain(istr,"e")) istr += "_"+to_string(dp);
        arb_clear(ri);
        ex n = numeric(cln::complex(cln::cl_R(rstr.c_str()),cln::cl_R(istr.c_str())));
        reset_precision();
        return n;
    }
    
    void mat_to_acb(acb_mat_t m, const matrix & mat, slong dp) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        acb_t z;
        acb_init(z);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            ex_to_acb(z, mat(r,c), dp);
            acb_set(acb_mat_entry(m,r,c), z);
        }
        acb_clear(z);
    }
    
    matrix acb_to_mat(acb_mat_t m, slong dp) {
        auto nr = acb_mat_nrows(m);
        auto nc = acb_mat_ncols(m);
        matrix mat(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            mat(r,c) = acb_to_ex(acb_mat_entry(m,r,c), dp);
        }
        return mat;
    }
    
    ex flint_to_ex(fmpz_poly_q_t f, const ex & x) {
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
    
    ex flint_to_ex(fmpz_t f, const ex & x) {
        auto str = fmpz_get_str(NULL,10,f);
        auto res = numeric(str);
        flint_free(str);
        return res;
    }
    
    ex flint_to_ex(fmpz_poly_t f, const ex & x) {
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
    
    MX::MX() { }
    MX::MX(const matrix & m) { init(m); }
    MX::MX(const MX & mx) { init(mx); }
    
    MX::~MX() { clear(); n=0; s=-1; }
        
    void MX::clear() {
        if(n<0) return;
        fmpz_poly_clear(Qx);
        if(M.size()>0) {
            for(int r=0; r<n; r++) {
                for(int c=0; c<n; c++) fmpz_poly_q_clear(M[r][c]);
            }
        }
        if(s>1 && _s_>0) {
            for(int i=0; i<=s; i++) {
                fmpz_mat_clear(QxM[i]);
                fmpz_clear(Qs[i]);
            }
        }
        n = -1;
        s = -1;
    }
    
    void MX::init(const matrix & m) {
        if(n>0) throw Error("MX already inited.");
        n = m.rows();
        if(m.cols()!=n) throw Error("MX only works for squared matrix.");
        M.resize(n);
        for(int r=0; r<n; r++) {
            M[r] = vector<fmpz_poly_q_t>(n);
            for(int c=0; c<n; c++) {
                fmpz_poly_q_init(M[r][c]);
                ex_to_flint(M[r][c],m(r,c));
            }
        }
        fmpz_poly_init(Qx);
    }
    
    void MX::init(const MX & mx) {
        if(n>0) throw Error("MX already inited.");
        if(mx.n<1) return;
        n = mx.n;
        fmpz_poly_init(Qx);
        fmpz_poly_set(Qx,mx.Qx);
        if(mx.M.size()>0) { // only copy Qm/QxM if s>0
            M.resize(n);
            for(int r=0; r<n; r++) {
                M[r] = vector<fmpz_poly_q_t>(n);
                for(int c=0; c<n; c++) {
                    fmpz_poly_q_init(M[r][c]);
                    fmpz_poly_q_set(M[r][c],mx.M[r][c]);
                }
            }
        } else if(mx.s>0) {
            s = mx.s;
            Qs = vector<fmpz_t>(s+1);
            QxM = vector<fmpz_mat_t>(s+1);
            for(int i=0; i<=s; i++) {
                fmpz_init(Qs[i]);
                fmpz_mat_init(QxM[i],n,n);
                fmpz_set(Qs[i], mx.Qs[i]);
                fmpz_mat_set(QxM[i], mx.QxM[i]);
            }
        }
    }
    
    MX & MX::add(const MX & mx) {
        for(int r=0; r<n; r++) {
            #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
            for(int c=0; c<n; c++) fmpz_poly_q_add(M[r][c],M[r][c],mx.M[r][c]);
        }
        return *this;
    }
    
    MX & MX::sub(const MX & mx) {
        for(int r=0; r<n; r++) {
            #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
            for(int c=0; c<n; c++) fmpz_poly_q_sub(M[r][c],M[r][c],mx.M[r][c]);
        }
        return *this;
    }
    
    MX & MX::mul(const MX & mx) { 
        if(this==&mx) {
            MX mx2(mx);
            return mul(mx2);
        }
        vector<fmpz_poly_q_t> row(n);
        for(int c=0; c<n; c++) fmpz_poly_q_init(row[c]);
        for(int r=0; r<n; r++) {
            for(int c=0; c<n; c++) fmpz_poly_q_set(row[c], M[r][c]);
            #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
            for(int c=0; c<n; c++) {
                fmpz_poly_q_zero(M[r][c]);
                for(int i=0; i<n; i++) fmpz_poly_q_addmul(M[r][c],row[i],mx.M[i][c]);
            }
        }
        for(int c=0; c<n; c++) fmpz_poly_q_clear(row[c]);
        return *this;
    }
    
    MX & MX::mul_left(const MX & mx) {
        if(this==&mx) {
            MX mx2(mx);
            return mul_left(mx2);
        }
        vector<fmpz_poly_q_t> col(n);
        for(int r=0; r<n; r++) fmpz_poly_q_init(col[r]);
        for(int c=0; c<n; c++) {
            for(int r=0; r<n; r++) fmpz_poly_q_set(col[r], M[r][c]);
            #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
            for(int r=0; r<n; r++) {
                fmpz_poly_q_zero(M[r][c]);
                for(int i=0; i<n; i++) fmpz_poly_q_addmul(M[r][c],mx.M[r][i],col[i]);
            }
        }
        for(int r=0; r<n; r++) fmpz_poly_q_clear(col[r]);
        return *this;
    }
    
    MX & MX::add(const matrix & m) { MX mx(m); return add(mx); }
    MX & MX::sub(const matrix & m) { MX mx(m); return sub(mx); }
    MX & MX::mul(const matrix & m) { MX mx(m); return mul(mx); }
    MX & MX::mul_left(const matrix & m) { MX mx(m); return mul_left(mx); }
    
    MX & MX::scale(const ex & s) {
        fmpz_poly_q_t f;
        fmpz_poly_q_init(f);
        ex_to_flint(f,s);
        for(int r=0; r<n; r++) {
            #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
            for(int c=0; c<n; c++) fmpz_poly_q_mul(M[r][c],M[r][c],f);
        }
        fmpz_poly_q_clear(f);
        return *this;
    }
    
    MX & MX::dx() {
        for(int r=0; r<n; r++) {
            #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
            for(int c=0; c<n; c++) fmpz_poly_q_derivative(M[r][c],M[r][c]);
        }
        return *this;
    }
        
    MX & MX::balance(const matrix & P) {
        static symbol x("x");
        MX mx1(P);
        MX mx2(mx1);
        mx1.scale(x);
        mx2.scale(1/x);
        auto coP = ex_to<matrix>(unit_matrix(n)).sub(P);
        MX mx3(coP);
        auto mx4(mx3);
        mx3.sub(mx1);
        mx4.sub(mx2);
        return mul_left(mx3).mul(mx4).add(mx2);
    }

    MX & MX::transform(const matrix & t, const matrix & ti) {
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
        ex_to_fmpz(qn,nd.op(0));
        ex_to_fmpz(qd,nd.op(1));
        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
        for(int r=0; r<n; r++) {
            fmpz_t z;
            fmpz_init(z);
            for(int c=0; c<n; c++) {
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
        }
        fmpz_clear(qn); 
        fmpz_clear(qd);
        return *this;
    }
    
    matrix MX::operator()(const ex & x) {
        if(n<1) return matrix();
        matrix rmat(n,n);
        symtab st;
        for(int r=0; r<n; r++) for(int c=0; c<n; c++) {
            rmat(r,c) = flint_to_ex(M[r][c],x);
        }
        return rmat;
    }
    
    int MX::prank() {
        int pr;
        for(int r=0; r<n; r++) for(int c=0; c<n; c++) {
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
        matrix m(n,n);
        for(int r=0; r<n; r++) for(int c=0; c<n; c++) {
            auto num = fmpz_poly_q_numref(M[r][c]);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto ipr_num = ldegree(num);
            auto ipr_den = ldegree(den);
            if(ipr_den-ipr_num<pr+1) m(r,c)=0;
            else if(ipr_num>0) {
                if(ipr_den!=0) throw Error("something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num);
                m(r,c) = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,den,0);
                m(r,c) /= flint_to_ex(z,1);
            } else {
                if(ipr_num!=0) throw Error("something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,0);
                m(r,c) = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,den,ipr_den);
                m(r,c) /= flint_to_ex(z,1);
            }
        }
        fmpz_clear(z);
        return m;
    }
    
    pair<matrix,matrix> MX::a01() {
        fmpz_t z;
        fmpz_init(z);
        auto pr = prank();
        matrix m0(n,n), m1(n,n);
        for(int r=0; r<n; r++) for(int c=0; c<n; c++) {
            auto num = fmpz_poly_q_numref(M[r][c]);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto ipr_num = ldegree(num);
            auto ipr_den = ldegree(den);
            if(ipr_den-ipr_num<pr) { m0(r,c)=0; m1(r,c)=0; continue; }
            ex n0, n1, d0, d1;
            if(ipr_num>0) {
                if(ipr_den!=0) throw Error("1: something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num);
                n0 = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num+1);
                n1 = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,den,0);
                d0 = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,den,1);
                d1 = flint_to_ex(z,1);
            } else {
                if(ipr_num!=0) throw Error("2: something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,0);
                n0 = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,num,1);
                n1 = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,den,ipr_den);
                d0 = flint_to_ex(z,1);
                fmpz_poly_get_coeff_fmpz(z,den,ipr_den+1);
                d1 = flint_to_ex(z,1);
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
        for(int r=0; r<n; r++) for(int c=0; c<n; c++) {
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            if(fmpz_poly_length(den)!=1) throw Error("degree: denominator is NOT constant");
            auto cdeg = fmpz_poly_degree(num);
            if(r==0 && c==0) deg = cdeg;
            else if(deg<cdeg) deg = cdeg;
        }
        return deg;
    }
    
    matrix MX::coeff(int i,bool ef) {
        fmpz_t z;
        fmpz_init(z);
        matrix rmat(n,n);
        for(int r=0; r<n; r++) for(int c=0; c<n; c++) {
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            if(fmpz_poly_length(den)!=1) throw Error("degree: denominator is NOT constant");
            fmpz_poly_get_coeff_fmpz(z,den,0);
            ex _den = flint_to_ex(z,1);
            if(ef) _den = evalf(_den);
            fmpz_poly_get_coeff_fmpz(z,num,(slong)i);
            ex _num = flint_to_ex(z,1);
            if(ef) _num = evalf(_num);
            rmat(r,c) = _num/_den;
        }
        fmpz_clear(z);
        return rmat;
    }
    
    void MX::lcm(bool init_QxM) { // Q is LCM of denominators in x*M, M is now x*M
        fmpz_poly_t x;
        fmpz_poly_init(x);
        vector<fmpz_poly_t> lcms(n);
        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
        for(int r=0; r<n; r++) {
            auto lcm = lcms[r];
            fmpz_poly_init(lcm);
            fmpz_poly_set_str(lcm, "1  1");
            for(int c=0; c<n; c++) {
                auto num = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_shift_left(num,num,1); // x*M
                fmpz_poly_q_canonicalise(M[r][c]);
                auto den = fmpz_poly_q_denref(M[r][c]);
                fmpz_poly_lcm(lcm, lcm, den);
            }
        }
        fmpz_poly_set_str(Qx, "1  1");
        for(int r=0; r<n; r++) {
            fmpz_poly_lcm(Qx,Qx,lcms[r]);
            fmpz_poly_clear(lcms[r]);
        }
        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
        for(int r=0; r<n; r++) {
            for(int c=0; c<n; c++) {
                auto num = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_mul(num, num, Qx);
                fmpz_poly_q_canonicalise(M[r][c]);
            }
        }
        
        // now initialize Qs & QxM
        s = fmpz_poly_degree(Qx);
        auto deg = degree();
        if(deg>s) s = deg;
        
        if(!init_QxM) return;
        
        _s_ = 1;
        Qs = vector<fmpz_t>(s+1);
        for(slong i=0; i<=s; i++) {
            fmpz_init(Qs[i]);
            fmpz_poly_get_coeff_fmpz(Qs[i],Qx,i);
        }
        
        QxM = vector<fmpz_mat_t>(s+1);
        for(int i=0; i<=s; i++) fmpz_mat_init(QxM[i],n,n); 
        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
        for(int r=0; r<n; r++) {
            fmpz_t z;
            fmpz_init(z);
            for(int c=0; c<n; c++) {
                auto den = fmpz_poly_q_denref(M[r][c]);
                auto num = fmpz_poly_q_numref(M[r][c]);
                if(fmpz_poly_length(den)!=1) throw Error("inconsistant: denominator is NOT constant.");
                fmpz_poly_get_coeff_fmpz(z,den,0);
                if(!fmpz_is_one(z)) throw Error("inconsistant: denominator is NOT one.");
                for(slong m=0; m<=s; m++) {
                    fmpz_poly_get_coeff_fmpz(acb_mat_entry(QxM[m],r,c),num,m);
                }
            }
            fmpz_clear(z);
        }
        
        if(true) { // remove M to save mem
            for(int r=0; r<n; r++) {
                for(int c=0; c<n; c++) fmpz_poly_q_clear(M[r][c]);
            }
            M.clear();
            M.shrink_to_fit();
        }
    }
    
    // note that we rescale with q0
    void MX::Q(arb_t Qm, slong m, slong fp) {
        arb_set_fmpz(Qm, Qs[m]);
        arb_div_fmpz(Qm, Qm, Qs[0], fp);
    }
    
    // QxM[m] - a*Q[m], note QxM/Q can be scaled by a factor, we choose q0
    void MX::B(arb_mat_t Bm, slong m, arb_t a, slong fp) {
        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
        for(int r=0; r<n; r++) {
            for(int c=0; c<n; c++) {   
                arb_set_fmpz(arb_mat_entry(Bm,r,c), fmpz_mat_entry(QxM[m],r,c));
                arb_div_fmpz(arb_mat_entry(Bm,r,c),arb_mat_entry(Bm,r,c),Qs[0],fp);
            }
            arb_t b;
            arb_init(b);
            arb_mul_fmpz(b,a,Qs[m],fp);
            arb_div_fmpz(b,b,Qs[0],fp);
            arb_sub(arb_mat_entry(Bm,r,r),arb_mat_entry(Bm,r,r),b,fp);
            arb_clear(b);
        }
    }
    
    // note that we rescale with q0
    void MX::Q(acb_t Qm, slong m, slong fp) {
        acb_set_fmpz(Qm, Qs[m]);
        acb_div_fmpz(Qm, Qm, Qs[0], fp);
    }
    
    // QxM[m] - a*Q[m], note QxM/Q can be scaled by a factor, we choose q0
    void MX::B(acb_mat_t Bm, slong m, acb_t a, slong fp) {
        if(fmpz_is_zero(Qs[0])) throw Error("Qs[0] is zero!");
        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
        for(int r=0; r<n; r++) {
            for(int c=0; c<n; c++) {
                acb_set_fmpz(acb_mat_entry(Bm,r,c), fmpz_mat_entry(QxM[m],r,c));
                acb_div_fmpz(arb_mat_entry(Bm,r,c),arb_mat_entry(Bm,r,c),Qs[0],fp);
            }
            acb_t b;
            acb_init(b);
            acb_mul_fmpz(b,a,Qs[m],fp);
            acb_div_fmpz(b,b,Qs[0],fp);
            acb_sub(acb_mat_entry(Bm,r,r),acb_mat_entry(Bm,r,r),b,fp);
            acb_clear(b);
        }
    }
    
    // note that we rescale with q0
    void MX::Q(fmpq_t Qm, slong m) { 
        fmpq_set_fmpz_frac(Qm, Qs[m], Qs[0]);
    }
    
    // QxM[m] - a*Q[m], note QxM/Q can be scaled by a factor, we choose q0
    void MX::B(fmpq_mat_t Bm, slong m, fmpq_t a) {
        if(fmpz_is_zero(Qs[0])) throw Error("Qs[0] is zero!");
        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
        for(int r=0; r<n; r++) {
            for(int c=0; c<n; c++) {
                fmpq_set_fmpz_frac(acb_mat_entry(Bm,r,c), fmpz_mat_entry(QxM[m],r,c), Qs[0]);
            }
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_fmpz_frac(q,Qs[m],Qs[0]);
            fmpq_mul(q,q,a);
            fmpq_sub(acb_mat_entry(Bm,r,r),acb_mat_entry(Bm,r,r),q);
            fmpq_clear(q);
        }
    }
    
    // other functions
    
    ex factor_flint(const ex & expr) {
        ex x = syms(expr);
        if(is_zero(x)) return expr;
        auto sx = ex2str(x).c_str();
        symtab st;
        st[sx] = x;
        fmpz_poly_q_t nd;
        fmpz_poly_q_init(nd);
        ex_to_flint(nd,expr); 
        
        auto num = fmpz_poly_q_numref(nd);
        auto den = fmpz_poly_q_denref(nd);
        if(fmpz_poly_length(den)!=1) throw Error("factor_flint: Denominator's degree is NOT 1.");
        ex res = flint_to_ex(den,x);
        
        fmpz_poly_factor_t fs;
        fmpz_poly_factor_init(fs);
        fmpz_poly_factor_zassenhaus(fs,num);
        
        fmpz_poly_t f;
        fmpz_poly_init(f);
        fmpz_poly_set_fmpz(f,&fs->c);
        res *= flint_to_ex(f,x);
        
        for (int i=0; i<fs->num; i++) {
            fmpz_poly_set(f,fs->p+i);
            ex fx = flint_to_ex(f,x);
            if(fs->exp[i]==1) res *= fx;
            else res *= pow(fx,fs->exp[i]);
        }
        
        fmpz_poly_clear(f);
        fmpz_poly_factor_clear(fs);
        fmpz_poly_q_clear(nd);
        return res;
    }
    
    ex den_lcm_flint(const ex & expr) {
        ex x = syms(expr);
        if(is_zero(x)) throw Error("no symbol found.");
        auto sx = ex2str(x).c_str();
        symtab st;
        st[sx] = x;
        
        fmpz_poly_t lcm;
        fmpz_poly_init(lcm);
        fmpz_poly_set_str(lcm, "1  1");
        
        fmpz_poly_q_t f;
        fmpz_poly_q_init(f);
        
        for(int i=0; i<expr.nops(); i++) {
            ex_to_flint(f, expr.op(i));
            fmpz_poly_lcm(lcm, lcm, fmpz_poly_q_denref(f));
        }
        
        ex res = flint_to_ex(lcm,x);
        fmpz_poly_q_clear(f);
        fmpz_poly_clear(lcm);
        return res;
    }
    
    ex normal_flint(const ex & expr) {
        ex x = syms(expr);
        if(is_zero(x)) return expr;
        auto sx = ex2str(x).c_str();
        symtab st;
        st[sx] = x;
        fmpz_poly_q_t f;
        fmpz_poly_q_init(f);
        ex_to_flint(f, expr);
        auto res = flint_to_ex(f,x);
        fmpz_poly_q_clear(f);
        return res;
    }
    
    
        
}
