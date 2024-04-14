
#include "exFlint.h"

namespace HepLib {

    namespace {
        
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
    
    MX::MX() { }
    MX::MX(const matrix & m) { init(m); }
    MX::MX(const MX & mx) { init(mx); }
    MX::MX(const vector<vector<fmpz_poly_q_struct*>> & m) { init(m); }
    
    MX::~MX() { clear(); }
        
    void MX::clear() {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_clear(M[r][c]); // fmpz_poly_q_init() @ init
                flint_free(M[r][c]); // flint_malloc @ init
            }
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
            M[r].resize(nc);
            for(int c=0; c<nc; c++) {
                M[r][c] = (fmpz_poly_q_struct *) flint_malloc(sizeof(fmpz_poly_q_struct)); // flint_free @ clear
                fmpz_poly_q_init(M[r][c]); // call fmpz_poly_q_clear() @ ~MX
                fmpz_poly_q_zero(M[r][c]);
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
            M[r].resize(nc);
            for(int c=0; c<nc; c++) {
                M[r][c] = (fmpz_poly_q_struct *) flint_malloc(sizeof(fmpz_poly_q_struct)); // flint_free @ clear
                fmpz_poly_q_init(M[r][c]); // call fmpz_poly_q_clear() @ ~MX
                fmpz_poly_q_set(M[r][c], mx.M[r][c]);
            }
        }
    }
    
    void MX::init(const vector<vector<fmpz_poly_q_struct*>> & m) {
        clear();
        nr = m.size();
        if(nr<1) return;
        nc = m[0].size();
        if(nc<1) return;
        M.resize(nr);
        for(int r=0; r<nr; r++) {
            M[r].resize(nc);
            for(int c=0; c<nc; c++) {
                M[r][c] = (fmpz_poly_q_struct *) flint_malloc(sizeof(fmpz_poly_q_struct)); // flint_free @ clear
                fmpz_poly_q_init(M[r][c]); // call fmpz_poly_q_clear() @ ~MX
                fmpz_poly_q_set(M[r][c],m[r][c]);
            }
        }
    }
    
    MX & MX::add(const MX & mx) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_add(M[r][c], M[r][c], mx.M[r][c]);
            }
        }
        return *this;
    }
    
    MX & MX::sub(const MX & mx) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_sub(M[r][c], M[r][c], mx.M[r][c]);
            }
        }
        return *this;
    }
    
    MX & MX::mul(const MX & mx) { 
        if(this!=&mx && mx.nc==nc) {
            vector<fmpz_poly_q_struct*> row(nc); // keep a row
            for(int c=0; c<nc; c++) {
                row[c] = (fmpz_poly_q_struct*) flint_malloc(sizeof(fmpz_poly_q_struct));
                fmpz_poly_q_init(row[c]);
            }
            for(int r=0; r<nr; r++) {
                for(int c=0; c<nc; c++) fmpz_poly_q_set(row[c], M[r][c]);
                for(int c=0; c<nc; c++) {
                    fmpz_poly_q_zero(M[r][c]);
                    for(int i=0; i<nc; i++) fmpz_poly_q_addmul(M[r][c], row[i], mx.M[i][c]);
                }
            }
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_clear(row[c]);
                flint_free(row[c]);
            }
        } else {
            vector<vector<fmpz_poly_q_struct*>> mat(nr);
            int nc2 = mx.nc;
            for(int r=0; r<nr; r++) {
                mat[r].resize(nc2);
                for(int c=0; c<nc2; c++) {
                    mat[r][c] = (fmpz_poly_q_struct*) flint_malloc(sizeof(fmpz_poly_q_struct));
                    fmpz_poly_q_init(mat[r][c]);
                    fmpz_poly_q_zero(mat[r][c]);
                    for(int i=0; i<nc; i++) fmpz_poly_q_addmul(mat[r][c], M[r][i], mx.M[i][c]);
                }
            }
            init(mat);
            for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) {
                fmpz_poly_q_clear(mat[r][c]);
                flint_free(mat[r][c]);
            }
        }
        return *this;
    }
    
    MX & MX::mul_left(const MX & mx) {
        if(this!=&mx && mx.nr==nr) {
            vector<fmpz_poly_q_struct*> col(nr);
            for(int r=0; r<nr; r++) {
                col[r] = (fmpz_poly_q_struct*) flint_malloc(sizeof(fmpz_poly_q_struct));
                fmpz_poly_q_init(col[r]);
            }
            for(int c=0; c<nc; c++) {
                for(int r=0; r<nr; r++) fmpz_poly_q_set(col[r], M[r][c]);
                for(int r=0; r<nr; r++) {
                    fmpz_poly_q_zero(M[r][c]);
                    for(int i=0; i<nr; i++) fmpz_poly_q_addmul(M[r][c], mx.M[r][i], col[i]);
                }
            }
            for(int r=0; r<nr; r++) {
                fmpz_poly_q_clear(col[r]);
                flint_free(col[r]);
            }
        } else {
            int nr2 = mx.nr;
            vector<vector<fmpz_poly_q_struct*>> mat(nr2);
            for(int r=0; r<nr2; r++) {
                mat[r].resize(nc);
                for(int c=0; c<nc; c++) {
                    mat[r][c] = (fmpz_poly_q_struct*) flint_malloc(sizeof(fmpz_poly_q_struct));
                    fmpz_poly_q_init(mat[r][c]);
                    fmpz_poly_q_zero(mat[r][c]);
                    for(int i=0; i<nr; i++) fmpz_poly_q_addmul(mat[r][c],mx.M[r][i],M[i][c]);
                }
            }
            init(mat);
            for(int r=0; r<nr2; r++) for(int c=0; c<nc; c++) {
                fmpz_poly_q_clear(mat[r][c]);
                flint_free(mat[r][c]);
            }
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
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_mul(M[r][c], M[r][c], f);
            }
        }
        fmpz_poly_q_clear(f);
        return *this;
    }
    
    MX & MX::scale(fmpz_poly_q_t f) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_mul(M[r][c],M[r][c],f);
            }
        }
        return *this;
    }
    
    MX & MX::scale(fmpz_poly_t f) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto & nref = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_mul(nref, nref, f);
                fmpz_poly_q_canonicalise(M[r][c]);
            }
        }
        return *this;
    }
    
    MX & MX::dx() {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_derivative(M[r][c], M[r][c]);
            }
        }
        return *this;
    }
        
    MX & MX::balance(const matrix & P) { // balance between 0 and infinity
        if(nr!=nc) throw Error("nr!=nc");
        static symbol x("xyz");
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
        fmpz_t z;
        fmpz_init(z);
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto num_ref = fmpz_poly_q_numref(M[r][c]);
                auto nn = fmpz_poly_length(num_ref);
                auto den_ref = fmpz_poly_q_denref(M[r][c]);
                auto dn = fmpz_poly_length(den_ref);
                slong nr2 = nn;
                if(nn>dn) {
                    nr2 = nn;
                    fmpz_pow_ui(z,qd,nn-dn);
                    fmpz_poly_scalar_mul_fmpz(den_ref, den_ref, z);
                } else if(dn>nn) {
                    nr2 = dn;
                    fmpz_pow_ui(z,qd,dn-nn);
                    fmpz_poly_scalar_mul_fmpz(num_ref, num_ref, z);
                }
                for(int i=0; i<nn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(num_ref, i);
                    if(ci==NULL) throw Error("ci == NULL");
                    fmpz_pow_ui(z,qd,nn-i);
                    fmpz_mul(ci, ci, z);
                }
                for(int i=0; i<dn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(den_ref, i);
                    fmpz_pow_ui(z,qd,dn-i);
                    fmpz_mul(ci, ci, z);
                }
                fmpz* rs;
                rs = _fmpz_vec_init(nr2);
                for(int i=0; i<nr2; i++) fmpz_set(rs+i,qn);
                _fmpz_poly_newton_to_monomial(num_ref->coeffs, rs, nn);
                _fmpz_poly_newton_to_monomial(den_ref->coeffs, rs, dn);
                for(int i=0; i<nn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(num_ref, i);
                    fmpz_pow_ui(z, qd, i);
                    fmpz_mul(ci, ci, z);
                }
                for(int i=0; i<dn; i++) {
                    auto ci = fmpz_poly_get_coeff_ptr(den_ref, i);
                    fmpz_pow_ui(z, qd, i);
                    fmpz_mul(ci, ci, z);
                }
                _fmpz_vec_clear(rs, nr2);
                fmpz_poly_q_canonicalise(M[r][c]);
            }
        }
        fmpz_clear(z);
        fmpz_clear(qn);
        fmpz_clear(qd);
        return *this;
    }
    
    matrix MX::operator()(const ex & x) {
        if(nr<1||nc<1) return matrix();
        matrix rmat(nr,nc);
        symtab st;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            rmat(r,c) = _to_(x, M[r][c]);
        }
        return rmat;
    }
    
    void MX::operator()(vector<vector<fmpz_poly_struct*>> & m) {
        int nr_m = m.size();
        if(nr_m<1) throw Error("MX: dimension < 1.");
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            fmpz_poly_set(m[r][c], fmpz_poly_q_numref(M[r][c]));
        }
    }
    
    void MX::operator()(vector<vector<fmpz_poly_q_struct*>> & m) {
        int nr_m = m.size();
        if(nr_m<1) throw Error("MX: dimension < 1.");
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            fmpz_poly_q_set(m[r][c], M[r][c]);
        }
    }
    
    void MX::operator()(vector<vector<acb_poly_struct*>> & m, slong fp) {
        int nr_m = m.size();
        if(nr_m<1) throw Error("MX: dimension < 1.");
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            acb_poly_set_fmpz_poly(m[r][c], fmpz_poly_q_numref(M[r][c]),fp);
        }
    }
    
    int MX::prank() {
        int pr;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto ipr = ldegree(fmpz_poly_q_numref(M[r][c]));
            if(ipr==0) ipr = -ldegree(fmpz_poly_q_denref(M[r][c]));
            if(r==0 && c==0) pr = ipr;
            else if(ipr<pr) pr = ipr;
        }
        return -1-pr;
    }
    
    matrix MX::a0() {
        auto pr = prank();
        matrix m(nr,nc);
        if(pr<0) return m;
        fmpz_t z;
        fmpz_init(z);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto num = fmpz_poly_q_numref(M[r][c]);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto ipr_num = ldegree(num);
            auto ipr_den = ldegree(den);
            if(ipr_den-ipr_num<pr+1) m(r,c)=0;
            else if(ipr_num>0) {
                if(ipr_den!=0) throw Error("a0: something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num);
                m(r,c) = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,0);
                m(r,c) /= _to_(z);
            } else {
                if(ipr_num!=0) throw Error("a0: something is wrong here.");
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
        auto pr = prank();
        if(pr<0) pr = 0;
        fmpz_t z;
        fmpz_init(z);
        matrix m0(nr,nc), m1(nr,nc);
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto num = fmpz_poly_q_numref(M[r][c]);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto ipr_num = ldegree(num);
            auto ipr_den = ldegree(den);
            if(ipr_den-ipr_num<pr) { m0(r,c)=0; m1(r,c)=0; continue; }
            ex n0, n1, d0, d1;
            if(ipr_num>0) {
                if(ipr_den!=0) throw Error("a01: something is wrong here.");
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num);
                n0 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,num,ipr_num+1);
                n1 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,0);
                d0 = _to_(z);
                fmpz_poly_get_coeff_fmpz(z,den,1);
                d1 = _to_(z);
            } else {
                if(ipr_num!=0) throw Error("a01: something is wrong here.");
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
        fmpz_t z;
        fmpz_init(z);
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
            fmpq_poly_get_denominator(z,num_q);
            if(fmpz_is_zero(z)) throw Error("MX::series, den is zero.");
            fmpz_poly_set_fmpz(den,z);
            if(dn>0) fmpz_poly_shift_left(den,den,dn);
            fmpq_poly_clear(den_q);
            fmpq_poly_clear(num_q);
            fmpz_poly_q_canonicalise(M[r][c]);
        }
        fmpz_clear(z);
    }
    
    int MX::denlcm(fmpz_poly_t dl) { // M will be updated
        fmpz_poly_one(dl);
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_lcm(dl, dl, fmpz_poly_q_denref(M[r][c]));
            }
        }
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto & nref = fmpz_poly_q_numref(M[r][c]);
                auto & dref = fmpz_poly_q_denref(M[r][c]);
                fmpz_poly_div(dref, dl, dref);
                fmpz_poly_mul(nref, nref, dref);
                fmpz_poly_one(dref);
                fmpz_poly_q_canonicalise(M[r][c]);
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
            fmpz_poly_get_coeff_fmpz(zd, den, ldeg);
            if(fmpz_is_zero(zd)) throw Error("denominator is zero.");
            if(i+ldeg>=0) {
                fmpz_poly_get_coeff_fmpz(zn, num, (slong)(i+ldeg));
                fmpq_set_fmpz_frac(item, zn, zd);
            } else fmpq_zero(item);
        }
        fmpz_clear(zn);
        fmpz_clear(zd);
    }
    
    void MX::coeff(acb_mat_t m, int i, slong fp) {
        fmpz_t zn,zd;
        fmpq_t q;
        fmpz_init(zn);
        fmpz_init(zd);
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
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_clear(M[r][c]); // fmpz_poly_q_init() @ init
                flint_free(M[r][c]); // flint_malloc @ init
            }
        }
        nr=-1; nc=-1;
    }
    
    void MQ::init(const vector<vector<fmpz_poly_q_struct*>> & m) {
        clear();
        nr = m.size();
        if(nr<1) throw Error("MQ::init, nr<1");
        nc = m[0].size();
        if(nc<1) throw Error("MQ::init, nc<1");
        M.resize(nr);
        for(int r=0; r<nr; r++) {
            M[r].resize(nc);
            for(int c=0; c<nc; c++) {
                M[r][c] = (fmpz_poly_q_struct*) flint_malloc(sizeof(fmpz_poly_q_struct)); // flint_free @ clear
                fmpz_poly_q_init(M[r][c]); // call fmpz_poly_q_clear() @ MQ
                fmpz_poly_q_set(M[r][c], m[r][c]);
            }
        }
    }
    
    void MQ::scale(fmpz_poly_q_t f) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_q_mul(M[r][c], M[r][c], f);
            }
        }
    }
    
    void MQ::scale(fmpz_poly_t f) {
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto & nref = fmpz_poly_q_numref(M[r][c]);
                fmpz_poly_mul(nref, nref, f);
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
                fmpz_poly_q_print_pretty(M[r][c], "x");
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
        fmpz_poly_one(dl);
        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                fmpz_poly_lcm(dl, dl, fmpz_poly_q_denref(M[r][c]));
            }
        }

        for(int r=0; r<nr; r++) {
            for(int c=0; c<nc; c++) {
                auto & nref = fmpz_poly_q_numref(M[r][c]);
                auto & dref = fmpz_poly_q_denref(M[r][c]);
                fmpz_poly_div(dref, dl, dref);
                fmpz_poly_mul(nref, nref, dref);
                fmpz_poly_one(dref);
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
        fmpz_t zn, zd;
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
    
    void MQ::coeff(gr_mat_t m, int i, gr_ctx_t ctx) {
        fmpz_t zn,zd;
        fmpz_init(zn);
        fmpz_init(zd);
        fmpq_t q;
        fmpq_init(q);
        int status = GR_SUCCESS;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            auto item = gr_mat_entry_ptr(m, r, c, ctx);
            auto den = fmpz_poly_q_denref(M[r][c]);
            auto num = fmpz_poly_q_numref(M[r][c]);
            slong ldeg = ldegree(den);
            if(fmpz_poly_degree(den)!=ldeg) throw Error("coeff: not a monomial in denominator.");
            fmpz_poly_get_coeff_fmpz(zd,den,ldeg);
            if(fmpz_is_zero(zd)) throw Error("denominator is zero.");
            if(i+ldeg>=0) {
                fmpz_poly_get_coeff_fmpz(zn, num, (slong)(i+ldeg));
                fmpq_set_fmpz_frac(q, zn, zd);
                status |= gr_set_fmpq(item, q, ctx);
            } else status |= gr_zero(item, ctx);
        }
        fmpz_clear(zn);
        fmpz_clear(zd);
        fmpq_clear(q);
        if(status != GR_SUCCESS) throw Error("coeff: status != GR_SUCCESS!");
    }
    
    void MQ::operator()(vector<vector<fmpz_poly_struct*>> & m) {
        int nr_m = m.size();
        if(nr_m<1) throw Error("MQ: nr_m<1.");
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            fmpz_poly_set(m[r][c], fmpz_poly_q_numref(M[r][c]));
        }
    }
    
    void MQ::operator()(vector<vector<fmpz_poly_q_struct*>> & m) {
        int nr_m = m.size();
        if(nr_m<1) throw Error("MQ: nr_m<1.");
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            fmpz_poly_q_set(m[r][c], M[r][c]);
        }
    }
    
    void MQ::operator()(vector<vector<acb_poly_struct*>> & m, slong fp) {
        int nr_m = m.size();
        if(nr_m<1) throw Error("MQ: nr_m<1.");
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            acb_poly_set_fmpz_poly(m[r][c],fmpz_poly_q_numref(M[r][c]),fp);
        }
    }
    
    void MQ::operator()(vector<vector<gr_poly_struct*>> & m, gr_ctx_t ctx) {
        int nr_m = m.size();
        if(nr_m<1) throw Error("MQ: nr_m<1.");
        int nc_m = m[0].size();
        if(nr!=nr_m || nc!=nc_m) throw Error("matrix dimension not match.");
        int status = GR_SUCCESS;
        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
            if(!fmpz_poly_is_one(fmpz_poly_q_denref(M[r][c]))) throw Error("the denominator is not 1.");
            status |= gr_poly_set_fmpz_poly(m[r][c], fmpz_poly_q_numref(M[r][c]), ctx);
        }
        if(status != GR_SUCCESS) throw Error("operator(): status != GR_SUCCESS");
    }

    
}
