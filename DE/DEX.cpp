/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"
#include "vspace.h"

namespace HepLib {

    //=*********************************************************************=
    
    DEX::DEX(const symbol & _x) : x(_x) { }
    DEX::~DEX() { clear(); }
    void DEX::clear() {
        int nbs = bs.size();
        if(nbs<1) return;
        for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
            for(auto & row : Mat[br][bc]) for(auto & item : row) fmpz_poly_q_clear(item);
            if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) fmpq_mat_clear(item[0]);
            if(taylor_inited) for(auto & row : QMat[br][bc]) for(auto & item : row) fmpz_poly_clear(item);
            if(ntaylor_inited) for(auto & row : TMat[br][bc]) for(auto & item : row) acb_poly_clear(item);
        }
        bs.clear();
        Mat.clear();
        if(taylor_inited) for(int br=0; br<nbs; br++) fmpz_poly_clear(QD[br][0]);
        if(ntaylor_inited) for(int br=0; br<nbs; br++) acb_poly_clear(TD[br][0]);
        if(fuchsified) { 
            U0.clear(); 
            for(auto & item : qlas) fmpq_clear(item[0]);
            qlas.clear();
        }
        if(taylor_inited) { QMat.clear(); QD.clear(); }
        if(ntaylor_inited) { TMat.clear(); TD.clear(); }
        fuchsified = false;
        taylor_inited = false;
        ntaylor_inited = false;
    }
    
    matrix DEX::T() {
        matrix t = ex_to<matrix>(unit_matrix(N));
        for(auto ti : Ts) t = t.mul(ti);
        return t;
    }
    
    matrix DEX::M() { 
        matrix mat(N,N);
        int nbs = bs.size();
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int bc=0; bc<=br; bc++) {
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    mat(r0+r,c0+c) = _to_(x,Mat[br][bc][r][c]);
                }
            }
        }
        return mat;
    }
        
    //=*********************************************************************=
    // U-Series - rational
    //=*********************************************************************=
    
    umat_t DEX::series(int xn, const lst & slas) { // RU[a][b][la][k][n]
        if(!fuchsified) fuchsify();
        vector<fmpq_t> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        
        block_umat_fmpq_mat_t U;
        series(U,xn,qslas);
        for(int i=0; i<qslas.size(); i++) fmpq_clear(qslas[i]);
    
        int nbs = bs.size();
        umat_t RU;
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int bc=0; bc<=br; bc++) { // cycle columns
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                for(auto & kv : U[br][bc]) {
                    auto ila = kv.first;
                    auto la = las[ila];
                    int kmax = kv.second.size();
                    if(RU.find(la)==RU.end()) {
                        int kmax = IK[nbs-1][ila];
                        RU[la].resize(kmax);
                        for(int k=0; k<kmax; k++) {
                            RU[la][k].resize(xn+1);
                            for(int n=0; n<=xn; n++) RU[la][k][n] = matrix(N,N);
                        }
                    }
                    for(int k=0; k<kmax; k++) {
                        for(int n=0; n<=xn; n++) {
                            auto mat = _to_(U[br][bc][ila][k][n]);
                            fmpq_mat_clear(U[br][bc][ila][k][n]);
                            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RU[la][k][n](r0+r,c0+c) = mat(r,c);
                        }
                    }
                }
            }
        }
        return RU;
    }
    
    
    //=*********************************************************************=
    // U-Series - acb
    //=*********************************************************************=
    
    umat_t DEX::series(int xn, slong dp, const lst & slas) { // RU[a][b][la][k][n]
        if(!fuchsified) fuchsify();
        vector<fmpq_t> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        block_umat_acb_mat_t U;
        series(U, xn,dp,qslas);
        for(int i=0; i<qslas.size(); i++) fmpq_clear(qslas[i]);
        
        int nbs = bs.size();
        auto fp = dp2fp(dp);
        umat_t RU;
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int bc=br; bc>=0; bc--) { // cycle columns
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                for(auto & kv : U[br][bc]) {
                    auto ila = kv.first;
                    auto la = las[ila];
                    int kmax = kv.second.size();
                    if(RU.find(la)==RU.end()) {
                        int kmax = IK[nbs-1][ila];
                        RU[la].resize(kmax);
                        for(int k=0; k<kmax; k++) {
                            RU[la][k].resize(xn+1);
                            for(int n=0; n<=xn; n++) RU[la][k][n] = matrix(N,N);
                        }
                    }
                    for(int k=0; k<kmax; k++) {
                        for(int n=0; n<=xn; n++) {
                            auto mat = _to_(U[br][bc][ila][k][n],fp);
                            acb_mat_clear(U[br][bc][ila][k][n]);
                            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RU[la][k][n](r0+r,c0+c) = mat(r,c);
                        }
                    }
                }
            }
        }
        return RU;
    }
    
    matrix DEX::series(int xn, slong dp, const ex & x0, const lst & slas) { // RU[a][b][la][k][n]
        if(!fuchsified) fuchsify();
        vector<fmpq_t> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        block_umat_acb_mat_t U;
        series(U,xn,dp,qslas);
        for(int i=0; i<qslas.size(); i++) fmpq_clear(qslas[i]);
        
        int nbs = bs.size();
        auto fp = dp2fp(dp);
        acb_t z0, lz0, z;
        acb_init(z0);
        acb_init(lz0);
        acb_init(z);
        _to_(z0,x0,fp);
        acb_log(lz0,z0,fp);
        vector<vector<acb_t>> zln(nlas); // cache z0^(la+n) zln[ila][n]
        for(int ila=0; ila<nlas; ila++) {
            zln[ila] = vector<acb_t>(xn+1);
            arb_t r;
            arb_init(r);
            auto qla = qlas[ila][0];
            arb_set_fmpq(r,qla,fp);
            acb_pow_arb(z,z0,r,fp); // z^la
            for(int n=0; n<=xn; n++) {
                acb_init(zln[ila][n]);
                acb_set(zln[ila][n],z);
                acb_mul(z,z,z0,fp); // next z^(la+n)
            }
            arb_clear(r);
        }
        vector<acb_t> lzk(kmmax); // cache ln^k z/k!
        acb_one(z);
        for(int k=0; k<kmmax; k++) {
            if(k>1) acb_div_si(z,z,k,fp); // div k!
            acb_init(lzk[k]);
            acb_set(lzk[k],z);
            acb_mul(z,z,lz0,fp); // next ln^k z
        }

        acb_mat_t cmat;
        acb_mat_init(cmat,N,N);
        acb_mat_zero(cmat);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int bc=0; bc<=br; bc++) { // cycle columns
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                for(auto & kv : U[br][bc]) {
                    auto ila = kv.first;
                    int kmax = kv.second.size();
                    for(int k=0; k<kmax; k++) {
                        auto nmax = kv.second[k].size();
                        for(int n=0; n<nmax; n++) {
                            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                                auto item = acb_mat_entry(cmat,r0+r,c0+c);
                                acb_t tz;
                                acb_init(tz);
                                acb_mul(tz,acb_mat_entry(kv.second[k][n],r,c),zln[ila][n],fp);
                                acb_mul(tz,tz,lzk[k],fp);
                                acb_add(item,item,tz,fp); 
                                acb_clear(tz);
                            }
                            acb_mat_clear(U[br][bc][ila][k][n]);
                        }
                    }
                }
            }
        }
        
        acb_clear(z0);
        acb_clear(lz0);
        acb_clear(z);
        for(int ila=0; ila<nlas; ila++) for(int n=0; n<=xn; n++) acb_clear(zln[ila][n]);
        for(int k=0; k<kmmax; k++) acb_clear(lzk[k]);
        auto mat =  _to_(cmat,fp);
        acb_mat_clear(cmat);
        return mat;
    }
    
    
    //=*********************************************************************=
    // I-Series - expansion - rational
    //=*********************************************************************=
    
    imat_t DEX::series(int xn, const matrix & m, const lst & slas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        int nc = m.cols();
        
        vector<fmpq_mat_t> qmat0(nbs);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            fmpq_mat_init(qmat0[br],nr,nc);
            _to_(qmat0[br],ex_to<matrix>(sub_matrix(m, r0, nr, 0, nc)));
        }
        
        block_imat_fmpq_mat_t In0(nbs);
        for(int br=0; br<nbs; br++) {
            int a = br;
            int nr = bs[br].second;
            fmpq_mat_t qmat;
            fmpq_mat_init(qmat,nr,nc);
            for(int bc=0; bc<=br; bc++) {
                int b = bc;
                for(auto & kv : U0[a][b]) {
                    auto ila = kv.first;
                    int kmax = kv.second.size();
                    if(In0[a].find(ila)==In0[a].end()) {
                        auto kmax2 = IK[a][ila];
                        In0[a][ila].resize(kmax2);
                        for(int k=0; k<kmax2; k++) {
                            In0[a][ila][k] = vector<fmpq_mat_t>(1); 
                            fmpq_mat_init(In0[a][ila][k][0],nr,nc);
                        }
                    }
                    for(int k=0; k<kmax; k++) {
                        fmpq_mat_mul(qmat,kv.second[k][0],qmat0[bc]);
                        fmpq_mat_add(In0[a][ila][k][0],In0[a][ila][k][0],qmat);
                    }
                }
            }
            fmpq_mat_clear(qmat);
        }
        for(int br=0; br<nbs; br++) fmpq_mat_clear(qmat0[br]);
        
        vector<fmpq_t> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        
        block_imat_fmpq_mat_t I;
        series(I,xn,In0,nc,qslas);
        for(int a=0; a<nbs; a++) for(auto & kv : In0[a]) for(auto & item : kv.second) fmpq_mat_clear(item[0]);
        for(int i=0; i<qslas.size(); i++) fmpq_clear(qslas[i]);
        
        imat_t RI;
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(auto & kv : I[br]) {
                auto ila = kv.first;
                auto la = las[ila];
                int kmax = kv.second.size();
                if(RI.find(la)==RI.end()) {
                    int kmax = IK[nbs-1][ila];
                    RI[la].resize(kmax);
                    for(int k=0; k<kmax; k++) {
                        RI[la][k].resize(xn+1);
                        for(int n=0; n<=xn; n++) RI[la][k][n] = matrix(N,nc);
                    }
                }
                for(int k=0; k<kmax; k++) {
                    for(int n=0; n<=xn; n++) {
                        auto mat = _to_(I[br][ila][k][n]);
                        fmpq_mat_clear(I[br][ila][k][n]);
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RI[la][k][n](r0+r,c) = mat(r,c);
                    }
                }
            }
        }
        return RI;
        
    }
    
    
    //=*********************************************************************=
    // I-Series - acb
    //=*********************************************************************=
    
    imat_t DEX::series(int xn, const matrix & m, slong dp, const lst & slas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        int nc = m.cols();
        vector<acb_mat_t> qmat0(nbs);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            acb_mat_init(qmat0[br],nr,nc);
            _to_(qmat0[br],ex_to<matrix>(sub_matrix(m, r0, nr, 0, nc)),fp);
        }
        
        block_imat_acb_mat_t In0(nbs);
        for(int br=0; br<nbs; br++) {
            int a = br;
            int nr = bs[br].second;
            acb_mat_t qmat;
            acb_mat_init(qmat,nr,nc);
            for(int bc=0; bc<=br; bc++) {
                int nr2 = bs[bc].second;
                acb_mat_t umat;
                acb_mat_init(umat,nr,nr2);
                int b = bc;
                for(auto & kv : U0[a][b]) {
                    auto ila = kv.first;
                    int kmax = kv.second.size();
                    if(In0[a].find(ila)==In0[a].end()) {
                        auto kmax2 = IK[a][ila];
                        In0[a][ila].resize(kmax2);
                        for(int k=0; k<kmax2; k++) {
                            In0[a][ila][k] = vector<acb_mat_t>(1);
                            acb_mat_init(In0[a][ila][k][0],nr,nc);
                        }
                    }
                    for(int k=0; k<kmax; k++) {
                        acb_mat_set_fmpq_mat(umat,kv.second[k][0],fp);
                        acb_mat_mul(qmat,umat,qmat0[bc],fp);
                        acb_mat_add(In0[a][ila][k][0],In0[a][ila][k][0],qmat,fp);
                    }
                }
                acb_mat_clear(umat);
            }
            acb_mat_clear(qmat);
        }
        for(int br=0; br<nbs; br++) acb_mat_clear(qmat0[br]);
        
        vector<fmpq_t> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        
        block_imat_acb_mat_t I;
        series(I,xn,In0,nc,dp,qslas);
        for(int a=0; a<nbs; a++) for(auto & kv : In0[a]) for(auto & item : kv.second) acb_mat_clear(item[0]);
        for(int i=0; i<qslas.size(); i++) fmpq_clear(qslas[i]);
        
        mag_t mag;
        mag_init(mag);
        imat_t RI;
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(auto & kv : I[br]) {
                auto ila = kv.first;
                auto la = las[ila];
                int kmax = kv.second.size();
                if(RI.find(la)==RI.end()) {
                    int kmax = IK[nbs-1][ila];
                    RI[la].resize(kmax);
                    for(int k=0; k<kmax; k++) {
                        RI[la][k].resize(xn+1);
                        for(int n=0; n<=xn; n++) RI[la][k][n] = matrix(N,nc);
                    }
                }
                for(int k=0; k<kmax; k++) {
                    for(int n=0; n<=xn; n++) {
                        // - error check
                        int nr = acb_mat_nrows(I[br][ila][k][n]);
                        int nc = acb_mat_ncols(I[br][ila][k][n]);
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                            auto item = acb_mat_entry(I[br][ila][k][n],r,c);
                            auto ri = acb_realref(item);
                            if(arb_rel_error_bits(ri)>-rel_fp) {
                                arb_get_mag(mag,ri);
                                if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                    cout << endl; arb_printd(ri,5); cout << endl;
                                }
                            }
                            ri = acb_imagref(item);
                            if(arb_rel_error_bits(ri)>-rel_fp) {
                                if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                    cout << endl; arb_printd(ri,5); cout << endl;  
                                }
                            }
                        }
                        // - error check end
                            
                        auto mat = _to_(I[br][ila][k][n],fp);
                        acb_mat_clear(I[br][ila][k][n]);
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RI[la][k][n](r0+r,c) = mat(r,c);
                    }
                }
            }
        }
        mag_clear(mag);
        return RI;
        
    }
    
    matrix DEX::series(int xn, const matrix & m, slong dp, const ex & x0, const lst & slas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        int nc = m.cols();
        vector<acb_mat_t> qmat0(nbs);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            acb_mat_init(qmat0[br],nr,nc);
            _to_(qmat0[br],ex_to<matrix>(sub_matrix(m, r0, nr, 0, nc)),fp);
        }
        
        block_imat_acb_mat_t In0(nbs);
        for(int br=0; br<nbs; br++) {
            int a = br;
            int nr = bs[br].second;
            acb_mat_t qmat;
            acb_mat_init(qmat,nr,nc);
            for(int bc=0; bc<=br; bc++) {
                int nr2 = bs[bc].second;
                acb_mat_t umat;
                acb_mat_init(umat,nr,nr2);
                int b = bc;
                for(auto & kv : U0[a][b]) {
                    auto ila = kv.first;
                    int kmax = kv.second.size();
                    if(In0[a].find(ila)==In0[a].end()) {
                        auto kmax2 = IK[a][ila];
                        In0[a][ila].resize(kmax2);
                        for(int k=0; k<kmax2; k++) {
                            In0[a][ila][k] = vector<acb_mat_t>(1);
                            acb_mat_init(In0[a][ila][k][0],nr,nc);
                        }
                    }
                    for(int k=0; k<kmax; k++) {
                        acb_mat_set_fmpq_mat(umat,kv.second[k][0],fp);
                        acb_mat_mul(qmat,umat,qmat0[bc],fp);
                        acb_mat_add(In0[a][ila][k][0],In0[a][ila][k][0],qmat,fp);
                    }
                }
                acb_mat_clear(umat);
            }
            acb_mat_clear(qmat);
        }
        for(int br=0; br<nbs; br++) acb_mat_clear(qmat0[br]);
        
        vector<fmpq_t> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        
        block_imat_acb_mat_t I;
        series(I,xn,In0,nc,dp,qslas);
        for(int a=0; a<nbs; a++) for(auto & kv : In0[a]) for(auto & item : kv.second) acb_mat_clear(item[0]);
        for(int i=0; i<qslas.size(); i++) fmpq_clear(qslas[i]);
        
        acb_t z0, lz0, z;
        acb_init(z0);
        acb_init(lz0);
        acb_init(z);
        _to_(z0,x0,fp);
        acb_log(lz0,z0,fp);
        vector<vector<acb_t>> zln(nlas); // cache z0^(la+n) zln[ila][n]
        for(int ila=0; ila<nlas; ila++) {
            zln[ila] = vector<acb_t>(xn+1);
            arb_t r;
            arb_init(r);
            auto qla = qlas[ila][0];
            arb_set_fmpq(r,qla,fp);
            acb_pow_arb(z,z0,r,fp); // z^la
            for(int n=0; n<=xn; n++) {
                acb_init(zln[ila][n]);
                acb_set(zln[ila][n],z);
                acb_mul(z,z,z0,fp); // next z^(la+n)
            }
            arb_clear(r);
        }
        vector<acb_t> lzk(kmmax); // cache ln^k z/k!
        acb_one(z);
        for(int k=0; k<kmmax; k++) {
            if(k>1) acb_div_si(z,z,k,fp); // div k!
            acb_init(lzk[k]);
            acb_set(lzk[k],z);
            acb_mul(z,z,lz0,fp); // next ln^k z
        }
        
        mag_t mag;
        mag_init(mag);
        acb_mat_t imat;
        acb_mat_init(imat,N,nc);
        acb_mat_zero(imat);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(auto & kv : I[br]) {
                auto ila = kv.first;
                int kmax = kv.second.size();
                for(int k=0; k<kmax; k++) {
                    for(int n=0; n<=xn; n++) {
                        // - error check
                        int nr = acb_mat_nrows(I[br][ila][k][n]);
                        int nc = acb_mat_ncols(I[br][ila][k][n]);
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                            auto item = acb_mat_entry(I[br][ila][k][n],r,c);
                            auto ri = acb_realref(item);
                            if(arb_rel_error_bits(ri)>-rel_fp) {
                                arb_get_mag(mag,ri);
                                if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                    cout << endl; arb_printd(ri,5); cout << endl;
                                }
                            }
                            ri = acb_imagref(item);
                            if(arb_rel_error_bits(ri)>-rel_fp) {
                                if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                    cout << endl; arb_printd(ri,5); cout << endl;  
                                }
                            }
                        }
                        // - error check end
                            
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                            auto item = acb_mat_entry(imat,r0+r,c);
                            acb_t tz;
                            acb_init(tz);
                            acb_mul(tz,acb_mat_entry(kv.second[k][n],r,c),zln[ila][n],fp);
                            acb_mul(tz,tz,lzk[k],fp);
                            acb_add(item,item,tz,fp); 
                            acb_clear(tz);
                        }
                        acb_mat_clear(I[br][ila][k][n]);
                    }
                }
            }
        }
        mag_clear(mag);
        acb_clear(z0);
        acb_clear(lz0);
        acb_clear(z);
        for(int ila=0; ila<nlas; ila++) for(int n=0; n<=xn; n++) acb_clear(zln[ila][n]);
        for(int k=0; k<kmmax; k++) acb_clear(lzk[k]);
        auto mat =  _to_(imat,fp);
        acb_mat_clear(imat);
        return mat;
    }
    
    //=*********************************************************************=
    // Taylor expansion
    //=*********************************************************************=
    
    vector<matrix> DEX::taylor(int xn, const matrix I0, const ex & x0) {  
        auto nbs = bs.size();
        int nc = I0.cols();
 
        vector<vector<fmpq_mat_t>> I; // I[a][n]
        taylor(I,xn,I0,x0);
        
        vector<matrix> RI(xn+1);
        for(int n=0; n<=xn; n++) RI[n] = matrix(N,nc);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int n=0; n<=xn; n++) {
                auto mat = _to_(I[br][n]);
                fmpq_mat_clear(I[br][n]);
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RI[n](r0+r,c) = mat(r,c);
            }
        }
        return RI;
    }
    
    vector<matrix> DEX::taylor(int xn, const matrix I0, const ex & x0, slong dp) {  
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        int nc = I0.cols();
        
        vector<vector<acb_mat_t>> I; // I[a][n]
        taylor(I,xn,I0,x0,dp);
        
        mag_t mag;
        mag_init(mag);
        vector<matrix> RI(xn+1);
        for(int n=0; n<=xn; n++) RI[n] = matrix(N,nc);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int n=0; n<=xn; n++) {
                
                // - error check
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    auto item = acb_mat_entry(I[br][n],r,c);
                    auto ri = acb_realref(item);
                    if(arb_rel_error_bits(ri)>-rel_fp) {
                        arb_get_mag(mag,ri);
                        if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                            cout << endl; arb_printd(ri,5); cout << endl;
                        }
                    }
                    ri = acb_imagref(item);
                    if(arb_rel_error_bits(ri)>-rel_fp) {
                        if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                            cout << endl; arb_printd(ri,5); cout << endl;  
                        }
                    }
                }
                // - error check end
                    
                auto mat = _to_(I[br][n],fp);
                acb_mat_clear(I[br][n]);
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RI[n](r0+r,c) = mat(r,c);
            }
        }
        mag_clear(mag);
        
        return RI;
    }
    
    matrix DEX::taylor(int xn, const matrix I0, const ex & x0, slong dp, const ex & dx) {  
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        int nc = I0.cols();
        
        vector<vector<acb_mat_t>> I; // I[a][n]
        taylor(I,xn,I0,x0,dp);
        
        acb_t z0, z;
        acb_init(z0);
        acb_init(z);
        _to_(z0,dx,fp);
        vector<acb_t> zn(xn+1);
        acb_one(z);
        for(int n=0; n<=xn; n++) {
            acb_init(zn[n]);
            acb_set(zn[n],z);
            acb_mul(z,z,z0,fp); // next z^n
        }
        
        mag_t mag;
        mag_init(mag);
        acb_mat_t imat;
        acb_mat_init(imat,N,nc);
        acb_mat_zero(imat);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int n=0; n<=xn; n++) {
                
                // - error check
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    auto item = acb_mat_entry(I[br][n],r,c);
                    auto ri = acb_realref(item);
                    if(arb_rel_error_bits(ri)>-rel_fp) {
                        arb_get_mag(mag,ri);
                        if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                            cout << endl; arb_printd(ri,5); cout << endl;
                        }
                    }
                    ri = acb_imagref(item);
                    if(arb_rel_error_bits(ri)>-rel_fp) {
                        if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                            cout << endl; arb_printd(ri,5); cout << endl;  
                        }
                    }
                }
                // - error check end
                    
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    auto item = acb_mat_entry(imat,r0+r,c);
                    acb_t tz;
                    acb_init(tz);
                    acb_mul(tz,acb_mat_entry(I[br][n],r,c),zn[n],fp);
                    acb_add(item,item,tz,fp); 
                    acb_clear(tz);
                }
                acb_mat_clear(I[br][n]);
            }
        }
        mag_clear(mag);
        acb_clear(z0);
        acb_clear(z);
        for(int n=0; n<=xn; n++) acb_clear(zn[n]);
        auto mat =  _to_(imat,fp);
        acb_mat_clear(imat);
        return mat;
    }
    
    //=*********************************************************************=
    
    matrix u2m(const umat_t & umat, const ex & x) {
        matrix mat;
        bool first = true;
        for(auto kv : umat) {
            auto item = kv.second;
            int kmax = item.size();
            ex lxk = 1;
            for(int k=0; k<kmax; k++) {
                int nmax = item[k].size();
                ex xln = pow(x,kv.first);
                for(int n=0; n<nmax; n++) {
                    auto m = item[k][n].mul_scalar(xln*lxk);
                    if(first) { first=false; mat = m; }
                    else mat = mat.add(m);
                    xln *= x;
                }
                lxk *= log(x)/(k+1);
            }
        }
        return mat;
    }
    
    matrix i2m(const imat_t & imat, const ex & x) { return u2m(imat,x); }
    
    matrix i2m(const vector<matrix> & tmat, const ex & x) {
        matrix mat;
        bool first = true;
        int nmax = tmat.size();
        ex xn = 1;
        for(int n=0; n<nmax; n++) {
            auto m = tmat[n].mul_scalar(xn);
            if(first) { first=false; mat = m; }
            else mat = mat.add(m);
            xn *= x;
        }
        return mat;
    }
    
    //=*********************************************************************=
        
}

