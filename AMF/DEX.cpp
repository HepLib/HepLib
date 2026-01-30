/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "AMF.h"
#include "vspace.h"

namespace HepLib {

    //=*********************************************************************=
    
    DEX::DEX(const symbol & _x, const string & _pre) : x(_x), pre(_pre) { }
    DEX::~DEX() { clear(); }
    void DEX::clear() {
        int nbs = bs.size();
        if(nbs<1) return;
        for(auto & i1 : Mat) for(auto & i2 : i1) for(auto & i3 : i2) for(auto & i4 : i3) {
            fmpz_poly_q_clear(i4); // fmpz_poly_q_init() in initialize
            flint_free(i4);
        }
        Mat.clear();
        for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
            if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) {
                fmpq_mat_clear(item[0]); // fmpq_mat_init() in fuchsify
                flint_free(item[0]);
            }
            if(taylor_inited) for(auto & row : QMat[br][bc]) for(auto & item : row) {
                fmpz_poly_clear(item); // fmpz_poly_init() in taylor
                flint_free(item);
            }
            if(gr_taylor_inited) for(auto & row : GMat[br][bc]) for(auto & item : row) {
                gr_poly_clear(item, _ctx_); // gr_poly_init() in gr-taylor
                flint_free(item);
            }
        }
        bs.clear();
        if(taylor_inited) for(int br=0; br<nbs; br++) {
            fmpz_poly_clear(QD[br]); // fmpz_poly_init() in taylor
            flint_free(QD[br]);
        }
        if(gr_taylor_inited) for(int br=0; br<nbs; br++) {
            gr_poly_clear(GD[br], _ctx_); // gr_poly_init() in gr-taylor
            flint_free(GD[br]);
        }
        if(fuchsified) {
            U0.clear(); 
            for(auto & item : qlas) {
                fmpq_clear(item); // fmpq_init() in fuchsify
                flint_free(item);
            }
            qlas.clear();
        }
        if(taylor_inited) { QMat.clear(); QD.clear(); }
        if(gr_taylor_inited) { GMat.clear(); GD.clear(); gr_ctx_clear(_ctx_); } // gr_ctx_init() in gr-taylor
        fuchsified = false;
        taylor_inited = false;
        gr_taylor_inited = false;
    }
    
    matrix DEX::T() {
        int n = Ts.size();
        int n2 = n/2;
        vector<matrix> Ts2(n2);
        for(int i=0; i<n2; i++) Ts2[i] = Ts[2*i].mul(Ts[2*i+1]);
        matrix t = ex_to<matrix>(unit_matrix(N));
        if(n%2==1) t = Ts[n-1];
        for(int i=n2-1; i>=0; i--) t = Ts2[i].mul(t);
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
        vector<fmpq*> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            qslas[i] = (fmpq*) flint_malloc(sizeof(fmpq));
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        int nbs = bs.size();
        
        abikn_fmpq_mat_t U;
        //if(nbs<50) all_series(U,xn,qslas);
        //else ab_series(U,xn,qslas);
        ab_series(U,xn,qslas,1); // non-parallel

        for(int i=0; i<qslas.size(); i++) {
            fmpq_clear(qslas[i]);
            flint_free(qslas[i]);
        }

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
                            fmpq_mat_clear(U[br][bc][ila][k][n]); // fmpq_mat_init() in _series
                            flint_free(U[br][bc][ila][k][n]);
                            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RU[la][k][n](r0+r,c0+c) = mat(r,c);
                        }
                    }
                }
            }
        }
        return RU;
    }
        
    matrix DEX::series(int xn, gr_ptr z0, gr_ctx_t ctx, const lst & slas) { // x -> z0
        if(!fuchsified) fuchsify();
        vector<fmpq*> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            qslas[i] = (fmpq*) flint_malloc(sizeof(fmpq));
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        int nbs = bs.size();
        
        abikn_gr_mat_t U;
        // TODO: add all_series version
        //ab_series(U, xn, ctx, qslas);
        ab_series(U, xn, ctx, qslas,1); // non-parallel
        for(int i=0; i<qslas.size(); i++) {
            fmpq_clear(qslas[i]);
            flint_free(qslas[i]);
        }
        
        int status = GR_SUCCESS;
        gr_ptr lz0 = gr_heap_init(ctx);
        gr_ptr z = gr_heap_init(ctx);
        status |= gr_log(lz0, z0, ctx);
        vector<vector<gr_ptr>> zln(nlas); // cache z0^(la+n) zln[ila][n]
        for(int ila=0; ila<nlas; ila++) {
            zln[ila].resize(xn+1);
            status |= gr_pow_fmpq(z, z0, qlas[ila], ctx); // z^la
            for(int n=0; n<=xn; n++) {
                zln[ila][n] = gr_heap_init(ctx);
                status |= gr_set(zln[ila][n], z, ctx);
                status |= gr_mul(z, z, z0, ctx); // next z^(la+n)
            }
        }
        vector<gr_ptr> lzk(kmmax); // cache ln^k z/k!
        status |= gr_one(z, ctx);
        for(int k=0; k<kmmax; k++) {
            if(k>1) status |= gr_div_si(z, z, k, ctx); // div k!
            lzk[k] = gr_heap_init(ctx);
            status |= gr_set(lzk[k], z, ctx);
            status |= gr_mul(z, z, lz0, ctx); // next ln^k z
        }

        acb_mat_t cmat;
        acb_mat_init(cmat, N, N);
        acb_mat_zero(cmat);
        acb_t zb;
        acb_init(zb);
        slong fp;
        status |= gr_ctx_get_real_prec(&fp, ctx);
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
                                auto item = acb_mat_entry(cmat, r0+r, c0+c);
                                status |= gr_mul(z, gr_mat_entry_ptr(kv.second[k][n],r,c,ctx), zln[ila][n], ctx);
                                status |= gr_mul(z, z, lzk[k], ctx);
                                _to_(zb, (acf_struct*)z);
                                acb_add(item, item, zb, fp);
                            }
                            gr_mat_clear(U[br][bc][ila][k][n], ctx); // call gr_mat_init() in _series
                            flint_free(U[br][bc][ila][k][n]);
                        }
                    }
                }
            }
        }
        acb_clear(zb);
        gr_heap_clear(lz0, ctx);
        gr_heap_clear(z, ctx);
        for(int ila=0; ila<nlas; ila++) for(int n=0; n<=xn; n++) gr_heap_clear(zln[ila][n], ctx);
        for(int k=0; k<kmmax; k++) gr_heap_clear(lzk[k], ctx);
        auto mat = _to_(cmat, fp);
        acb_mat_clear(cmat);
        if(status != GR_SUCCESS) throw Error("u-series: status != GR_SUCCESS");
        return mat;
    }
    
    
    //=*********************************************************************=
    // I-Series - expansion - rational
    //=*********************************************************************=
    
    imat_t DEX::series(int xn, const matrix & m, const lst & slas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        int nc = m.cols();
        
        vector<fmpq_mat_struct*> qmat0(nbs);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            qmat0[br] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
            fmpq_mat_init(qmat0[br],nr,nc);
            _to_(qmat0[br],ex_to<matrix>(sub_matrix(m, r0, nr, 0, nc)));
        }
        
        aikn_fmpq_mat_t In0(nbs);
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
                            In0[a][ila][k].resize(1);
                            In0[a][ila][k][0] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
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
        for(int br=0; br<nbs; br++) {
            fmpq_mat_clear(qmat0[br]);
            flint_free(qmat0[br]);
        }
        
        vector<fmpq*> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            qslas[i] = (fmpq*) flint_malloc(sizeof(fmpq));
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
        
        aikn_fmpq_mat_t I;
        series(I,xn,In0,nc,qslas);
        //pseries(I,xn,In0,nc,qslas); // TODO: to be added
        for(int a=0; a<nbs; a++) for(auto & kv : In0[a]) for(auto & item : kv.second) {
            fmpq_mat_clear(item[0]);
            flint_free(item[0]);
        }
        for(int i=0; i<qslas.size(); i++) {
            fmpq_clear(qslas[i]);
            flint_free(qslas[i]);
        }
        
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
                        fmpq_mat_clear(I[br][ila][k][n]); // call fmpq_mat_init() in _series
                        flint_free(I[br][ila][k][n]);
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RI[la][k][n](r0+r,c) = mat(r,c);
                    }
                }
            }
        }
        return RI;
        
    }
    
    
    //=*********************************************************************=
    // I-Series
    //=*********************************************************************=
    
    matrix DEX::series(int xn, const matrix & m, gr_ptr z0, gr_ctx_t ctx, const lst & slas) { // x -> z0
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        int nc = m.cols();
        vector<gr_mat_struct*> vmat(nbs);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            vmat[br] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
            gr_mat_init(vmat[br], nr, nc, ctx);
            _to_(vmat[br], ex_to<matrix>(sub_matrix(m, r0, nr, 0, nc)), ctx);
        }
        
        int status = GR_SUCCESS;
        aikn_gr_mat_t In0(nbs);
        for(int br=0; br<nbs; br++) {
            int a = br;
            int nr = bs[br].second;
            gr_mat_t qmat;
            gr_mat_init(qmat, nr, nc, ctx);
            for(int bc=0; bc<=br; bc++) {
                int nc2 = bs[bc].second;
                gr_mat_t umat;
                gr_mat_init(umat, nr, nc2, ctx);
                int b = bc;
                for(auto & kv : U0[a][b]) {
                    auto ila = kv.first;
                    int kmax = kv.second.size();
                    if(In0[a].find(ila)==In0[a].end()) {
                        auto kmax2 = IK[a][ila];
                        In0[a][ila].resize(kmax2);
                        for(int k=0; k<kmax2; k++) {
                            In0[a][ila][k].resize(1);
                            In0[a][ila][k][0] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                            gr_mat_init(In0[a][ila][k][0], nr, nc, ctx);
                            status |= gr_mat_zero(In0[a][ila][k][0], ctx);
                        }
                    }
                    for(int k=0; k<kmax; k++) {
                        status |= gr_mat_set_fmpq_mat(umat, kv.second[k][0], ctx);
                        status |= gr_mat_mul(qmat, umat, vmat[bc], ctx);
                        status |= gr_mat_add(In0[a][ila][k][0], In0[a][ila][k][0], qmat, ctx);
                    }
                }
                gr_mat_clear(umat, ctx);
            }
            gr_mat_clear(qmat, ctx);
        }
        for(int br=0; br<nbs; br++) {
            gr_mat_clear(vmat[br], ctx);
            flint_free(vmat[br]);
        }
        
        vector<fmpq*> qslas(slas.nops());
        for(int i=0; i<qslas.size(); i++) {
            qslas[i] = (fmpq*) flint_malloc(sizeof(fmpq));
            fmpq_init(qslas[i]);
            _to_(qslas[i], slas.op(i));
        }
                
        aikn_gr_mat_t I;
        //if(nbs<50) all_series(I, xn, In0, nc, ctx, qslas);
        //else a_series(I, xn, In0, nc, ctx, qslas);
        a_series(I, xn, In0, nc, ctx, qslas, 1); // non-parallel
        
        for(int a=0; a<nbs; a++) for(auto & kv : In0[a]) for(auto & item : kv.second) {
            gr_mat_clear(item[0], ctx);
            flint_free(item[0]);
        }
        for(int i=0; i<qslas.size(); i++) {
            fmpq_clear(qslas[i]);
            flint_free(qslas[i]);
        }
        
        gr_ptr lz0 = gr_heap_init(ctx);
        gr_ptr z = gr_heap_init(ctx);
        status |= gr_log(lz0, z0, ctx);
        vector<vector<gr_ptr>> zln(nlas); // cache z0^(la+n) zln[ila][n]
        for(int ila=0; ila<nlas; ila++) {
            zln[ila].resize(xn+1);
            status |= gr_pow_fmpq(z, z0, qlas[ila], ctx); // z^la
            for(int n=0; n<=xn; n++) {
                zln[ila][n] = gr_heap_init(ctx);
                status |= gr_set(zln[ila][n], z, ctx);
                status |= gr_mul(z, z, z0, ctx); // next z^(la+n)
            }
        }
        vector<gr_ptr> lzk(kmmax); // cache ln^k z/k!
        status |= gr_one(z, ctx);
        for(int k=0; k<kmmax; k++) {
            if(k>1) status |= gr_div_si(z, z, k, ctx); // div k!
            lzk[k] = gr_heap_init(ctx);
            status |= gr_set(lzk[k], z, ctx);
            status |= gr_mul(z, z, lz0, ctx); // next ln^k z
        }

        acb_mat_t imat;
        acb_mat_init(imat,N,nc);
        acb_mat_zero(imat);
        acb_t zb;
        acb_init(zb);
        slong fp;
        status |= gr_ctx_get_real_prec(&fp, ctx);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(auto & kv : I[br]) {
                auto ila = kv.first;
                int kmax = kv.second.size();
                for(int k=0; k<kmax; k++) {
                    for(int n=0; n<=xn; n++) {
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                            auto item = acb_mat_entry(imat,r0+r,c);
                            status |= gr_mul(z, gr_mat_entry_ptr(kv.second[k][n],r,c,ctx), zln[ila][n], ctx);
                            status |= gr_mul(z, z, lzk[k], ctx);
                            _to_(zb, (acf_struct*)z);
                            acb_add(item, item, zb, fp);
                        }
                        gr_mat_clear(I[br][ila][k][n], ctx); // call gr_mat_init() in _series
                        flint_free(I[br][ila][k][n]);
                    }
                }
            }
        }
        acb_clear(zb);
        gr_heap_clear(lz0, ctx);
        gr_heap_clear(z, ctx);
        for(int ila=0; ila<nlas; ila++) for(int n=0; n<=xn; n++) gr_heap_clear(zln[ila][n], ctx);
        for(int k=0; k<kmmax; k++) gr_heap_clear(lzk[k], ctx);
        auto mat =  _to_(imat,fp);
        acb_mat_clear(imat);
        if(status != GR_SUCCESS) throw Error("i-series: status != GR_SUCCESS");
        return mat;
    }
    
    //=*********************************************************************=
    // Taylor expansion
    //=*********************************************************************=
    
    vector<matrix> DEX::taylor(int xn, const matrix I0, const ex & x0) {  
        auto nbs = bs.size();
        int nc = I0.cols();
 
        vector<vector<fmpq_mat_struct*>> I; // I[a][n]
        taylor(I,xn,I0,x0);
        
        vector<matrix> RI(xn+1);
        for(int n=0; n<=xn; n++) RI[n] = matrix(N,nc);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int n=0; n<=xn; n++) {
                auto mat = _to_(I[br][n]);
                fmpq_mat_clear(I[br][n]); // call fmpq_mat_init() in taylor
                flint_free(I[br][n]);
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) RI[n](r0+r,c) = mat(r,c);
            }
        }
        return RI;
    }
    
    // imat as input and output
    void DEX::taylor(int xn, gr_mat_t imat, gr_ptr z0, gr_ptr dz, gr_ctx_t ctx, const string & es) {
        auto nbs = bs.size();
        int nc = gr_mat_ncols(imat, ctx);
        
        vector<vector<gr_mat_struct*>> I; // I[a][n]
        //if(nbs<50) an_taylor(I, xn, imat, z0, ctx, es);
        //else a_taylor(I, xn, imat, z0, ctx, es);
        a_taylor(I, xn, imat, z0, ctx, es, 1); // non - parallel
        
        vector<gr_ptr> zn(xn+1);
        int status = GR_SUCCESS;
        gr_ptr z = gr_heap_init(ctx);
        status |= gr_one(z, ctx);
        for(int n=0; n<=xn; n++) {
            zn[n] = gr_heap_init(ctx);
            status |= gr_set(zn[n], z, ctx);
            status |= gr_mul(z, z, dz, ctx); // next dz^n
        }
        
        status |= gr_mat_zero(imat, ctx);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int n=0; n<=xn; n++) {
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    auto item = gr_mat_entry_ptr(imat, r0+r, c, ctx);
                    status |= gr_mul(z, gr_mat_entry_ptr(I[br][n], r, c, ctx), zn[n], ctx);
                    status |= gr_add(item, item, z, ctx);
                }
                gr_mat_clear(I[br][n], ctx); // call gr_mat_init() in _taylor
                flint_free(I[br][n]);
            }
        }
        gr_heap_clear(z, ctx);
        for(int n=0; n<=xn; n++) gr_heap_clear(zn[n], ctx);
        if(status != GR_SUCCESS) throw Error("taylor: status != GR_SUCCESS");
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

