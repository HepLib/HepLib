/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"
#include "vspace.h"

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
    }

    //=*********************************************************************=
    
    DEX::DEX(const symbol & _x) : x(_x) { }
    DEX::~DEX() { clear(); }
    void DEX::clear() {
        int nbs = bs.size();
        if(nbs<1) return;
        for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
            for(auto & row : Mat[br][bc]) for(auto & item : row) fmpz_poly_q_clear(item);
            if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) fmpq_mat_clear(item[0]);
            if(ntaylor_inited) for(auto & row : TMat[br][bc]) for(auto & item : row) acb_poly_clear(item);
        }
        bs.clear();
        Mat.clear();
        if(ntaylor_inited) for(int br=0; br<nbs; br++) acb_poly_clear(TD[br][0]);
        if(fuchsified) { 
            U0.clear(); 
            for(auto & item : qlas) fmpq_clear(item[0]);
            qlas.clear();
        }
        if(ntaylor_inited) { TMat.clear(); TD.clear(); }
        fuchsified = false;
        ntaylor_inited = false;
    }
    
    // init to lower triangle block matrix
    void DEX::init(const matrix & m) {
        auto n = m.rows();
        bool_mat_t bm;
        bool_mat_init(bm,n,n);
        slong vs[n];
        for(int r=0; r<n; r++) for(int c=0; c<n; c++) bool_mat_set_entry(bm,r,c,!is_zero(m(r,c)));
        auto nb = bool_mat_get_strongly_connected_components(vs,bm);
        bool_mat_clear(bm);
        set<slong> sets[nb];
        bs.clear();
        for(int i=0; i<n; i++) sets[vs[i]].insert(i);
        int ci = 0;
        for(auto s : sets) {
            bs.push_back(make_pair(ci,s.size()));
            for(auto si : s) {
                vs[ci] = si;
                ci++;
            }
        }
        matrix pmat(n,n);
        for(int r=0; r<n; r++) pmat(vs[r],r) = 1;
        auto mat = pmat.inverse().mul(m).mul(pmat);
        Ts.push_back(pmat);
        N = n;
        // init Mat[br][bc][r][c]
        int nbs = bs.size();
        Mat.resize(nbs);
        for(int br=0; br<nbs; br++) {
            Mat[br].resize(br+1);
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int bc=0; bc<=br; bc++) {
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                Mat[br][bc].resize(nr);
                for(int r=0; r<nr; r++) {
                    Mat[br][bc][r] = vector<fmpz_poly_q_t>(nc);
                    for(int c=0; c<nc; c++) {
                        fmpz_poly_q_init(Mat[br][bc][r][c]);
                        _to_(Mat[br][bc][r][c], mat(r0+r,c0+c));
                    }
                }
            }
        }
        fuchsified = false;
    }
    
    vector<vector<matrix>> DEX::m2b(const matrix & m) { // m is low triangular block w.r.t. bs
        int nbs = bs.size();
        vector<vector<matrix>> bm(nbs);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            bm[br].resize(br+1);
            for(int bc=0; bc<=br; bc++) {
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                bm[br][bc] = ex_to<matrix>(sub_matrix(m, r0, nr, c0, nc));
            }
        }
        return bm;
    }
    
    matrix DEX::b2m(vector<vector<matrix>> & bm) {
        int nbs = bs.size();
        matrix m(N,N);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int bc=0; bc<=br; bc++) {
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) m(r0+r,c0+c) = bm[br][bc](r,c);
            }
        }
        return m;
    }
    
    vector<matrix> DEX::c2b(const matrix & m) {
        int nbs = bs.size();
        int nc = m.cols();
        vector<matrix> bc(nbs);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            bc[br] = ex_to<matrix>(sub_matrix(m, r0, nr, 0, nc));
        }
        return bc;
    }
    
    matrix DEX::b2c(vector<matrix> & bc) {
        int nbs = bs.size();
        int nc = bc[0].cols();
        matrix m(N,nc);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) m(r0+r,c) = bc[br](r,c);
        }
        return m;
    }
    
    map<ex,vector<vector<matrix>>,ex_is_less> DEX::b2m(const block_umat_t & bu) { // U[la][k][n]
        map<ex,vector<vector<matrix>>,ex_is_less> umat;
        int nbs = bs.size();
        for(int br=nbs-1; br>=0; br--) for(int bc=0; bc<=br; bc++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            int c0 = bs[bc].first;
            int nc = bs[bc].second;
            for(auto kv : bu[br][bc]) {
                auto la = kv.first;
                auto kmax = kv.second.size();
                if(umat.find(la)==umat.end()) umat[la].resize(kmax);
                for(int k=0; k<kmax; k++) {
                    auto nmax = kv.second[k].size();
                    if(umat[la][k].size()==0) {
                        umat[la][k].resize(nmax);
                        for(int n=0; n<nmax; n++) umat[la][k][n] = matrix(N,N);
                    }
                    for(int n=0; n<nmax; n++) {
                        auto & mat = umat[la][k][n];
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                            mat(r0+r,c0+c) += kv.second[k][n](r,c);
                        }
                    }
                }
            }
        }
        return umat;
    }
    
    matrix DEX::b2m(const block_umat_t & bu, const ex & x0) {
        int nbs = bs.size();
        matrix mat(N, N);
        for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            int c0 = bs[bc].first;
            int nc = bs[bc].second;
            for(auto kv : bu[br][bc]) {
                auto la = kv.first;
                auto kmax = kv.second.size();
                for(int k=0; k<kmax; k++) {
                    auto nmax = kv.second[k].size();
                    for(int n=0; n<nmax; n++) {
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                            if(is_zero(x0-1)) { if(k==0) mat(r0+r,c0+c) += kv.second[k][n](r,c); }
                            else mat(r0+r,c0+c) += kv.second[k][n](r,c) * pow(x0,la+n) * pow(log(x0),k)/factorial(k);
                        }
                    }
                }
            }
        }
        return mat;
    }
    
    matrix DEX::b2m(const block_imat_t & bi, const ex & x0, const int nc) {
        int nbs = bs.size();
        matrix mat(N, nc);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(auto kv : bi[br]) {
                auto la = kv.first;
                auto kmax = kv.second.size();
                for(int k=0; k<kmax; k++) {
                    auto nmax = kv.second[k].size();
                    for(int n=0; n<nmax; n++) {
                        for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                            mat(r0+r,c) += kv.second[k][n](r,c) * pow(x0,la+n) * pow(log(x0),k)/factorial(k);
                        }
                    }
                }
            }
        }
        return mat;
    }
    
    matrix DEX::b2m(const vector<vector<matrix>> & bi, const ex & x0, const int nc) {
        int nbs = bs.size();
        matrix mat(N, nc);
        for(int br=0; br<nbs; br++) {
            int r0 = bs[br].first;
            int nr = bs[br].second;
            int nmax = bi[br].size();
            for(int n=0; n<nmax; n++) {
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    mat(r0+r,c) += bi[br][n](r,c) * pow(x0,n);
                }
            }
        }
        return mat;
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
    
    void DEX::fuchsify() { // first on diagonal blocks, then off diagnoal ones
        if(fuchsified) return;
        matrix t(N,N);
        auto nbs = bs.size();
        vector<MX> mxt(nbs),mxti(nbs);
        for(int bi=0; bi<nbs; bi++) { // diagonal blocks
            if(!In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                    \r" << flush;
                cout << "  \\--reducing diagonal blocks: " << nbs << "|" << bi+1 << flush;
            }
            auto n0 = bs[bi].first;
            auto n = bs[bi].second;
            matrix ti = ex_to<matrix>(unit_matrix(n,n));
            MX mx(Mat[bi][bi]);
            while(true) { // fuchsify diagonal blocks
                auto pr = mx.prank();
                if(pr<1) break;
                auto p = proj_mat(mx.a01());
                mx.balance(p);
                matrix cop = ex_to<matrix>(unit_matrix(n)).sub(p);
                auto px = p.mul_scalar(1/x);
                ti = ti.mul(cop.sub(px));
                ti = normal_flint(ti);
            }
            while(true) { // shearing transformation on diagonal blocks
                auto pr = mx.prank();
                if(pr<0) break;
                auto a0 = mx.a0();
                auto qj = jordan(a0);
                auto tm = qj.first;
                
                bool shearing = false;
                matrix smat = ex_to<matrix>(unit_matrix(n));
                int cpos = 0;
                int sm = 0;
                for(auto kv : qj.second) {
                    auto ev = kv.first.subs(lst{d==4,ep==0});
                    ex fx = 1;
                    if(sm<=0 && ev<-1/ex(2)) { sm=-1; fx = 1/x; }
                    else if(sm>=0 && ev>=1/ex(2)) { sm=1; fx = x; }
                    for(int j=0; j<kv.second; j++) smat(cpos+j, cpos+j) = fx;
                    if(sm!=0) shearing = true;
                    cpos += kv.second;
                }
                if(shearing) tm = tm.mul(smat);
                mx.transform(tm,tm.inverse());
                ti = ti.mul(tm);
                if(!shearing) break;
            }
            mx(Mat[bi][bi]); // back to Mat[bi][bi]
            mx.clear();
            for(int r=0; r<n; r++) for(int c=0; c<n; c++) t(n0+r,n0+c) = ti(r,c);
            mxt[bi].init(ti);
            mxti[bi].init(ti.inverse());
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        #pragma omp parallel for schedule(runtime)
        for(int br=0; br<nbs; br++) for(int bc=0; bc<br; bc++) {
            MX mx(Mat[br][bc]);
            mx.mul_left(mxti[br]).mul(mxt[bc])(Mat[br][bc]);
            mx.clear();
            flint_cleanup();
        }
        t = normal_flint(t);
        Ts.push_back(t);
        for(int bi=0; bi<nbs; bi++) { mxt[bi].clear(); mxti[bi].clear(); }
        
        for(int br=0; br<nbs; br++) { // off-diagonal blocks
            for(int bc=br-1; bc>=0; bc--) {
                auto nr=bs[br].second, nc=bs[bc].second;
                if(!In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                    \r" << flush;
                    cout << "  \\--reducing off-diagonal blocks: " << nbs << "|" << br+1 << "|" << (br-bc+1) << flush; 
                } 
                matrix sdmat = ex_to<matrix>(symbolic_matrix(nr,nc,"x"));
                matrix trc(nr,nc);
                MX mrr(Mat[br][br]);
                MX mcc(Mat[bc][bc]);
                MX mrc(Mat[br][bc]);
                auto a0_rr = mrr.a0();
                auto a0_cc = mcc.a0();
                while(true) {
                    auto pr = mrc.prank();
                    if(pr<1) break;
                    auto a0_rc = mrc.a0();
                    auto eq_mat = sdmat.mul_scalar(pr).add(a0_rr.mul(sdmat)).sub(sdmat.mul(a0_cc));
                    lst eqs, xs;
                    for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                        xs.append(sdmat(r,c));
                        eqs.append(eq_mat(r,c)+a0_rc(r,c)==0);
                    }
                    auto dmat = ex_to<matrix>(subs(sdmat,lsolve(eqs,xs))).mul_scalar(pow(x,-pr));
                    MX dm1(dmat);
                    MX dm2(dm1), dm3(dm1);
                    dm1.scale(pr/x); // dmat.mul_scalar(pr/x)
                    dm2.mul_left(mrr); // mrr.mul(dmat)
                    dm3.mul(mcc); // dmat.mul(mcc)
                    mrc.add(dm1).add(dm2).sub(dm3);
                    dm1.clear();
                    dm2.clear();
                    dm3.clear();
                    trc = trc.add(dmat);
                }
                if(!trc.is_zero_matrix()) {
                    t = ex_to<matrix>(unit_matrix(N,N));
                    for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                        t(bs[br].first+r, bs[bc].first+c) = trc(r,c);
                    }
                    Ts.push_back(t);
                    mrc(Mat[br][bc]); // back to Mat[br][bc]
                    for(int b=0; b<bc; b++) {
                        MX m1(Mat[br][b]), m2(Mat[bc][b]);
                        m2.mul_left(trc);
                        m1.sub(m2)(Mat[br][b]); // back to Mat[br][b]
                    }
                    for(int a=br+1; a<nbs; a++) {
                        MX m1(Mat[a][bc]), m2(Mat[a][br]);
                        m2.mul(trc);
                        m1.add(m2)(Mat[a][bc]); // back to Mat[a][bc]
                    }
                }
                mrc.clear();
                mrr.clear();
                mcc.clear();
            }
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << "  \\--initializing x^Ao ..." << flush;

        // matrix exponetial, note ln^k x --> ln^k x/k!
        matrix a0(N,N);
        for(int br=0; br<nbs; br++) {
            auto r0 = bs[br].first;
            auto nr = bs[br].second;
            for(int bc=0; bc<=br; bc++) {
                auto c0 = bs[bc].first;
                auto nc = bs[bc].second;
                MX mx(Mat[br][bc]);
                auto m = mx.a0();
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) a0(r0+r,c0+c) = m(r,c);
            }
        }
        auto qj = jordan(a0);
        symbol lnx("lnx");
        exmap la2s; // la 2 symbol
        matrix m(N,N);
        int cur_pos = 0;
        exmap s2x;
        for(auto kv : qj.second) {
            if(la2s.find(kv.first)==la2s.end()) {
                symbol s;
                la2s[kv.first] = s;
                s2x[s] = pow(x,kv.first);
            }
            for(int ir=0;ir<kv.second;ir++) {
                for(int ic=0;ic<kv.second;ic++) {
                    if(ic>=ir) {
                        m(cur_pos+ir, cur_pos+ic)=la2s[kv.first]*pow(lnx,ic-ir); //  removed (ic-ir)! from denominator
                    }
                }
            }
            cur_pos += kv.second;
        }
        m = qj.first.mul(m).mul(qj.first.inverse());
        //cout << subs(m,s2x) << endl;
        
        // initialize UK[a][b][ila] & U0[a][b][ila][k]
        UK.resize(nbs);
        IK.resize(nbs);
        U0.resize(nbs);
        la2i.clear();
        las.clear();
        kmmax = -1;
        nlas = -1;
        int la_i = 0;
        for(int br=0; br<nbs; br++) {
            UK[br].resize(br+1);
            U0[br].resize(br+1);
            int r0 = bs[br].first;
            int nr = bs[br].second;
            for(int bc=br; bc>=0; bc--) {
                int c0 = bs[bc].first;
                int nc = bs[bc].second;
                for(auto kv: la2s) {
                    auto la = kv.first;
                    auto sla = kv.second;
                    int kmax = 0;
                    if(la2i.find(la)==la2i.end()) {
                        la2i[la] = la_i;
                        las.push_back(la);
                        la_i++;
                    }
                    auto ila = la2i[la];
                    if(br>0 && br-1>=bc && UK[br-1][bc].find(ila)!=UK[br-1][bc].end()) kmax = UK[br-1][bc][ila];
                    if(br>=bc+1 && UK[br][bc+1].find(ila)!=UK[br][bc+1].end() && UK[br][bc+1][ila]>kmax) kmax = UK[br][bc+1][ila];
                    for(int r=0; r<nr; r++) {    
                        for(int c=0; c<nc; c++) {
                            auto item = expand_ex(m(r0+r,c0+c),sla).coeff(sla);
                            if(item.is_zero()) continue;
                            auto kd = expand_ex(item,lnx).degree(lnx)+1;
                            if(kd>kmax) kmax = kd;
                        }
                    }
                    if(kmax>0) {
                        if(kmmax<kmax) kmmax = kmax;
                        UK[br][bc][ila] = kmax;
                        U0[br][bc][ila].resize(kmax);
                        for(int k=0; k<kmax; k++) {
                            U0[br][bc][ila][k] = vector<fmpq_mat_t>(1);
                            fmpq_mat_init(U0[br][bc][ila][k][0],nr,nc);
                            matrix mat(nr,nc);
                            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                                mat(r,c) = expand_ex(m(r0+r,c0+c),lst{sla,lnx}).coeff(sla).coeff(lnx,k);
                            }
                            _to_(U0[br][bc][ila][k][0],mat);
                        }
                    }
                }
            }
            for(auto kv : UK[br][0]) IK[br][kv.first] = kv.second;
        }
        nlas = las.size();
        qlas.clear();
        qlas.resize(nlas);
        for(int i=0; i<nlas; i++) {
            qlas[i] = vector<fmpq_t>(1);
            fmpq_init(qlas[i][0]);
            _to_(qlas[i][0],las[i]);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        fuchsified = true;
    }
    
    namespace {
        inline bool is_resonant(const vector<fmpq_t> & qs, fmpq_t q) {
            fmpq_t dq;
            fmpq_init(dq);
            mpz_t nz,dz;
            mpz_init(nz);
            mpz_init(dz);
            fmpz_t fdz;
            fmpz_init(fdz);
            bool res = false;
            for(auto & qi : qs) {
                fmpq_sub(dq,qi,q);
                fmpq_get_mpz_frac(nz,dz,dq);
                fmpz_set_mpz(fdz,dz);
                if(fmpz_is_pm1(fdz)) {
                    res = true;
                    break;
                }
            }
            fmpq_clear(dq);
            mpz_clear(nz);
            mpz_clear(dz);
            fmpz_clear(fdz);
            return res;
        }
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
    
    // U[a][b][la][k][n], block_umat_fmpq_mat_t U; no need to initialize U
    void DEX::series(block_umat_fmpq_mat_t & U, int xn, const vector<fmpq_t> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        U.resize(nbs);
        for(int br=0; br<nbs; br++) { // cycle rows
            int nr = bs[br].second;
            U[br].resize(br+1);
            vector<MX> A(br+1); // M=A/x
            fmpz_poly_t lcm, ilcm, olcm, rlcm; // D=lcm
            fmpz_poly_init(lcm);
            fmpz_poly_init(ilcm);
            fmpz_poly_init(olcm);
            fmpz_poly_init(rlcm);
            fmpz_poly_set_str(lcm, "1  1");
            for(int bc=br; bc>=0; bc--) { // cycle columns
                int nc = bs[bc].second;

                // lcm each block A[br][bc] -> A[bc]
                A[bc].init(Mat[br][bc]);
                A[bc].scale(x); // A=x*M

                A[bc].denlcm(ilcm);
                fmpz_poly_set(olcm,lcm);
                fmpz_poly_lcm(lcm, olcm, ilcm);
                fmpz_poly_div(rlcm, lcm, olcm);
                for(int bc2=br; bc2>bc; bc2--) A[bc2].scale(rlcm);
                fmpz_poly_div(rlcm, lcm, ilcm);
                A[bc].scale(rlcm);
                
                // now we use a=br, b=bc, note that M=A/x
                // [(ila+n)D0-A0aa].Uab(ila,k,n) = -Uab(ila,k+1,n)D0
                //   - sum_{0<m<=n} [(ila+n)Uab(ila,k,n-m)+Uab(ila,k+1,n-m)]Dm 
                //   + sum_{b<=c<=a,0<=m<=n,NO(a=b,m=0)} Amac.Ucb(ila,k,n-m)
                int a = br, b = bc;
                fmpq_mat_t A0aa;
                fmpq_mat_init(A0aa,nr,nr);
                A[a].coeff(A0aa,0);
                fmpz_t D0;
                fmpz_init(D0);
                fmpz_poly_get_coeff_fmpz(D0,lcm,0);
                
                int nla = UK[a][b].size();
                if(!In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--series: " << nbs << "|" << br+1 << "|" << (br-bc+1);
                    cout << " [" << nr << "\u2A09" << nc << "]";
                    cout << " \u03BB" << nla << " n" << xn << flush;
                }
                for(auto kv : UK[a][b]) U[a][b][kv.first].resize(kv.second); // insert kv first and make sure thread-safe
                bool la_parallel = (qslas.size()>0 ? qslas.size() : nla) >= omp_get_max_threads();
                #pragma omp parallel for schedule(runtime) if(la_parallel)
                for(int cla=0; cla<nla; cla++) { // 1-cycle over lambda
                    fmpq_mat_t invA;
                    fmpq_mat_init(invA,nr,nr);
                    fmpq_t q;
                    fmpq_init(q);
                    fmpq_mat_t smat;
                    fmpq_mat_init(smat,nr,nc);
                
                    auto kv = UK[a][b].begin();
                    advance(kv,cla);
                    auto ila = kv->first;
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) {
                        if(!In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--series: " << nbs << "|" << br+1 << "|" << (br-bc+1);
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << " n" << xn << flush;
                        }
                        continue;
                    }
                    auto kmax = kv->second;
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][ila][k] = vector<fmpq_mat_t>(xn+1);
                        fmpq_mat_init(U[a][b][ila][k][0],nr,nc);
                        fmpq_mat_set(U[a][b][ila][k][0],U0[a][b][ila][k][0]);
                        
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            if(omp_get_num_threads()==1 && !In_GiNaC_Parallel && Verbose>5) {
                                cout << "\r                                                              \r" << flush;
                                cout << "  \\--series: " << nbs << "|" << br+1 << "|" << (br-bc+1);
                                cout << " [" << nr << "\u2A09" << nc << "]";
                                cout << " \u03BB" << nla << "|" << cla+1;
                                cout << " k" << kmax << "|" << (kmax-k);
                                cout << " n" << xn << "|" << n << flush;
                            }
                            
                            fmpq_mat_zero(smat);
                            fmpq_add_si(q,qlas[ila][0],n); // q = la+n

                            if(true) { 
                                slong s = fmpz_poly_degree(lcm);
                                if(s>n) s=n;
                                vector<fmpq_mat_t> mat_vec(s+1);
                                for(int m=0; m<=s; m++) fmpq_mat_init(mat_vec[m],nr,nc);
                                #pragma omp parallel for schedule(runtime)
                                for(int m=0; m<=s; m++) {
                                    auto mat = mat_vec[m];
                                    fmpz_t Dm;
                                    fmpz_init(Dm);
                                    fmpq_t qm;
                                    fmpq_init(qm);
                                    if(m==0 && k+1<kmax) {
                                        fmpz_neg(Dm,D0);
                                        fmpq_mat_scalar_mul_fmpz(mat,U[a][b][ila][k+1][n],Dm);
                                    } else if(m>0) {
                                        fmpq_sub_si(qm,q,m);
                                        fmpq_mat_scalar_mul_fmpq(mat,U[a][b][ila][k][n-m],qm);
                                        if(k+1<kmax) fmpq_mat_add(mat,mat,U[a][b][ila][k+1][n-m]);
                                        fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                                        fmpz_neg(Dm,Dm);
                                        fmpq_mat_scalar_mul_fmpz(mat,mat,Dm);
                                    }
                                    fmpz_clear(Dm);
                                    fmpq_clear(qm);
                                    flint_cleanup();
                                }
                                for(int m=0; m<=s; m++) {
                                    fmpq_mat_add(smat,smat,mat_vec[m]);
                                    fmpq_mat_clear(mat_vec[m]);
                                }
                            }

                            if(true) {
                                int ab1 = a-b+1;
                                vector<fmpq_mat_t> smat_vec(ab1);
                                for(int i=0; i<ab1; i++) fmpq_mat_init(smat_vec[i],nr,nc);
                                #pragma omp parallel for schedule(runtime) if(a>omp_get_max_threads())
                                for(int c=b; c<=a; c++) {
                                    if(U[c][b].find(ila)!=U[c][b].end() && U[c][b][ila].size()>k) {
                                        slong s = A[c].degree();
                                        if(s>n) s = n;
                                        vector<fmpq_mat_t> mat_vec(s+1);
                                        for(int i=0; i<=s; i++) fmpq_mat_init(mat_vec[i],nr,nc);
                                        #pragma omp parallel for schedule(runtime)
                                        for(int m=0; m<=s; m++) {
                                            if(c!=a || m!=0) {
                                                int nc2 = bs[c].second;
                                                fmpq_mat_t Amac;
                                                fmpq_mat_init(Amac,nr,nc2);
                                                A[c].coeff(Amac,m);
                                                fmpq_mat_mul(mat_vec[m],Amac,U[c][b][ila][k][n-m]);
                                                fmpq_mat_clear(Amac);
                                            }
                                            flint_cleanup();
                                        }
                                        fmpq_mat_zero(smat_vec[c-b]);
                                        for(int i=0; i<=s; i++) {
                                            fmpq_mat_add(smat_vec[c-b],smat_vec[c-b],mat_vec[i]);
                                            fmpq_mat_clear(mat_vec[i]);
                                        }
                                    }
                                    flint_cleanup();
                                }
                                for(int i=0; i<ab1; i++) {
                                    fmpq_mat_add(smat,smat,smat_vec[i]);
                                    fmpq_mat_clear(smat_vec[i]);
                                }
                            }
                            
                            fmpq_mul_fmpz(q,q,D0); // q = (la+n)D0
                            fmpq_mat_one(invA);
                            fmpq_mat_scalar_mul_fmpq(invA,invA,q);
                            fmpq_mat_sub(invA,invA,A0aa);
                            fmpq_mat_inv(invA,invA);

                            fmpq_mat_init(U[a][b][ila][k][n],nr,nc);
                            fmpq_mat_mul(U[a][b][ila][k][n],invA,smat);
                        }
                    }
                    fmpq_mat_clear(invA);
                    fmpq_clear(q);
                    fmpq_mat_clear(smat);
                    flint_cleanup();
                }
                fmpq_mat_clear(A0aa);
                fmpz_clear(D0);
            }

            for(int bc=br; bc>=0; bc--) A[bc].clear();
            fmpz_poly_clear(lcm);
            fmpz_poly_clear(ilcm);
            fmpz_poly_clear(olcm);
            fmpz_poly_clear(rlcm);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
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
    
    // U[a][b][la][k][n], block_umat_acb_mat_t U; no need to initialize U
    void DEX::series(block_umat_acb_mat_t & U, int xn, slong dp, const vector<fmpq_t> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        U.resize(nbs);
        for(int br=0; br<nbs; br++) { // cycle rows
            int nr = bs[br].second;
            U[br].resize(br+1);
            vector<MX> A(br+1); // M=A/x
            fmpz_poly_t lcm, ilcm, olcm, rlcm; // D=lcm
            fmpz_poly_init(lcm);
            fmpz_poly_init(ilcm);
            fmpz_poly_init(olcm);
            fmpz_poly_init(rlcm);
            fmpz_poly_set_str(lcm, "1  1");
            for(int bc=br; bc>=0; bc--) { // cycle columns
                int nc = bs[bc].second;

                // lcm each block A[br][bc] -> A[bc]
                A[bc].init(Mat[br][bc]);
                A[bc].scale(x); // A=x*M

                A[bc].denlcm(ilcm);
                fmpz_poly_set(olcm,lcm);
                fmpz_poly_lcm(lcm, olcm, ilcm);
                fmpz_poly_div(rlcm, lcm, olcm);
                for(int bc2=br; bc2>bc; bc2--) A[bc2].scale(rlcm);
                fmpz_poly_div(rlcm, lcm, ilcm);
                A[bc].scale(rlcm);
                
                // now we use a=br, b=bc, note that M=A/x
                // [(ila+n)D0-A0aa].Uab(ila,k,n) = -Uab(ila,k+1,n)D0
                //   - sum_{0<m<=n} [(ila+n)Uab(ila,k,n-m)+Uab(ila,k+1,n-m)]Dm 
                //   + sum_{b<=c<=a,0<=m<=n,NO(a=b,m=0)} Amac.Ucb(ila,k,n-m)
                int a = br, b = bc;
                fmpq_mat_t A0aa;
                fmpq_mat_init(A0aa,nr,nr);
                A[a].coeff(A0aa,0);
                fmpz_t D0;
                fmpz_init(D0);
                fmpz_poly_get_coeff_fmpz(D0,lcm,0);
                
                int nla = UK[a][b].size();
                if(!In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--nseries: " << nbs << "|" << br+1 << "|" << (br-bc+1);
                    cout << " [" << nr << "\u2A09" << nc << "]";
                    cout << " \u03BB" << nla << " n" << xn << flush;
                }
                for(auto kv : UK[a][b]) U[a][b][kv.first].resize(kv.second); // insert kv first and make sure thread-safe
                bool la_parallel = (qslas.size()>0 ? qslas.size() : nla) >= omp_get_max_threads();
                #pragma omp parallel for schedule(runtime) if(la_parallel)
                for(int cla=0; cla<nla; cla++) { // 1-cycle over lambda
                    fmpq_mat_t invA;
                    fmpq_mat_init(invA,nr,nr);
                    acb_mat_t invAf;
                    acb_mat_init(invAf,nr,nr);
                    fmpq_t q;
                    fmpq_init(q);
                    acb_mat_t smat;
                    acb_mat_init(smat,nr,nc);
                
                    auto kv = UK[a][b].begin();
                    advance(kv,cla);
                    auto ila = kv->first;
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) {
                        if(!In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--nseries: " << nbs << "|" << br+1 << "|" << (br-bc+1);
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << " n" << xn << flush;
                        }
                        continue;
                    }
                    auto kmax = kv->second;
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][ila][k] = vector<acb_mat_t>(xn+1);
                        acb_mat_init(U[a][b][ila][k][0],nr,nc);
                        acb_mat_set_fmpq_mat(U[a][b][ila][k][0],U0[a][b][ila][k][0],fp);
                        
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            if(omp_get_num_threads()==1 && !In_GiNaC_Parallel && Verbose>5) {
                                cout << "\r                                                              \r" << flush;
                                cout << "  \\--nseries: " << nbs << "|" << br+1 << "|" << (br-bc+1);
                                cout << " [" << nr << "\u2A09" << nc << "]";
                                cout << " \u03BB" << nla << "|" << cla+1;
                                cout << " k" << kmax << "|" << (kmax-k);
                                cout << " n" << xn << "|" << n << flush;
                            }
                            
                            acb_mat_zero(smat);
                            fmpq_add_si(q,qlas[ila][0],n); // q = la+n

                            if(true) { 
                                slong s = fmpz_poly_degree(lcm);
                                if(s>n) s=n;
                                vector<acb_mat_t> mat_vec(s+1);
                                for(int m=0; m<=s; m++) acb_mat_init(mat_vec[m],nr,nc);
                                #pragma omp parallel for schedule(runtime)
                                for(int m=0; m<=s; m++) {
                                    auto mat = mat_vec[m];
                                    fmpz_t Dm;
                                    fmpz_init(Dm);
                                    fmpq_t qm;
                                    fmpq_init(qm);
                                    if(m==0 && k+1<kmax) {
                                        fmpz_neg(Dm,D0);
                                        acb_mat_scalar_mul_fmpz(mat,U[a][b][ila][k+1][n],Dm,fp);
                                    } else if(m>0) {
                                        fmpq_sub_si(qm,q,m);
                                        acb_t fqm; acb_init(fqm);
                                        acb_set_fmpq(fqm,qm,fp);
                                        acb_mat_scalar_mul_acb(mat,U[a][b][ila][k][n-m],fqm,fp);
                                        acb_clear(fqm);
                                        if(k+1<kmax) acb_mat_add(mat,mat,U[a][b][ila][k+1][n-m],fp);
                                        fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                                        fmpz_neg(Dm,Dm);
                                        acb_mat_scalar_mul_fmpz(mat,mat,Dm,fp);
                                    }
                                    fmpz_clear(Dm);
                                    fmpq_clear(qm);
                                    flint_cleanup();
                                }
                                for(int m=0; m<=s; m++) {
                                    acb_mat_add(smat,smat,mat_vec[m],fp);
                                    acb_mat_clear(mat_vec[m]);
                                }
                            }

                            if(true) {
                                int ab1 = a-b+1;
                                vector<acb_mat_t> smat_vec(ab1);
                                for(int i=0; i<ab1; i++) acb_mat_init(smat_vec[i],nr,nc);
                                #pragma omp parallel for schedule(runtime) if(a>omp_get_max_threads())
                                for(int c=b; c<=a; c++) {
                                    if(U[c][b].find(ila)!=U[c][b].end() && U[c][b][ila].size()>k) {
                                        slong s = A[c].degree();
                                        if(s>n) s = n;
                                        vector<acb_mat_t> mat_vec(s+1);
                                        for(int i=0; i<=s; i++) acb_mat_init(mat_vec[i],nr,nc);
                                        #pragma omp parallel for schedule(runtime)
                                        for(int m=0; m<=s; m++) {
                                            if(c!=a || m!=0) {
                                                int nc2 = bs[c].second;
                                                acb_mat_t Amac;
                                                acb_mat_init(Amac,nr,nc2);
                                                A[c].coeff(Amac,m,fp);
                                                acb_mat_mul(mat_vec[m],Amac,U[c][b][ila][k][n-m],fp);
                                                acb_mat_clear(Amac);
                                            }
                                            flint_cleanup();
                                        }
                                        acb_mat_zero(smat_vec[c-b]);
                                        for(int i=0; i<=s; i++) {
                                            acb_mat_add(smat_vec[c-b],smat_vec[c-b],mat_vec[i],fp);
                                            acb_mat_clear(mat_vec[i]);
                                        }
                                    }
                                    flint_cleanup();
                                }
                                for(int i=0; i<ab1; i++) {
                                    acb_mat_add(smat,smat,smat_vec[i],fp);
                                    acb_mat_clear(smat_vec[i]);
                                }
                            }
                            
                            fmpq_mul_fmpz(q,q,D0); // q = (la+n)D0
                            fmpq_mat_one(invA);
                            fmpq_mat_scalar_mul_fmpq(invA,invA,q);
                            fmpq_mat_sub(invA,invA,A0aa);
                            fmpq_mat_inv(invA,invA);
                            acb_mat_set_fmpq_mat(invAf,invA,fp);

                            acb_mat_init(U[a][b][ila][k][n],nr,nc);
                            acb_mat_mul(U[a][b][ila][k][n],invAf,smat,fp);
                        }
                    }
                    fmpq_mat_clear(invA);
                    acb_mat_clear(invAf);
                    fmpq_clear(q);
                    acb_mat_clear(smat);
                    flint_cleanup();
                }
                fmpq_mat_clear(A0aa);
                fmpz_clear(D0);
            }

            for(int bc=br; bc>=0; bc--) A[bc].clear();
            fmpz_poly_clear(lcm);
            fmpz_poly_clear(ilcm);
            fmpz_poly_clear(olcm);
            fmpz_poly_clear(rlcm);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
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
    
    // I[a][la][k][n] & In0[a][ila][k], no need to initialize I
    void DEX::series(block_imat_fmpq_mat_t & I, int xn, block_imat_fmpq_mat_t & In0, int nc, const vector<fmpq_t> & qslas) { 
        
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        I.resize(nbs);
        for(int br=0; br<nbs; br++) { // cycle rows
            int nr = bs[br].second;
            vector<MX> A(br+1); // M=A/x
            fmpz_poly_t lcm, rlcm; // D=lcm
            vector<fmpz_poly_t> lcm_vec(br+1);
            fmpz_poly_init(lcm);
            fmpz_poly_init(rlcm);
            fmpz_poly_set_str(lcm, "1  1");
            for(int bc=br; bc>=0; bc--) { 
                // lcm for each block A[br][bc] -> A[bc] 
                A[bc].init(Mat[br][bc]);
                A[bc].scale(x); // A=x*M
                fmpz_poly_init(lcm_vec[bc]);
                A[bc].denlcm(lcm_vec[bc]);
                fmpz_poly_lcm(lcm, lcm, lcm_vec[bc]);
            }
            for(int bc=br; bc>=0; bc--) {
                fmpz_poly_div(rlcm, lcm, lcm_vec[bc]);
                A[bc].scale(rlcm);
                fmpz_poly_clear(lcm_vec[bc]);
            }
            fmpz_poly_clear(rlcm);
            
            // now we use a=br
            // [(la+n)D0-A0aa].Ia(la,k,n) = -D0 Ia(la,k+1,n)
            //   - sum_{0<m<=n} Dm [(la+n-m).Ia(la,k,n-m)+Ia(la,k+1,n-m)] 
            //   + sum_{not(b=a|m=0)} Amab.Ib(la,k,n-m)
            int a = br;
            fmpq_mat_t A0aa;
            fmpq_mat_init(A0aa,nr,nr);
            A[a].coeff(A0aa,0);
            fmpz_t D0;
            fmpz_init(D0);
            fmpz_poly_get_coeff_fmpz(D0,lcm,0);
            
            int nla = IK[a].size();
            if(!In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                              \r" << flush;
                cout << "  \\--series: " << nbs << "|" << br+1;
                cout << " [" << nr << "\u2A09" << nc << "]";
                cout << " \u03BB" << nla << " n" << xn << flush;
            }
            for(auto kv : IK[a]) I[a][kv.first].resize(kv.second); // insert kv first and make sure thread-safe
            bool la_parallel = (qslas.size()>0 ? qslas.size() : nla) >= omp_get_max_threads();
            #pragma omp parallel for schedule(runtime) if(la_parallel)
            for(int cla=0; cla<nla; cla++) { // 1-cycle over lambda
                fmpq_mat_t invA;
                fmpq_mat_init(invA,nr,nr);
                fmpq_t q;
                fmpq_init(q);
                fmpq_mat_t smat;
                fmpq_mat_init(smat,nr,nc);
                
                auto kv = IK[a].begin();
                advance(kv,cla);
                auto ila = kv->first;
                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) {
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                                              \r" << flush;
                        cout << "  \\--series: " << nbs << "|" << br+1;
                        cout << " [" << nr << "\u2A09" << nc << "]";
                        cout << " \u03BB" << nla << " n" << xn << flush;
                    }
                    continue;
                }
                auto kmax = kv->second;  
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    I[a][ila][k] = vector<fmpq_mat_t>(xn+1);
                    fmpq_mat_init(I[a][ila][k][0],nr,nc);
                    fmpq_mat_set(I[a][ila][k][0],In0[a][ila][k][0]);
                    
                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                    
                        if(omp_get_num_threads()==1 && !In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--series: " << nbs << "|" << br+1;
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << "|" << cla+1;
                            cout << " k" << kmax << "|" << (kmax-k);
                            cout << " n" << xn << "|" << n << flush;
                        }
                  
                        fmpq_mat_zero(smat);
                        fmpq_add_si(q,qlas[ila][0],n); // q = la+n

                        if(true) { 
                            slong s = fmpz_poly_degree(lcm);
                            if(s>n) s=n;
                            vector<fmpq_mat_t> mat_vec(s+1);
                            for(int m=0; m<=s; m++) fmpq_mat_init(mat_vec[m],nr,nc);
                            #pragma omp parallel for schedule(runtime)
                            for(int m=0; m<=s; m++) {
                                auto mat = mat_vec[m];
                                fmpz_t Dm;
                                fmpz_init(Dm);
                                fmpq_t qm;
                                fmpq_init(qm);
                                if(m==0 && k+1<kmax) {
                                    fmpz_neg(Dm,D0);
                                    fmpq_mat_scalar_mul_fmpz(mat,I[a][ila][k+1][n],Dm);
                                } else if(m>0) {
                                    fmpq_sub_si(qm,q,m);
                                    fmpq_mat_scalar_mul_fmpq(mat,I[a][ila][k][n-m],qm);
                                    if(k+1<kmax) fmpq_mat_add(mat,mat,I[a][ila][k+1][n-m]);
                                    fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                                    fmpz_neg(Dm,Dm);
                                    fmpq_mat_scalar_mul_fmpz(mat,mat,Dm);
                                }
                                fmpz_clear(Dm);
                                fmpq_clear(qm);
                                flint_cleanup();
                            }
                            for(int m=0; m<=s; m++) {
                                fmpq_mat_add(smat,smat,mat_vec[m]);
                                fmpq_mat_clear(mat_vec[m]);
                            }
                        }
                        
                        if(true) {
                            vector<fmpq_mat_t> smat_vec(a+1);
                            for(int i=0; i<=a; i++) fmpq_mat_init(smat_vec[i],nr,nc);
                            #pragma omp parallel for schedule(runtime) if(a>omp_get_max_threads())
                            for(int b=0; b<=a; b++) {
                                if(I[b].find(ila)!=I[b].end() && I[b][ila].size()>k) {
                                    slong s = A[b].degree();
                                    if(s>n) s = n;
                                    vector<fmpq_mat_t> mat_vec(s+1);
                                    for(int i=0; i<=s; i++) fmpq_mat_init(mat_vec[i],nr,nc);
                                    #pragma omp parallel for schedule(runtime)
                                    for(int m=0; m<=s; m++) {
                                        if(b!=a || m!=0) {
                                            int nc2 = bs[b].second;
                                            fmpq_mat_t Amab;
                                            fmpq_mat_init(Amab,nr,nc2);
                                            A[b].coeff(Amab,m);
                                            fmpq_mat_mul(mat_vec[m],Amab,I[b][ila][k][n-m]);
                                            fmpq_mat_clear(Amab);
                                        }
                                        flint_cleanup();
                                    }
                                    fmpq_mat_zero(smat_vec[b]);
                                    for(int i=0; i<=s; i++) {
                                        fmpq_mat_add(smat_vec[b],smat_vec[b],mat_vec[i]);
                                        fmpq_mat_clear(mat_vec[i]);
                                    }
                                }
                                flint_cleanup();
                            }

                            for(int i=0; i<=a; i++) {
                                fmpq_mat_add(smat,smat,smat_vec[i]);
                                fmpq_mat_clear(smat_vec[i]);
                            }
                        }
                        
                        fmpq_mul_fmpz(q,q,D0); // q = (la+n)D0
                        fmpq_mat_one(invA);
                        fmpq_mat_scalar_mul_fmpq(invA,invA,q);
                        fmpq_mat_sub(invA,invA,A0aa);
                        fmpq_mat_inv(invA,invA);
                        
                        fmpq_mat_init(I[a][ila][k][n],nr,nc);
                        fmpq_mat_mul(I[a][ila][k][n],invA,smat);
                    }
                }
                fmpq_mat_clear(invA);
                fmpq_clear(q);
                fmpq_mat_clear(smat);
                flint_cleanup();
            }
            for(int bc=br; bc>=0; bc--) A[bc].clear();
            fmpz_poly_clear(lcm);
            fmpq_mat_clear(A0aa);
            fmpz_clear(D0);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        
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
        int nc = m.cols();
        auto cb = c2b(m);
        
        auto fp = dp2fp(dp);
        block_imat_acb_mat_t In0(nbs);
        for(int br=0; br<nbs; br++) {
            int a = br;
            int nr = bs[br].second;
            acb_mat_t qmat;
            acb_mat_init(qmat,nr,nc);
            for(int bc=0; bc<=br; bc++) {
                int nr2 = bs[bc].second;
                acb_mat_t qmat2, umat;
                acb_mat_init(qmat2,nr2,nc);
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
                        _to_(qmat2,cb[bc],fp);
                        acb_mat_set_fmpq_mat(umat,kv.second[k][0],fp);
                        acb_mat_mul(qmat,umat,qmat2,fp);
                        acb_mat_add(In0[a][ila][k][0],In0[a][ila][k][0],qmat,fp);
                    }
                }
                acb_mat_clear(qmat2);
                acb_mat_clear(umat);
            }
            acb_mat_clear(qmat);
        }
        
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
    
    // I[a][la][k][n] & In0[a][la][k], no need to initialize I
    void DEX::series(block_imat_acb_mat_t & I, int xn, block_imat_acb_mat_t & In0, int nc, slong dp, const vector<fmpq_t> & qslas) {  
    
        if(!fuchsified) fuchsify();
        auto fp = dp2fp(dp);
        auto nbs = bs.size();
        I.resize(nbs);
        for(int br=0; br<nbs; br++) { // cycle rows
            int nr = bs[br].second;
            vector<MX> A(br+1); // M=A/x
            fmpz_poly_t zlcm, rlcm; // D=lcm
            vector<fmpz_poly_t> lcm_vec(br+1);
            fmpz_poly_init(zlcm);
            fmpz_poly_init(rlcm);
            fmpz_poly_set_str(zlcm, "1  1");
            for(int bc=br; bc>=0; bc--) { 
                // lcm for each block A[br][bc] -> A[bc] 
                A[bc].init(Mat[br][bc]);
                A[bc].scale(x); // A=x*M
                fmpz_poly_init(lcm_vec[bc]);
                A[bc].denlcm(lcm_vec[bc]);
                fmpz_poly_lcm(zlcm, zlcm, lcm_vec[bc]);
            }
            for(int bc=br; bc>=0; bc--) {
                fmpz_poly_div(rlcm, zlcm, lcm_vec[bc]);
                A[bc].scale(rlcm);
                fmpz_poly_clear(lcm_vec[bc]);
            }
            acb_poly_t lcm;
            acb_poly_init(lcm);
            acb_poly_set_fmpz_poly(lcm,zlcm,fp);
            fmpz_poly_clear(zlcm);
            fmpz_poly_clear(rlcm);
            
            // now we use a=br
            // [(la+n)D0-A0aa].Ia(la,k,n) = -D0 Ia(la,k+1,n)
            //   - sum_{0<m<=n} Dm [(la+n-m).Ia(la,k,n-m)+Ia(la,k+1,n-m)] 
            //   + sum_{not(b=a|m=0)} Amab.Ib(la,k,n-m)
            int a = br;
            acb_mat_t A0aa;
            acb_mat_init(A0aa,nr,nr);
            A[a].coeff(A0aa,0,fp);
            acb_t D0;
            acb_init(D0);
            acb_poly_get_coeff_acb(D0,lcm,0);
            
            int nla = IK[a].size();
            if(!In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                              \r" << flush;
                cout << "  \\--nseries: " << nbs << "|" << br+1;
                cout << " [" << nr << "\u2A09" << nc << "]";
                cout << " \u03BB" << nla << " n" << xn << flush;
            }
            for(auto kv : IK[a]) I[a][kv.first].resize(kv.second); // insert kv first and make sure thread-safe
            bool la_parallel = (qslas.size()>0 ? qslas.size() : nla) >= omp_get_max_threads();
            #pragma omp parallel for schedule(runtime) if(la_parallel)
            for(int cla=0; cla<nla; cla++) { // 1-cycle over lambda
                acb_mat_t invA;
                acb_mat_init(invA,nr,nr);
                acb_t q;
                acb_init(q);
                acb_mat_t smat;
                acb_mat_init(smat,nr,nc);
            
                auto kv = IK[a].begin();
                advance(kv,cla);
                auto ila = kv->first;
                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) {
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                                              \r" << flush;
                        cout << "  \\--nseries: " << nbs << "|" << br+1;
                        cout << " [" << nr << "\u2A09" << nc << "]";
                        cout << " \u03BB" << nla << " n" << xn << flush;
                    }
                    continue;
                }
                auto kmax = kv->second;
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    I[a][ila][k] = vector<acb_mat_t>(xn+1);
                    acb_mat_init(I[a][ila][k][0],nr,nc);
                    acb_mat_set(I[a][ila][k][0],In0[a][ila][k][0]);
                    
                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                    
                        if(omp_get_num_threads()==1 && !In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--nseries: " << nbs << "|" << br+1;
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << "|" << cla+1;
                            cout << " k" << kmax << "|" << (kmax-k);
                            cout << " n" << xn << "|" << n << flush;
                        }
                  
                        acb_mat_zero(smat);
                        acb_set_fmpq(q,qlas[ila][0],fp);
                        acb_add_si(q,q,n,fp); // q = la+n

                        if(true) {
                            slong s = acb_poly_degree(lcm);
                            if(s>n) s=n;
                            vector<acb_mat_t> mat_vec(s+1);
                            for(int m=0; m<=s; m++) acb_mat_init(mat_vec[m],nr,nc);
                            #pragma omp parallel for schedule(runtime)
                            for(int m=0; m<=s; m++) {
                                auto mat = mat_vec[m];
                                acb_t Dm;
                                acb_init(Dm);
                                acb_t qm;
                                acb_init(qm);
                                if(m==0 && k+1<kmax) {
                                    acb_neg(Dm,D0);
                                    acb_mat_scalar_mul_acb(mat,I[a][ila][k+1][n],Dm,fp);
                                } else if(m>0) {
                                    acb_sub_si(qm,q,m,fp);
                                    acb_mat_scalar_mul_acb(mat,I[a][ila][k][n-m],qm,fp);
                                    if(k+1<kmax) acb_mat_add(mat,mat,I[a][ila][k+1][n-m],fp);
                                    acb_poly_get_coeff_acb(Dm,lcm,m);
                                    acb_neg(Dm,Dm);
                                    acb_mat_scalar_mul_acb(mat,mat,Dm,fp);
                                }
                                acb_clear(Dm);
                                acb_clear(qm);
                                flint_cleanup();
                            }
                            for(int m=0; m<=s; m++) {
                                acb_mat_add(smat,smat,mat_vec[m],fp);
                                acb_mat_clear(mat_vec[m]);
                            }
                        }

                        if(true) {
                            vector<acb_mat_t> smat_vec(a+1);
                            for(int i=0; i<=a; i++) acb_mat_init(smat_vec[i],nr,nc);
                            #pragma omp parallel for schedule(runtime) if(a>omp_get_max_threads())
                            for(int b=0; b<=a; b++) {
                                if(I[b].find(ila)!=I[b].end() && I[b][ila].size()>k) {
                                    slong s = A[b].degree();
                                    if(s>n) s = n;
                                    vector<acb_mat_t> mat_vec(s+1);
                                    for(int i=0; i<=s; i++) acb_mat_init(mat_vec[i],nr,nc);
                                    #pragma omp parallel for schedule(runtime)
                                    for(int m=0; m<=s; m++) {
                                        if(b!=a || m!=0) {
                                            int nc2 = bs[b].second;                                            
                                            acb_mat_t Amab;
                                            acb_mat_init(Amab,nr,nc2);
                                            A[b].coeff(Amab,m,fp);
                                            acb_mat_mul(mat_vec[m],Amab,I[b][ila][k][n-m],fp);
                                            acb_mat_clear(Amab);
                                        }
                                        flint_cleanup();
                                    }
                                    acb_mat_zero(smat_vec[b]);
                                    for(int i=0; i<=s; i++) {
                                        acb_mat_add(smat_vec[b],smat_vec[b],mat_vec[i],fp);
                                        acb_mat_clear(mat_vec[i]);
                                    }
                                }
                                flint_cleanup();
                            }

                            for(int i=0; i<=a; i++) {
                                acb_mat_add(smat,smat,smat_vec[i],fp);
                                acb_mat_clear(smat_vec[i]);
                            }
                        }

                        acb_mul(q,q,D0,fp); // q = (la+n)D0
                        acb_mat_one(invA);
                        acb_mat_scalar_mul_acb(invA,invA,q,fp);
                        acb_mat_sub(invA,invA,A0aa,fp);
                        acb_mat_inv(invA,invA,fp);
                        
                        acb_mat_init(I[a][ila][k][n],nr,nc);
                        acb_mat_mul(I[a][ila][k][n],invA,smat,fp);
                    }
                }
                acb_mat_clear(invA);
                acb_clear(q);
                acb_mat_clear(smat);
                flint_cleanup();
            }
            for(int bc=br; bc>=0; bc--) A[bc].clear();
            acb_poly_clear(lcm);
            acb_mat_clear(A0aa);
            acb_clear(D0);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        
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
    
    void DEX::taylor(vector<vector<fmpq_mat_t>> & I, int xn, const matrix I0, const ex & x0) {  
        auto nbs = bs.size();
        int nc = I0.cols();
 
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            
            vector<MX> A(br+1);
            vector<fmpz_poly_t> lcm_vec(br+1);
            fmpz_poly_t lcm, rlcm;
            fmpz_poly_init(lcm);
            fmpz_poly_init(rlcm);
            fmpz_poly_set_str(lcm, "1  1");
            for(int bc=br; bc>=0; bc--) { 
                // lcm for each block A[br][bc] -> A[bc] 
                A[bc].init(Mat[br][bc]);
                A[bc].shift(x0); // x -> x+x0
                fmpz_poly_init(lcm_vec[bc]);
                A[bc].denlcm(lcm_vec[bc]);
                fmpz_poly_lcm(lcm, lcm, lcm_vec[bc]);
            }
            for(int bc=br; bc>=0; bc--) {
                fmpz_poly_div(rlcm, lcm, lcm_vec[bc]);
                A[bc].scale(rlcm);
                fmpz_poly_clear(lcm_vec[bc]);
            }
            
            // now we use a=br
            // n D0 Ia(n) = - sum_{0<m<n} (n-m) Dm Ia(n-m) + sum_{b<=a} sum_{0<=m<n} Amab.Ib(n-1-m)
            int a = br;
            fmpz_t D0;
            fmpz_init(D0);
            fmpz_poly_get_coeff_fmpz(D0,lcm,0);
            if(fmpz_is_zero(D0)) throw Error("taylor: D0 is zero.");
            
            I[a] = vector<fmpq_mat_t>(xn+1);
            fmpq_mat_init(I[a][0],nr,nc);
            _to_(I[a][0],ex_to<matrix>(sub_matrix(I0, r0, nr, 0, nc)));
            
            fmpq_t q;
            fmpq_init(q);
            fmpq_mat_t smat;
            fmpq_mat_init(smat,nr,nc);
                
            for(int n=1; n<=xn; n++) { // 3-cycle over n
            
                if(omp_get_num_threads()==1 && !In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--taylor: " << nbs << "|" << br+1;
                    cout << " [" << nr << "\u2A09" << nc << "] n" << xn << "|" << n << flush;
                }
                
                fmpq_mat_zero(smat);
                
                if(true) {
                    slong s = fmpz_poly_degree(lcm);
                    if(s>n-1) s=n-1;
                    vector<fmpq_mat_t> mat_vec(s+1);
                    for(int m=1; m<=s; m++) fmpq_mat_init(mat_vec[m],nr,nc);
                    #pragma omp parallel for schedule(runtime)
                    for(int m=0; m<=s; m++) {
                        auto mat = mat_vec[m];
                        fmpz_t Dm;
                        fmpz_init(Dm);
                        fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                        fmpz_mul_si(Dm,Dm,n-m);
                        fmpz_neg(Dm,Dm);
                        fmpq_mat_scalar_mul_fmpz(mat,I[a][n-m],Dm);
                        fmpz_clear(Dm);
                        flint_cleanup();
                    }
                    for(int m=1; m<=s; m++) {
                        fmpq_mat_add(smat,smat,mat_vec[m]);
                        fmpq_mat_clear(mat_vec[m]);
                    }
                }
          
                if(true) {
                    vector<fmpq_mat_t> smat_vec(a+1);
                    for(int i=0; i<=a; i++) fmpq_mat_init(smat_vec[i],nr,nc);
                    #pragma omp parallel for schedule(runtime) if(a>omp_get_max_threads())
                    for(int b=0; b<=a; b++) {
                        slong s = A[b].degree();
                        if(s>n-1) s = n-1;
                        int nc2 = bs[b].second;
                        
                        vector<fmpq_mat_t> mat_vec(s+1);
                        for(int i=0; i<=s; i++) fmpq_mat_init(mat_vec[i],nr,nc);
                        #pragma omp parallel for schedule(runtime)
                        for(int m=0; m<=s; m++) {
                            fmpq_mat_t Amab;
                            fmpq_mat_init(Amab,nr,nc2); 
                            A[b].coeff(Amab,m);
                            fmpq_mat_mul(mat_vec[m],Amab,I[b][n-1-m]);
                            fmpq_mat_clear(Amab);
                            flint_cleanup();
                        }
                        fmpq_mat_zero(smat_vec[b]);
                        for(int i=0; i<=s; i++) {
                            fmpq_mat_add(smat_vec[b],smat_vec[b],mat_vec[i]);
                            fmpq_mat_clear(mat_vec[i]);
                        }                        
                        flint_cleanup();
                    }

                    for(int i=0; i<=a; i++) {
                        fmpq_mat_add(smat,smat,smat_vec[i]);
                        fmpq_mat_clear(smat_vec[i]);
                    }
                }

                fmpq_set_si(q,n,1); 
                fmpq_mul_fmpz(q,q,D0); // q = n D0
                fmpq_inv(q,q); // q = 1 / (n D0)
                fmpq_mat_init(I[a][n],nr,nc);
                fmpq_mat_scalar_mul_fmpq(I[a][n],smat,q);
            }
            
            for(int bc=br; bc>=0; bc--) A[bc].clear();
            fmpz_poly_clear(lcm);
            fmpq_clear(q);
            fmpz_clear(D0);
            fmpq_mat_clear(smat);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
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
    
    void DEX::taylor(vector<vector<acb_mat_t>> & I, int xn, const matrix I0, const ex & x0, slong dp) {  
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        int nc = I0.cols();
        
        if(!ntaylor_inited) {
            TMat.resize(nbs);
            TD.resize(nbs);
            fmpz_poly_t zlcm, rlcm; // D=lcm
            for(int br=0; br<nbs; br++) { // cycle rows
                TMat[br].resize(br+1);
                TD[br] = vector<acb_poly_t>(1);
                int nr = bs[br].second;
                vector<MX> AX(br+1);
                vector<fmpz_poly_t> lcm_vec(br+1);
                fmpz_poly_init(zlcm);
                fmpz_poly_init(rlcm);
                fmpz_poly_set_str(zlcm, "1  1");
                for(int bc=br; bc>=0; bc--) { 
                    // lcm for each block A[br][bc] -> A[bc] 
                    AX[bc].init(Mat[br][bc]);
                    fmpz_poly_init(lcm_vec[bc]);
                    AX[bc].denlcm(lcm_vec[bc]);
                    fmpz_poly_lcm(zlcm, zlcm, lcm_vec[bc]);
                }
                for(int bc=br; bc>=0; bc--) {
                    int nc2 = bs[bc].second;
                    TMat[br][bc].resize(nr);
                    for(int r=0; r<nr; r++) TMat[br][bc][r] = vector<acb_poly_t>(nc2);
                    fmpz_poly_div(rlcm, zlcm, lcm_vec[bc]);
                    AX[bc].scale(rlcm);
                    fmpz_poly_clear(lcm_vec[bc]);
                    AX[bc](TMat[br][bc],fp); // back to TMat
                    AX[bc].clear();
                }
                acb_poly_init(TD[br][0]);
                acb_poly_set_fmpz_poly(TD[br][0],zlcm,fp);
            }
            fmpz_poly_clear(zlcm);
            fmpz_poly_clear(rlcm);
            ntaylor_inited = true;
        }
        
        acb_t z0;
        acb_init(z0); 
        _to_(z0,x0,fp);
        I.resize(nbs);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            
            acb_poly_t lcm;
            acb_poly_init(lcm);
            acb_poly_set(lcm,TD[br][0]);
            acb_poly_taylor_shift(lcm,lcm,z0,fp); // shift: x -> x+x0
            
            // now we use a=br
            // n D0 Ia(n) = - sum_{0<m<n} (n-m) Dm Ia(n-m) + sum_{b<=a} sum_{0<=m<n} Amab.Ib(n-1-m)
            int a = br;
            acb_t q,D0;
            acb_init(q);
            acb_init(D0);
            acb_poly_get_coeff_acb(D0,lcm,0);
            if(acb_is_zero(D0)) throw Error("taylor: D0 is zero.");
            
            I[a] = vector<acb_mat_t>(xn+1);
            acb_mat_init(I[a][0],nr,nc);
            _to_(I[a][0],ex_to<matrix>(sub_matrix(I0, r0, nr, 0, nc)),fp);
            
            acb_mat_t smat;
            acb_mat_init(smat,nr,nc);
                
            for(int n=1; n<=xn; n++) { // 3-cycle over n
            
                if(omp_get_num_threads()==1 && !In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--ntaylor: " << nbs << "|" << br+1;
                    cout << " [" << nr << "\u2A09" << nc << "] n" << xn << "|" << n << flush;
                }
                
                acb_mat_zero(smat);
                
                if(true) {
                    slong s = acb_poly_degree(lcm);
                    if(s>n-1) s=n-1;
                    vector<acb_mat_t> mat_vec(s+1);
                    for(int m=1; m<=s; m++) acb_mat_init(mat_vec[m],nr,nc);
                    #pragma omp parallel for schedule(runtime)
                    for(int m=0; m<=s; m++) {
                        auto mat = mat_vec[m];
                        acb_t Dm;
                        acb_init(Dm);
                        acb_poly_get_coeff_acb(Dm,lcm,m);
                        acb_mul_si(Dm,Dm,n-m,fp);
                        acb_neg(Dm,Dm);
                        acb_mat_scalar_mul_acb(mat,I[a][n-m],Dm,fp);
                        acb_clear(Dm);
                        flint_cleanup();
                    }
                    for(int m=1; m<=s; m++) {
                        acb_mat_add(smat,smat,mat_vec[m],fp);
                        acb_mat_clear(mat_vec[m]);
                    }
                }
          
                if(true) {
                    vector<acb_mat_t> smat_vec(a+1);
                    for(int i=0; i<=a; i++) acb_mat_init(smat_vec[i],nr,nc);
                    #pragma omp parallel for schedule(runtime) if(a>omp_get_max_threads())
                    for(int b=0; b<=a; b++) {
                        int nc2 = bs[b].second; 
                        acb_poly_t TM[nr][nc2];
                        slong s = 0;
                        for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) { // shift
                            acb_poly_init(TM[r][c]);
                            acb_poly_set(TM[r][c],TMat[a][b][r][c]);
                            acb_poly_taylor_shift(TM[r][c],TM[r][c],z0,fp); // shift: x -> x+x0
                            int ss = acb_poly_degree(TM[r][c]);
                            if(ss>s) s = ss;
                        }
                        if(s>n-1) s = n-1;
                        
                        vector<acb_mat_t> mat_vec(s+1);
                        for(int i=0; i<=s; i++) acb_mat_init(mat_vec[i],nr,nc);
                        #pragma omp parallel for schedule(runtime)
                        for(int m=0; m<=s; m++) {
                            acb_mat_t Amab;
                            acb_mat_init(Amab,nr,nc2); // coefficients
                            for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) { 
                                acb_poly_get_coeff_acb(acb_mat_entry(Amab,r,c),TM[r][c],m);
                            }
                            acb_mat_mul(mat_vec[m],Amab,I[b][n-1-m],fp);
                            acb_mat_clear(Amab);
                            flint_cleanup();
                        }
                        acb_mat_zero(smat_vec[b]);
                        for(int i=0; i<=s; i++) {
                            acb_mat_add(smat_vec[b],smat_vec[b],mat_vec[i],fp);
                            acb_mat_clear(mat_vec[i]);
                        }
                        
                        for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) acb_poly_clear(TM[r][c]);
                        flint_cleanup();
                    }

                    for(int i=0; i<=a; i++) {
                        acb_mat_add(smat,smat,smat_vec[i],fp);
                        acb_mat_clear(smat_vec[i]);
                    }
                }

                acb_mul_si(q,D0,n,fp); // q = n D0
                acb_inv(q,q,fp); // q = 1 / (n D0)
                acb_mat_init(I[a][n],nr,nc);
                acb_mat_scalar_mul_acb(I[a][n],smat,q,fp);
            }
                
            acb_poly_clear(lcm);
            acb_clear(q);
            acb_clear(D0);
            acb_mat_clear(smat);
        }
        acb_clear(z0); 
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
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
        
}

