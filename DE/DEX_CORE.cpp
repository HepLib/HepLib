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
        
        inline void fmpq_poly_taylor_shift(fmpq_poly_t p1, fmpq_poly_t p2, fmpq_t q) {
            fmpq_poly_t p;
            fmpq_poly_init(p);
            fmpq_poly_set(p,p2); // p = p2
            fmpq_poly_rescale(p,p,q);
            fmpz_t z1;
            fmpz_init(z1);
            fmpz_one(z1);
            fmpz_poly_t np;
            fmpz_poly_init(np);
            fmpq_poly_get_numerator(np,p);
            fmpz_poly_taylor_shift(np,np,z1);
            fmpz_clear(z1);
            fmpq_poly_set_fmpz_poly(p1,np);
            auto & dp = fmpq_poly_denref(p);
            auto & dp1 = fmpq_poly_denref(p1);
            fmpz_mul(dp1,dp1,dp);
            fmpq_poly_canonicalise(p1);
            fmpq_poly_clear(p);
            fmpz_poly_clear(np);
            fmpq_t qi;
            fmpq_init(qi);
            fmpq_inv(qi,q);
            fmpq_poly_rescale(p1,p1,qi);
            fmpq_clear(qi);
            fmpq_poly_canonicalise(p1);
        }
        
    }

    //=*********************************************************************=
    
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
        
    //=*********************************************************************=
    // U-Series - rational
    //=*********************************************************************=
    
    // U[a][b][la][k][n], block_umat_fmpq_mat_t U; no need to initialize U
    void DEX::series(block_umat_fmpq_mat_t & U, int xn, const vector<fmpq_t> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        U.resize(nbs);
        for(int br=0; br<nbs; br++) U[br].resize(br+1);
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        fmpz_poly_set_str(x1, "2  0 1");
        fmpz_poly_set_str(x0, "1  1");
        
        int tot = 1;
        #pragma omp parallel for schedule(dynamic,1)
        for(int b=nbs-1; b>=0; b--) { // cycle columns
            if(!In_GiNaC_Parallel && Verbose>5) {
                #pragma omp critical
                { 
                    cout << "\r                                 \r" << flush;
                    cout << "  \\--series: x^" << xn << " U" << nbs << "|" << (tot++) << flush; 
                }
            }
            for(int a=b; a<nbs; a++) { // cycle rows
                // [(ila+n)D0-A0aa].Uab(ila,k,n) = -Uab(ila,k+1,n)D0
                //   - sum_{0<m<=n} [(ila+n)Uab(ila,k,n-m)+Uab(ila,k+1,n-m)]Dm 
                //   + sum_{b<=c<=a,0<=m<=n,NO(a=b,m=0)} Amac.Ucb(ila,k,n-m)
                
                int ab1 = a-b+1;
                int nr = bs[a].second;
                int nc = bs[b].second;
                
                vector<MQ> A(ab1); // M=A/x
                vector<int> sdeg(a+1); // A[bi].degree()
                fmpz_poly_t lcm, olcm, rlcm; // D=lcm
                fmpz_poly_init(lcm);
                fmpz_poly_init(olcm);
                fmpz_poly_init(rlcm);
                fmpz_poly_set(lcm, x0);
                vector<fmpz_poly_t> lcms(ab1);
                for(int c=b; c<=a; c++) {
                    int bi = c-b;
                    A[bi].init(Mat[a][c]);
                    A[bi].scale(x1); // A=x*M
                    fmpz_poly_init(lcms[bi]);
                    A[bi].denlcm(lcms[bi]);
                    fmpz_poly_set(olcm,lcm);
                    fmpz_poly_lcm(lcm, olcm, lcms[bi]);
                }
                for(int c=b; c<=a; c++) {
                    int bi = c-b;
                    fmpz_poly_div(rlcm, lcm, lcms[bi]);
                    A[bi].scale(rlcm);
                    fmpz_poly_clear(lcms[bi]);
                    sdeg[bi] = A[bi].degree();
                }
                fmpz_poly_clear(olcm);
                fmpz_poly_clear(rlcm);
                
                
                for(auto kv : UK[a][b]) { // insert kv first and make sure thread-safe
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[kv.first][0])) continue;
                    U[a][b][kv.first].resize(kv.second); 
                }
                int nla = UK[a][b].size();
                
                #pragma omp parallel for schedule(dynamic,1) if(b<4)
                for(int cla=0; cla<nla; cla++) { // 1-cycle over lambda
                
                    fmpq_mat_t A0aa,invA,smat,tmat;
                    fmpq_mat_init(A0aa,nr,nr);
                    fmpq_mat_init(invA,nr,nr);
                    fmpq_mat_init(smat,nr,nc);
                    fmpq_mat_init(tmat,nr,nc);
                    fmpz_t D0,Dm;
                    fmpz_init(D0);                
                    fmpz_init(Dm);
                    fmpq_t q,qm;
                    fmpq_init(q);
                    fmpq_init(qm);
                    
                    A[a-b].coeff(A0aa,0);
                    fmpz_poly_get_coeff_fmpz(D0,lcm,0);
                    
                    auto kv = UK[a][b].begin();
                    advance(kv,cla);
                    auto ila = kv->first;
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) continue;
                    auto kmax = kv->second;
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][ila][k] = vector<fmpq_mat_t>(xn+1);
                        fmpq_mat_init(U[a][b][ila][k][0],nr,nc);
                        fmpq_mat_set(U[a][b][ila][k][0],U0[a][b][ila][k][0]);

                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            fmpq_mat_zero(smat);
                            fmpq_add_si(q,qlas[ila][0],n); // q = la+n

                            if(true) { 
                                slong s = fmpz_poly_degree(lcm);
                                if(s>n) s=n;
                                for(int m=0; m<=s; m++) {
                                    if(m==0 && k+1<kmax) {
                                        fmpz_neg(Dm,D0);
                                        fmpq_mat_scalar_mul_fmpz(tmat,U[a][b][ila][k+1][n],Dm);
                                        fmpq_mat_add(smat,smat,tmat);
                                    } else if(m>0) {
                                        fmpq_sub_si(qm,q,m);
                                        fmpq_mat_scalar_mul_fmpq(tmat,U[a][b][ila][k][n-m],qm);
                                        if(k+1<kmax) fmpq_mat_add(tmat,tmat,U[a][b][ila][k+1][n-m]);
                                        fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                                        fmpz_neg(Dm,Dm);
                                        fmpq_mat_scalar_mul_fmpz(tmat,tmat,Dm);
                                        fmpq_mat_add(smat,smat,tmat);
                                    }
                                }
                            }

                            if(true) {
                                for(int c=b; c<=a; c++) {
                                    fmpq_mat_t Amac;
                                    int nc2 = bs[c].second;
                                    fmpq_mat_init(Amac,nr,nc2);
                                    int ci = c-b;
                                    auto kv = U[c][b].find(ila);
                                    if(kv!=U[c][b].end() && kv->second.size()>k) {
                                        slong s = sdeg[ci];
                                        if(s>n) s = n;
                                        for(int m=0; m<=s; m++) {
                                            if(c!=a || m!=0) {
                                                A[ci].coeff(Amac,m);
                                                fmpq_mat_mul(tmat,Amac,kv->second[k][n-m]);
                                                fmpq_mat_add(smat,smat,tmat);
                                            }
                                        }
                                    }
                                    fmpq_mat_clear(Amac);
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
                    fmpq_mat_clear(A0aa);
                    fmpq_mat_clear(invA);
                    fmpq_mat_clear(smat);
                    fmpq_mat_clear(tmat);
                    fmpz_clear(D0);
                    fmpz_clear(Dm);
                    fmpq_clear(q);
                    fmpq_clear(qm);
                    flint_cleanup();   
                }
                for(int c=b; c<=a; c++) A[c-b].clear();
                fmpz_poly_clear(lcm);
            }
            flint_cleanup();
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
    }
    
    //=*********************************************************************=
    // U-Series - acb
    //=*********************************************************************=
        
    // U[a][b][la][k][n], block_umat_acb_mat_t U; no need to initialize U
    void DEX::series(block_umat_acb_mat_t & U, int xn, slong dp, const vector<fmpq_t> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        U.resize(nbs);
        for(int br=0; br<nbs; br++) U[br].resize(br+1);
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        fmpz_poly_set_str(x1, "2  0 1");
        fmpz_poly_set_str(x0, "1  1");
        
        int tot = 1;
        #pragma omp parallel for schedule(dynamic,1)
        for(int b=nbs-1; b>=0; b--) { // cycle columns
            if(!In_GiNaC_Parallel && Verbose>5) {
                #pragma omp critical
                { 
                    cout << "\r                                 \r" << flush;
                    cout << "  \\--nseries: x^" << xn << " U" << nbs << "|" << (tot++) << flush; 
                }
            }
            for(int a=b; a<nbs; a++) { // cycle rows
                // [(ila+n)D0-A0aa].Uab(ila,k,n) = -Uab(ila,k+1,n)D0
                //   - sum_{0<m<=n} [(ila+n)Uab(ila,k,n-m)+Uab(ila,k+1,n-m)]Dm 
                //   + sum_{b<=c<=a,0<=m<=n,NO(a=b,m=0)} Amac.Ucb(ila,k,n-m)
                
                int ab1 = a-b+1;
                int nr = bs[a].second;
                int nc = bs[b].second;
                
                vector<MQ> A(ab1); // M=A/x
                vector<int> sdeg(a+1); // A[bi].degree();
                fmpz_poly_t lcm, olcm, rlcm; // D=lcm
                fmpz_poly_init(lcm);
                fmpz_poly_init(olcm);
                fmpz_poly_init(rlcm);
                fmpz_poly_set(lcm, x0);
                vector<fmpz_poly_t> lcms(ab1);
                for(int c=b; c<=a; c++) {
                    int bi = c-b;
                    A[bi].init(Mat[a][c]);
                    A[bi].scale(x1); // A=x*M
                    fmpz_poly_init(lcms[bi]);
                    A[bi].denlcm(lcms[bi]);
                    fmpz_poly_set(olcm,lcm);
                    fmpz_poly_lcm(lcm, olcm, lcms[bi]);
                }
                for(int c=b; c<=a; c++) {
                    int bi = c-b;
                    fmpz_poly_div(rlcm, lcm, lcms[bi]);
                    A[bi].scale(rlcm);
                    fmpz_poly_clear(lcms[bi]);
                    sdeg[bi] = A[bi].degree();
                }
                fmpz_poly_clear(olcm);
                fmpz_poly_clear(rlcm);
            
                for(auto kv : UK[a][b]) { // insert kv first and make sure thread-safe
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[kv.first][0])) continue;
                    U[a][b][kv.first].resize(kv.second); 
                }
                int nla = UK[a][b].size();
                
                #pragma omp parallel for schedule(runtime) if(b<4)
                for(int cla=0; cla<nla; cla++) { // 1-cycle over lambda
                
                    fmpq_mat_t A0aa,invA;
                    fmpq_mat_init(A0aa,nr,nr);
                    fmpq_mat_init(invA,nr,nr);
                    fmpz_t D0,Dm;
                    fmpz_init(D0);
                    fmpz_init(Dm);
                    fmpq_t q,qm;
                    fmpq_init(q);
                    fmpq_init(qm);
                
                    acb_mat_t invAf,smat,tmat;
                    acb_mat_init(invAf,nr,nr);
                    acb_mat_init(smat,nr,nc);
                    acb_mat_init(tmat,nr,nc);
                    
                    A[a-b].coeff(A0aa,0);
                    fmpz_poly_get_coeff_fmpz(D0,lcm,0);
                
                    auto kv = UK[a][b].begin();
                    advance(kv,cla);
                    auto ila = kv->first;
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) continue;
                    auto kmax = kv->second;
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][ila][k] = vector<acb_mat_t>(xn+1);
                        acb_mat_init(U[a][b][ila][k][0],nr,nc);
                        acb_mat_set_fmpq_mat(U[a][b][ila][k][0],U0[a][b][ila][k][0],fp);
                        
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            acb_mat_zero(smat);
                            fmpq_add_si(q,qlas[ila][0],n); // q = la+n

                            if(true) { 
                                slong s = fmpz_poly_degree(lcm);
                                if(s>n) s=n;
                                for(int m=0; m<=s; m++) {
                                    if(m==0 && k+1<kmax) {
                                        fmpz_neg(Dm,D0);
                                        acb_mat_scalar_mul_fmpz(tmat,U[a][b][ila][k+1][n],Dm,fp);
                                        acb_mat_add(smat,smat,tmat,fp);
                                    } else if(m>0) {
                                        fmpq_sub_si(qm,q,m);
                                        acb_t fqm; acb_init(fqm);
                                        acb_set_fmpq(fqm,qm,fp);
                                        acb_mat_scalar_mul_acb(tmat,U[a][b][ila][k][n-m],fqm,fp);
                                        acb_clear(fqm);
                                        if(k+1<kmax) acb_mat_add(tmat,tmat,U[a][b][ila][k+1][n-m],fp);
                                        fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                                        fmpz_neg(Dm,Dm);
                                        acb_mat_scalar_mul_fmpz(tmat,tmat,Dm,fp);
                                        acb_mat_add(smat,smat,tmat,fp);
                                    }
                                }
                            }

                            if(true) {
                                for(int c=b; c<=a; c++) {
                                    acb_mat_t Amac;
                                    int nc2 = bs[c].second;            
                                    acb_mat_init(Amac,nr,nc2);
                                    int ci = c-b;
                                    auto kv = U[c][b].find(ila);
                                    if(kv!=U[c][b].end() && kv->second.size()>k) {
                                        slong s = sdeg[ci];
                                        if(s>n) s = n;
                                        for(int m=0; m<=s; m++) {
                                            if(c!=a || m!=0) {
                                                A[ci].coeff(Amac,m,fp);
                                                acb_mat_mul(tmat,Amac,kv->second[k][n-m],fp);
                                                acb_mat_add(smat,smat,tmat,fp);
                                            }
                                        }
                                    }
                                    acb_mat_clear(Amac);
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
                    
                    fmpq_mat_clear(A0aa);
                    fmpq_mat_clear(invA);
                    fmpz_clear(D0);
                    fmpz_clear(Dm);
                    fmpq_clear(q);
                    fmpq_clear(qm);
                    acb_mat_clear(invAf);
                    acb_mat_clear(smat);
                    acb_mat_clear(tmat);
                    flint_cleanup();
                }
                for(int c=b; c<=a; c++) A[c-b].clear();
                fmpz_poly_clear(lcm);
            }
            flint_cleanup();
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
    }
    
    //=*********************************************************************=
    // I-Series - expansion - rational
    //=*********************************************************************=
        
    // I[a][la][k][n] & In0[a][ila][k], no need to initialize I
    void DEX::series(block_imat_fmpq_mat_t & I, int xn, block_imat_fmpq_mat_t & In0, int nc, const vector<fmpq_t> & qslas) { 
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        I.resize(nbs);
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        fmpz_poly_set_str(x1, "2  0 1");
        fmpz_poly_set_str(x0, "1  1");
        
        for(int a=0; a<nbs; a++) { // cycle rows
            // [(la+n)D0-A0aa].Ia(la,k,n) = -D0 Ia(la,k+1,n)
            //   - sum_{0<m<=n} Dm [(la+n-m).Ia(la,k,n-m)+Ia(la,k+1,n-m)] 
            //   + sum_{not(b=a|m=0)} Amab.Ib(la,k,n-m)
            
            int nr = bs[a].second;
            vector<MQ> A(a+1); // M=A/x
            vector<int> sdeg(a+1); // A[b].degree();
            fmpz_poly_t lcm, rlcm; // D=lcm
            vector<fmpz_poly_t> lcm_vec(a+1);
            fmpz_poly_init(lcm);
            fmpz_poly_init(rlcm);
            fmpz_poly_set(lcm,x0);
            for(int b=0; b<=a; b++) { 
                // lcm for each block A[b] in each row
                A[b].init(Mat[a][b]);
                A[b].scale(x1); // A=x*M
                fmpz_poly_init(lcm_vec[b]);
                A[b].denlcm(lcm_vec[b]);
                fmpz_poly_lcm(lcm, lcm, lcm_vec[b]);
            }
            for(int b=0; b<=a; b++) {
                fmpz_poly_div(rlcm, lcm, lcm_vec[b]);
                A[b].scale(rlcm);
                fmpz_poly_clear(lcm_vec[b]);
                sdeg[b] = A[b].degree();
            }
            fmpz_poly_clear(rlcm);
            
            fmpq_mat_t A0aa;
            fmpq_mat_init(A0aa,nr,nr);
            A[a].coeff(A0aa,0);
            fmpz_t D0;
            fmpz_init(D0);
            fmpz_poly_get_coeff_fmpz(D0,lcm,0);
            
            for(auto kv : IK[a]) { // insert kv first and make sure thread-safe
                if(qslas.size()>0 && !is_resonant(qslas,qlas[kv.first][0])) continue;
                I[a][kv.first].resize(kv.second); 
            }
            
            int nla = IK[a].size();
            bool la_parallel = (qslas.size()>0 ? qslas.size() : nla) >= omp_get_max_threads();
            if(la_parallel) { // to ultilize the thread_pool, we use following trick
                if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--series: " << nbs << "|" << a+1;
                    cout << " [" << nr << "\u2A09" << nc << "]";
                    cout << " \u03BB" << nla << " n" << xn << flush;
                }
                //-------------------------------------------------//
                // THE SAME AS BELOW, EXCEPT OpenMP 
                //-------------------------------------------------//
                #pragma omp parallel for schedule(runtime)
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
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) { // seleted lambda set
                        if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--series: " << nbs << "|" << a+1;
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << "|" << cla+1 << " n" << xn << flush;
                        }
                        continue;
                    }
                    auto kmax = kv->second;  
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        I[a][ila][k] = vector<fmpq_mat_t>(xn+1);
                        fmpq_mat_init(I[a][ila][k][0],nr,nc);
                        fmpq_mat_set(I[a][ila][k][0],In0[a][ila][k][0]);
                        
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                                cout << "\r                                                              \r" << flush;
                                cout << "  \\--series: " << nbs << "|" << a+1;
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
                                int tmax = omp_get_max_threads();
                                vector<fmpq_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_init(mat_sum[i],nr,nc);
                                    fmpq_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax)
                                for(int m=0; m<=s; m++) {
                                    auto tid = omp_get_thread_num();
                                    auto & mat = mat_tmp[tid];
                                    fmpz_t Dm;
                                    fmpz_init(Dm);
                                    fmpq_t qm;
                                    fmpq_init(qm);
                                    if(m==0 && k+1<kmax) {
                                        fmpz_neg(Dm,D0);
                                        fmpq_mat_scalar_mul_fmpz(mat,I[a][ila][k+1][n],Dm);
                                        fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat);
                                    } else if(m>0) {
                                        fmpq_sub_si(qm,q,m);
                                        fmpq_mat_scalar_mul_fmpq(mat,I[a][ila][k][n-m],qm);
                                        if(k+1<kmax) fmpq_mat_add(mat,mat,I[a][ila][k+1][n-m]);
                                        fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                                        fmpz_neg(Dm,Dm);
                                        fmpq_mat_scalar_mul_fmpz(mat,mat,Dm);
                                        fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat);
                                    }
                                    fmpz_clear(Dm);
                                    fmpq_clear(qm);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_add(smat,smat,mat_sum[i]);
                                    fmpq_mat_clear(mat_sum[i]);
                                    fmpq_mat_clear(mat_tmp[i]);
                                }
                            }
                            
                            if(true) {
                                int tmax = omp_get_max_threads();
                                vector<fmpq_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_init(mat_sum[i],nr,nc);
                                    fmpq_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax)
                                for(int b=0; b<=a; b++) {
                                    auto tid = omp_get_thread_num();
                                    int nc2 = bs[b].second;
                                    fmpq_mat_t Amab;
                                    fmpq_mat_init(Amab,nr,nc2);
                                    slong s = sdeg[b];
                                    if(s>n) s=n;
                                    for(int m=0; m<=s; m++) {
                                        auto kv = I[b].find(ila);
                                        if(kv!=I[b].end() && kv->second.size()>k) {
                                            if(b==a && m==0) continue;
                                            A[b].coeff(Amab,m);
                                            fmpq_mat_mul(mat_tmp[tid],Amab,kv->second[k][n-m]);
                                            fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat_tmp[tid]);
                                        }
                                    }
                                    fmpq_mat_clear(Amab);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_add(smat,smat,mat_sum[i]);
                                    fmpq_mat_clear(mat_sum[i]);
                                    fmpq_mat_clear(mat_tmp[i]);
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
                }
                //-------------------------------------------------//
            } else {
                //-------------------------------------------------//
                // THE SAME AS ABOVE, EXCEPT OpenMP 
                //-------------------------------------------------//
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
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) { // seleted lambda set
                        if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--series: " << nbs << "|" << a+1;
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << "|" << cla+1 << " n" << xn << flush;
                        }
                        continue;
                    }
                    auto kmax = kv->second;  
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        I[a][ila][k] = vector<fmpq_mat_t>(xn+1);
                        fmpq_mat_init(I[a][ila][k][0],nr,nc);
                        fmpq_mat_set(I[a][ila][k][0],In0[a][ila][k][0]);
                        
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                                cout << "\r                                                              \r" << flush;
                                cout << "  \\--series: " << nbs << "|" << a+1;
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
                                int tmax = omp_get_max_threads();
                                vector<fmpq_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_init(mat_sum[i],nr,nc);
                                    fmpq_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax)
                                for(int m=0; m<=s; m++) {
                                    auto tid = omp_get_thread_num();
                                    auto & mat = mat_tmp[tid];
                                    fmpz_t Dm;
                                    fmpz_init(Dm);
                                    fmpq_t qm;
                                    fmpq_init(qm);
                                    if(m==0 && k+1<kmax) {
                                        fmpz_neg(Dm,D0);
                                        fmpq_mat_scalar_mul_fmpz(mat,I[a][ila][k+1][n],Dm);
                                        fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat);
                                    } else if(m>0) {
                                        fmpq_sub_si(qm,q,m);
                                        fmpq_mat_scalar_mul_fmpq(mat,I[a][ila][k][n-m],qm);
                                        if(k+1<kmax) fmpq_mat_add(mat,mat,I[a][ila][k+1][n-m]);
                                        fmpz_poly_get_coeff_fmpz(Dm,lcm,m);
                                        fmpz_neg(Dm,Dm);
                                        fmpq_mat_scalar_mul_fmpz(mat,mat,Dm);
                                        fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat);
                                    }
                                    fmpz_clear(Dm);
                                    fmpq_clear(qm);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_add(smat,smat,mat_sum[i]);
                                    fmpq_mat_clear(mat_sum[i]);
                                    fmpq_mat_clear(mat_tmp[i]);
                                }
                            }
                            
                            if(true) {
                                int tmax = omp_get_max_threads();
                                vector<fmpq_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_init(mat_sum[i],nr,nc);
                                    fmpq_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax)
                                for(int b=0; b<=a; b++) {
                                    auto tid = omp_get_thread_num();
                                    int nc2 = bs[b].second;
                                    fmpq_mat_t Amab;
                                    fmpq_mat_init(Amab,nr,nc2);
                                    slong s = sdeg[b];
                                    if(s>n) s=n;
                                    for(int m=0; m<=s; m++) {
                                        auto kv = I[b].find(ila);
                                        if(kv!=I[b].end() && kv->second.size()>k) {
                                            if(b==a && m==0) continue;
                                            A[b].coeff(Amab,m);
                                            fmpq_mat_mul(mat_tmp[tid],Amab,kv->second[k][n-m]);
                                            fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat_tmp[tid]);
                                        }
                                    }
                                    fmpq_mat_clear(Amab);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    fmpq_mat_add(smat,smat,mat_sum[i]);
                                    fmpq_mat_clear(mat_sum[i]);
                                    fmpq_mat_clear(mat_tmp[i]);
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
                }
                //-------------------------------------------------//
            }
            for(int b=0; b<=a; b++) A[b].clear();
            fmpz_poly_clear(lcm);
            fmpq_mat_clear(A0aa);
            fmpz_clear(D0);
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
    }
    
    
    //=*********************************************************************=
    // I-Series - acb
    //=*********************************************************************=
    
    // I[a][la][k][n] & In0[a][la][k], no need to initialize I
    void DEX::series(block_imat_acb_mat_t & I, int xn, block_imat_acb_mat_t & In0, int nc, slong dp, const vector<fmpq_t> & qslas) {  
        if(!fuchsified) fuchsify();
        auto fp = dp2fp(dp);
        auto nbs = bs.size();
        I.resize(nbs);
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        fmpz_poly_set_str(x1, "2  0 1");
        fmpz_poly_set_str(x0, "1  1");
        
        for(int a=0; a<nbs; a++) { // cycle rows
            // [(la+n)D0-A0aa].Ia(la,k,n) = -D0 Ia(la,k+1,n)
            //   - sum_{0<m<=n} Dm [(la+n-m).Ia(la,k,n-m)+Ia(la,k+1,n-m)] 
            //   + sum_{not(b=a|m=0)} Amab.Ib(la,k,n-m)
            
            int nr = bs[a].second;
            vector<MQ> A(a+1); // M=A/x
            vector<int> sdeg(a+1); // A[a].degree()
            fmpz_poly_t zlcm, rlcm; // D=lcm
            vector<fmpz_poly_t> lcm_vec(a+1);
            fmpz_poly_init(zlcm);
            fmpz_poly_init(rlcm);
            fmpz_poly_set(zlcm, x0);
            for(int b=0; b<=a; b++) { 
                // lcm for each block A[b] in each row 
                A[b].init(Mat[a][b]);
                A[b].scale(x1); // A=x*M
                fmpz_poly_init(lcm_vec[b]);
                A[b].denlcm(lcm_vec[b]);
                fmpz_poly_lcm(zlcm, zlcm, lcm_vec[b]);
            }
            for(int b=0; b<=a; b++) {
                fmpz_poly_div(rlcm, zlcm, lcm_vec[b]);
                A[b].scale(rlcm);
                fmpz_poly_clear(lcm_vec[b]);
                sdeg[b] = A[b].degree();
            }
            acb_poly_t lcm;
            acb_poly_init(lcm);
            acb_poly_set_fmpz_poly(lcm,zlcm,fp);
            fmpz_poly_clear(zlcm);
            fmpz_poly_clear(rlcm);
            
            acb_mat_t A0aa;
            acb_mat_init(A0aa,nr,nr);
            A[a].coeff(A0aa,0,fp);
            acb_t D0;
            acb_init(D0);
            acb_poly_get_coeff_acb(D0,lcm,0);
            
            for(auto kv : IK[a]) { // insert kv first and make sure thread-safe
                if(qslas.size()>0 && !is_resonant(qslas,qlas[kv.first][0])) continue;
                I[a][kv.first].resize(kv.second); 
            }
            
            int nla = IK[a].size();
            bool la_parallel = (qslas.size()>0 ? qslas.size() : nla) >= omp_get_max_threads();
            if(la_parallel) { // to ultilize the thread_pool, we use following trick
                if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--nseries: " << nbs << "|" << a+1;
                    cout << " [" << nr << "\u2A09" << nc << "]";
                    cout << " \u03BB" << nla << " n" << xn << flush;
                }
                //-------------------------------------------------//
                // THE SAME AS BELOW, EXCEPT OpenMP 
                //-------------------------------------------------//
                #pragma omp parallel for schedule(runtime)
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
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) { // selected lambda set
                        if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--nseries: " << nbs << "|" << a+1;
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << "|" << cla+1 << " n" << xn << flush;
                        }
                        continue;
                    }
                    auto kmax = kv->second;
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        I[a][ila][k] = vector<acb_mat_t>(xn+1);
                        acb_mat_init(I[a][ila][k][0],nr,nc);
                        acb_mat_set(I[a][ila][k][0],In0[a][ila][k][0]);
                        
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                                cout << "\r                                                              \r" << flush;
                                cout << "  \\--nseries: " << nbs << "|" << a+1;
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
                                int tmax = omp_get_max_threads();
                                vector<acb_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_init(mat_sum[i],nr,nc);
                                    acb_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax)
                                for(int m=0; m<=s; m++) {
                                    auto tid = omp_get_thread_num();
                                    auto & mat = mat_tmp[tid];
                                    acb_t Dm;
                                    acb_init(Dm);
                                    acb_t qm;
                                    acb_init(qm);
                                    if(m==0 && k+1<kmax) {
                                        acb_neg(Dm,D0);
                                        acb_mat_scalar_mul_acb(mat,I[a][ila][k+1][n],Dm,fp);
                                        acb_mat_add(mat_sum[tid],mat_sum[tid],mat,fp);
                                    } else if(m>0) {
                                        acb_sub_si(qm,q,m,fp);
                                        acb_mat_scalar_mul_acb(mat,I[a][ila][k][n-m],qm,fp);
                                        if(k+1<kmax) acb_mat_add(mat,mat,I[a][ila][k+1][n-m],fp);
                                        acb_poly_get_coeff_acb(Dm,lcm,m);
                                        acb_neg(Dm,Dm);
                                        acb_mat_scalar_mul_acb(mat,mat,Dm,fp);
                                        acb_mat_add(mat_sum[tid],mat_sum[tid],mat,fp);
                                    }
                                    acb_clear(Dm);
                                    acb_clear(qm);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_add(smat,smat,mat_sum[i],fp);
                                    acb_mat_clear(mat_sum[i]);
                                    acb_mat_clear(mat_tmp[i]);
                                }
                            }

                            if(true) {
                                int tmax = omp_get_max_threads();
                                vector<acb_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_init(mat_sum[i],nr,nc);
                                    acb_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax) 
                                for(int b=0; b<=a; b++) {
                                    auto tid = omp_get_thread_num();
                                    int nc2 = bs[b].second;                                            
                                    acb_mat_t Amab;
                                    acb_mat_init(Amab,nr,nc2);
                                    slong s = sdeg[b];
                                    if(s>n) s=n;
                                    for(int m=0; m<=s; m++) {
                                        if(b==a && m==0) continue;
                                        auto kv = I[b].find(ila);
                                        if(kv!=I[b].end() && kv->second.size()>k) {
                                            A[b].coeff(Amab,m,fp);
                                            acb_mat_mul(mat_tmp[tid],Amab,kv->second[k][n-m],fp);
                                            acb_mat_add(mat_sum[tid],mat_sum[tid],mat_tmp[tid],fp);
                                        }
                                    }
                                    acb_mat_clear(Amab);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_add(smat,smat,mat_sum[i],fp);
                                    acb_mat_clear(mat_sum[i]);
                                    acb_mat_clear(mat_tmp[i]);
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
                }
                //-------------------------------------------------//
            } else {
                //-------------------------------------------------//
                // THE SAME AS ABOVE, EXCEPT OpenMP 
                //-------------------------------------------------//
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
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila][0])) { // selected lambda set
                        if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--nseries: " << nbs << "|" << a+1;
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << "|" << cla+1 << " n" << xn << flush;
                        }
                        continue;
                    }
                    auto kmax = kv->second;
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        I[a][ila][k] = vector<acb_mat_t>(xn+1);
                        acb_mat_init(I[a][ila][k][0],nr,nc);
                        acb_mat_set(I[a][ila][k][0],In0[a][ila][k][0]);
                        
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                                cout << "\r                                                              \r" << flush;
                                cout << "  \\--nseries: " << nbs << "|" << a+1;
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
                                int tmax = omp_get_max_threads();
                                vector<acb_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_init(mat_sum[i],nr,nc);
                                    acb_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax)
                                for(int m=0; m<=s; m++) {
                                    auto tid = omp_get_thread_num();
                                    auto & mat = mat_tmp[tid];
                                    acb_t Dm;
                                    acb_init(Dm);
                                    acb_t qm;
                                    acb_init(qm);
                                    if(m==0 && k+1<kmax) {
                                        acb_neg(Dm,D0);
                                        acb_mat_scalar_mul_acb(mat,I[a][ila][k+1][n],Dm,fp);
                                        acb_mat_add(mat_sum[tid],mat_sum[tid],mat,fp);
                                    } else if(m>0) {
                                        acb_sub_si(qm,q,m,fp);
                                        acb_mat_scalar_mul_acb(mat,I[a][ila][k][n-m],qm,fp);
                                        if(k+1<kmax) acb_mat_add(mat,mat,I[a][ila][k+1][n-m],fp);
                                        acb_poly_get_coeff_acb(Dm,lcm,m);
                                        acb_neg(Dm,Dm);
                                        acb_mat_scalar_mul_acb(mat,mat,Dm,fp);
                                        acb_mat_add(mat_sum[tid],mat_sum[tid],mat,fp);
                                    }
                                    acb_clear(Dm);
                                    acb_clear(qm);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_add(smat,smat,mat_sum[i],fp);
                                    acb_mat_clear(mat_sum[i]);
                                    acb_mat_clear(mat_tmp[i]);
                                }
                            }

                            if(true) {
                                int tmax = omp_get_max_threads();
                                vector<acb_mat_t> mat_sum(tmax), mat_tmp(tmax);
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_init(mat_sum[i],nr,nc);
                                    acb_mat_init(mat_tmp[i],nr,nc); 
                                }
                                #pragma omp parallel for schedule(runtime) num_threads(tmax) 
                                for(int b=0; b<=a; b++) {
                                    auto tid = omp_get_thread_num();
                                    int nc2 = bs[b].second;                                            
                                    acb_mat_t Amab;
                                    acb_mat_init(Amab,nr,nc2);
                                    slong s = sdeg[b];
                                    if(s>n) s=n;
                                    for(int m=0; m<=s; m++) {
                                        if(b==a && m==0) continue;
                                        auto kv = I[b].find(ila);
                                        if(kv!=I[b].end() && kv->second.size()>k) {
                                            A[b].coeff(Amab,m,fp);
                                            acb_mat_mul(mat_tmp[tid],Amab,kv->second[k][n-m],fp);
                                            acb_mat_add(mat_sum[tid],mat_sum[tid],mat_tmp[tid],fp);
                                        }
                                    }
                                    acb_mat_clear(Amab);
                                    flint_cleanup();
                                }
                                for(int i=0; i<tmax; i++) {
                                    acb_mat_add(smat,smat,mat_sum[i],fp);
                                    acb_mat_clear(mat_sum[i]);
                                    acb_mat_clear(mat_tmp[i]);
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
                }
                //-------------------------------------------------//
            }
            for(int b=0; b<=a; b++) A[b].clear();
            acb_poly_clear(lcm);
            acb_mat_clear(A0aa);
            acb_clear(D0);
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
    }
    
    //=*********************************************************************=
    // Taylor expansion - rational
    //=*********************************************************************=
    
    void DEX::taylor(vector<vector<fmpq_mat_t>> & I, int xn, const matrix I0, const ex & x0) {  
        if(ntaylor_inited) throw Error("ntaylor_inited = true");
        if(!x0.info(info_flags::rational)) throw Error("x0 is NOT rational.");
        auto nbs = bs.size();
        int nc = I0.cols();
        
        if(!taylor_inited) { // initialize TMat/TD, cache for later taylor call
            QMat.resize(nbs);
            QD.resize(nbs);
            fmpz_poly_t zlcm, rlcm; // D=lcm
            for(int br=0; br<nbs; br++) { // cycle rows
                QMat[br].resize(br+1);
                QD[br] = vector<fmpz_poly_t>(1);
                int nr = bs[br].second;
                vector<MQ> AX(br+1);
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
                    QMat[br][bc].resize(nr);
                    for(int r=0; r<nr; r++) QMat[br][bc][r] = vector<fmpz_poly_t>(nc2);
                    fmpz_poly_div(rlcm, zlcm, lcm_vec[bc]);
                    AX[bc].scale(rlcm);
                    fmpz_poly_clear(lcm_vec[bc]);
                    AX[bc](QMat[br][bc]); // back to TMat
                    AX[bc].clear();
                }
                fmpz_poly_init(QD[br][0]);
                fmpz_poly_set(QD[br][0],zlcm);
            }
            fmpz_poly_clear(zlcm);
            fmpz_poly_clear(rlcm);
            taylor_inited = true;
            
            if(taylor_clearq) { // clearq or not
                for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
                    for(auto & row : Mat[br][bc]) for(auto & item : row) fmpz_poly_q_clear(item);
                    if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) fmpq_mat_clear(item[0]);
                }
                Mat.clear();
                if(fuchsified) { 
                    U0.clear(); 
                    for(auto & item : qlas) fmpq_clear(item[0]);
                    qlas.clear();
                }
                fuchsified = false;
            }
        }
 
        fmpq_t q0;
        fmpq_init(q0); 
        _to_(q0,x0);  
        I.resize(nbs);
        for(int a=0; a<nbs; a++) { // cycle rows
            // n D0 Ia(n) = - sum_{0<m<n} (n-m) Dm Ia(n-m) + sum_{b<=a} sum_{0<=m<n} Amab.Ib(n-1-m)
            
            int r0 = bs[a].first;
            int nr = bs[a].second;
            
            // taylor shift: x -> x+x0
            vector<int> sdeg(a+1); 
            vector<vector<vector<fmpq_poly_t>>> TM(a+1); // TM[b][r][c] for each row a
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                TM[b].resize(nr);
                int s = -1;
                for(int r=0; r<nr; r++) {
                    TM[b][r] = vector<fmpq_poly_t>(nc);
                    for(int c=0; c<nc; c++) {
                        auto & item = TM[b][r][c];
                        fmpq_poly_init(item);
                        fmpq_poly_set_fmpz_poly(item,QMat[a][b][r][c]);
                        fmpq_poly_taylor_shift(item,item,q0); // shift: x -> x+x0
                        int ss = fmpq_poly_degree(item);
                        if(ss>s) s = ss;
                    }
                }
                sdeg[b] = s;
            }
            
            fmpq_poly_t lcm;
            fmpq_poly_init(lcm);
            fmpq_poly_set_fmpz_poly(lcm,QD[a][0]);
            fmpq_poly_taylor_shift(lcm,lcm,q0); // shift: x -> x+x0
            
            fmpq_t q,D0;
            fmpq_init(q);
            fmpq_init(D0);
            fmpq_poly_get_coeff_fmpq(D0,lcm,0);
            if(fmpq_is_zero(D0)) throw Error("taylor: D0 is zero.");
            
            I[a] = vector<fmpq_mat_t>(xn+1);
            fmpq_mat_init(I[a][0],nr,nc);
            _to_(I[a][0],ex_to<matrix>(sub_matrix(I0, r0, nr, 0, nc)));
            
            fmpq_mat_t smat;
            fmpq_mat_init(smat,nr,nc);
                
            for(int n=1; n<=xn; n++) { // 3-cycle over n
            
                if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--taylor: " << nbs << "|" << a+1;
                    cout << " [" << nr << "\u2A09" << nc << "] n" << xn << "|" << n << flush;
                }
                
                fmpq_mat_zero(smat);
                
                if(true) {
                    slong s = fmpq_poly_degree(lcm);
                    if(s>n-1) s=n-1;
                    int tmax = omp_get_max_threads();
                    vector<fmpq_mat_t> mat_sum(tmax), mat_tmp(tmax);
                    for(int i=0; i<tmax; i++) {
                        fmpq_mat_init(mat_sum[i],nr,nc);
                        fmpq_mat_init(mat_tmp[i],nr,nc); 
                    }
                    #pragma omp parallel for schedule(runtime) num_threads(tmax)
                    for(int m=0; m<=s; m++) {
                        auto tid = omp_get_thread_num();
                        auto & mat = mat_tmp[tid];
                        fmpq_t Dm;
                        fmpq_init(Dm);
                        fmpq_poly_get_coeff_fmpq(Dm,lcm,m);
                        fmpq_mul_si(Dm,Dm,n-m);
                        fmpq_neg(Dm,Dm);
                        fmpq_mat_scalar_mul_fmpq(mat,I[a][n-m],Dm);
                        fmpq_clear(Dm);
                        fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat);
                        flint_cleanup();
                    }
                    for(int i=0; i<tmax; i++) {
                        fmpq_mat_add(smat,smat,mat_sum[i]);
                        fmpq_mat_clear(mat_sum[i]);
                        fmpq_mat_clear(mat_tmp[i]);
                    }
                }
          
                if(true) {
                    int tmax = omp_get_max_threads();
                    vector<fmpq_mat_t> mat_sum(tmax), mat_tmp(tmax);
                    for(int i=0; i<tmax; i++) {
                        fmpq_mat_init(mat_sum[i],nr,nc);
                        fmpq_mat_init(mat_tmp[i],nr,nc); 
                    }
                    #pragma omp parallel for schedule(runtime) num_threads(tmax)
                    for(int b=0; b<=a; b++) {
                        slong s = sdeg[b];
                        if(s>n-1) s=n-1;
                        auto tid = omp_get_thread_num();
                        int nc2 = bs[b].second; 
                        fmpq_mat_t Amab;
                        fmpq_mat_init(Amab,nr,nc2); 
                        for(int m=0; m<=s; m++) {
                            for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) { // coefficients
                                fmpq_poly_get_coeff_fmpq(fmpq_mat_entry(Amab,r,c),TM[b][r][c],m);
                            }
                            fmpq_mat_mul(mat_tmp[tid],Amab,I[b][n-1-m]);
                            fmpq_mat_add(mat_sum[tid],mat_sum[tid],mat_tmp[tid]);
                        }
                        fmpq_mat_clear(Amab);
                        flint_cleanup();
                    }
                    for(int i=0; i<tmax; i++) {
                        fmpq_mat_add(smat,smat,mat_sum[i]);
                        fmpq_mat_clear(mat_sum[i]);
                        fmpq_mat_clear(mat_tmp[i]); 
                    }
                }
                
                fmpq_mul_si(q,D0,n); // q = n D0
                fmpq_inv(q,q); // q = 1 / (n D0)
                fmpq_mat_init(I[a][n],nr,nc);
                fmpq_mat_scalar_mul_fmpq(I[a][n],smat,q);
            }
                
            fmpq_poly_clear(lcm);
            fmpq_clear(q);
            fmpq_clear(D0);
            fmpq_mat_clear(smat);
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) fmpq_poly_clear(TM[b][r][c]);
            }
        }
        fmpq_clear(q0); 
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
    }
    
    //=*********************************************************************=
    // Taylor expansion - acb
    //=*********************************************************************=
    
    void DEX::taylor(vector<vector<acb_mat_t>> & I, int xn, const matrix I0, const ex & x0, slong dp) {  
        if(taylor_inited) throw Error("taylor_inited = true");
        auto nbs = bs.size();
        auto fp = dp2fp(dp);
        int nc = I0.cols();
        
        if(!ntaylor_inited) { // initialize TMat/TD, cache for later taylor call
            TMat.resize(nbs);
            TD.resize(nbs);
            fmpz_poly_t zlcm, rlcm; // D=lcm
            for(int br=0; br<nbs; br++) { // cycle rows
                TMat[br].resize(br+1);
                TD[br] = vector<acb_poly_t>(1);
                int nr = bs[br].second;
                vector<MQ> AX(br+1);
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
            
            if(taylor_clearq) { // clearq or not
                for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
                    for(auto & row : Mat[br][bc]) for(auto & item : row) fmpz_poly_q_clear(item);
                    if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) fmpq_mat_clear(item[0]);
                }
                Mat.clear();
                if(fuchsified) { 
                    U0.clear(); 
                    for(auto & item : qlas) fmpq_clear(item[0]);
                    qlas.clear();
                }
                fuchsified = false;
            }
        }
        
        acb_t z0;
        acb_init(z0); 
        _to_(z0,x0,fp);        
        I.resize(nbs);
        for(int a=0; a<nbs; a++) { // cycle rows
            // n D0 Ia(n) = - sum_{0<m<n} (n-m) Dm Ia(n-m) + sum_{b<=a} sum_{0<=m<n} Amab.Ib(n-1-m)
            
            int r0 = bs[a].first;
            int nr = bs[a].second;
            
            // taylor shift: x -> x+x0
            vector<int> sdeg(a+1);
            vector<vector<vector<acb_poly_t>>> TM(a+1); // TM[b][r][c] for each row a
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                TM[b].resize(nr);
                int s = -1;
                for(int r=0; r<nr; r++) {
                    TM[b][r] = vector<acb_poly_t>(nc);
                    for(int c=0; c<nc; c++) {
                        auto & item = TM[b][r][c];
                        acb_poly_init(item);
                        acb_poly_set(item,TMat[a][b][r][c]);
                        acb_poly_taylor_shift(item,item,z0,fp); // shift: x -> x+x0
                        int ss = acb_poly_degree(item);
                        if(ss>s) s = ss;
                    }
                }
                sdeg[b] = s;
            }
            
            acb_poly_t lcm;
            acb_poly_init(lcm);
            acb_poly_set(lcm,TD[a][0]);
            acb_poly_taylor_shift(lcm,lcm,z0,fp); // shift: x -> x+x0
            
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
            
                if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                              \r" << flush;
                    cout << "  \\--ntaylor: " << nbs << "|" << a+1;
                    cout << " [" << nr << "\u2A09" << nc << "] n" << xn << "|" << n << flush;
                }
                
                acb_mat_zero(smat);
                
                if(true) {
                    slong s = acb_poly_degree(lcm);
                    if(s>n-1) s=n-1;
                    int tmax = omp_get_max_threads();
                    vector<acb_mat_t> mat_sum(tmax), mat_tmp(tmax);
                    for(int i=0; i<tmax; i++) {
                        acb_mat_init(mat_sum[i],nr,nc);
                        acb_mat_init(mat_tmp[i],nr,nc); 
                    }
                    #pragma omp parallel for schedule(runtime) num_threads(tmax)
                    for(int m=0; m<=s; m++) {
                        auto tid = omp_get_thread_num();
                        auto & mat = mat_tmp[tid];
                        acb_t Dm;
                        acb_init(Dm);
                        acb_poly_get_coeff_acb(Dm,lcm,m);
                        acb_mul_si(Dm,Dm,n-m,fp);
                        acb_neg(Dm,Dm);
                        acb_mat_scalar_mul_acb(mat,I[a][n-m],Dm,fp);
                        acb_clear(Dm);
                        acb_mat_add(mat_sum[tid],mat_sum[tid],mat,fp);
                        flint_cleanup();
                    }
                    for(int i=0; i<tmax; i++) {
                        acb_mat_add(smat,smat,mat_sum[i],fp);
                        acb_mat_clear(mat_sum[i]);
                        acb_mat_clear(mat_tmp[i]);
                    }
                }
          
                if(true) {
                    int tmax = omp_get_max_threads();
                    vector<acb_mat_t> mat_sum(tmax), mat_tmp(tmax);
                    for(int i=0; i<tmax; i++) {
                        acb_mat_init(mat_sum[i],nr,nc);
                        acb_mat_init(mat_tmp[i],nr,nc); 
                    }
                    #pragma omp parallel for schedule(runtime) num_threads(tmax)
                    for(int b=0; b<=a; b++) {
                        slong s = sdeg[b];
                        if(s>n-1) s=n-1;
                        auto tid = omp_get_thread_num();
                        int nc2 = bs[b].second; 
                        acb_mat_t Amab;
                        acb_mat_init(Amab,nr,nc2); 
                        for(int m=0; m<=s; m++) { 
                            for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) { // coefficients
                                acb_poly_get_coeff_acb(acb_mat_entry(Amab,r,c),TM[b][r][c],m);
                            } 
                            acb_mat_mul(mat_tmp[tid],Amab,I[b][n-1-m],fp);
                            acb_mat_add(mat_sum[tid],mat_sum[tid],mat_tmp[tid],fp);
                        }
                        acb_mat_clear(Amab);
                        flint_cleanup();
                    }
                    for(int i=0; i<tmax; i++) {
                        acb_mat_add(smat,smat,mat_sum[i],fp);
                        acb_mat_clear(mat_sum[i]);
                        acb_mat_clear(mat_tmp[i]); 
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
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) acb_poly_clear(TM[b][r][c]);
            }
        }
        acb_clear(z0); 
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
    }
    
    //=*********************************************************************=
        
}

