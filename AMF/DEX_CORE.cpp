/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "AMF.h"
#include "vspace.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>

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
        
        inline bool is_resonant(const vector<fmpq*> & qs, fmpq_t q) {
            fmpq_t dq;
            fmpq_init(dq);
            bool res = false;
            for(auto & qi : qs) {
                fmpq_sub(dq, qi, q);
                if(fmpz_is_pm1(fmpq_denref(dq))) { // check dq is integer or not
                    res = true;
                    break;
                }
            }
            fmpq_clear(dq);
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
        matrix pmat(n,n), pmat_inv(n,n);
        for(int r=0; r<n; r++) {
            pmat(vs[r],r) = 1;
            pmat_inv(r,vs[r]) = 1;
        }
        auto mat = pmat_inv.mul(m).mul(pmat);
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
                    Mat[br][bc][r].resize(nc);
                    for(int c=0; c<nc; c++) {
                        Mat[br][bc][r][c] = (fmpz_poly_q_struct*) flint_malloc(sizeof(fmpz_poly_q_struct));
                        fmpz_poly_q_init(Mat[br][bc][r][c]); // call fmpz_poly_q_clear() in ~DEX
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
        exset evs;
        evs.insert(0); // since a0 may be zero
        ev_am_t ev_am_blk;
        for(int bi=0; bi<nbs; bi++) { // diagonal blocks
            if(!In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                    \r" << flush;
                cout << pre << "\\--reducing diagonal blocks: " << nbs << "|" << bi+1 << flush;
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
                ti = exnormal(ti);
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
                    bool any_found = false;
                    ex diff = 0;
                    for(auto evi : evs) {
                        diff = exnormal(kv.first-evi);
                        if(diff.info(info_flags::integer)) {
                            any_found = true;
                            break;
                        }
                    }
                    if(!any_found) {
                        evs.insert(kv.first);
                        diff = 0;
                    }
                    ex fx = 1;
                    if(sm<=0 && diff<0) { sm=-1; fx = 1/x; shearing = true; }
                    else if(sm>=0 && diff>0) { sm=1; fx = x; shearing = true; }
                    for(int j=0; j<kv.second; j++) smat(cpos+j, cpos+j) = fx;
                    cpos += kv.second;
                }
                if(shearing) {
                    tm = tm.mul(smat);
                    mx.transform(tm, tm.inverse(solve_algo::gauss));
                    ti = ti.mul(tm);
                } else break;
            }
            auto ev_am_i = ev_am(mx.a0());
            for(auto kv : ev_am_i) ev_am_blk[kv.first] += kv.second;
            mx(Mat[bi][bi]); // back to Mat[bi][bi]
            mx.clear();
            for(int r=0; r<n; r++) for(int c=0; c<n; c++) t(n0+r,n0+c) = ti(r,c);
            mxt[bi].init(ti);
            mxti[bi].init(ti.inverse(solve_algo::gauss));
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int br=0; br<nbs; br++) for(int bc=0; bc<nbs-1; bc++) {
            if(bc>=br) continue;
            MX mx(Mat[br][bc]);
            mx.mul_left(mxti[br]).mul(mxt[bc])(Mat[br][bc]);
            mx.clear();
        }
        t = exnormal(t);
        Ts.push_back(t);
        for(int bi=0; bi<nbs; bi++) { mxt[bi].clear(); mxti[bi].clear(); }
        
        for(int br=0; br<nbs; br++) { // off-diagonal blocks
            vector<matrix> ts_vec;
            for(int bc=br-1; bc>=0; bc--) {
                auto nr=bs[br].second, nc=bs[bc].second;
                if(!In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                    \r" << flush;
                    cout << pre << "\\--reducing off-diagonal blocks: " << nbs << "|" << br+1 << "|" << (br-bc+1) << flush; 
                } 
                matrix sdmat = ex_to<matrix>(symbolic_matrix(nr,nc,"xyz"));
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
                    auto sol = lsolve(eqs,xs);
                    if(sol.nops()<1) throw Error("No solution found!");
                    auto dmat = ex_to<matrix>(subs(sdmat,sol)).mul_scalar(pow(x,-pr));
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
                    auto tt = ex_to<matrix>(unit_matrix(N,N));
                    for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                        tt(bs[br].first+r, bs[bc].first+c) = trc(r,c);
                    }
                    ts_vec.push_back(tt);
                    mrc(Mat[br][bc]); // back to Mat[br][bc]
                    MX trc_mx(trc);
                    for(int b=0; b<bc; b++) {
                        MX m1(Mat[br][b]), m2(Mat[bc][b]);
                        m2.mul_left(trc_mx);
                        m1.sub(m2)(Mat[br][b]); // back to Mat[br][b]
                        m2.clear();
                    }
                    for(int a=br+1; a<nbs; a++) {
                        MX m1(Mat[a][bc]), m2(Mat[a][br]);
                        m2.mul(trc_mx);
                        m1.add(m2)(Mat[a][bc]); // back to Mat[a][bc]
                        m1.clear();
                        m2.clear();
                    }
                }
                mrc.clear();
                mrr.clear();
                mcc.clear();
            }
            if(ts_vec.size()>0) { // mul ts_vec
                int n = ts_vec.size();
                int n2 = n/2;
                vector<matrix> ts2(n2);
                for(int i=0; i<n2; i++) ts2[i] = ts_vec[2*i].mul(ts_vec[2*i+1]);
                matrix t = ex_to<matrix>(unit_matrix(N));
                if(n%2==1) t = ts_vec[n-1];
                for(int i=n2-1; i>=0; i--) t = ts2[i].mul(t);
                if(!unit_matrix(N,N).is_equal(t)) Ts.push_back(t);
            }
        }
        if(nbs>1 && !In_GiNaC_Parallel && Verbose>5) cout << endl;
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--initializing x^Ao" << flush;

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

        auto qj = jordan(a0, ev_am_blk);
        
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
        m = qj.first.mul(m).mul(qj.first.inverse(solve_algo::gauss));
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
                            U0[br][bc][ila][k].resize(1);
                            U0[br][bc][ila][k][0] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                            fmpq_mat_init(U0[br][bc][ila][k][0],nr,nc); // call fmpq_mat_clear() in ~DEX
                            matrix mat(nr,nc);
                            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                                mat(r,c) = expand_ex(m(r0+r,c0+c),lst{sla,lnx}).coeff(sla).coeff(lnx,k);
                            }
                            _to_(U0[br][bc][ila][k][0], mat);
                        }
                    }
                }
            }
            for(auto kv : UK[br][0]) IK[br][kv.first] = kv.second;
        }
        nlas = las.size();
        qlas.resize(nlas);
        for(int i=0; i<nlas; i++) {
            qlas[i] = (fmpq*) flint_malloc(sizeof(fmpq));
            fmpq_init(qlas[i]); // call fmpq_clear() in ~DEX
            _to_(qlas[i], las[i]);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
        fuchsified = true;
        
    }
        
    //=*********************************************************************=
    // U-Series - rational
    //=*********************************************************************=
    
    // U[a][b][la][k][n], abikn_fmpq_mat_t U; no need to initialize U
    void DEX::series(abikn_fmpq_mat_t & U, int xn, const vector<fmpq*> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        U.resize(nbs);
        for(int br=0; br<nbs; br++) U[br].resize(br+1);
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        if(fmpz_poly_set_str(x1, "2  0 1")) throw Error("fmpz_poly_set_str failed.");
        fmpz_poly_one(x0);
        
        int tot = 1;
        #pragma omp parallel for schedule(dynamic,1)
        for(int b=nbs-1; b>=0; b--) { // cycle columns
            if(!In_GiNaC_Parallel && Verbose>5) {
                #pragma omp critical
                { 
                    cout << "\r                                 \r" << flush;
                    cout << pre << "\\--series: x^" << xn << " U" << nbs << "|" << (tot++) << flush;
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
                fmpz_poly_t lcm; // D=lcm
                fmpz_poly_init(lcm);
                fmpz_poly_set(lcm, x0);
                vector<fmpz_poly_struct*> lcms(ab1);
                for(int c=b; c<=a; c++) {
                    int bi = c-b;
                    A[bi].init(Mat[a][c]);
                    A[bi].scale(x1); // A=x*M
                    lcms[bi] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                    fmpz_poly_init(lcms[bi]);
                    A[bi].denlcm(lcms[bi]);
                    fmpz_poly_lcm(lcm, lcm, lcms[bi]);
                }
                for(int c=b; c<=a; c++) {
                    int bi = c-b;
                    fmpz_poly_div(lcms[bi], lcm, lcms[bi]);
                    A[bi].scale(lcms[bi]);
                    fmpz_poly_clear(lcms[bi]);
                    flint_free(lcms[bi]);
                    sdeg[bi] = A[bi].degree();
                }
                
                for(auto kv : UK[a][b]) { // insert kv first and make sure thread-safe
                    if(qslas.size()>0 && !is_resonant(qslas, qlas[kv.first])) continue;
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
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue;
                    auto kmax = kv->second;
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][ila][k].resize(xn+1);
                        U[a][b][ila][k][0] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(U[a][b][ila][k][0],nr,nc);
                        fmpq_mat_set(U[a][b][ila][k][0],U0[a][b][ila][k][0]);

                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            fmpq_mat_zero(smat);
                            fmpq_add_si(q,qlas[ila],n); // q = la+n

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
                            int nok = fmpq_mat_inv(invA,invA);
                            if(!nok) {
                                cout << "fmpq_mat_inv failed." << endl;
                                abort();
                            }
                            
                            U[a][b][ila][k][n] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
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
                    //flint_cleanup();
                }
                for(int c=b; c<=a; c++) A[c-b].clear();
                fmpz_poly_clear(lcm);
            }
            //flint_cleanup();
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
    }
    
    struct u_item_type {
        int a, b, ila, k, n, kmax;
        bool to_notify = false;
        list<int> ms;
        list<pair<int,int>> cms;
        u_item_type(int a_, int b_, int ila_, int k_, int n_, int kmax_, bool to_notify_, list<int> && ms_, list<pair<int,int>> && cms_) : a(a_), b(b_), ila(ila_), k(k_), n(n_), kmax(kmax_), to_notify(to_notify_), ms(std::move(ms_)), cms(std::move(cms_)) { }
        u_item_type() { }
    };
    
    // U[a][b][la][k][n], Parallel fmpq version
    void DEX::all_series(abikn_fmpq_mat_t & U, int xn, const vector<fmpq*> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        
        // Parallel Data
        vector<vector<vector<MQ>>> A(nbs); // M=A/x for A[a][b]
        vector<vector<fmpz_poly_struct*>> D(nbs); // Denominator, i.e., lcm
        vector<vector<map<int,vector<vector<bool>>>>> Done(nbs); // [a][b][ila][k][n] : ready or not
        vector<vector<map<int,vector<vector<bool>>>>> Used(nbs); // [a][b][ila][k][n] : used or not
        vector<vector<map<int,vector<vector<list<int>>>>>> S1(nbs); // [a][b][ila][k][n] : { m, ... }
        vector<vector<map<int,vector<vector<list<pair<int,int>>>>>>> S2(nbs); // [a][b][ila][k][n] -> { <c,m>, ... }
        unsigned long long worker_left = 0;
                
        U.resize(nbs);
        for(int a=0; a<nbs; a++) {
            A[a].resize(a+1);
            D[a].resize(a+1);
            U[a].resize(a+1);
            Done[a].resize(a+1);
            Used[a].resize(a+1);
            S1[a].resize(a+1);
            S2[a].resize(a+1);
            for(int b=0; b<=a; b++) {
                int nr = bs[a].second;
                int nc = bs[b].second;
                for(auto const & kv : UK[a][b]) { // 1-cycle over lambda
                    auto ila = kv.first;
                    auto kmax = kv.second;
                    if(qslas.size()>0 && !is_resonant(qslas, qlas[ila])) continue; // selected lambda set
                    U[a][b][ila].resize(kmax);
                    Done[a][b][ila].resize(kmax);
                    Used[a][b][ila].resize(kmax);
                    S1[a][b][ila].resize(kmax);
                    S2[a][b][ila].resize(kmax);
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][ila][k].resize(xn+1);
                        U[a][b][ila][k][0] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(U[a][b][ila][k][0],nr,nc); // call fmpq_mat_clear() after this call
                        fmpq_mat_set(U[a][b][ila][k][0],U0[a][b][ila][k][0]);
                    }
                }
            }
        }
        if(xn==0) {
            if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--Series U(x^" << xn << ") @ " << now(false) << endl;
            return;
        }
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--initializing Series(U)" << flush;
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        if(fmpz_poly_set_str(x1, "2  0 1")) throw Error("fmpz_poly_set_str failed!");
        fmpz_poly_one(x0);
        #pragma omp parallel for schedule(dynamic,1) collapse(2) reduction(+: worker_left)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b>a) continue;
            int ab1 = a-b+1;
            int nr = bs[a].second;
            int nc = bs[b].second;
            A[a][b].resize(ab1); // M=A/x
            vector<int> sdeg(a+1); // A[a][b][bi].degree()
            D[a][b] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(D[a][b]);
            fmpz_poly_set(D[a][b], x0);
            vector<fmpz_poly_struct*> lcms(ab1);
            for(int c=b; c<=a; c++) {
                int ci = c-b;
                A[a][b][ci].init(Mat[a][c]);
                A[a][b][ci].scale(x1); // A=x*M
                lcms[ci] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                fmpz_poly_init(lcms[ci]);
                A[a][b][ci].denlcm(lcms[ci]);
                fmpz_poly_lcm(D[a][b], D[a][b], lcms[ci]);
            }
            for(int c=b; c<=a; c++) {
                int ci = c-b;
                fmpz_poly_div(lcms[ci], D[a][b], lcms[ci]);
                A[a][b][ci].scale(lcms[ci]);
                fmpz_poly_clear(lcms[ci]);
                flint_free(lcms[ci]);
                sdeg[ci] = A[a][b][ci].degree();
            }
            
            slong s0 = fmpz_poly_degree(D[a][b]);
            for(auto const & kv : UK[a][b]) { // 1-cycle over lambda
                auto ila = kv.first;
                auto kmax = kv.second;
                if(qslas.size()>0 && !is_resonant(qslas, qlas[ila])) continue; // selected lambda set
                
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    
                    Done[a][b][ila][k].resize(xn+1);
                    Done[a][b][ila][k][0] = true;
                    Used[a][b][ila][k].resize(xn+1);
                    Used[a][b][ila][k][0] = false;
                    S1[a][b][ila][k].resize(xn+1);
                    S2[a][b][ila][k].resize(xn+1);
                    
                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                        Done[a][b][ila][k][n] = true;
                        Used[a][b][ila][k][n] = false;
                        U[a][b][ila][k][n] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(U[a][b][ila][k][n],nr,nc); // call fmpq_mat_clear() after this call
                        fmpq_mat_zero(U[a][b][ila][k][n]);
                        
                        slong s = s0;
                        if(s>n) s=n;
                        for(int m=0; m<=s; m++) {
                            if((m==0 && k+1<kmax) || (m>0)) {
                                S1[a][b][ila][k][n].push_back(m);
                                worker_left++;
                            }
                        }

                        for(int c=b; c<=a; c++) {
                            int nc2 = bs[c].second;
                            slong ss = sdeg[c-b];
                            if(ss>n) ss=n;
                            for(int m=0; m<=ss; m++) {
                                if(c==a && m==0) continue;
                                auto itr = U[c][b].find(ila);
                                if(itr!=U[c][b].end() && itr->second.size()>k) {
                                    S2[a][b][ila][k][n].push_back(make_pair(c,m));
                                    worker_left++;
                                }
                            }
                        }
                        
                        if(!S1[a][b][ila][k][n].empty() || !S2[a][b][ila][k][n].empty()) Done[a][b][ila][k][n] = false;
                    }
                }
            }
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
                
        mutex worker_mutex;
        condition_variable worker_cond;
        
        auto worker_total = worker_left;
        list<u_item_type> worker_item_list;
        auto worker_main = [&]() {
            while(true) {
                list<u_item_type> wi_list;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(!worker_left) return true;
                        start:
                        int cnt = 0;
                        auto tot = worker_item_list.size();
                        while(worker_item_list.size()>0) {
                            auto itr = worker_item_list.begin();
                            wi_list.emplace_back(std::move(*itr));
                            worker_item_list.erase(itr);
                            cnt++;
                            if(cnt>10) break;
                        }
                        if(cnt) return true;
                        for(int a=0; a<nbs; a++) for(int b=0; b<=a; b++) { // cycle rows & cols
                            bool ab_Done = true;
                            for(auto const & kv : UK[a][b]) { // 1-cycle over lambda
                                int ila = kv.first;
                                int kmax = kv.second;
                                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                                    auto & SS1 = S1[a][b][ila][k];
                                    auto & SS2 = S2[a][b][ila][k];
                                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                                        if(ab_Done && !Done[a][b][ila][k][n]) ab_Done = false;
                                        if(Used[a][b][ila][k][n] || Done[a][b][ila][k][n]) continue;
                                        
                                        list<int> ms;
                                        list<pair<int,int>> cms;
                                        
                                        list<int> & s1 = SS1[n];
                                        for(auto itr = s1.begin(); itr != s1.end(); ) {
                                            int m = *itr;
                                            if(m==0) {
                                                if(Done[a][b][ila][k+1][n]) {
                                                    itr = s1.erase(itr);
                                                    ms.push_back(m);
                                                } else ++itr;
                                            } else {
                                                if(Done[a][b][ila][k][n-m] && (k+1>=kmax || Done[a][b][ila][k+1][n-m])) {
                                                    itr = s1.erase(itr);
                                                    ms.push_back(m);
                                                } else ++itr;
                                            }
                                        }
                                        
                                        list<pair<int,int>> & s2 = SS2[n];
                                        for(auto itr = s2.begin(); itr != s2.end(); ) {
                                            int c = itr->first;
                                            int m = itr->second;
                                            if(Done[c][b][ila][k][n-m]) {
                                                itr = s2.erase(itr);
                                                cms.push_back(make_pair(c,m));
                                            } else ++itr;
                                        }
                                        
                                        if(!ms.empty() || !cms.empty()) {
                                            Used[a][b][ila][k][n] = true;
                                            bool to_notify = false;
                                            if(s1.empty() && s2.empty()) to_notify = true;
                                            worker_item_list.emplace_back(a,b,ila,k,n,kmax,to_notify,std::move(ms),std::move(cms));
                                        }
                                    }
                                }
                            }
                            if(ab_Done) {
                                for(int c=b; c<=a; c++) A[a][b][c-b].clear();
                            }
                        }
                        if(worker_item_list.size()==0) return false;
                        goto start; // for return to check the start label
                    });
                    if(!worker_left) {
                        flint_cleanup();
                        return;
                    }
                }
                
                unsigned long long worker_nn = 0;
                for(auto const & worker_item : wi_list) {
                    int a = worker_item.a;
                    int b = worker_item.b;
                    int ila = worker_item.ila;
                    int k = worker_item.k;
                    int n = worker_item.n;
                    int kmax = worker_item.kmax;
                    bool to_notify = worker_item.to_notify;
                    auto & ms = worker_item.ms;
                    auto & cms = worker_item.cms;
                    
                    int nr = bs[a].second;
                    int nc = bs[b].second;
                    fmpq_mat_t mat;
                    fmpq_mat_init(mat,nr,nc);
                    
                    fmpz_t Dm;
                    fmpz_init(Dm);
                    fmpq_t q;
                    fmpq_init(q);
                    for(int m : ms) {
                        if(m==0 && k+1<kmax) {
                            fmpz_poly_get_coeff_fmpz(Dm,D[a][b],0);
                            fmpz_neg(Dm,Dm);
                            fmpq_mat_scalar_mul_fmpz(mat,U[a][b][ila][k+1][n],Dm);
                        } else if(m>0) {
                            fmpq_add_si(q, qlas[ila], n-m); // q = la+n-m
                            fmpq_mat_scalar_mul_fmpq(mat,U[a][b][ila][k][n-m],q);
                            if(k+1<kmax) fmpq_mat_add(mat,mat,U[a][b][ila][k+1][n-m]);
                            fmpz_poly_get_coeff_fmpz(Dm,D[a][b],m);
                            fmpz_neg(Dm,Dm);
                            fmpq_mat_scalar_mul_fmpz(mat,mat,Dm);
                        } else abort();
                        fmpq_mat_add(U[a][b][ila][k][n],U[a][b][ila][k][n],mat);
                        ++worker_nn;
                    }
                    fmpz_clear(Dm);
                    fmpq_clear(q);
                    for(auto kv : cms) {
                        auto c = kv.first;
                        auto m = kv.second;
                        int nc2 = bs[c].second;
                        fmpq_mat_t Amac;
                        fmpq_mat_init(Amac,nr,nc2);
                        A[a][b][c-b].coeff(Amac,m);
                        fmpq_mat_mul(mat,Amac,U[c][b][ila][k][n-m]);
                        fmpq_mat_clear(Amac);
                        fmpq_mat_add(U[a][b][ila][k][n],U[a][b][ila][k][n],mat);
                        ++worker_nn;
                    }
                    fmpq_mat_clear(mat);
                    if(to_notify) {
                        fmpz_t D0;
                        fmpz_init(D0);
                        fmpz_poly_get_coeff_fmpz(D0,D[a][b],0);
                        fmpq_t q;
                        fmpq_init(q);
                        fmpq_set(q,qlas[ila]);
                        fmpq_add_si(q,q,n); // q = la+n
                        fmpq_mul_fmpz(q,q,D0); // q = (la+n)D0
                        fmpq_mat_t invA;
                        fmpq_mat_init(invA,nr,nr);
                        fmpq_mat_one(invA);
                        fmpq_mat_scalar_mul_fmpq(invA,invA,q);
                        fmpq_mat_t A0aa;
                        fmpq_mat_init(A0aa,nr,nr);
                        A[a][b][a-b].coeff(A0aa,0);
                        fmpq_mat_sub(invA,invA,A0aa);
                        int nok = fmpq_mat_inv(invA,invA);
                        if(!nok) {
                            cout << "fmpq_mat_inv failed." << endl;
                            abort();
                        }
                        fmpq_mat_mul(U[a][b][ila][k][n],invA,U[a][b][ila][k][n]);
                        fmpz_clear(D0);
                        fmpq_clear(q);
                        fmpq_mat_clear(invA);
                        fmpq_mat_clear(A0aa);
                    }
                }
                
                bool any_notify = false;
                {
                    lock_guard<mutex> guard(worker_mutex);
                    for(auto const & worker_item : wi_list) {
                        int a = worker_item.a;
                        int b = worker_item.b;
                        int ila = worker_item.ila;
                        int k = worker_item.k;
                        int n = worker_item.n;
                        bool to_notify = worker_item.to_notify;
                        if(!any_notify && to_notify) any_notify = true;
                        Used[a][b][ila][k][n] = false;
                        if(to_notify) Done[a][b][ila][k][n] = true;
                    }
                    worker_left -= worker_nn;
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                \r" << flush;
                        cout << pre << "\\--Series U(x^" << xn << ") " << worker_total << "|" << (worker_total-worker_left) << flush;
                    }
                }
                if(any_notify) worker_cond.notify_all();
            }
        };

        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--Series U(x^" << xn << ") " << worker_total << "|0" << flush;
        int total_threads = Threads;
        if(total_threads<=0) total_threads = CpuCores();
        thread worker[total_threads];
        for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        for(int i=0; i<total_threads; ++i) worker[i].join();
        for(int a=0; a<nbs; a++) for(int b=0; b<=a; b++) {
            fmpz_poly_clear(D[a][b]);
            flint_free(D[a][b]);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;

    }
    
    // U[a][b][la][k][n], Parallel fmpq version
    void DEX::ab_series(abikn_fmpq_mat_t & U, int xn, const vector<fmpq*> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        
        // Parallel Data
        vector<vector<MQ>> A(nbs); // M=A/x for A[a][b]
        vector<vector<fmpz_poly_struct*>> DA(nbs);
        vector<vector<slong>> sA(nbs); // degree for A[a][b]
        vector<fmpz_poly_struct*> D(nbs); // Denominator, i.e., lcm
        vector<vector<map<int,vector<vector<bool>>>>> Done(nbs); // [a][b][ila][k][n] : ready or not
        vector<vector<map<int,vector<vector<bool>>>>> Used(nbs); // [a][b][ila][k][n] : used or not
        vector<vector<map<int,vector<vector<list<int>>>>>> S1(nbs); // [a][b][ila][k][n] : { m, ... }
        vector<vector<map<int,vector<vector<list<pair<int,int>>>>>>> S2(nbs); // [a][b][ila][k][n] -> { <c,m>, ... }
        unsigned long long worker_left = 0;
                
        U.resize(nbs);
        for(int a=0; a<nbs; a++) {
            A[a].resize(a+1); // M=A/x
            DA[a].resize(a+1); // M=A/x
            sA[a].resize(a+1);
            U[a].resize(a+1);
            Done[a].resize(a+1);
            Used[a].resize(a+1);
            S1[a].resize(a+1);
            S2[a].resize(a+1);
            for(int b=0; b<=a; b++) {
                int nr = bs[a].second;
                int nc = bs[b].second;
                for(auto const & kv : UK[a][b]) { // 1-cycle over lambda
                    auto ila = kv.first;
                    auto kmax = kv.second;
                    if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                    U[a][b][ila].resize(kmax);
                    Done[a][b][ila].resize(kmax);
                    Used[a][b][ila].resize(kmax);
                    S1[a][b][ila].resize(kmax);
                    S2[a][b][ila].resize(kmax);
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][ila][k].resize(xn+1);
                        U[a][b][ila][k][0] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(U[a][b][ila][k][0],nr,nc); // call fmpq_mat_clear() after this call
                        fmpq_mat_set(U[a][b][ila][k][0],U0[a][b][ila][k][0]);
                    }
                }
            }
        }
        if(xn==0) {
            if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--Series U(x^" << xn << ") @ " << now(false) << endl;
            return;
        }
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--initializing Series(U)" << flush;
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        if(fmpz_poly_set_str(x1, "2  0 1")) throw Error("fmpz_poly_set_str failed.");
        fmpz_poly_one(x0);

        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b>a) continue;
            A[a][b].init(Mat[a][b]);
            A[a][b].scale(x1); // A=x*M
            DA[a][b] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(DA[a][b]);
            A[a][b].denlcm(DA[a][b]);
        }
        #pragma omp parallel for schedule(dynamic,1)
        for(int a=0; a<nbs; a++) { // cycle rows
            D[a] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(D[a]);
            fmpz_poly_set(D[a], x0);
            for(int b=0; b<=a; b++) {
                fmpz_poly_lcm(D[a], D[a], DA[a][b]);
            }
        }
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b>a) continue;
            fmpz_poly_div(DA[a][b], D[a], DA[a][b]);
            A[a][b].scale(DA[a][b]);
            fmpz_poly_clear(DA[a][b]);
            flint_free(DA[a][b]);
            sA[a][b] = A[a][b].degree();
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);        
        #pragma omp parallel for schedule(dynamic,1) collapse(2) reduction(+: worker_left)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b>a) continue;
            int ab1 = a-b+1;
            int nr = bs[a].second;
            int nc = bs[b].second;
            slong s0 = fmpz_poly_degree(D[a]);
            
            for(auto const & kv : UK[a][b]) { // 1-cycle over lambda
                auto ila = kv.first;
                auto kmax = kv.second;
                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    
                    Done[a][b][ila][k].resize(xn+1);
                    Done[a][b][ila][k][0] = true;
                    Used[a][b][ila][k].resize(xn+1);
                    Used[a][b][ila][k][0] = false;
                    S1[a][b][ila][k].resize(xn+1);
                    S2[a][b][ila][k].resize(xn+1);
                    
                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                        Done[a][b][ila][k][n] = true;
                        Used[a][b][ila][k][n] = false;
                        U[a][b][ila][k][n] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(U[a][b][ila][k][n],nr,nc); // call fmpq_mat_clear() after this call
                        fmpq_mat_zero(U[a][b][ila][k][n]);
                        
                        slong s = s0;
                        if(s>n) s=n;
                        for(int m=0; m<=s; m++) {
                            if((m==0 && k+1<kmax) || (m>0)) {
                                S1[a][b][ila][k][n].push_back(m);
                                worker_left++;
                            }
                        }

                        for(int c=b; c<=a; c++) {
                            int nc2 = bs[c].second;
                            slong ss = sA[a][c];
                            if(ss>n) ss=n;
                            for(int m=0; m<=ss; m++) {
                                if(c==a && m==0) continue;
                                auto itr = U[c][b].find(ila);
                                if(itr!=U[c][b].end() && itr->second.size()>k) {
                                    S2[a][b][ila][k][n].push_back(make_pair(c,m));
                                    worker_left++;
                                }
                            }
                        }
                        
                        if(!S1[a][b][ila][k][n].empty() || !S2[a][b][ila][k][n].empty()) Done[a][b][ila][k][n] = false;
                    }
                }
            }
        }

        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
                
        mutex worker_mutex;
        condition_variable worker_cond;
        
        auto worker_total = worker_left;
        list<u_item_type> worker_item_list;
        auto worker_main = [&]() {
            while(true) {
                list<u_item_type> wi_list;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(!worker_left) return true;
                        start:
                        int cnt = 0;
                        auto tot = worker_item_list.size();
                        while(worker_item_list.size()>0) {
                            auto itr = worker_item_list.begin();
                            wi_list.emplace_back(std::move(*itr));
                            worker_item_list.erase(itr);
                            cnt++;
                            if(cnt>10) break;
                        }
                        if(cnt) return true;
                        for(int a=0; a<nbs; a++) for(int b=0; b<=a; b++) { // cycle rows & cols
                            for(auto const & kv : UK[a][b]) { // 1-cycle over lambda
                                int ila = kv.first;
                                int kmax = kv.second;
                                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                                    auto & SS1 = S1[a][b][ila][k];
                                    auto & SS2 = S2[a][b][ila][k];
                                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                                        if(Used[a][b][ila][k][n] || Done[a][b][ila][k][n]) continue;
                                        
                                        list<int> ms;
                                        list<pair<int,int>> cms;
                                        
                                        list<int> & s1 = SS1[n];
                                        for(auto itr = s1.begin(); itr != s1.end(); ) {
                                            int m = *itr;
                                            if(m==0) {
                                                if(Done[a][b][ila][k+1][n]) {
                                                    itr = s1.erase(itr);
                                                    ms.push_back(m);
                                                } else ++itr;
                                            } else {
                                                if(Done[a][b][ila][k][n-m] && (k+1>=kmax || Done[a][b][ila][k+1][n-m])) {
                                                    itr = s1.erase(itr);
                                                    ms.push_back(m);
                                                } else ++itr;
                                            }
                                        }
                                        
                                        list<pair<int,int>> & s2 = SS2[n];
                                        for(auto itr = s2.begin(); itr != s2.end(); ) {
                                            int c = itr->first;
                                            int m = itr->second;
                                            if(Done[c][b][ila][k][n-m]) {
                                                itr = s2.erase(itr);
                                                cms.push_back(make_pair(c,m));
                                            } else ++itr;
                                        }
                                        
                                        if(!ms.empty() || !cms.empty()) {
                                            Used[a][b][ila][k][n] = true;
                                            bool to_notify = false;
                                            if(s1.empty() && s2.empty()) to_notify = true;
                                            worker_item_list.emplace_back(a,b,ila,k,n,kmax,to_notify,std::move(ms),std::move(cms));
                                        }
                                    }
                                }
                            }
                        }
                        if(worker_item_list.size()==0) return false;
                        goto start; // for return to check the start label
                    });
                    if(!worker_left) {
                        flint_cleanup();
                        return;
                    }
                }
                
                unsigned long long worker_nn = 0;
                for(auto const & worker_item : wi_list) {
                    int a = worker_item.a;
                    int b = worker_item.b;
                    int ila = worker_item.ila;
                    int k = worker_item.k;
                    int n = worker_item.n;
                    int kmax = worker_item.kmax;
                    bool to_notify = worker_item.to_notify;
                    auto & ms = worker_item.ms;
                    auto & cms = worker_item.cms;
                    
                    int nr = bs[a].second;
                    int nc = bs[b].second;
                    fmpq_mat_t mat;
                    fmpq_mat_init(mat,nr,nc);
                    
                    fmpz_t Dm;
                    fmpz_init(Dm);
                    fmpq_t q;
                    fmpq_init(q);
                    for(int m : ms) {
                        if(m==0 && k+1<kmax) {
                            fmpz_poly_get_coeff_fmpz(Dm, D[a], 0);
                            fmpz_neg(Dm,Dm);
                            fmpq_mat_scalar_mul_fmpz(mat, U[a][b][ila][k+1][n], Dm);
                        } else if(m>0) {
                            fmpq_add_si(q, qlas[ila], n-m); // q = la+n-m
                            fmpq_mat_scalar_mul_fmpq(mat, U[a][b][ila][k][n-m], q);
                            if(k+1<kmax) fmpq_mat_add(mat, mat, U[a][b][ila][k+1][n-m]);
                            fmpz_poly_get_coeff_fmpz(Dm, D[a], m);
                            fmpz_neg(Dm,Dm);
                            fmpq_mat_scalar_mul_fmpz(mat, mat, Dm);
                        } else abort();
                        fmpq_mat_add(U[a][b][ila][k][n], U[a][b][ila][k][n], mat);
                        ++worker_nn;
                    }
                    fmpz_clear(Dm);
                    fmpq_clear(q);
                    for(auto kv : cms) {
                        auto c = kv.first;
                        auto m = kv.second;
                        int nc2 = bs[c].second;
                        fmpq_mat_t Amac;
                        fmpq_mat_init(Amac, nr, nc2);
                        A[a][c].coeff(Amac, m);
                        fmpq_mat_mul(mat, Amac, U[c][b][ila][k][n-m]);
                        fmpq_mat_clear(Amac);
                        fmpq_mat_add(U[a][b][ila][k][n], U[a][b][ila][k][n], mat);
                        ++worker_nn;
                    }
                    fmpq_mat_clear(mat);
                    if(to_notify) {
                        fmpz_t D0;
                        fmpz_init(D0);
                        fmpz_poly_get_coeff_fmpz(D0, D[a], 0);
                        fmpq_t q;
                        fmpq_init(q);
                        fmpq_set(q, qlas[ila]);
                        fmpq_add_si(q, q, n); // q = la+n
                        fmpq_mul_fmpz(q, q, D0); // q = (la+n)D0
                        fmpq_mat_t invA;
                        fmpq_mat_init(invA, nr, nr);
                        fmpq_mat_one(invA);
                        fmpq_mat_scalar_mul_fmpq(invA, invA, q);
                        fmpq_mat_t A0aa;
                        fmpq_mat_init(A0aa, nr, nr);
                        A[a][a].coeff(A0aa,0);
                        fmpq_mat_sub(invA, invA, A0aa);
                        int nok = fmpq_mat_inv(invA, invA);
                        if(!nok) {
                            cout << "fmpq_mat_inv failed." << endl;
                            abort();
                        }
                        fmpq_mat_mul(U[a][b][ila][k][n], invA, U[a][b][ila][k][n]);
                        fmpz_clear(D0);
                        fmpq_clear(q);
                        fmpq_mat_clear(invA);
                        fmpq_mat_clear(A0aa);
                    }
                }
                
                bool any_notify = false;
                {
                    lock_guard<mutex> guard(worker_mutex);
                    for(auto const & worker_item : wi_list) {
                        int a = worker_item.a;
                        int b = worker_item.b;
                        int ila = worker_item.ila;
                        int k = worker_item.k;
                        int n = worker_item.n;
                        bool to_notify = worker_item.to_notify;
                        if(!any_notify && to_notify) any_notify = true;
                        Used[a][b][ila][k][n] = false;
                        if(to_notify) Done[a][b][ila][k][n] = true;
                    }
                    worker_left -= worker_nn;
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                \r" << flush;
                        cout << pre << "\\--Series U(x^" << xn << ") " << worker_total << "|" << (worker_total-worker_left) << flush;
                    }
                }
                if(any_notify) worker_cond.notify_all();
            }
        };

        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--Series U(x^" << xn << ") " << worker_total << "|0" << flush;
        int total_threads = Threads;
        if(total_threads<=0) total_threads = CpuCores();
        thread worker[total_threads];
        for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        for(int i=0; i<total_threads; ++i) worker[i].join();
        for(int a=0; a<nbs; a++) {
            fmpz_poly_clear(D[a]);
            flint_free(D[a]);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;

    }
    
    //=*********************************************************************=
    // U-Series
    //=*********************************************************************=
    
    // U[a][b][la][k][n] - Parallel acf version - only on a & b, note acf version, aka gr version
    void DEX::ab_series(abikn_gr_mat_t & U, int xn, gr_ctx_t ctx, const vector<fmpq*> & qslas) {
        if(!fuchsified) fuchsify();
        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--initializing nSeries(U)" << flush;
        
        auto nbs = bs.size();
        
        // Parallel Data
        vector<vector<MQ>> A(nbs); // M=A/x for A[a][b]
        vector<vector<fmpz_poly_struct*>> DA(nbs);
        vector<vector<slong>> sA(nbs); // degree for A[a][b]
        vector<fmpz_poly_struct*> D(nbs); // Denominator
        vector<vector<bool>> Done(nbs); // [a][b] : ready or not
        vector<vector<bool>> Used(nbs); // [a][b] : used or not
        vector<vector<list<int>>> Sc(nbs); // [a][b] : { c, ... }
        unsigned long long worker_left = 0;
        
        U.resize(nbs);
        for(int a=0; a<nbs; a++) {
            A[a].resize(a+1); // M=A/x
            DA[a].resize(a+1); // M=A/x
            sA[a].resize(a+1);
            U[a].resize(a+1);
            Done[a].resize(a+1);
            Used[a].resize(a+1);
            Sc[a].resize(a+1);
        }
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        if(fmpz_poly_set_str(x1, "2  0 1")) throw Error("fmpz_poly_set_str failed.");
        fmpz_poly_one(x0);
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b>a) continue;
            A[a][b].init(Mat[a][b]);
            A[a][b].scale(x1); // A=x*M
            DA[a][b] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(DA[a][b]);
            A[a][b].denlcm(DA[a][b]);
        }
        #pragma omp parallel for schedule(dynamic,1)
        for(int a=0; a<nbs; a++) { // cycle rows
            D[a] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(D[a]);
            fmpz_poly_set(D[a], x0);
            for(int b=0; b<=a; b++) {
                fmpz_poly_lcm(D[a], D[a], DA[a][b]);
            }
        }
        #pragma omp parallel for schedule(dynamic,1) collapse(2) reduction(+: worker_left)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b>a) continue;
            int nr = bs[a].second;
            int nc = bs[b].second;
            fmpz_poly_div(DA[a][b], D[a], DA[a][b]);
            A[a][b].scale(DA[a][b]);
            fmpz_poly_clear(DA[a][b]);
            flint_free(DA[a][b]);
            sA[a][b] = A[a][b].degree();
            Done[a][b] = false;
            Used[a][b] = false;
            for(int c=b; c<a; c++) Sc[a][b].push_back(c);
            worker_left += a-b+1;
            for(auto const & kv : UK[a][b]) {
                auto ila = kv.first;
                auto kmax = kv.second;
                if(qslas.size()>0 && !is_resonant(qslas, qlas[ila])) continue; // selected lambda set
                U[a][b][ila].resize(kmax);
                for(int k=kmax-1; k>=0; k--) {
                    U[a][b][ila][k].resize(xn+1);
                    U[a][b][ila][k][0] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                    gr_mat_init(U[a][b][ila][k][0], nr, nc, ctx); // call fmpq_mat_clear() after this call
                    gr_mat_set_fmpq_mat(U[a][b][ila][k][0], U0[a][b][ila][k][0], ctx);
                    for(int n=1; n<=xn; n++) {
                        U[a][b][ila][k][n] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                        gr_mat_init(U[a][b][ila][k][n], nr, nc, ctx); // call fmpq_mat_clear() after this call
                        gr_mat_zero(U[a][b][ila][k][n], ctx);
                    }
                }
            }
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;

        mutex worker_mutex;
        condition_variable worker_cond;
        
        int total_threads = Threads;
        if(total_threads<=0) total_threads = CpuCores();
        auto worker_total = worker_left;
        auto worker_main = [&]() {
            while(true) {
                int a, b;
                list<int> cset;
                bool to_notify = false;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(!worker_left) return true;
                        for(a=0; a<nbs; a++) for(b=a; b>=0; b--) { // cycle rows & cols
                            if(Used[a][b] || Done[a][b]) continue;
                            list<int> & sc = Sc[a][b];
                            for(auto itr = sc.begin(); itr != sc.end(); ) {
                                int c = *itr;
                                if(Done[c][b]) {
                                    itr = sc.erase(itr);
                                    cset.push_back(c);
                                } else ++itr;
                            }
                            
                            if(!cset.empty() || a==b) { // note a==b
                                Used[a][b] = true;
                                if(sc.empty()) to_notify = true;
                                return true;
                            }
                        }
                        return false;
                    });
                    if(!worker_left) {
                        flint_cleanup();
                        return;
                    }
                }
                
                int status = GR_SUCCESS;
                int nr = bs[a].second;
                int nc = bs[b].second;
                gr_mat_t mat;
                gr_mat_init(mat, nr, nc, ctx);
                
                for(auto c : cset) {
                    int nc2 = bs[c].second;
                    gr_mat_t Amac;
                    gr_mat_init(Amac, nr, nc2, ctx);
                    for(auto const & kv : UK[a][b]) { // cycle over lambda
                        int ila = kv.first;
                        int kmax = kv.second;
                        if(qslas.size()>0 && !is_resonant(qslas, qlas[ila])) continue; // selected lambda set
                        for(int k=kmax-1; k>=0; k--) { // cycle over k
                            auto itr = UK[c][b].find(ila);
                            if(itr!=UK[c][b].end() && itr->second>k) {
                                for(int n=1; n<=xn; n++) { // cycle over n
                                    int s = sA[a][c];
                                    if(s>n) s = n;
                                    for(int m=0; m<=s; m++) {
                                        A[a][c].coeff(Amac, m, ctx);
                                        status |= gr_mat_mul(mat, Amac, U[c][b][ila][k][n-m], ctx);
                                        status |= gr_mat_add(U[a][b][ila][k][n], U[a][b][ila][k][n], mat, ctx);
                                    }
                                }
                            }
                        }
                    }
                    gr_mat_clear(Amac, ctx);
                }
                unsigned long long worker_nn = cset.size();
                if(to_notify) {
                    worker_nn++;
                    gr_mat_t Amaa;
                    gr_mat_init(Amaa, nr, nr, ctx);
                    fmpz_t Dm;
                    fmpz_init(Dm);
                    gr_ptr z = gr_heap_init(ctx);
                    
                    fmpz_t D0;
                    fmpz_init(D0);
                    fmpz_poly_get_coeff_fmpz(D0, D[a], 0);
                    fmpq_t q;
                    fmpq_init(q);
                    fmpq_mat_t invA;
                    fmpq_mat_init(invA, nr, nr);
                    fmpq_mat_t A0aa;
                    fmpq_mat_init(A0aa, nr, nr);
                    A[a][a].coeff(A0aa, 0);
                    gr_mat_t invA_acb;
                    gr_mat_init(invA_acb, nr, nr, ctx);
                    
                    auto saa = sA[a][a];
                    auto sa = fmpz_poly_degree(D[a]);
                    
                    for(auto const & kv : UK[a][b]) { // cycle over lambda
                        int ila = kv.first;
                        int kmax = kv.second;
                        if(qslas.size()>0 && !is_resonant(qslas, qlas[ila])) continue; // selected lambda set
                        for(int k=kmax-1; k>=0; k--) {
                            for(int n=1; n<=xn; n++) {
                                auto s = saa;
                                if(s>n) s = n;
                                for(int m=1; m<=s; m++) {
                                    A[a][a].coeff(Amaa, m, ctx);
                                    status |= gr_mat_mul(mat, Amaa, U[a][b][ila][k][n-m], ctx);
                                    status |= gr_mat_add(U[a][b][ila][k][n], U[a][b][ila][k][n], mat, ctx);
                                }
                                if(k+1<kmax) {
                                    status |= gr_set_fmpz(z, D0, ctx);
                                    status |= gr_mat_mul_scalar(mat, U[a][b][ila][k+1][n], z, ctx);
                                    status |= gr_mat_sub(U[a][b][ila][k][n], U[a][b][ila][k][n], mat, ctx);
                                }
                                s = sa;
                                if(s>n) s = n;
                                for(int m=1; m<=s; m++) {
                                    status |= gr_set_fmpq(z, qlas[ila], ctx);
                                    status |= gr_add_si(z, z, n-m, ctx); // la+n-m
                                    status |= gr_mat_mul_scalar(mat, U[a][b][ila][k][n-m], z, ctx);
                                    if(k+1<kmax) status |= gr_mat_add(mat, mat, U[a][b][ila][k+1][n-m], ctx);
                                    fmpz_poly_get_coeff_fmpz(Dm, D[a], m);
                                    fmpz_neg(Dm, Dm);
                                    status |= gr_set_fmpz(z, Dm, ctx);
                                    status |= gr_mat_addmul_scalar(U[a][b][ila][k][n], mat, z, ctx);
                                }
                                fmpq_add_si(q, qlas[ila], n); // q = la+n
                                fmpq_mul_fmpz(q, q, D0); // q = (la+n)D0
                                
                                fmpq_mat_one(invA);
                                fmpq_mat_scalar_mul_fmpq(invA, invA, q);

                                fmpq_mat_sub(invA, invA, A0aa);
                                int nok = fmpq_mat_inv(invA, invA);
                                if(!nok) {
                                    cout << "fmpq_mat_inv failed." << endl;
                                    abort();
                                }
                                status |= gr_mat_set_fmpq_mat(invA_acb, invA, ctx);
                                status |= gr_mat_mul(U[a][b][ila][k][n], invA_acb, U[a][b][ila][k][n], ctx);
                            }
                        }
                    }
                    fmpz_clear(Dm);
                    gr_heap_clear(z, ctx);
                    gr_mat_clear(Amaa, ctx);
                    fmpz_clear(D0);
                    fmpq_clear(q);
                    fmpq_mat_clear(A0aa);
                    fmpq_mat_clear(invA);
                    gr_mat_clear(invA_acb, ctx);
                }
                gr_mat_clear(mat, ctx);
                {
                    lock_guard<mutex> guard(worker_mutex);
                    worker_left -= worker_nn;
                    Used[a][b] = false;
                    if(to_notify) Done[a][b] = true;
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                \r" << flush;
                        cout << pre << "\\--nSeries U(x^" << xn << ") " << worker_total << "|" << (worker_total-worker_left) << flush;
                    }
                }
                if(to_notify) worker_cond.notify_all();

            }
        };

        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--nSeries U(x^" << xn << ") " << worker_total << "|0" << flush;
        thread worker[total_threads];
        for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        for(int i=0; i<total_threads; ++i) worker[i].join();
        for(int a=0; a<nbs; a++) {
            fmpz_poly_clear(D[a]);
            flint_free(D[a]);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
        
    }
    
    //=*********************************************************************=
    // I-Series - expansion - rational
    //=*********************************************************************=
    
    struct i_item_type {
        int a, ila, k, n, kmax;
        bool to_notify = false;
        list<int> ms;
        list<pair<int,int>> bms;
        i_item_type(int a_, int ila_, int k_, int n_, int kmax_, bool to_notify_, list<int> && ms_, list<pair<int,int>> && bms_) : a(a_), ila(ila_), k(k_), n(n_), kmax(kmax_), to_notify(to_notify_), ms(std::move(ms_)), bms(std::move(bms_)) { }
        i_item_type() { }
    };
    
    // I[a][la][k][n] & In0[a][ila][k], no need to initialize I
    void DEX::series(aikn_fmpq_mat_t & I, int xn, aikn_fmpq_mat_t & In0, int nc, const vector<fmpq*> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        I.resize(nbs);
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        if(fmpz_poly_set_str(x1, "2  0 1")) throw Error("fmpz_poly_set_str failed.");
        fmpz_poly_one(x0);
        
        for(int a=0; a<nbs; a++) { // cycle rows
            // [(la+n)D0-A0aa].Ia(la,k,n) = -D0 Ia(la,k+1,n)
            //   - sum_{0<m<=n} Dm [(la+n-m).Ia(la,k,n-m)+Ia(la,k+1,n-m)] 
            //   + sum_{not(b=a|m=0)} Amab.Ib(la,k,n-m)
            
            int nr = bs[a].second;
            vector<MQ> A(a+1); // M=A/x
            vector<int> sdeg(a+1); // A[b].degree();
            fmpz_poly_t lcm; // D=lcm
            vector<fmpz_poly_struct*> lcm_vec(a+1);
            fmpz_poly_init(lcm);
            fmpz_poly_set(lcm,x0);
            for(int b=0; b<=a; b++) { 
                // lcm for each block A[b] in each row
                A[b].init(Mat[a][b]);
                A[b].scale(x1); // A=x*M
                lcm_vec[b] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                fmpz_poly_init(lcm_vec[b]);
                A[b].denlcm(lcm_vec[b]);
                fmpz_poly_lcm(lcm, lcm, lcm_vec[b]);
            }
            for(int b=0; b<=a; b++) {
                fmpz_poly_div(lcm_vec[b], lcm, lcm_vec[b]);
                A[b].scale(lcm_vec[b]);
                fmpz_poly_clear(lcm_vec[b]);
                flint_free(lcm_vec[b]);
                sdeg[b] = A[b].degree();
            }
            
            fmpq_mat_t A0aa;
            fmpq_mat_init(A0aa,nr,nr);
            A[a].coeff(A0aa,0);
            fmpz_t D0;
            fmpz_init(D0);
            fmpz_poly_get_coeff_fmpz(D0,lcm,0);
            
            for(auto kv : IK[a]) { // insert kv first and make sure thread-safe
                if(qslas.size()>0 && !is_resonant(qslas,qlas[kv.first])) continue;
                I[a][kv.first].resize(kv.second); 
            }
            
            int nla = IK[a].size();
            if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                              \r" << flush;
                cout << pre << "\\--series: " << nbs << "|" << a+1;
                cout << " [" << nr << "\u2A09" << nc << "]";
                cout << " \u03BB" << nla << " x^" << xn << flush;
            }
            #pragma omp parallel for schedule(dynamic,1)
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
                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // seleted lambda set
                auto kmax = kv->second;  
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    I[a][ila][k].resize(xn+1);
                    I[a][ila][k][0] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                    fmpq_mat_init(I[a][ila][k][0],nr,nc); // call fmpq_mat_clear() after this call
                    fmpq_mat_set(I[a][ila][k][0],In0[a][ila][k][0]);
                    
                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                  
                        fmpq_mat_zero(smat);
                        fmpq_add_si(q,qlas[ila],n); // q = la+n

                        if(true) { 
                            slong s = fmpz_poly_degree(lcm);
                            if(s>n) s=n;
                            int tmax = omp_get_max_threads();
                            vector<fmpq_mat_struct*> mat_sum(tmax), mat_tmp(tmax);
                            for(int i=0; i<tmax; i++) {
                                mat_sum[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                                fmpq_mat_init(mat_sum[i],nr,nc);
                                mat_tmp[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                                fmpq_mat_init(mat_tmp[i],nr,nc);
                            }
                            #pragma omp parallel for schedule(runtime)
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
                            }
                            for(int i=0; i<tmax; i++) {
                                fmpq_mat_add(smat,smat,mat_sum[i]);
                                fmpq_mat_clear(mat_sum[i]);
                                flint_free(mat_sum[i]);
                                fmpq_mat_clear(mat_tmp[i]);
                                flint_free(mat_tmp[i]);
                            }
                        }
                        
                        if(true) {
                            int tmax = omp_get_max_threads();
                            vector<fmpq_mat_struct*> mat_sum(tmax), mat_tmp(tmax);
                            for(int i=0; i<tmax; i++) {
                                mat_sum[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                                fmpq_mat_init(mat_sum[i],nr,nc);
                                mat_tmp[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                                fmpq_mat_init(mat_tmp[i],nr,nc);
                            }
                            #pragma omp parallel for schedule(runtime) 
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
                            }
                            for(int i=0; i<tmax; i++) {
                                fmpq_mat_add(smat,smat,mat_sum[i]);
                                fmpq_mat_clear(mat_sum[i]);
                                flint_free(mat_sum[i]);
                                fmpq_mat_clear(mat_tmp[i]);
                                flint_free(mat_tmp[i]);
                            }
                        }
                        
                        fmpq_mul_fmpz(q,q,D0); // q = (la+n)D0
                        fmpq_mat_one(invA);
                        fmpq_mat_scalar_mul_fmpq(invA,invA,q);
                        fmpq_mat_sub(invA,invA,A0aa);
                        int nok = fmpq_mat_inv(invA,invA);
                        if(!nok) {
                            cout << "fmpq_mat_inv failed." << endl;
                            abort();
                        }
                        
                        I[a][ila][k][n] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(I[a][ila][k][n],nr,nc); // call fmpq_mat_clear() after this call
                        fmpq_mat_mul(I[a][ila][k][n],invA,smat);
                    }
                }
                fmpq_mat_clear(invA);
                fmpq_clear(q);
                fmpq_mat_clear(smat);
            }
            for(int b=0; b<=a; b++) A[b].clear();
            fmpz_poly_clear(lcm);
            fmpq_mat_clear(A0aa);
            fmpz_clear(D0);
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
    }
    
    
    //=*********************************************************************=
    // I-Series
    //=*********************************************************************=
    
    // I[a][la][k][n] - Parallel gr/acf version - on all a, ila, k, n
    void DEX::all_series(aikn_gr_mat_t & I, int xn, aikn_gr_mat_t & In0, int nc, gr_ctx_t ctx, const vector<fmpq*> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        I.resize(nbs);
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--initializing Series(I)" << flush;
                
        // Parallel Data
        vector<vector<MQ>> A(nbs); // M=A/x
        vector<vector<fmpz_poly_struct*>> DA(nbs);
        vector<vector<slong>> sA(nbs); // degree for A[a][b]
        vector<fmpz_poly_struct*> D(nbs); // Denominator, i.e., lcm
        vector<map<int,vector<vector<bool>>>> Done(nbs); // [a][ila][k][n] : ready or not
        vector<map<int,vector<vector<bool>>>> Used(nbs); // [a][ila][k][n] : used or not
        vector<map<int,vector<vector<list<int>>>>> S1(nbs); // [a][ila][k][n] : { m, ... }
        vector<map<int,vector<vector<list<pair<int,int>>>>>> S2(nbs); // [a][ila][k][n] -> { <b,m>, ... }
        unsigned long long worker_left = 0;
        
        int status = GR_SUCCESS;
        for(int a=0; a<nbs; a++) {
            int nr = bs[a].second;
            A[a].resize(a+1); // M=A/x
            DA[a].resize(a+1);
            sA[a].resize(a+1);
            for(auto kv : IK[a]) { // 1-cycle over lambda
                auto ila = kv.first;
                auto kmax = kv.second;
                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                I[a][ila].resize(kmax);
                Done[a][ila].resize(kmax);
                Used[a][ila].resize(kmax);
                S1[a][ila].resize(kmax);
                S2[a][ila].resize(kmax);
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    I[a][ila][k].resize(xn+1);
                    I[a][ila][k][0] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                    gr_mat_init(I[a][ila][k][0], nr, nc, ctx); // call fmpq_mat_clear() after this call
                    status |= gr_mat_set(I[a][ila][k][0], In0[a][ila][k][0], ctx);
                    Done[a][ila][k].resize(xn+1);
                    Done[a][ila][k][0] = true;
                    Used[a][ila][k].resize(xn+1);
                    Used[a][ila][k][0] = false;
                    S1[a][ila][k].resize(xn+1);
                    S2[a][ila][k].resize(xn+1);
                }
            }
        }
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        if(fmpz_poly_set_str(x1, "2  0 1")) throw Error("fmpz_poly_set_str failed.");
        fmpz_poly_one(x0);
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b<=a) {
                A[a][b].init(Mat[a][b]);
                A[a][b].scale(x1); // A=x*M
                DA[a][b] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                fmpz_poly_init(DA[a][b]);
                A[a][b].denlcm(DA[a][b]);
            }
        }
        #pragma omp parallel for schedule(dynamic,1)
        for(int a=0; a<nbs; a++) { // cycle rows
            D[a] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(D[a]);
            fmpz_poly_set(D[a], x0);
            for(int b=0; b<=a; b++) {
                fmpz_poly_lcm(D[a], D[a], DA[a][b]);
            }
        }
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b<=a) {
                fmpz_poly_div(DA[a][b], D[a], DA[a][b]);
                A[a][b].scale(DA[a][b]);
                fmpz_poly_clear(DA[a][b]);
                flint_free(DA[a][b]);
                sA[a][b] = A[a][b].degree();
            }
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        
        #pragma omp parallel for schedule(dynamic,1) reduction(+: worker_left)
        for(int a=0; a<nbs; a++) { // cycle rows
            int nr = bs[a].second;
            slong s0 = fmpz_poly_degree(D[a]);
            
            for(auto kv : IK[a]) { // 1-cycle over lambda
                auto ila = kv.first;
                auto kmax = kv.second;
                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                        Done[a][ila][k][n] = true;
                        Used[a][ila][k][n] = false;
                        I[a][ila][k][n] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                        gr_mat_init(I[a][ila][k][n], nr, nc, ctx); // call gr_mat_clear() after this call
                        status |= gr_mat_zero(I[a][ila][k][n], ctx);
                        
                        slong s = s0;
                        if(s>n) s=n;
                        for(int m=0; m<=s; m++) {
                            if((m==0 && k+1<kmax) || (m>0)) {
                                S1[a][ila][k][n].push_back(m);
                                worker_left++;
                            }
                        }

                        for(int b=0; b<=a; b++) {
                            int nc2 = bs[b].second;
                            slong ss = sA[a][b];
                            if(ss>n) ss=n;
                            for(int m=0; m<=ss; m++) {
                                if(b==a && m==0) continue;
                                auto itr = I[b].find(ila);
                                if(itr!=I[b].end() && itr->second.size()>k) {
                                    S2[a][ila][k][n].push_back(make_pair(b,m));
                                    worker_left++;
                                }
                            }
                        }
                        
                        if(!S1[a][ila][k][n].empty() || !S2[a][ila][k][n].empty()) Done[a][ila][k][n] = false;
                    }
                }
            }
        }
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
        
        mutex worker_mutex;
        condition_variable worker_cond;
        
        int total_threads = Threads;
        if(total_threads<=0) CpuCores();
        auto worker_total = worker_left;
        list<i_item_type> worker_item_list;
        auto worker_main = [&]() {
            while(true) {
                list<i_item_type> wi_list;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(!worker_left) return true;
                        start:
                        int cnt = 0;
                        auto tot = worker_item_list.size();
                        while(worker_item_list.size()>0) {
                            auto itr = worker_item_list.begin();
                            wi_list.emplace_back(std::move(*itr));
                            worker_item_list.erase(itr);
                            cnt++;
                            if(cnt>10) break;
                        }
                        if(cnt) return true;
                        for(int a=0; a<nbs; a++) { // cycle rows
                            for(auto const & kv : IK[a]) { // 1-cycle over lambda
                                int ila = kv.first;
                                int kmax = kv.second;
                                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                                    auto & SS1 = S1[a][ila][k];
                                    auto & SS2 = S2[a][ila][k];
                                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                                        if(Used[a][ila][k][n] || Done[a][ila][k][n]) continue;
                                        
                                        list<int> ms;
                                        list<pair<int,int>> bms;
                                        
                                        list<int> & s1 = SS1[n];
                                        for(auto itr = s1.begin(); itr != s1.end(); ) {
                                            int m = *itr;
                                            if(m==0) {
                                                if(Done[a][ila][k+1][n]) {
                                                    itr = s1.erase(itr);
                                                    ms.push_back(m);
                                                } else ++itr;
                                            } else {
                                                if(Done[a][ila][k][n-m] && (k+1>=kmax || Done[a][ila][k+1][n-m])) {
                                                    itr = s1.erase(itr);
                                                    ms.push_back(m);
                                                } else ++itr;
                                            }
                                        }
                                        
                                        list<pair<int,int>> & s2 = SS2[n];
                                        for(auto itr = s2.begin(); itr != s2.end(); ) {
                                            int b = itr->first;
                                            int m = itr->second;
                                            if(Done[b][ila][k][n-m]) {
                                                itr = s2.erase(itr);
                                                bms.push_back(make_pair(b,m));
                                            } else ++itr;
                                        }
                                        
                                        if(!ms.empty() || !bms.empty()) {
                                            Used[a][ila][k][n] = true;
                                            bool to_notify = false;
                                            if(s1.empty() && s2.empty()) to_notify = true;
                                            worker_item_list.emplace_back(a,ila,k,n,kmax,to_notify,std::move(ms),std::move(bms));
                                        }
                                    }
                                }
                            }
                        }
                        if(worker_item_list.size()==0) return false;
                        goto start; // for return to check the start label
                    });
                    if(!worker_left) {
                        flint_cleanup();
                        return;
                    }
                }
                
                unsigned long long worker_nn = 0;
                for(auto const & worker_item : wi_list) {
                    int a = worker_item.a;
                    int ila = worker_item.ila;
                    int k = worker_item.k;
                    int n = worker_item.n;
                    int kmax = worker_item.kmax;
                    bool to_notify = worker_item.to_notify;
                    auto & ms = worker_item.ms;
                    auto & bms = worker_item.bms;

                    int nr = bs[a].second;
                    gr_mat_t mat;
                    gr_mat_init(mat, nr, nc, ctx);
                    fmpz_t Dm;
                    fmpz_init(Dm);
                    gr_ptr qr = gr_heap_init(ctx);
                    for(int m : ms) {
                        if(m==0 && k+1<kmax) {
                            fmpz_poly_get_coeff_fmpz(Dm, D[a], 0);
                            fmpz_neg(Dm, Dm);
                            status |= gr_set_fmpz(qr, Dm, ctx);
                            status |= gr_mat_mul_scalar(mat, I[a][ila][k+1][n], qr, ctx);
                        } else if(m>0) {
                            status |= gr_set_fmpq(qr, qlas[ila], ctx);
                            status |= gr_add_si(qr, qr, n-m, ctx); // la+n-m
                            status |= gr_mat_mul_scalar(mat, I[a][ila][k][n-m], qr, ctx);
                            if(k+1<kmax) status |= gr_mat_add(mat, mat, I[a][ila][k+1][n-m], ctx);
                            fmpz_poly_get_coeff_fmpz(Dm, D[a], m);
                            fmpz_neg(Dm, Dm);
                            status |= gr_set_fmpz(qr, Dm, ctx);
                            status |= gr_mat_mul_scalar(mat, mat, qr, ctx);
                        } else abort();
                        status |= gr_mat_add(I[a][ila][k][n], I[a][ila][k][n], mat, ctx);
                        ++worker_nn;
                    }
                    gr_heap_clear(qr, ctx);
                    for(auto kv : bms) {
                        auto b = kv.first;
                        auto m = kv.second;
                        int nc2 = bs[b].second;
                        gr_mat_t Amab;
                        gr_mat_init(Amab, nr, nc2, ctx);
                        A[a][b].coeff(Amab, m, ctx);
                        status |= gr_mat_mul(mat, Amab, I[b][ila][k][n-m], ctx);
                        gr_mat_clear(Amab, ctx);
                        status |= gr_mat_add(I[a][ila][k][n], I[a][ila][k][n], mat, ctx);
                        ++worker_nn;
                    }
                    gr_mat_clear(mat, ctx);
                    if(to_notify) {
                        fmpz_poly_get_coeff_fmpz(Dm, D[a], 0); // D0
                        fmpq_t q;
                        fmpq_init(q);
                        fmpq_add_si(q, qlas[ila], n); // q = la+n
                        fmpq_mul_fmpz(q, q, Dm); // q = (la+n)D0
                        fmpq_mat_t invA;
                        fmpq_mat_init(invA, nr, nr);
                        fmpq_mat_one(invA);
                        fmpq_mat_scalar_mul_fmpq(invA, invA, q);
                        fmpq_mat_t A0aa;
                        fmpq_mat_init(A0aa, nr, nr);
                        A[a][a].coeff(A0aa, 0);
                        fmpq_mat_sub(invA, invA, A0aa);
                        int nok = fmpq_mat_inv(invA, invA);
                        if(!nok) {
                            cout << "fmpq_mat_inv failed." << endl;
                            abort();
                        }
                        fmpq_clear(q);
                        fmpq_mat_clear(A0aa);
                        gr_mat_t invA_acb;
                        gr_mat_init(invA_acb, nr, nr, ctx);
                        status |= gr_mat_set_fmpq_mat(invA_acb, invA, ctx);
                        status |= gr_mat_mul(I[a][ila][k][n], invA_acb, I[a][ila][k][n], ctx);
                        fmpq_mat_clear(invA);
                        gr_mat_clear(invA_acb, ctx);
                    }
                    fmpz_clear(Dm);
                }
                bool any_notify = false;
                {
                    lock_guard<mutex> guard(worker_mutex);
                    for(auto const & worker_item : wi_list) {
                        int a = worker_item.a;
                        int ila = worker_item.ila;
                        int k = worker_item.k;
                        int n = worker_item.n;
                        bool to_notify = worker_item.to_notify;
                        if(!any_notify && to_notify) any_notify = true;
                        Used[a][ila][k][n] = false;
                        if(to_notify) Done[a][ila][k][n] = true;
                    }
                    worker_left -= worker_nn;
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                \r" << flush;
                        cout << pre << "\\--nSeries I(x^" << xn << ") " << worker_total << "|" << (worker_total-worker_left) << flush;
                    }
                }
                if(any_notify) worker_cond.notify_all();
            }
        };

        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--nSeries I(x^" << xn << ") " << worker_total << "|0" << flush;
        thread worker[total_threads];
        for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        for(int i=0; i<total_threads; ++i) worker[i].join();
        for(int a=0; a<nbs; a++) {
            fmpz_poly_clear(D[a]);
            flint_free(D[a]);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;

    }
    
    // I[a][la][k][n] - Parallel gr/acf version - only on a
    void DEX::a_series(aikn_gr_mat_t & I, int xn, aikn_gr_mat_t & In0, int nc, gr_ctx_t ctx, const vector<fmpq*> & qslas) {
        if(!fuchsified) fuchsify();
        auto nbs = bs.size();
        I.resize(nbs);
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--initializing Series(I)" << flush;
                
        // Parallel Data
        vector<vector<MQ>> A(nbs); // M=A/x
        vector<vector<fmpz_poly_struct*>> DA(nbs);
        vector<vector<slong>> sA(nbs); // degree for A[a][b]
        vector<fmpz_poly_struct*> D(nbs); // Denominator, i.e., lcm
        vector<bool> Done(nbs); // [a]: ready or not
        vector<bool> Used(nbs); // [a]: used or not
        vector<list<int>> Sb(nbs); // [a] : { b, ... }
        unsigned long long worker_left = 0;
        
        int status = GR_SUCCESS;
        for(int a=0; a<nbs; a++) {
            int nr = bs[a].second;
            A[a].resize(a+1); // M=A/x
            DA[a].resize(a+1);
            sA[a].resize(a+1);
            for(auto kv : IK[a]) { // 1-cycle over lambda
                auto ila = kv.first;
                auto kmax = kv.second;
                if(qslas.size()>0 && !is_resonant(qslas,qlas[ila])) continue; // selected lambda set
                I[a][ila].resize(kmax);
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    I[a][ila][k].resize(xn+1);
                    I[a][ila][k][0] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                    gr_mat_init(I[a][ila][k][0], nr, nc, ctx); // call gr_mat_clear() after this call
                    status |= gr_mat_set(I[a][ila][k][0], In0[a][ila][k][0], ctx);
                    for(int n=1; n<=xn; n++) {
                        I[a][ila][k][n] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                        gr_mat_init(I[a][ila][k][n], nr, nc, ctx); // call gr_mat_clear() after this call
                        status |= gr_mat_zero(I[a][ila][k][n], ctx);
                    }
                }
            }
            Done[a] = false;
            Used[a] = false;
            for(int b=0; b<a; b++) Sb[a].push_back(b);
            worker_left += a+1;
        }
        
        fmpz_poly_t x1, x0;
        fmpz_poly_init(x1);
        fmpz_poly_init(x0);
        if(fmpz_poly_set_str(x1, "2  0 1")) throw Error("fmpz_poly_set_str failed.");
        fmpz_poly_one(x0);
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b<=a) {
                A[a][b].init(Mat[a][b]);
                A[a][b].scale(x1); // A=x*M
                DA[a][b] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                fmpz_poly_init(DA[a][b]);
                A[a][b].denlcm(DA[a][b]);
            }
        }
        #pragma omp parallel for schedule(dynamic,1)
        for(int a=0; a<nbs; a++) { // cycle rows
            D[a] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
            fmpz_poly_init(D[a]);
            fmpz_poly_set(D[a], x0);
            for(int b=0; b<=a; b++) {
                fmpz_poly_lcm(D[a], D[a], DA[a][b]);
            }
        }
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) for(int b=0; b<nbs; b++) { // cycle rows & cols
            if(b<=a) {
                fmpz_poly_div(DA[a][b], D[a], DA[a][b]);
                A[a][b].scale(DA[a][b]);
                fmpz_poly_clear(DA[a][b]);
                flint_free(DA[a][b]);
                sA[a][b] = A[a][b].degree();
            }
        }
        fmpz_poly_clear(x1);
        fmpz_poly_clear(x0);
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
        
        mutex worker_mutex;
        condition_variable worker_cond;
        
        auto worker_total = worker_left;
        auto worker_main = [&]() {
            while(true) {
                int a;
                list<int> bset;
                bool to_notify = false;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(!worker_left) return true;
                        for(a=0; a<nbs; a++) {
                            if(Used[a] || Done[a]) continue;
                            list<int> & sb = Sb[a];
                            for(auto itr = sb.begin(); itr != sb.end(); ) {
                                int b = *itr;
                                if(Done[b]) {
                                    itr = sb.erase(itr);
                                    bset.push_back(b);
                                } else ++itr;
                            }
                            
                            if(!bset.empty() || a==0) { // note a==b
                                Used[a] = true;
                                if(sb.empty()) to_notify = true;
                                return true;
                            }
                        }
                        return false;
                    });
                    if(!worker_left) {
                        flint_cleanup();
                        return;
                    }
                }
                
                int nr = bs[a].second;
                gr_mat_t mat;
                gr_mat_init(mat, nr, nc, ctx);
                
                for(auto b : bset) {
                    int nc2 = bs[b].second;
                    gr_mat_t Amab;
                    gr_mat_init(Amab, nr, nc2, ctx);
                    for(auto const & kv : IK[a]) { // cycle over lambda
                        int ila = kv.first;
                        int kmax = kv.second;
                        if(qslas.size()>0 && !is_resonant(qslas, qlas[ila])) continue; // selected lambda set
                        for(int k=kmax-1; k>=0; k--) { // cycle over k
                            auto itr = IK[b].find(ila);
                            if(itr!=IK[b].end() && itr->second>k) {
                                for(int n=1; n<=xn; n++) { // cycle over n
                                    int s = sA[a][b];
                                    if(s>n) s = n;
                                    for(int m=0; m<=s; m++) {
                                        A[a][b].coeff(Amab, m, ctx);
                                        status |= gr_mat_mul(mat, Amab, I[b][ila][k][n-m], ctx);
                                        status |= gr_mat_add(I[a][ila][k][n], I[a][ila][k][n], mat, ctx);
                                    }
                                }
                            }
                        }
                    }
                    gr_mat_clear(Amab, ctx);
                }
                unsigned long long worker_nn = bset.size();
                if(to_notify) {
                    worker_nn++;
                    gr_mat_t Amaa;
                    gr_mat_init(Amaa, nr, nr, ctx);
                    fmpz_t Dm;
                    fmpz_init(Dm);
                    gr_ptr z = gr_heap_init(ctx);
                    
                    fmpz_t D0;
                    fmpz_init(D0);
                    fmpz_poly_get_coeff_fmpz(D0, D[a], 0);
                    fmpq_t q;
                    fmpq_init(q);
                    fmpq_mat_t invA;
                    fmpq_mat_init(invA, nr, nr);
                    fmpq_mat_t A0aa;
                    fmpq_mat_init(A0aa, nr, nr);
                    A[a][a].coeff(A0aa, 0);
                    gr_mat_t invA_acb;
                    gr_mat_init(invA_acb, nr, nr, ctx);
                    
                    auto saa = sA[a][a];
                    auto sa = fmpz_poly_degree(D[a]);
                    
                    for(auto const & kv : IK[a]) { // cycle over lambda
                        int ila = kv.first;
                        int kmax = kv.second;
                        if(qslas.size()>0 && !is_resonant(qslas, qlas[ila])) continue; // selected lambda set
                        for(int k=kmax-1; k>=0; k--) {
                            for(int n=1; n<=xn; n++) {
                                auto s = saa;
                                if(s>n) s = n;
                                for(int m=1; m<=s; m++) {
                                    A[a][a].coeff(Amaa, m, ctx);
                                    status |= gr_mat_mul(mat, Amaa, I[a][ila][k][n-m], ctx);
                                    status |= gr_mat_add(I[a][ila][k][n], I[a][ila][k][n], mat, ctx);
                                }
                                if(k+1<kmax) {
                                    status |= gr_set_fmpz(z, D0, ctx);
                                    status |= gr_mat_mul_scalar(mat, I[a][ila][k+1][n], z, ctx);
                                    status |= gr_mat_sub(I[a][ila][k][n], I[a][ila][k][n], mat, ctx);
                                }
                                s = sa;
                                if(s>n) s = n;
                                for(int m=1; m<=s; m++) {
                                    status |= gr_set_fmpq(z, qlas[ila], ctx);
                                    status |= gr_add_si(z, z, n-m, ctx); // la+n-m
                                    status |= gr_mat_mul_scalar(mat, I[a][ila][k][n-m], z, ctx);
                                    if(k+1<kmax) status |= gr_mat_add(mat, mat, I[a][ila][k+1][n-m], ctx);
                                    fmpz_poly_get_coeff_fmpz(Dm, D[a], m);
                                    fmpz_neg(Dm, Dm);
                                    status |= gr_set_fmpz(z, Dm, ctx);
                                    status |= gr_mat_addmul_scalar(I[a][ila][k][n], mat, z, ctx);
                                }
                                fmpq_add_si(q, qlas[ila], n); // q = la+n
                                fmpq_mul_fmpz(q, q, D0); // q = (la+n)D0
                                
                                fmpq_mat_one(invA);
                                fmpq_mat_scalar_mul_fmpq(invA, invA, q);

                                fmpq_mat_sub(invA, invA, A0aa);
                                int nok = fmpq_mat_inv(invA, invA);
                                if(!nok) {
                                    cout << "fmpq_mat_inv failed." << endl;
                                    abort();
                                }
                                status |= gr_mat_set_fmpq_mat(invA_acb, invA, ctx);
                                status |= gr_mat_mul(I[a][ila][k][n], invA_acb, I[a][ila][k][n], ctx);
                            }
                        }
                    }
                    fmpz_clear(Dm);
                    gr_heap_clear(z ,ctx);
                    gr_mat_clear(Amaa ,ctx);
                    fmpz_clear(D0);
                    fmpq_clear(q);
                    fmpq_mat_clear(A0aa);
                    fmpq_mat_clear(invA);
                    gr_mat_clear(invA_acb ,ctx);
                }
                gr_mat_clear(mat ,ctx);
                {
                    lock_guard<mutex> guard(worker_mutex);
                    worker_left -= worker_nn;
                    Used[a] = false;
                    if(to_notify) Done[a] = true;
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                \r" << flush;
                        cout << pre << "\\--nSeries I(x^" << xn << ") " << worker_total << "|" << (worker_total-worker_left) << flush;
                    }
                }
                if(to_notify) worker_cond.notify_all();
            }
        };

        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--nSeries I(x^" << xn << ") " << worker_total << "|0" << flush;
        int total_threads = Threads;
        if(total_threads<=0) total_threads = CpuCores();
        thread worker[total_threads];
        for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        for(int i=0; i<total_threads; ++i) worker[i].join();
        for(int a=0; a<nbs; a++) {
            fmpz_poly_clear(D[a]);
            flint_free(D[a]);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;

    }
    
    //=*********************************************************************=
    // Taylor expansion - rational
    //=*********************************************************************=
    
    void DEX::taylor(vector<vector<fmpq_mat_struct*>> & I, int xn, const matrix I0, const ex & x0) {
        if(gr_taylor_inited) throw Error("gr_taylor_inited = true");
        if(!x0.info(info_flags::rational)) throw Error("x0 is NOT rational.");
        auto nbs = bs.size();
        int nc = I0.cols();
        
        if(!taylor_inited) { // initialize QMat/QD, cache for later taylor call
            QMat.resize(nbs);
            QD.resize(nbs);
            #pragma omp parallel for schedule(dynamic,1)
            for(int br=0; br<nbs; br++) { // cycle rows
                QMat[br].resize(br+1);
                int nr = bs[br].second;
                vector<MQ> AX(br+1);
                vector<fmpz_poly_struct*> lcm_vec(br+1);
                fmpz_poly_t lcm; // D=lcm
                fmpz_poly_init(lcm);
                fmpz_poly_set_str(lcm, "1  1");
                for(int bc=br; bc>=0; bc--) { 
                    // lcm for each block A[br][bc] -> A[bc] 
                    AX[bc].init(Mat[br][bc]);
                    lcm_vec[bc] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                    fmpz_poly_init(lcm_vec[bc]);
                    AX[bc].denlcm(lcm_vec[bc]);
                    fmpz_poly_lcm(lcm, lcm, lcm_vec[bc]);
                }
                for(int bc=br; bc>=0; bc--) {
                    int nc = bs[bc].second;
                    QMat[br][bc].resize(nr);
                    for(int r=0; r<nr; r++) {
                        QMat[br][bc][r].resize(nc);
                        for(int c=0; c<nc; c++) {
                            QMat[br][bc][r][c] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                            fmpz_poly_init(QMat[br][bc][r][c]);
                        }
                    }
                    fmpz_poly_div(lcm_vec[bc], lcm, lcm_vec[bc]);
                    AX[bc].scale(lcm_vec[bc]);
                    fmpz_poly_clear(lcm_vec[bc]);
                    flint_free(lcm_vec[bc]);
                    AX[bc](QMat[br][bc]); // back to QMat
                    AX[bc].clear();
                }
                QD[br] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                fmpz_poly_init(QD[br]);
                fmpz_poly_set(QD[br],lcm);
                fmpz_poly_clear(lcm);
            }
            taylor_inited = true;
            
            if(taylor_clearq) { // clearq or not
                for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
                    for(auto & row : Mat[br][bc]) for(auto & item : row) {
                        fmpz_poly_q_clear(item);
                        flint_free(item);
                    }
                    if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) {
                        fmpq_mat_clear(item[0]);
                        flint_free(item[0]);
                    }
                }
                Mat.clear();
                if(fuchsified) { 
                    U0.clear(); 
                    for(auto & item : qlas) {
                        fmpq_clear(item);
                        flint_free(item);
                    }
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
            vector<vector<vector<fmpq_poly_struct*>>> TM(a+1); // TM[b][r][c] for each row a
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                TM[b].resize(nr);
                int s = -1;
                for(int r=0; r<nr; r++) {
                    TM[b][r].resize(nc);
                    for(int c=0; c<nc; c++) {
                        auto & item = TM[b][r][c];
                        item = (fmpq_poly_struct*) flint_malloc(sizeof(fmpq_poly_struct));
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
            fmpq_poly_set_fmpz_poly(lcm,QD[a]);
            fmpq_poly_taylor_shift(lcm,lcm,q0); // shift: x -> x+x0
            
            fmpq_t q,D0;
            fmpq_init(q);
            fmpq_init(D0);
            fmpq_poly_get_coeff_fmpq(D0,lcm,0);
            if(fmpq_is_zero(D0)) throw Error("taylor: D0 is zero.");
            
            I[a].resize(xn+1);
            I[a][0] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
            fmpq_mat_init(I[a][0],nr,nc); // call fmpq_mat_clear() after this call
            _to_(I[a][0],ex_to<matrix>(sub_matrix(I0, r0, nr, 0, nc)));
            
            fmpq_mat_t smat;
            fmpq_mat_init(smat,nr,nc);
            
            if(omp_get_active_level()==0 && !In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                              \r" << flush;
                cout << pre << "\\--taylor: " << nbs << "|" << a+1;
                cout << " [" << nr << "\u2A09" << nc << "] x^" << xn << flush;
            }
                
            for(int n=1; n<=xn; n++) { // 3-cycle over n
            
                fmpq_mat_zero(smat);
                
                if(true) {
                    slong s = fmpq_poly_degree(lcm);
                    if(s>n-1) s=n-1;
                    int tmax = omp_get_max_threads();
                    vector<fmpq_mat_struct*> mat_sum(tmax), mat_tmp(tmax);
                    for(int i=0; i<tmax; i++) {
                        mat_sum[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(mat_sum[i],nr,nc);
                        mat_tmp[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(mat_tmp[i],nr,nc);
                    }
                    #pragma omp parallel for schedule(runtime) num_threads(tmax)
                    for(int m=1; m<=s; m++) {
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
                    }
                    for(int i=0; i<tmax; i++) {
                        fmpq_mat_add(smat,smat,mat_sum[i]);
                        fmpq_mat_clear(mat_sum[i]);
                        flint_free(mat_sum[i]);
                        fmpq_mat_clear(mat_tmp[i]);
                        flint_free(mat_tmp[i]);
                    }
                }
          
                if(true) {
                    int tmax = omp_get_max_threads();
                    vector<fmpq_mat_struct*> mat_sum(tmax), mat_tmp(tmax);
                    for(int i=0; i<tmax; i++) {
                        mat_sum[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                        fmpq_mat_init(mat_sum[i],nr,nc);
                        mat_tmp[i] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
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
                    }
                    for(int i=0; i<tmax; i++) {
                        fmpq_mat_add(smat,smat,mat_sum[i]);
                        fmpq_mat_clear(mat_sum[i]);
                        flint_free(mat_sum[i]);
                        fmpq_mat_clear(mat_tmp[i]);
                        flint_free(mat_tmp[i]);
                    }
                }
                
                fmpq_mul_si(q,D0,n); // q = n D0
                fmpq_inv(q,q); // q = 1 / (n D0)
                I[a][n] = (fmpq_mat_struct*) flint_malloc(sizeof(fmpq_mat_struct));
                fmpq_mat_init(I[a][n],nr,nc); // call fmpq_mat_clear() after this call
                fmpq_mat_scalar_mul_fmpq(I[a][n],smat,q);
            }
                
            fmpq_poly_clear(lcm);
            fmpq_clear(q);
            fmpq_clear(D0);
            fmpq_mat_clear(smat);
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    fmpq_poly_clear(TM[b][r][c]);
                    flint_free(TM[b][r][c]);
                }
            }
        }
        fmpq_clear(q0); 
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
    }
    
    
    struct t_item_type {
        int a, n;
        bool to_notify = false;
        list<int> ms;
        list<pair<int,int>> bms;
        t_item_type(int a_, int n_, bool to_notify_, list<int> && ms_, list<pair<int,int>> && bms_) : a(a_), n(n_), to_notify(to_notify_), ms(std::move(ms_)), bms(std::move(bms_)) { }
        t_item_type() { }
    };
    
    
    //=*********************************************************************=
    // Taylor expansion - Parallel all: a & n, acf version
    //=*********************************************************************=
    
    void DEX::an_taylor(vector<vector<gr_mat_struct*>> & I, int xn, gr_mat_t imat, gr_ptr z0, gr_ctx_t ctx, const string & es) {

        if(taylor_inited) throw Error("taylor_inited = true");
        auto nbs = bs.size();
        int nc = gr_mat_ncols(imat, ctx);

        if(!gr_taylor_inited) { // initialize GMat/GD, cache for later taylor call
            slong fp;
            int status = GR_SUCCESS;
            status |= gr_ctx_get_real_prec(&fp, ctx);
            if(fp>1000) fp += 500;
            else fp += fp/2;
            gr_ctx_init_complex_float_acf(_ctx_, fp);
            GMat.resize(nbs);
            GD.resize(nbs);
            #pragma omp parallel for schedule(dynamic,1)
            for(int br=0; br<nbs; br++) { // cycle rows
                fmpz_poly_t lcm; // D=lcm
                GMat[br].resize(br+1);
                int nr = bs[br].second;
                vector<MQ> AX(br+1);
                vector<fmpz_poly_struct*> lcm_vec(br+1);
                fmpz_poly_init(lcm);
                fmpz_poly_set_str(lcm, "1  1");
                for(int bc=br; bc>=0; bc--) {
                    // lcm for each block A[br][bc] -> A[bc]
                    AX[bc].init(Mat[br][bc]);
                    lcm_vec[bc] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                    fmpz_poly_init(lcm_vec[bc]);
                    AX[bc].denlcm(lcm_vec[bc]);
                    fmpz_poly_lcm(lcm, lcm, lcm_vec[bc]);
                }
                for(int bc=br; bc>=0; bc--) {
                    int nc = bs[bc].second;
                    GMat[br][bc].resize(nr);
                    for(int r=0; r<nr; r++) {
                        GMat[br][bc][r].resize(nc);
                        for(int c=0; c<nc; c++) {
                            GMat[br][bc][r][c] = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
                            gr_poly_init(GMat[br][bc][r][c], ctx);
                        }
                    }
                    fmpz_poly_div(lcm_vec[bc], lcm, lcm_vec[bc]);
                    AX[bc].scale(lcm_vec[bc]);
                    fmpz_poly_clear(lcm_vec[bc]);
                    flint_free(lcm_vec[bc]);
                    AX[bc](GMat[br][bc], _ctx_); // back to GMat
                    AX[bc].clear();
                }
                GD[br] = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
                gr_poly_init(GD[br], _ctx_);
                status |= gr_poly_set_fmpz_poly(GD[br], lcm, _ctx_);
                fmpz_poly_clear(lcm);
            }
            if(status != GR_SUCCESS) {
                cout << "an_taylor: status is Not GR_SUCCESS!" << endl;
                abort();
            }
            gr_taylor_inited = true;
            
            if(taylor_clearq) { // clearq or not
                for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
                    for(auto & row : Mat[br][bc]) for(auto & item : row) fmpz_poly_q_clear(item);
                    if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) {
                        fmpq_mat_clear(item[0]);
                        flint_free(item[0]);
                    }
                }
                Mat.clear();
                if(fuchsified) {
                    U0.clear();
                    for(auto & item : qlas) {
                        fmpq_clear(item);
                        flint_free(item);
                    }
                    qlas.clear();
                }
                fuchsified = false;
            }
        }
        
        I.resize(nbs);
        
        // Parallel Data
        vector<vector<vector<vector<gr_poly_struct*>>>> GM(nbs); // GM[a][b][r][c]
        vector<gr_poly_struct*> D(nbs); // Denominator, i.e., lcm
        vector<vector<bool>> Done(nbs); // [a][n] : ready or not
        vector<vector<bool>> Used(nbs); // [a][n] : used or not
        vector<vector<list<int>>> S1(nbs); // [a][n] : { m, ... }
        vector<vector<list<pair<int,int>>>> S2(nbs); // [a][n] : { <b,m>, ... }
        unsigned long long worker_left = 0;
        
        // init data
        for(int a=0; a<nbs; a++) {
            int r0 = bs[a].first;
            int nr = bs[a].second;
            GM[a].resize(a+1);
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                GM[a][b].resize(nr);
                for(int r=0; r<nr; r++) {
                    GM[a][b][r].resize(nc);
                }
            }
        }
        
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) {
            for(int b=0; b<=nbs; b++) {
                if(b>a) continue;
                int nr = bs[a].second;
                int nc = bs[b].second;
                for(int r=0; r<nr; r++) {
                    for(int c=0; c<nc; c++) {
                        auto & item = GM[a][b][r][c];
                        item = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
                        gr_poly_init(item, ctx);
                        gr_poly_taylor_shift(item, GMat[a][b][r][c], z0, ctx); // shift: x -> x+x0
                    }
                }
            }
        }
        
        #pragma omp parallel for schedule(dynamic,1) reduction(+: worker_left)
        for(int a=0; a<nbs; a++) { // cycle rows
            int r0 = bs[a].first;
            int nr = bs[a].second;
            // taylor shift: x -> x+x0
            vector<int> sdeg(a+1);
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                int s = -1;
                for(int r=0; r<nr; r++) {
                    for(int c=0; c<nc; c++) {
                        int ss = gr_poly_length(GM[a][b][r][c], ctx)-1;
                        if(ss>s) s = ss;
                    }
                }
                sdeg[b] = s;
            }
            
            D[a] = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
            gr_poly_init(D[a], ctx);
            gr_poly_taylor_shift(D[a], GD[a], z0, ctx); // shift: x -> x+x0
            
            I[a].resize(xn+1);
            I[a][0] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
            gr_mat_init(I[a][0], nr, nc, ctx); // call gr_mat_clear() after this call
            
            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                gr_set(gr_mat_entry_ptr(I[a][0], r, c, ctx), gr_mat_entry_ptr(imat, r0+r, c, ctx), ctx);
            }
                    
            Done[a].resize(xn+1);
            Done[a][0] = true;
            Used[a].resize(xn+1);
            Used[a][0] = false;
            S1[a].resize(xn+1);
            S2[a].resize(xn+1);

            for(int n=1; n<=xn; n++) { // cycle over n
                Done[a][n] = true;
                Used[a][n] = false;
                I[a][n] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                gr_mat_init(I[a][n], nr, nc, ctx); // call gr_mat_clear() after this call
                
                slong s = gr_poly_length(D[a], ctx) - 1;
                if(s>n-1) s=n-1;
                for(int m=1; m<=s; m++) {
                    S1[a][n].push_back(m);
                    worker_left++;
                }
                for(int b=0; b<=a; b++) {
                    slong s = sdeg[b];
                    if(s>n-1) s=n-1;
                    for(int m=0; m<=s; m++) {
                        S2[a][n].push_back(make_pair(b,m));
                        worker_left++;
                    }
                }
                if(!S1[a][n].empty() || !S2[a][n].empty()) Done[a][n] = false;
            }
        }
        
        mutex worker_mutex;
        condition_variable worker_cond;
        
        int total_threads = Threads;
        if(total_threads<=0) total_threads = CpuCores();
        auto worker_total = worker_left;
        list<t_item_type> worker_item_list;
        auto worker_main = [&]() {
            while(true) {
                list<t_item_type> wi_list;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(!worker_left) return true;
                        start:
                        int cnt = 0;
                        auto tot = worker_item_list.size();
                        while(tot>0) {
                            auto itr = worker_item_list.begin();
                            wi_list.emplace_back(std::move(*itr));
                            worker_item_list.erase(itr);
                            cnt++;
                            tot--;
                            if(cnt>10) break;
                        }
                        if(cnt) return true;
                        for(int a=0; a<nbs; a++) { // cycle over row a
                            for(int n=1; n<=xn; n++) { // cycle over n
                                if(Used[a][n] || Done[a][n]) continue;
                                
                                list<int> ms;
                                list<pair<int,int>> bms;
                                
                                list<int> & s1 = S1[a][n];
                                for(auto itr = s1.begin(); itr != s1.end(); ) {
                                    int m = *itr;
                                    if(Done[a][n-m]) {
                                        itr = s1.erase(itr);
                                        ms.push_back(m);
                                    } else ++itr;
                                }
                                
                                list<pair<int,int>> & s2 = S2[a][n];
                                for(auto itr = s2.begin(); itr != s2.end(); ) {
                                    int b = itr->first;
                                    int m = itr->second;
                                    if(Done[b][n-1-m]) {
                                        itr = s2.erase(itr);
                                        bms.push_back(make_pair(b,m));
                                    } else ++itr;
                                }
                                
                                if(!ms.empty() || !bms.empty()) {
                                    Used[a][n] = true;
                                    bool to_notify = false;
                                    if(s1.empty() && s2.empty()) to_notify = true;
                                    worker_item_list.emplace_back(a,n,to_notify,std::move(ms),std::move(bms));
                                }
                            }
                        }
                        if(worker_item_list.size()==0) return false;
                        goto start; // for return to check the start label
                    });
                    if(!worker_left) {
                        flint_cleanup();
                        return;
                    }
                }
                
                int status = GR_SUCCESS;
                unsigned long long worker_nn = 0;
                for(auto const & worker_item : wi_list) {
                    int a = worker_item.a;
                    int n = worker_item.n;
                    bool to_notify = worker_item.to_notify;
                    auto & ms = worker_item.ms;
                    auto & bms = worker_item.bms;
                    
                    gr_mat_t mat;
                    int nr = bs[a].second;
                    gr_mat_init(mat, nr, nc, ctx);
                    gr_ptr Dm = gr_heap_init(ctx);
                    for(int m : ms) {
                        status |= gr_poly_get_coeff_scalar(Dm, D[a], m, ctx);
                        status |= gr_mul_si(Dm, Dm, n-m, ctx);
                        status |= gr_neg(Dm, Dm, ctx);
                        status |= gr_mat_mul_scalar(mat, I[a][n-m], Dm, ctx);
                        status |= gr_mat_add(I[a][n], I[a][n], mat, ctx);
                        ++worker_nn;
                    }
                    gr_heap_clear(Dm, ctx);
                    for(auto kv : bms) {
                        auto b = kv.first;
                        auto m = kv.second;
                        int nc2 = bs[b].second;
                        gr_mat_t Amab;
                        gr_mat_init(Amab, nr, nc2, ctx);
                        for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) { // coefficients
                            status |= gr_poly_get_coeff_scalar(gr_mat_entry_ptr(Amab,r,c,ctx), GM[a][b][r][c], m, ctx);
                        }
                        status |= gr_mat_mul(mat, Amab, I[b][n-1-m], ctx);
                        gr_mat_clear(Amab, ctx);
                        status |= gr_mat_add(I[a][n], I[a][n], mat, ctx);
                        ++worker_nn;
                    }
                    gr_mat_clear(mat, ctx);
                    if(to_notify) {
                        gr_ptr z = gr_heap_init(ctx);
                        status |= gr_poly_get_coeff_scalar(z, D[a], 0, ctx);
                        status |= gr_mul_si(z, z, n, ctx); // n D0
                        status |= gr_inv(z, z, ctx); // 1 / (n D0)
                        status |= gr_mat_mul_scalar(I[a][n], I[a][n], z, ctx);
                        gr_heap_clear(z, ctx);
                    }
                }
                if(status != GR_SUCCESS) {
                    cout << "an_taylor: status is Not GR_SUCCESS!" << endl;
                    abort();
                }
                bool any_notify = false;
                {
                    lock_guard<mutex> guard(worker_mutex);
                    for(auto const & worker_item : wi_list) {
                        int a = worker_item.a;
                        int n = worker_item.n;
                        bool to_notify = worker_item.to_notify;
                        if(!any_notify && to_notify) any_notify = true;
                        Used[a][n] = false;
                        if(to_notify) Done[a][n] = true;
                    }
                    worker_left -= worker_nn;
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                \r" << flush;
                        cout << pre << "\\--nTaylor" << es << " I(x^" << xn << ") " << worker_total << "|" << (worker_total-worker_left) << flush;
                    }
                }
                if(any_notify) worker_cond.notify_all();
            }
        };

        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--nTaylor" << es << " I(x^" << xn << ") " << worker_total << "|0" << flush;
        thread worker[total_threads];
        for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        for(int i=0; i<total_threads; ++i) worker[i].join();
        for(int a=0; a<nbs; a++) {
            gr_poly_clear(D[a], ctx);
            flint_free(D[a]);
        }
        
        for(auto & i1 : GM) for(auto & i2 : i1) for(auto & i3 : i2) for(auto & i4 : i3) {
            gr_poly_clear(i4, ctx);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
        
    }
    
    //=*********************************************************************=
    // Taylor expansion - acf - Parallel - only on a, note acf version
    //=*********************************************************************=
    
    void DEX::a_taylor(vector<vector<gr_mat_struct*>> & I, int xn, gr_mat_t imat, gr_ptr z0, gr_ctx_t ctx, const string & es) {

        if(taylor_inited) throw Error("taylor_inited = true");
        auto nbs = bs.size();
        int nc = gr_mat_ncols(imat, ctx);
        
        if(!gr_taylor_inited) { // initialize GMat/GD, cache for later taylor call
            slong fp;
            int status = GR_SUCCESS;
            status |= gr_ctx_get_real_prec(&fp, ctx);
            if(fp>1000) fp += 500;
            else fp += fp/2;
            gr_ctx_init_complex_float_acf(_ctx_, fp);
            GMat.resize(nbs);
            GD.resize(nbs);
            #pragma omp parallel for schedule(dynamic,1)
            for(int br=0; br<nbs; br++) { // cycle rows
                fmpz_poly_t lcm; // D=lcm
                GMat[br].resize(br+1);
                int nr = bs[br].second;
                vector<MQ> AX(br+1);
                vector<fmpz_poly_struct*> lcm_vec(br+1);
                fmpz_poly_init(lcm);
                fmpz_poly_set_str(lcm, "1  1");
                for(int bc=br; bc>=0; bc--) {
                    // lcm for each block A[br][bc] -> A[bc]
                    AX[bc].init(Mat[br][bc]);
                    lcm_vec[bc] = (fmpz_poly_struct*) flint_malloc(sizeof(fmpz_poly_struct));
                    fmpz_poly_init(lcm_vec[bc]);
                    AX[bc].denlcm(lcm_vec[bc]);
                    fmpz_poly_lcm(lcm, lcm, lcm_vec[bc]);
                }
                for(int bc=br; bc>=0; bc--) {
                    int nc = bs[bc].second;
                    GMat[br][bc].resize(nr);
                    for(int r=0; r<nr; r++) {
                        GMat[br][bc][r].resize(nc);
                        for(int c=0; c<nc; c++) {
                            GMat[br][bc][r][c] = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
                            gr_poly_init(GMat[br][bc][r][c], ctx);
                        }
                    }
                    fmpz_poly_div(lcm_vec[bc], lcm, lcm_vec[bc]);
                    AX[bc].scale(lcm_vec[bc]);
                    fmpz_poly_clear(lcm_vec[bc]);
                    flint_free(lcm_vec[bc]);
                    AX[bc](GMat[br][bc], _ctx_); // back to GMat
                    AX[bc].clear();
                }
                GD[br] = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
                gr_poly_init(GD[br], _ctx_);
                status |= gr_poly_set_fmpz_poly(GD[br], lcm, _ctx_);
                fmpz_poly_clear(lcm);
            }
            if(status != GR_SUCCESS) {
                cout << "a_taylor: status is Not GR_SUCCESS!" << endl;
                abort();
            }
            gr_taylor_inited = true;
            
            if(taylor_clearq) { // clearq or not
                for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
                    for(auto & row : Mat[br][bc]) for(auto & item : row) {
                        fmpz_poly_q_clear(item); // fmpz_poly_q_init() in initialize
                        flint_free(item);
                    }
                    if(fuchsified) for(auto & kv : U0[br][bc]) for(auto & item : kv.second) {
                        fmpq_mat_clear(item[0]); // fmpq_mat_init() in initialize
                        flint_free(item[0]);
                    }
                }
                Mat.clear();
                if(fuchsified) {
                    U0.clear();
                    for(auto & item : qlas) {
                        fmpq_clear(item); // fmpq_init() in initialize
                        flint_free(item);
                    }
                    qlas.clear();
                }
                fuchsified = false;
            }
        }
        
        // Parallel Data
        vector<vector<vector<vector<gr_poly_struct*>>>> GM(nbs); // GM[a][b]
        vector<vector<slong>> sM(nbs); // degree for A[a][b]
        vector<gr_poly_struct*> D(nbs); // Denominator, i.e., lcm
        vector<bool> Done(nbs); // [a] : ready or not
        vector<bool> Used(nbs); // [a] : used or not
        vector<list<int>> Sb(nbs); // [a] : { b, ... }
        unsigned long long worker_left = 0;
        
        // init data
        I.resize(nbs);
        for(int a=0; a<nbs; a++) {
            int r0 = bs[a].first;
            int nr = bs[a].second;
            GM[a].resize(a+1);
            sM[a].resize(a+1);
            for(int b=0; b<=a; b++) {
                int nc = bs[b].second;
                GM[a][b].resize(nr);
                for(int r=0; r<nr; r++) {
                    GM[a][b][r].resize(nc);
                }
            }
        }
        
        #pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(int a=0; a<nbs; a++) {
            for(int b=0; b<=nbs; b++) {
                if(b>a) continue;
                int nr = bs[a].second;
                int nc = bs[b].second;
                slong s = -1;
                for(int r=0; r<nr; r++) {
                    for(int c=0; c<nc; c++) {
                        auto & item = GM[a][b][r][c];
                        item = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
                        gr_poly_init(item, ctx);
                        gr_poly_set(item, GMat[a][b][r][c], ctx);
                        gr_poly_taylor_shift(item, item, z0, ctx); // shift: x -> x+x0
                        auto ss = gr_poly_length(item, ctx) - 1;
                        if(ss>s) s = ss;
                    }
                }
                sM[a][b] = s;
            }
        }
        
        #pragma omp parallel for schedule(dynamic,1) reduction(+: worker_left)
        for(int a=0; a<nbs; a++) { // cycle rows
            int r0 = bs[a].first;
            int nr = bs[a].second;
            D[a] = (gr_poly_struct*) flint_malloc(sizeof(gr_poly_struct));
            gr_poly_init(D[a], ctx);
            gr_poly_set(D[a], GD[a], ctx);
            gr_poly_taylor_shift(D[a], D[a], z0, ctx); // shift: x -> x+x0
            I[a].resize(xn+1);
            I[a][0] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
            gr_mat_init(I[a][0], nr, nc, ctx);
            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                gr_set(gr_mat_entry_ptr(I[a][0],r,c,ctx), gr_mat_entry_ptr(imat,r0+r,c,ctx), ctx);
            }
            for(int n=1; n<=xn; n++) {
                I[a][n] = (gr_mat_struct*) flint_malloc(sizeof(gr_mat_struct));
                gr_mat_init(I[a][n], nr, nc, ctx);
                gr_mat_zero(I[a][n], ctx);
            }
            Done[a] = false;
            Used[a] = false;
            for(int b=0; b<a; b++) Sb[a].push_back(b); // note b<a here
            worker_left += a+1;
        }
        
        mutex worker_mutex;
        condition_variable worker_cond;
        
        auto worker_total = worker_left;
        auto worker_main = [&]() {
            while(true) {
                int a;
                list<int> bset;
                bool to_notify = false;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(!worker_left) return true;
                        for(a=0; a<nbs; a++) {
                            if(Used[a] || Done[a]) continue;
                                
                            list<int> & sb = Sb[a];
                            for(auto itr = sb.begin(); itr != sb.end(); ) {
                                int b = *itr;
                                if(Done[b]) {
                                    itr = sb.erase(itr);
                                    bset.push_back(b);
                                } else ++itr;
                            }
                            
                            if(!bset.empty() || a==0) { // note a==0
                                Used[a] = true;
                                if(sb.empty()) to_notify = true;
                                return true;
                            }
                        }
                        return false;
                    });
                    if(!worker_left) {
                        flint_cleanup();
                        return;
                    }
                }
                
                int status = GR_SUCCESS;
                int nr = bs[a].second;
                gr_mat_t mat;
                gr_mat_init(mat, nr, nc, ctx);
                for(auto b : bset) {
                    int nc2 = bs[b].second;
                    gr_mat_t Amab;
                    gr_mat_init(Amab, nr, nc2, ctx);
                    int s = sM[a][b];
                    if(s>xn-1) s = xn-1;
                    for(int m=0; m<=s; m++) {
                        for(int r=0; r<nr; r++) for(int c=0; c<nc2; c++) { // coefficients
                            status |= gr_poly_get_coeff_scalar(gr_mat_entry_ptr(Amab,r,c,ctx), GM[a][b][r][c], m, ctx);
                        }
                        for(int n=xn; n>0 && n-1-m>=0; n--) {
                            status |= gr_mat_mul(mat, Amab, I[b][n-1-m], ctx);
                            status |= gr_mat_add(I[a][n], I[a][n], mat, ctx);
                        }
                    }
                    gr_mat_clear(Amab, ctx);
                }
                unsigned long long worker_nn = bset.size();
                if(to_notify) {
                    worker_nn++;
                    gr_mat_t Amaa;
                    gr_mat_init(Amaa, nr, nr, ctx);
                    gr_ptr Dm = gr_heap_init(ctx);
                    for(int n=1; n<=xn; n++) {
                        auto s = sM[a][a];
                        if(s>n-1) s = n-1;
                        for(int m=0; m<=s; m++) {
                            for(int r=0; r<nr; r++) for(int c=0; c<nr; c++) { // coefficients
                                status |= gr_poly_get_coeff_scalar(gr_mat_entry_ptr(Amaa,r,c,ctx), GM[a][a][r][c], m, ctx);
                            }
                            status |= gr_mat_mul(mat, Amaa, I[a][n-1-m], ctx);
                            status |= gr_mat_add(I[a][n], I[a][n], mat, ctx);
                        }
                        
                        s = gr_poly_length(D[a], ctx) - 1;
                        if(s>n-1) s = n-1;
                        for(int m=1; m<=s; m++) {
                            status |= gr_poly_get_coeff_scalar(Dm, D[a], m, ctx);
                            status |= gr_mul_si(Dm, Dm, n-m, ctx);
                            status |= gr_neg(Dm, Dm, ctx);
                            status |= gr_mat_addmul_scalar(I[a][n], I[a][n-m], Dm, ctx);
                        }
                        
                        status |= gr_poly_get_coeff_scalar(Dm, D[a], 0, ctx);
                        status |= gr_mul_si(Dm, Dm, n, ctx); // n D0
                        status |= gr_inv(Dm, Dm, ctx); // 1 / (n D0)
                        status |= gr_mat_mul_scalar(I[a][n], I[a][n], Dm, ctx);
                    }
                    gr_heap_clear(Dm, ctx);
                    gr_mat_clear(Amaa, ctx);
                }
                gr_mat_clear(mat, ctx);
                {
                    lock_guard<mutex> guard(worker_mutex);
                    worker_left -= worker_nn;
                    Used[a] = false;
                    if(to_notify) Done[a] = true;
                    if(!In_GiNaC_Parallel && Verbose>5) {
                        cout << "\r                                \r" << flush;
                        cout << pre << "\\--nTaylor" << es << " I(x^" << xn << ") " << worker_total << "|" << (worker_total-worker_left) << flush;
                    }
                }
                if(to_notify) worker_cond.notify_all();
                                
            }
        };

        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << "\\--nTaylor" << es << " I(x^" << xn << ") " << worker_total << "|0" << flush;
        int total_threads = Threads;
        if(total_threads<=0) total_threads = CpuCores();
        thread worker[total_threads];
        for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        for(int i=0; i<total_threads; ++i) worker[i].join();
        for(int a=0; a<nbs; a++) {
            gr_poly_clear(D[a], ctx);
            flint_free(D[a]);
        }

        for(auto & i1 : GM) for(auto & i2 : i1)
        for(auto & i3 : i2) for(auto & i4 : i3) {
            gr_poly_clear(i4, ctx);
            flint_free(i4);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << " @ " << now(false) << endl;
    }
    
    //=*********************************************************************=
        
}

