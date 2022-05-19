/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"
#include "vspace.h"

namespace HepLib {

    ev_am_t ev_am(const matrix & mat) { return ED::ev_am(mat); }
    pair<matrix,jcf_t> jordan(const matrix &m) { return ED::jordan(m); }
    pair<matrix,jcf_t> jordan(const matrix &m, ev_am_t & ev2am) { return ED::jordan(m,ev2am); }
    matrix proj_mat(const pair<matrix,matrix> & a01) { return ED::proj_mat(a01.first,a01.second); }
    matrix normal(const matrix & mat) { return ED::normal(mat); }
    matrix exnormal(const matrix & mat) { 
        auto m = mat;
        GiNaC_Parallel_Verb["NM"] = 0;
        auto res = GiNaC_Parallel(m.nops(), [&m](int idx) {
            return exnormal(m.op(idx));
        }, "NM");
        for(unsigned i=0; i<m.nops(); i++) m.let_op(i) = res[i];
        return m;
    }
    
    int prank(const matrix & mat, const symbol &x) { // each item has been normalized -> N/D
        if(mat.is_zero_matrix()) return -100000;
        int pr = -100000;
        for(int i=0; i<mat.nops(); i++) {
            auto nd = mat.op(i).numer_denom();
            auto p = nd.op(1).ldegree(x)-nd.op(0).ldegree(x)-1;
            if(pr<p) pr = p;
        }
        return pr;
    }
    
    ex xpow(const ex & e, const ex & x) {
        auto cvs = collect_lst(e,x);
        ex res = 0;
        for(auto cv : cvs) {
            auto pat = cv.op(1);
            if(!is_a<mul>(pat)) pat = lst{pat};
            ex cp = 1;
            ex xn = 0;
            for(auto pi : pat) {
                if(pi.match(pow(x,w))) xn += pi.op(1);
                else cp *= pi;
            }
            res += cv.op(0) * cp * pow(x,xn);
        }
        return res;
    }
    
    void xpow(matrix & mat, const ex & x) {
        auto n = mat.nops();
        for(int i=0; i<n; i++) mat.let_op(i) = xpow(mat.op(i),x);
    }
    
    void subs(matrix & mat, const ex & s, unsigned opt) {
        auto n = mat.nops();
        for(int i=0; i<n; i++) mat.let_op(i) = mat.op(i).subs(s, opt);
    }
    
    matrix a0_mat(const matrix &mat, const symbol &x, int pr) {
        int p = (pr == 19790923 ? prank(mat, x) : pr);
        int nr=mat.rows(), nc=mat.cols();
        matrix a0(nr,nc);
        for(int i=0; i<mat.nops(); i++) {
            auto nd = mat.op(i).numer_denom();
            auto nl = nd.op(0).ldegree(x);
            auto dl = nd.op(1).ldegree(x);
            if(dl-nl-1>pr) throw Error("input pr > prank!");
            if(dl-nl-1<pr) a0.let_op(i) = 0;
            else a0.let_op(i) = nd.op(0).coeff(x,nl) / nd.op(1).coeff(x,dl);
        }
        return a0;
    }
    
    pair<matrix,matrix> a01_mat(const matrix &mat, const symbol &x, int pr) {
        int p = (pr == 19790923 ? prank(mat, x) : pr);
        int nr=mat.rows(), nc=mat.cols();
        matrix a0(nr,nc), a1(nr,nc);
        for(int i=0; i<mat.nops(); i++) {
            auto nd = mat.op(i).numer_denom();
            auto nl = nd.op(0).ldegree(x);
            auto dl = nd.op(1).ldegree(x);
            if(dl-nl-1>pr) throw Error("input pr > prank!");
            if(dl-nl<pr) {
                a0.let_op(i) = 0;
                a1.let_op(i) = 0;
            } else if(dl-nl-1<pr) { //dl-nl=pr
                a0.let_op(i) = 0;
                a1.let_op(i) = nd.op(0).coeff(x,nl) / nd.op(1).coeff(x,dl);
            } else { //dl-nl-1=pr
                auto n0 = nd.op(0).coeff(x,nl);
                auto n1 = nd.op(0).coeff(x,nl+1);
                auto d0 = nd.op(1).coeff(x,dl);
                auto d1 = nd.op(1).coeff(x,dl+1);
                a0.let_op(i) = n0/d0;
                a1.let_op(i) = (d0*n1-d1*n0)/(d0*d0);
            }
        }
        return make_pair(a0, a1);
    }
    
    matrix matrix_diff(const matrix & m, const symbol &x) {
        auto mat = m;
        for(int i=0; i<mat.nops(); i++) mat.let_op(i) = mat.op(i).diff(x);
        return mat;
    }

    // Dx J = M.J ---> J = T.J' & Dx J' = M'.J' with M' = Ti.M.T - Ti.Dx T
    matrix transform(const matrix &m, const matrix &t, const symbol &x) {
        return t.inverse().mul(m.mul(t).sub(matrix_diff(t,x)));
    }
    
    matrix x2y(const matrix & mat, const ex & y, const symbol & x) {
        static symbol xx("xx");
        ex det = diff_ex(y,x);
        auto m = mat;
        for(int i=0; i<m.nops(); i++) {
            m.let_op(i) = det * m.op(i).subs(x==xx,nopat).subs(xx==y,nopat);
        }
        return m;
    }
    
    //=*********************************************************************=
    
    DE::DE(const symbol & _x) : x(_x) { }
    
    // init to lower triangle block matrix
    void DE::init(const matrix & m) {
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
        Mat = pmat.inverse().mul(m).mul(pmat);
        Mat = exnormal(Mat); // make sure M is normalized
        Ts.push_back(pmat);
        N = n;
    }
    
    vector<vector<matrix>> DE::m2b(const matrix & m) { // m is low triangular block w.r.t. bs
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
    
    matrix DE::b2m(vector<vector<matrix>> & bm) {
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
    
    vector<matrix> DE::c2b(const matrix & m) {
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
    
    matrix DE::b2c(vector<matrix> & bc) {
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
    
    map<ex,vector<vector<matrix>>,ex_is_less> DE::b2m(const block_umat_t & bu) { // U[la][k][n]
        map<ex,vector<vector<matrix>>,ex_is_less> umat;
        int nbs = bs.size();
        for(int br=0; br<nbs; br++) for(int bc=0; bc<=br; bc++) {
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
    
    matrix DE::b2m(block_umat_t & bu, const ex & x) {
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
                            mat(r0+r,c0+c) += kv.second[k][n](r,c) * pow(x,la+n) * pow(log(x),k)/factorial(k);
                        }
                    }
                }
            }
        }
        return mat;
    }
    
    matrix DE::b2m(block_imat_t & bi, const ex & x) {
        int nbs = bs.size();
        int nc = bi[0].begin()->second[0][0].cols();
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
                            mat(r0+r,c) += kv.second[k][n](r,c) * pow(x,la+n) * pow(log(x),k)/factorial(k);
                        }
                    }
                }
            }
        }
        return mat;
    }
    
    matrix DE::T() {
        matrix t = ex_to<matrix>(unit_matrix(N));
        for(auto ti : Ts) t = t.mul(ti);
        return t;
    }
    
    matrix DE::M() { return Mat; }
    
    void DE::fuchsify() { // first on diagonal blocks, then off diagnoal ones
        matrix t(N,N);
        auto nbs = bs.size();
        for(int bi=0; bi<nbs; bi++) { // diagonal blocks
            if(!In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                    \r" << flush;
                cout << "  \\--reducing diagonal blocks: " << nbs << "|" << bi+1 << flush;
            }
            auto n0 = bs[bi].first;
            auto n = bs[bi].second;
            matrix m = ex_to<matrix>(sub_matrix(Mat, n0, n, n0, n));
            m = exnormal(m);
            matrix ti = ex_to<matrix>(unit_matrix(n,n));
            while(true) { // fuchsify diagonal blocks
                auto pr = prank(m,x);
                if(pr<1) break;
                auto p = proj_mat(a01_mat(m,x,pr));
                matrix cop = ex_to<matrix>(unit_matrix(n)).sub(p);
                auto px = p.mul_scalar(1/x);
                m = cop.sub(p.mul_scalar(x)).mul(m).mul(cop.sub(px)).add(px);
                m = exnormal(m);
                ti = ti.mul(cop.sub(px));
                ti = exnormal(ti);
            }
            while(true) { // shearing transformation on diagonal blocks
                auto a0 = a0_mat(m,x,0);
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
                m = transform(m,tm,x);
                m = exnormal(m);
                ti = ti.mul(tm);
                if(!shearing) break;
            }
            for(int r=0; r<n; r++) for(int c=0; c<n; c++) t(n0+r,n0+c) = ti(r,c);
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        t = exnormal(t);
        Ts.push_back(t);
        Mat = transform(Mat,t,x);
        Mat = exnormal(Mat);
        
        for(int br=0; br<nbs; br++) { // off-diagonal blocks
            for(int bc=br-1; bc>=0; bc--) {
                auto nr=bs[br].second, nc=bs[bc].second;
                if(!In_GiNaC_Parallel && Verbose>5) {
                    cout << "\r                                                    \r" << flush;
                    cout << "  \\--reducing off-diagonal blocks: " << nbs << "|" << br+1 << flush;
                    cout << " " << br << "|" << (br-bc) << flush; 
                } 
                matrix mrc = ex_to<matrix>(sub_matrix(Mat, bs[br].first, nr, bs[bc].first, nc));
                matrix mrr = ex_to<matrix>(sub_matrix(Mat, bs[br].first, nr, bs[br].first, nr));
                matrix mcc = ex_to<matrix>(sub_matrix(Mat, bs[bc].first, nc, bs[bc].first, nc));
                auto a0_rr = a0_mat(mrr,x,0);
                auto a0_cc = a0_mat(mcc,x,0);
                matrix sdmat = ex_to<matrix>(symbolic_matrix(nr,nc,"x"));
                matrix trc(nr,nc);
                while(true) {
                    auto pr = prank(mrc,x);
                    if(pr<1) break;
                    auto a0_rc = a0_mat(mrc,x,pr);
                    auto eq_mat = sdmat.mul_scalar(pr).add(a0_rr.mul(sdmat)).sub(sdmat.mul(a0_cc));
                    lst eqs, xs;
                    for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                        xs.append(sdmat(r,c));
                        eqs.append(eq_mat(r,c)+a0_rc(r,c)==0);
                    }
                    auto dmat = ex_to<matrix>(subs(sdmat,lsolve(eqs,xs))).mul_scalar(pow(x,-pr));
                    mrc = mrc.add(dmat.mul_scalar(pr/x)).add(mrr.mul(dmat)).sub(dmat.mul(mcc));
                    mrc = exnormal(mrc);
                    trc = trc.add(dmat);
                }
                if(!trc.is_zero_matrix()) {
                    t = ex_to<matrix>(unit_matrix(N,N));
                    for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                        t(bs[br].first+r, bs[bc].first+c) = trc(r,c);
                    }
                    Ts.push_back(t);
                    Mat = transform(Mat,t,x);
                    Mat = exnormal(Mat);
                }
            }
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;

        // matrix exponetial, note ln^k x --> ln^k x/k!
        auto a0 = a0_mat(Mat,x,0);
        auto qj = jordan(a0);
        symbol lnx("lnx");
        exmap la2s; // la 2 symbol
        matrix m(N,N);
        int cur_pos = 0;
        for(auto kv : qj.second) {
            if(la2s.find(kv.first)==la2s.end()) {
                symbol s;
                la2s[kv.first] = s;
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
        m = exnormal(m);
        
        // initialize UK[a][b][la] & U0[a][b][la][k]
        UK.resize(nbs);
        IK.resize(nbs);
        U0.resize(nbs);
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
                    if(br>0 && br-1>=bc && UK[br-1][bc].find(la)!=UK[br-1][bc].end()) kmax = UK[br-1][bc][la];
                    if(br>=bc+1 && UK[br][bc+1].find(la)!=UK[br][bc+1].end() && UK[br][bc+1][la]>kmax) kmax = UK[br][bc+1][la];
                    for(int r=0; r<nr; r++) {    
                        for(int c=0; c<nc; c++) {
                            auto item = expand_ex(m(r0+r,c0+c),sla).coeff(sla);
                            if(item.is_zero()) continue;
                            auto kd = expand_ex(item,lnx).degree(lnx)+1;
                            if(kd>kmax) kmax = kd;
                        }
                    }
                    if(kmax>0) {
                        UK[br][bc][la] = kmax;
                        U0[br][bc][la].resize(kmax);
                        for(int k=0; k<kmax; k++) {
                            matrix mat(nr,nc);
                            for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                                mat(r,c) = expand_ex(m(r0+r,c0+c),lst{sla,lnx}).coeff(sla).coeff(lnx,k);
                            }
                            U0[br][bc][la][k] = mat;
                        }
                    }
                }
            }
            for(auto kv : UK[br][0]) IK[br][kv.first] = kv.second;
        }
    }
    
    block_umat_t DE::series(int xn) { // U[a][b][la][k][n]
        auto nbs = bs.size();
        block_umat_t U(nbs);
        vector<vector<matrix>> A(nbs);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            U[br].resize(br+1);
            A[br].resize(br+1); // M=A/x
            for(int bc=br; bc>=0; bc--) { // cycle columns
                int c0 = bs[bc].first;
                int nc = bs[bc].second;

                // rational series for each block A[br][bc] 
                matrix mat = ex_to<matrix>(sub_matrix(Mat, r0, nr, c0, nc));
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    auto item = mat(r,c);
                    auto nd = item.numer_denom();
                    auto num = nd.op(0);
                    auto den = nd.op(1);
                    auto ld = den.ldegree(x);
                    int exn = xn;
                    if(ld<1) { // ld=0, otherwise not fuchsified yet
                        exn--; num *= x;
                    }
                    num /= den.coeff(x,ld);
                    den = den/den.coeff(x,ld)/pow(x,ld);
                    den = normal(den);
                    mat(r,c) = 0;
                    for(int n=0; n<=exn; n++) {
                        if(is_zero(den-1)) {
                            auto cvs = collect_lst(num,x);
                            for(auto cv : cvs) {
                                if(cv.op(1).degree(x)>xn) continue;
                                mat(r,c) += cv.op(0) * cv.op(1);
                            }
                            break;
                        } else {
                            auto cvs = collect_lst(num*pow(1-den,n),x);
                            for(auto cv : cvs) {
                                if(cv.op(1).degree(x)>xn) continue;
                                mat(r,c) += cv.op(0) * cv.op(1);
                            }
                        }
                    }
                }
                A[br][bc] = mat;
                
                // now we use a=br, b=bc
                // (la+n-A0aa).Uab(la,k,n) = -Uab(la,k+1,n)
                //   + sum_{1<=m<=n} Amaa(la,k,m).Uab(la,k,n-m) + sum_{b<=c<a,0<=m<=n} Amac.Ucb(la,k,n-m)
                int a = br, b = bc;
                matrix A0aa(nr,nr);
                for(int r=0; r<nr; r++) for(int c=0; c<nr; c++) A0aa(r,c) = A[a][a](r,c).coeff(x,0);
                
                int nla = UK[a][b].size();
                int cla = 0;
                for(auto kv : UK[a][b]) { // 1-cycle over lambda
                    cla++;
                    auto la = kv.first;
                    auto kmax = kv.second;
                    U[a][b][la].resize(kmax);
                    for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                        U[a][b][la][k].resize(xn+1);
                        U[a][b][la][k][0] = U0[a][b][la][k];
                        
                        bool abk_chk = (U[a][b].find(la)!=U[a][b].end() && U[a][b][la].size()>k);
                        for(int n=1; n<=xn; n++) { // 3-cycle over n
                        
                            if(!In_GiNaC_Parallel && Verbose>5) {
                                cout << "\r                                                              \r" << flush;
                                cout << "  \\--series: r" << nbs << "|" << br+1;
                                cout << " c" << br << "|" << (br-bc);
                                cout << " [" << nr << "\u2A09" << nc << "]";
                                cout << " \u03BB" << nla << "|" << cla;
                                cout << " k" << kmax << "|" << (kmax-k);
                                cout << " n" << xn << "|" << n << flush;
                            }

                            matrix smat(nr,nc);
                            if(k+1<kmax) smat = smat.sub(U[a][b][la][k+1][n]);
                            
                            if(abk_chk)
                            for(int m=1; m<=n; m++) {
                                matrix Amaa(nr,nr);
                                for(int r=0; r<nr; r++) for(int c=0; c<nr; c++) Amaa(r,c) = A[a][a](r,c).coeff(x,m);
                                smat = smat.add(Amaa.mul(U[a][b][la][k][n-m]));
                            }

                            for(int c=b; c<a; c++) for(int m=0; m<=n; m++) {
                                if(U[c][b].find(la)!=U[c][b].end() && U[c][b][la].size()>k) {
                                    int nc2 = bs[c].second;
                                    matrix Amac(nr,nc2);
                                    for(int i=0; i<nr; i++) for(int j=0; j<nc2; j++) Amac(i,j) = A[a][c](i,j).coeff(x,m);
                                    smat = smat.add(Amac.mul(U[c][b][la][k][n-m]));
                                }
                            }
                            
                            matrix invA = ex_to<matrix>(unit_matrix(nr));
                            invA = invA.mul_scalar(la+n).sub(A0aa);
                            invA = invA.inverse();
                            U[a][b][la][k][n] = invA.mul(smat);
                        }
                    }
                }
            }
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        
        return U;
    }
    
    block_imat_t DE::series(int xn, const vector<matrix> & cb) {
        auto nbs = bs.size();
        block_imat_n0_t In0(nbs);
        for(int br=0; br<nbs; br++) {
            int a = br;
            for(int bc=0; bc<=br; bc++) {
                int b = bc;
                for(auto & kv : U0[a][b]) {
                    auto la = kv.first;
                    int kmax = kv.second.size();
                    if(In0[a].find(la)==In0[a].end()) In0[a][la].resize(IK[a][la]);
                    for(int k=0; k<kmax; k++) {
                        if(bc==0) In0[a][la][k] = kv.second[k].mul(cb[bc]);
                        else In0[a][la][k] = In0[a][la][k].add(kv.second[k].mul(cb[bc]));
                    }
                }
            }
        }
        return series(xn,In0);
    }
    
    block_imat_t DE::series(int xn, block_imat_n0_t & In0) { // I[a][la][k][n] 
        auto nbs = bs.size();
        block_imat_t I(nbs);
        vector<vector<matrix>> A(nbs);
        for(int br=0; br<nbs; br++) { // cycle rows
            int r0 = bs[br].first;
            int nr = bs[br].second;
            A[br].resize(br+1); // M=A/x
            for(int bc=br; bc>=0; bc--) {
                int c0 = bs[bc].first;
                int nc = bs[bc].second;

                // rational series for each block A[br][bc] 
                matrix mat = ex_to<matrix>(sub_matrix(Mat, r0, nr, c0, nc));
                for(int r=0; r<nr; r++) for(int c=0; c<nc; c++) {
                    auto item = mat(r,c);
                    auto nd = item.numer_denom();
                    auto num = nd.op(0);
                    auto den = nd.op(1);
                    auto ld = den.ldegree(x);
                    int exn = xn;
                    if(ld<1) { // ld=0, otherwise not fuchsified yet
                        exn--; num *= x;
                    }
                    num /= den.coeff(x,ld);
                    den = den/den.coeff(x,ld)/pow(x,ld);
                    den = normal(den);
                    mat(r,c) = 0;
                    for(int n=0; n<=exn; n++) {
                        if(is_zero(den-1)) {
                            auto cvs = collect_lst(num,x);
                            for(auto cv : cvs) {
                                if(cv.op(1).degree(x)>xn) continue;
                                mat(r,c) += cv.op(0) * cv.op(1);
                            }
                            break;
                        } else {
                            auto cvs = collect_lst(num*pow(1-den,n),x);
                            for(auto cv : cvs) {
                                if(cv.op(1).degree(x)>xn) continue;
                                mat(r,c) += cv.op(0) * cv.op(1);
                            }
                        }
                    }
                }
                A[br][bc] = mat;
            }
                
            // now we use a=br
            // (la+n-A0aa).Ia(la,k,n) = -Ia(la,k+1,n)
            //   + sum_{1<=m<=n} Amaa(la,k,m).Ia(la,k,n-m) + sum_{b<a,0<=m<=n} Amab.Ib(la,k,n-m)
            int a = br;
            matrix A0aa(nr,nr);
            for(int r=0; r<nr; r++) for(int c=0; c<nr; c++) A0aa(r,c) = A[a][a](r,c).coeff(x,0);
                
            int nla = IK[a].size();
            int cla = 0;
            for(auto kv : IK[a]) { // 1-cycle over lambda
                cla++;
                auto la = kv.first;
                auto kmax = kv.second;
                I[a][la].resize(kmax);
                for(int k=kmax-1; k>=0; k--) { // 2-cycle over k
                    I[a][la][k].resize(xn+1);
                    I[a][la][k][0] = In0[a][la][k];
                    int nc = I[a][la][k][0].cols();
                    
                    bool abk_chk = (I[a].find(la)!=I[a].end() && I[a][la].size()>k);
                    for(int n=1; n<=xn; n++) { // 3-cycle over n
                    
                        if(!In_GiNaC_Parallel && Verbose>5) {
                            cout << "\r                                                              \r" << flush;
                            cout << "  \\--series: b" << nbs << "|" << br+1;
                            cout << " [" << nr << "\u2A09" << nc << "]";
                            cout << " \u03BB" << nla << "|" << cla;
                            cout << " k" << kmax << "|" << (kmax-k);
                            cout << " n" << xn << "|" << n << flush;
                        }

                        matrix smat(nr,nc);
                        if(k+1<kmax) smat = smat.sub(I[a][la][k+1][n]);
                        
                        if(abk_chk)
                        for(int m=1; m<=n; m++) {
                            matrix Amaa(nr,nr);
                            for(int r=0; r<nr; r++) for(int c=0; c<nr; c++) Amaa(r,c) = A[a][a](r,c).coeff(x,m);
                            smat = smat.add(Amaa.mul(I[a][la][k][n-m]));
                        }

                        for(int b=0; b<a; b++) for(int m=0; m<=n; m++) {
                            if(I[b].find(la)!=I[b].end() && I[b][la].size()>k) {
                                int nc2 = bs[b].second;
                                matrix Amac(nr,nc2);
                                for(int i=0; i<nr; i++) for(int j=0; j<nc2; j++) Amac(i,j) = A[a][b](i,j).coeff(x,m);
                                smat = smat.add(Amac.mul(I[b][la][k][n-m]));
                            }
                        }
                        
                        matrix invA = ex_to<matrix>(unit_matrix(nr));
                        invA = invA.mul_scalar(la+n).sub(A0aa);
                        invA = invA.inverse();
                        I[a][la][k][n] = invA.mul(smat);
                    }
                }
            }
            
        }
        if(!In_GiNaC_Parallel && Verbose>5) cout << endl;
        
        return I;
    }
    
    //=*********************************************************************=
    
    
    
    //=*********************************************************************=
}

