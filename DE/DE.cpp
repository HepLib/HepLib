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
    void DE::init(matrix m) {
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
    
    matrix DE::T() {
        matrix t = ex_to<matrix>(unit_matrix(N));
        for(auto ti : Ts) t = t.mul(ti);
        return t;
    }
    
    void DE::fuchsify() { // first on diagonal blocks, then off diagnoal ones
        matrix t(N,N);
        auto nbs = bs.size();
        for(int bi=0; bi<nbs; bi++) { // diagonal blocks
            if(!In_GiNaC_Parallel && Verbose>5) {
                cout << "\r                                                    \r" << flush;
                cout << "  \\--reducing diagonal blocks: " << bi+1 << "/" << nbs << flush;
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
                    cout << "  \\--reducing off-diagonal blocks: " << bc << "-" << br+1 << "/" << nbs << flush;
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
    }
    
    //=*********************************************************************=
    
}

