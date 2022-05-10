/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */
    
#include "DE.h"

namespace HepLib {

    using namespace EoD;
    
    NDEH::NDEH(const symbol & _x) : x(_x), a("a") { }
    NDEH::NDEH(const matrix & _mat, const symbol & _x) : Mat(_mat), x(_x), a("a") { }
    NDEH::NDEH(const symbol & _x, const matrix & _mat) : Mat(_mat), x(_x), a("a") { }
    NDEH::NDEH(const NDEH & b) : Mat(b.Mat), x(b.x), Ts(b.Ts), a("a") { }
        
    // Dx J = M.J ---> J = T.J' & Dx J' = M'.J' with M' = Ti.M.T - Ti.Dx T
    void NDEH::Apply(const matrix & t, bool st) {
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        if(st) Ts.push_back(t);
    }
    
    void NDEH::Apply(const lst & diag, bool st) {
        auto matN = diag.nops();
        matrix t(matN, matN);
        for(int i=0; i<matN; i++) t(i,i) = diag.op(i);
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        if(st) Ts.push_back(t);
    }
    
    void NDEH::x2y(const ex & y) {
        static symbol xx("xx");
        ex det = diff_ex(y,x);
        int matN = Mat.rows();
        for(int r=0; r<matN; r++) {
            for(int c=0; c<matN; c++) {
                Mat(r,c) = det * Mat(r,c).subs(x==xx,nopat).subs(xx==y,nopat);
            }
        }
        auto n = Ts.size();
        for(int i=0; i<n; i++) {
            for(int r=0; r<matN; r++) {
                for(int c=0; c<matN; c++) {
                    Ts[i](r,c) = Ts[i](r,c).subs(x==xx,nopat).subs(xx==y,nopat);
                }
            }
        }
    }
    
    void NDEH::xpow() {
        int n = Mat.nops();
        for(int i=0; i<n; i++) Mat.let_op(i) = HepLib::xpow(Mat.op(i),x);
    }
    
    matrix NDEH::MatT() {
        int matN = Mat.rows();
        matrix t = ex_to<matrix>(unit_matrix(matN));
        for(auto ti : Ts) t = t.mul(ti);
        return t;
    }
        
    void NDEH::subs(const ex & sub, unsigned opt) {
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
    }
    
    void NDEH::subs(const exmap & sub, unsigned opt) {
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
    }
    
    void NDEH::subs(const lst & l, const lst & r, unsigned opt) {
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(l, r, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(l, r, opt);
        }
    }
    
    // fuchsify at x=0
    void NDEH::Fuchsify() {
        int matN = Mat.rows();
        mx.clear();
        mx.init(Mat);
        auto a01 = mx.a01();
        a0 = a01.first;
        matrix a1 = a01.second;
        symbol lambda("L");
        matrix t= ex_to<matrix>(unit_matrix(matN));
        
        int pr = mx.prank();
        if(pr<1) {
            if(!In_GiNaC_Parallel && Verbose>1) 
                cout << "  \\--Fuchsified: Poincare rank is alreay " << pr << endl;
            if(pr<0) a0 = matrix(matN,matN); // a0 = 0
            Ts.push_back(t);
            fuchsified = true;
            return;
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsifying @ " << now() << endl;
        
        while(true) {
            if(!In_GiNaC_Parallel && Verbose>1) 
                cout << "     \\--Poincare rank: " << pr << ",  Ao rank: " << a0.rank() << endl;
            
            // n = (a0 a1-l)
            //     (0  a0  )
            matrix n(2*a0.rows(), 2*a0.cols());
            for (unsigned i = 0; i < a0.rows(); i++) {
                for (unsigned j = 0; j < a0.cols(); j++) {
                    n(i, j) = a0(i, j);
                    n(i + a0.rows(), j + a0.cols()) = a0(i, j);
                }
            }
            for (unsigned i = 0; i < a1.rows(); i++) {
                for (unsigned j = 0; j < a1.cols(); j++) {
                    n(i, j + a0.cols()) = (i == j) ? a1(i, j) - lambda : a1(i, j);
                }
            }
            
            vspace ns = nullspace(n);
            vspace ws(a0.cols());
            
            // Find the span of coefficients of the second half of ns as a poly in lambda.
            for (unsigned i = 0; i < ns.dim(); i++) {
                int maxdeg = 0;
                for (unsigned j = 0; j < a0.cols(); j++) {
                    auto &&e = ns.basis_rows()(i, j + a0.cols());
                    maxdeg = max(maxdeg, e.degree(lambda));
                }
                matrix s(maxdeg + 1, a0.cols());
                for (unsigned j = 0; j < a0.cols(); j++) {
                    // It would be great to only expand by lambda here.
                    ex e = expand(ns.basis_rows()(i, j + a0.cols()));
                    for (int deg = 0; deg <= maxdeg; deg++) {
                        s(deg, j) = exfactor(normal(e.coeff(lambda, deg)));
                    }
                }
                ws.add_rows(s);
            }
            
            if(ws.basis_rows().has(lambda)) throw Error("Error: ws.basis_rows().has(lambda)");
            
            ws.normalize();
            matrix wscols = ws.basis_cols();
            if (ws.dim() == 0) {
                throw Error("Fuchsify: matrix is Moser-irreducible");
            }
            
            auto dualbs = dual_basis(wscols);
            
            if(dualbs.size()<=0) throw Error("Error: dualbs.size()<=0");
            
            auto dualb = dualbs[0];
            matrix p = normal(wscols.mul(dualb));
            
            auto tk = imatrix(matN).sub(p).sub(p.mul_scalar(1/x));
            //auto tki = fermat_inv(tk);
            //mx.transform(tk,tki);
            mx.balance(p);
            t = t.mul(tk);
            
            pr = mx.prank();
            if(pr<1) break;
            a01 = mx.a01();
            a0 = a01.first;
            a1 = a01.second;
        }

        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsified @ " << now() << endl;
        
        if(pr<0) a0 = matrix(matN,matN); // a0 = 0
        else a0 = mx.a0();
        Ts.push_back(t);
        fuchsified = true;
    }
    
    // shearing transformation at x=0
    void NDEH::Shear() {
        if(!fuchsified) Fuchsify();
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: start @ " << now() << endl;

        int matN = Mat.rows();
        matrix t = ex_to<matrix>(unit_matrix(matN));
        matrix ti = t;
        
        while(true) {
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--jordan decomposition ..." << endl;
            auto qj = jordan(a0);
            auto tm = qj.first;
            auto tmi = qj.first.inverse();
            
            bool shearing = false;
            matrix smat = ex_to<matrix>(unit_matrix(matN));
            int cpos = 0;
            int sm = 0;
            for(auto kv : qj.second) {
                ex fx = 1;
                if(sm<=0 && kv.first<-1/ex(2)) { sm=-1; fx = 1/x; }
                else if(sm>=0 && kv.first>=1/ex(2)) { sm=1; fx = x; }
                for(int j=0; j<kv.second; j++) smat(cpos+j, cpos+j) = fx;
                if(sm!=0) shearing = true;
                cpos += kv.second;
            }
            if(shearing) {
                if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--shearing transformation ..." << endl;
                tm = qj.first.mul(smat);
                tmi = smat.inverse().mul(qj.first.inverse());
            }
            mx.transform(tm,tmi);
            t = t.mul(tm);
            a0 = mx.a0();
            if(a0.has(x)) throw Error("a0 still has x.");
            if(!shearing) break;
        }
        Ts.push_back(t);
        
        // ---------------------------------------------------------------
        // parse a0 information
        // ---------------------------------------------------------------
        lst _lass;
        vector<pair<ex,int>> js; // {la, size}
        int lastN = 0;
        for(int r=0; r<matN; r++) {
            lastN++;
            if(r==matN-1 || a0(r,r+1)==0) { // r==matN-1 : the last row
                auto la = a0(r,r);
                _lass.append(la);
                js.push_back(make_pair(la, lastN));
                if(lastN-1 > ST.K[la]) ST.K[la] = lastN-1;
                lastN = 0;
            }
        }
        _lass.sort().unique();
        sort_lst(_lass);
        for(auto la : _lass) {
            las.push_back(la);
            ST.las.push_back(la);
        }
        
        ST.kmax = -1;
        ST.C0.clear();
        for(int idx=0; idx<ST.las.size(); idx++) {
            auto la = ST.las[idx];
            int kla = ST.K[la];
            if(ST.kmax<kla) ST.kmax = kla;
            vector<matrix> vm;
            matrix zm(matN,matN);
            for(int li=0; li<=kla; li++) vm.push_back(zm);
            ST.C0[la] = vm;
        }
        
        //init ST.C0
        int cur_pos = 0;
        for(int ji=0; ji<js.size(); ji++) {
            auto la = js[ji].first;
            int size = js[ji].second;
            for(int ir=0;ir<size;ir++) 
                for(int ic=0;ic<size;ic++)
                    if(ic-ir>=0) ST.C0[la][ic-ir](cur_pos+ir, cur_pos+ic) = 1; 
            cur_pos += size;
        }
                        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: fininshed @ " << now() << endl;
        sheared = true;
    }
    
    void NDEH::Shear0() {
        if(!fuchsified) Fuchsify();
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: start @ " << now() << endl;

        int matN = Mat.rows();
        matrix t = ex_to<matrix>(unit_matrix(matN));
        matrix ti = t;
        
        while(!is_sheared_form(a0)) {
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--eigen values of Ao ..." << endl;
            auto ev_map = eigenvalues(a0);
            vector<exvector> ev_groups;
            for(auto &kv : ev_map) {
                auto ev = kv.first;
                for(auto &gi : ev_groups) {
                    for(auto &ev_i : gi) {
                        auto diff_ev = normal(ev_i-ev);
                        if(diff_ev.info(info_flags::integer)) {
                            gi.push_back(ev);
                            goto gi_done;
                        }
                    }
                }
                if(true) { // no found
                    exvector gi_new;
                    gi_new.push_back(ev);
                    ev_groups.push_back(gi_new);
                }
                gi_done: ;
            }
            
            bool is_normalized = true;
            for(auto &gi : ev_groups) {
                if(gi.size()>1) { is_normalized = false; break; }
            }
            
            for(auto &gi : ev_groups) {
                sort(gi.begin(),gi.end(),[&](const auto &a, const auto &b){
                    return normal(a-b).info(info_flags::positive);
                });
            }
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--jordan decomposition ..." << endl;
            auto qj = jordan(a0,ev_map);
            auto tm = qj.first;
            auto tmi = qj.first.inverse();
                        
            if(!is_normalized) {
                if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--shearing transformation ..." << endl;
                matrix smat = imatrix(matN);
                int cpos = 0;
                for(auto kv : qj.second) {
                    bool is_last = false;
                    for(auto &gi : ev_groups) {
                        if(normal(gi[gi.size()-1] - kv.first).is_zero()) {
                            is_last = true;
                            break;
                        }
                    }
                    if(!is_last) for(int j=0; j<kv.second; j++) smat(cpos+j, cpos+j) = x;
                    cpos += kv.second;
                }
                tm = qj.first.mul(smat);
                tmi = smat.inverse().mul(qj.first.inverse());
            }
            
            mx.transform(tm,tmi);
            t = t.mul(tm);
            a0 = mx.a0();
            if(a0.has(x)) throw Error("a0 still has x.");
        }
        Ts.push_back(t);
        
        // ---------------------------------------------------------------
        // parse a0 information
        // ---------------------------------------------------------------
        lst _lass;
        vector<pair<ex,int>> js; // {la, size}
        int lastN = 0;
        for(int r=0; r<matN; r++) {
            lastN++;
            if(r==matN-1 || a0(r,r+1)==0) { // r==matN-1 : the last row
                auto la = a0(r,r);
                _lass.append(la);
                js.push_back(make_pair(la, lastN));
                if(lastN-1 > ST.K[la]) ST.K[la] = lastN-1;
                lastN = 0;
            }
        }
        _lass.sort().unique();
        sort_lst(_lass);
        for(auto la : _lass) {
            las.push_back(la);
            ST.las.push_back(la);
        }
        
        ST.kmax = -1;
        ST.C0.clear();
        for(int idx=0; idx<ST.las.size(); idx++) {
            auto la = ST.las[idx];
            int kla = ST.K[la];
            if(ST.kmax<kla) ST.kmax = kla;
            vector<matrix> vm;
            matrix zm(matN,matN);
            for(int li=0; li<=kla; li++) vm.push_back(zm);
            ST.C0[la] = vm;
        }
        
        //init ST.C0
        int cur_pos = 0;
        for(int ji=0; ji<js.size(); ji++) {
            auto la = js[ji].first;
            int size = js[ji].second;
            for(int ir=0;ir<size;ir++) 
                for(int ic=0;ic<size;ic++)
                    if(ic-ir>=0) ST.C0[la][ic-ir](cur_pos+ir, cur_pos+ic) = 1; 
            cur_pos += size;
        }
                        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: fininshed @ " << now() << endl;
        sheared = true;
    }
        
    matrix NDEH::Series(const ex & x0, const int xN, const lst & _las) { 
    
        if(!is_a<numeric>(x0)) {
            auto CMats = Series(xN,_las); //  map<ex,vector<vector<matrix>>,ex_is_less>
            auto matN = mx.n;
            matrix MatF(matN,matN);
            for(const auto & kv : CMats) {
                ex la = kv.first;
                auto xla = pow(x0,la);
                auto nk = kv.second.size();
                auto lx = log(x0);
                ex lxn = 1;
                for(int k=0; k<nk; k++) {
                    auto cmats = kv.second[k];
                    auto nn = cmats.size();
                    ex xn = 1;
                    for(int n=0; n<nn; n++) {
                        auto xterm = xla*lxn*xn/factorial(k);
                        MatF = MatF.add(cmats[n].mul_scalar(xterm));
                        xn *= x0;
                    }
                    lxn *= lx;
                }
            }
            return MatF;
        }
        
        if(!fuchsified) Fuchsify();
        if(!sheared) Shear();
        if(!is_sheared_form(a0)) throw Error("a0 is NOT sheared form.");

        if(mx.s<0) {    
            if(!In_GiNaC_Parallel && Verbose>10) cout << "     \\--LCM " << flush; 
            mx.lcm();
            if(!In_GiNaC_Parallel && Verbose>10) cout << "\r                 \r     \\--LCM s = " << mx.s << endl;
        }
        int matN = mx.n;        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << NN(x0,2) << endl;
                
        int s = mx.s;
        if(s>xN) s = xN;
        if(s<1) s = 1;
        
        auto _fp_ = fp+10000;
        auto _dp = dp+5000;
        acb_poly_t Qx, M[matN][matN];
        acb_poly_init(Qx);
        acb_poly_set_fmpz_poly(Qx,mx.Qx,_fp_);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
            acb_poly_init(M[r][c]);
            auto num = fmpz_poly_q_numref(mx.M[r][c]);
            acb_poly_set_fmpz_poly(M[r][c],num,_fp_); 
        }
        if(mx_clear) mx.clear();
        
        acb_t z0; acb_init(z0);
        acb_t z0n; acb_init(z0n);
        to_acb(z0,x0,_dp);
        
        acb_t Qs[s+1], q0;
        acb_mat_t QxM[s+1];
        acb_init(q0);
        acb_poly_get_coeff_acb(q0, Qx, 0);
        acb_one(z0n);
        for(int i=0; i<=s; i++) {
            acb_init(Qs[i]);
            acb_mat_init(QxM[i],matN,matN);
            acb_poly_get_coeff_acb(Qs[i], Qx, i);
            acb_mul(Qs[i],Qs[i],z0n,_fp_);
            acb_div(Qs[i],Qs[i],q0,_fp_);
            for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                auto item = acb_mat_entry(QxM[i],r,c);
                acb_poly_get_coeff_acb(item, M[r][c], i);
                acb_mul(item,item,z0n,_fp_);
                acb_div(item,item,q0,_fp_);
            }
            acb_mul(z0n,z0n,z0,_fp_);
        }
        
        // clear Qx & M
        acb_clear(z0);
        acb_clear(z0n);
        acb_clear(q0);
        acb_poly_clear(Qx);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) acb_poly_clear(M[r][c]);
        flint_cleanup();
        
        flint_set_num_threads(omp_get_num_procs()-1);
        
        // now cycle w.r.t lambda set
        acb_mat_t MatF, _MatF;
        acb_mat_init(MatF,matN,matN);
        acb_mat_init(_MatF,matN,matN);
        acb_mat_zero(MatF);
        
        int las_size = ST.las.size();
        if(_las.nops()>0) las_size = _las.nops();
        int _idx = 0;
        for(int idx=0; idx<ST.las.size(); idx++) { // cycle w.r.t lambda set
            auto la = ST.las[idx];
            if(_las.nops()>0) { // pick up the seleted la
                bool ok = false;
                for(auto _la : _las) {
                    ex diff = la-_la;
                    if(diff.info(info_flags::integer)) {
                        ok = true;
                        break;
                    }
                }
                if(!ok) continue;
            }
            _idx++;
            int kla = ST.K[la];
            
            vector<acb_struct*> z_cls;
            vector<acb_mat_struct*> m_cls;
            
            // FLINT variables            
            acb_t z, la_;
            acb_init(z);
            z_cls.push_back(z);
            acb_init(la_); 
            z_cls.push_back(la_);
            to_acb(la_,la,_dp);
            
            // setup B0 & its inverse
            acb_mat_t B0;
            acb_mat_init(B0,matN,matN); 
            m_cls.push_back(B0);
            acb_mat_t iB0[kla+1]; 
            for(int i=0; i<=kla; i++) {
                acb_mat_init(iB0[i],matN,matN);
                m_cls.push_back(iB0[i]);
            }
            
            // setup CMats & init CMats[0][k]
            acb_t xla_;
            acb_init(xla_);
            to_acb(xla_,pow(x0,la),_dp);
            acb_mat_t CMats[s][kla+1], _CMats[kla+1];
            acb_mat_set(_MatF,MatF); // get last MatF
            for(int k=0; k<=kla; k++) {
                for(int i=0; i<s; i++) {
                    acb_mat_init(CMats[i][k],matN,matN);
                    m_cls.push_back(CMats[i][k]);
                }
                acb_mat_init(_CMats[k],matN,matN);
                m_cls.push_back(_CMats[k]);
                
                to_acb(_CMats[k],ST.C0[la][k],_dp);
                acb_mat_scalar_mul_acb(_CMats[k],_CMats[k],xla_,_fp_);
                if(k==0) acb_mat_add(_MatF,_MatF,_CMats[k],_fp_);
            }
            acb_clear(xla_);

            auto _fp = fp + 500;
            restart: ; // precision restart here
            for(int k=0; k<=kla; k++) acb_mat_set(CMats[0][k],_CMats[k]);
            acb_mat_set(MatF,_MatF);
            
            for(int cn=1; cn<=xN; cn++) {
                if(!In_GiNaC_Parallel && Verbose>10) {
                    cout << "\r                                  \r" << flush;
                    cout << "     \\--C Matrix [" << _idx << "/" << las_size << "][" << cn << "/" << xN << "]" << flush;
                }
                
                // initialize iB0, depend only on la+n
                acb_add_si(z,la_,cn,_fp_); // z=la+n, alpha in B
                acb_mat_one(B0);
                acb_mat_scalar_mul_acb(B0,B0,z,_fp_); // q0 is set to 1
                acb_mat_sub(B0,QxM[0],B0,_fp_);
                if(!acb_mat_inv(iB0[0],B0,_fp_)) throw Error("acb_mat_inv gets wrong.");
                for(int i=1; i<=kla; i++) acb_mat_mul_threaded(iB0[i], iB0[i-1], iB0[0], _fp_);
                
                for(int k=0; k<=kla; k++) { // k-th row
                    int tot = (cn<s ? cn : s);
                    acb_mat_t mat_vec[tot];
                    for(int i=0; i<tot; i++) acb_mat_init(mat_vec[i],matN,matN);
                    #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
                    for(int cm=1; cm<=tot; cm++) { // sum m from 1 to s in DESS
                        acb_mat_zero(mat_vec[cm-1]);
                        acb_t z;
                        acb_init(z);
                        acb_mat_t Bm, T;
                        acb_mat_init(Bm,matN,matN);
                        acb_mat_init(T,matN,matN);
                        acb_add_si(z,la_,cn-cm,_fp_); // z = la+n-m
                        acb_mul(z,z,Qs[cm],_fp_);
                        acb_mat_set(Bm,QxM[cm]);
                        #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
                        for(int r=0; r<matN; r++) {
                            acb_sub(acb_mat_entry(Bm,r,r),acb_mat_entry(Bm,r,r),z,_fp_);
                            flint_cleanup();
                        }
                        for(int j=0; j<=kla; j++) { // T is T[k,j] in DESS
                            // T[k,j] = -iB0[j-k].Bm + iB0[j-k-1].Qm, q0 set 1
                            if(j<k) continue; // zero
                            acb_mat_mul_threaded(T,iB0[j-k],Bm,_fp);
                            acb_mat_neg(T,T);
                            if(j>k) acb_mat_scalar_addmul_acb(T,iB0[j-k-1],Qs[cm],_fp);
                            // get Tkj.Cj
                            acb_mat_mul_threaded(T,T,CMats[(cn-cm)%s][j],_fp);
                            acb_mat_add(mat_vec[cm-1],mat_vec[cm-1],T,_fp_);
                        }
                        acb_clear(z);
                        acb_mat_clear(Bm);
                        acb_mat_clear(T);
                        flint_cleanup();
                    } 
                    
                    acb_mat_t CMat;
                    acb_mat_init(CMat,matN,matN);
                    acb_mat_zero(CMat);
                    for(int i=0; i<tot; i++) {
                        acb_mat_add(CMat,CMat,mat_vec[i],_fp_);
                        acb_mat_clear(mat_vec[i]);
                    }
                    
                    acb_mat_set(CMats[cn%s][k],CMat); // before here, CMats should not be modified
                    if(k==0) acb_mat_add(MatF,MatF,CMat,_fp_);
                    acb_mat_clear(CMat);
                    
                    if(k==0) { // precision check
                        mag_t mag;
                        mag_init(mag);
                        bool ok = true;
                        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                            auto item = acb_mat_entry(MatF,r,c);
                            auto ri = acb_realref(item);
                            if(arb_rel_error_bits(ri)<-rel_fp) continue;
                            arb_get_mag(mag,ri);
                            if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                if(Verbose>500) { cout << endl; arb_printd(ri,5); cout << endl; }
                                ok = false;
                                goto done; 
                            }
                            if(arb_contains_si(ri,0)) arb_zero(ri);
                            ri = acb_imagref(item);
                            if(arb_rel_error_bits(ri)<-rel_fp) continue;
                            if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                if(Verbose>500) { cout << endl; arb_printd(ri,5); cout << endl; }
                                ok = false;
                                goto done; 
                            }
                        }
                        done: ;
                        mag_clear(mag); 
                        if(!ok) {
                            _fp += _fp_/10;
                            flint_cleanup();
                            if(_fp>_fp_) throw Error("still too low precision with high fp!");
                            goto restart;
                        }
                    }
                }
            }
            
            for(auto z : z_cls) acb_clear(z);
            for(auto m : m_cls) acb_mat_clear(m);
            if(!In_GiNaC_Parallel && Verbose>10 && xN>0) cout << endl;
            fp = _fp;
        }
        
        for(int i=0; i<=s; i++) {
            acb_clear(Qs[i]);
            acb_mat_clear(QxM[i]);
        }
        auto mat = acb_to_mat(MatF,_dp);
        acb_mat_clear(MatF);
        acb_mat_clear(_MatF);
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << now() << endl;
        return mat;
    }
    
    CMatrix NDEH::Series(const int xN, const lst & _las) { 
        if(!fuchsified) Fuchsify();
        if(!sheared) Shear();
        if(!is_sheared_form(a0)) throw Error("a0 is NOT sheared form.");

        if(mx.s<0) {    
            if(!In_GiNaC_Parallel && Verbose>10) cout << "     \\--LCM ..." << flush; 
            mx.lcm();
            if(!In_GiNaC_Parallel && Verbose>10) cout << "\r                 \r     \\--LCM s = " << mx.s << endl; 
        }
        int matN = Mat.rows();
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ xN=" << xN << endl;
                
        int s = mx.s;
        if(s>xN) s = xN;
        if(s<1) s = 1;
        
        fmpq_poly_t Qx, M[matN][matN];
        fmpq_poly_init(Qx);
        fmpq_poly_set_fmpz_poly(Qx,mx.Qx);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
            fmpq_poly_init(M[r][c]);
            auto num = fmpz_poly_q_numref(mx.M[r][c]);
            fmpq_poly_set_fmpz_poly(M[r][c],num); 
        }
        if(mx_clear) mx.clear();
        
        fmpq_t Qs[s+1], q0;
        fmpq_mat_t QxM[s+1];
        fmpq_init(q0);
        fmpq_poly_get_coeff_fmpq(q0, Qx, 0);
//fmpq_print(q0);cout << endl;
        for(int i=0; i<=s; i++) {
            fmpq_init(Qs[i]);
            fmpq_mat_init(QxM[i],matN,matN);
            fmpq_poly_get_coeff_fmpq(Qs[i], Qx, i);
            fmpq_div(Qs[i],Qs[i],q0);
            for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                auto item = acb_mat_entry(QxM[i],r,c);
                fmpq_poly_get_coeff_fmpq(item, M[r][c], i);
                fmpq_div(item,item,q0);
            }
        }
        
        // clear Qx & M
        fmpq_poly_clear(Qx);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) fmpq_poly_clear(M[r][c]);
        
        CMatrix CMatF; // CMatF[la][k][n] : coefficient of x^la*(log(x)^k/k!)*x^n
        int las_size = ST.las.size();
        if(_las.nops()>0) las_size = _las.nops();
        int _idx = 0;
        for(int idx=0; idx<ST.las.size(); idx++) {
            auto la = ST.las[idx];
            if(_las.nops()>0) { // pick up the seleted la
                bool ok = false;
                for(auto _la : _las) {
                    ex diff = la-_la;
                    if(diff.info(info_flags::integer)) {
                        ok = true;
                        break;
                    }
                }
                if(!ok) continue;
            }
            _idx++;
            int kla = ST.K[la];
            
            CMatF[la].resize(kla+1);
            for(int k=0; k<=kla; k++) {
                auto t = ST.C0[la][k];
                CMatF[la][k].resize(xN+1);
                CMatF[la][k][0] = t;
            }
            
            vector<fmpq_mat_struct*> mat_to_clear; // fmpq_mat_t to be cleared

            fmpq_mat_t B0, Bm;
            fmpq_mat_init(B0,matN,matN); mat_to_clear.push_back(B0);
            fmpq_mat_init(Bm,matN,matN); mat_to_clear.push_back(Bm);
            fmpq_t q, la_;
            fmpq_init(q); fmpq_init(la_);
            ex_to_fmpq(la_,la);
            
            // setup the inverse of BJF
            fmpq_mat_t iB0[kla+1];
            for(int i=0; i<=kla; i++) {
                fmpq_mat_init(iB0[i],matN,matN);
                mat_to_clear.push_back(iB0[i]);
            }
            
            fmpq_mat_t CMat, T, m;
            fmpq_mat_init(CMat,matN,matN); // k-th row of C in DESS
            mat_to_clear.push_back(CMat);
            fmpq_mat_init(T,matN,matN);
            mat_to_clear.push_back(T);
            fmpq_mat_init(m,matN,matN);
            mat_to_clear.push_back(m);
            
            // setup CMats & init CMats[0][k]
            fmpq_mat_t CMats[s][kla+1];
            for(int k=0; k<=kla; k++) {
                for(int i=0; i<s; i++) {
                    fmpq_mat_init(CMats[i][k],matN,matN);
                    mat_to_clear.push_back(CMats[i][k]);
                }
                mat_to_fmpq(CMats[0][k],ST.C0[la][k]);
            }
            
            for(int cn=1; cn<=xN; cn++) {
                if(!In_GiNaC_Parallel && Verbose>10) {
                    cout << "\r                                  \r" << flush;
                    cout << "     \\--C Matrix [" << _idx << "/" << las_size << "][" << cn << "/" << xN << "]" << flush;
                }
                
                // initialize iBJF, depend only on la+n
                fmpq_add_si(q,la_,cn); // q=la+n, alpha in B
                fmpq_mat_one(B0);
                fmpq_mat_scalar_mul_fmpq(B0,B0,q);
                fmpq_mat_sub(B0,QxM[0],B0);
                if(!fmpq_mat_inv(iB0[0],B0)) throw Error("fmpq_mat_inv gets wrong.");;
                for(int i=1; i<=kla; i++) fmpq_mat_mul(iB0[i], iB0[i-1], iB0[0]);

                for(int k=0; k<=kla; k++) { // k-th row
                    int tot = (cn<s ? cn : s);
                    fmpq_mat_t mat_vec[tot];
                    for(int i=0; i<tot; i++) fmpq_mat_init(mat_vec[i],matN,matN);
                    #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
                    for(int cm=1; cm<=tot; cm++) { // sum m from 1 to s in DESS
                        fmpq_mat_zero(mat_vec[cm-1]);
                        fmpq_t q;
                        fmpq_init(q);
                        fmpq_mat_t Bm,T;
                        fmpq_mat_init(Bm,matN,matN);
                        fmpq_mat_init(T,matN,matN);
                        fmpq_add_si(q,la_,cn-cm); // q = la+n-m
                        fmpq_mul(q,q,Qs[cm]);
                        fmpq_mat_one(Bm);
                        fmpq_mat_scalar_mul_fmpq(Bm,Bm,q);
                        fmpq_mat_sub(Bm,QxM[cm],Bm);                        
                        for(int j=0; j<=kla; j++) { // T is T[k,j] in DESS
                            // T[k,j] = -iB0[j-k].Bm + iB0[j-k-1].Qm, q0 set 1
                            if(j<k) continue; // zero
                            fmpq_mat_mul(T,iB0[j-k],Bm);
                            fmpq_mat_neg(T,T);
                            if(j>k) {
                                fmpq_mat_scalar_mul_fmpq(m,iB0[j-k-1],Qs[cm]);
                                fmpq_mat_add(T,T,m);
                            }
                            // get Tkj.Cj
                            fmpq_mat_mul(T,T,CMats[(cn-cm)%s][j]);
                            fmpq_mat_add(mat_vec[cm-1],mat_vec[cm-1],T);
                        }
                        fmpq_clear(q);
                        fmpq_mat_clear(Bm);
                        fmpq_mat_clear(T);
                        flint_cleanup();
                    }
                    
                    fmpq_mat_zero(CMat);
                    for(int i=0; i<tot; i++) {
                        fmpq_mat_add(CMat,CMat,mat_vec[i]);
                        fmpq_mat_clear(mat_vec[i]);
                    }
                     
                    fmpq_mat_set(CMats[cn%s][k],CMat); // before here, CMats should not be modified
                    CMatF[la][k][cn] = fmpq_to_mat(CMat);
                }
            }

            for(auto m : mat_to_clear) fmpq_mat_clear(m);
            fmpq_clear(q); fmpq_clear(la_);
            if(!In_GiNaC_Parallel && Verbose>10 && xN>0) cout << endl;
        }
        
        for(int i=0; i<=s; i++) {
            fmpq_clear(Qs[i]);
            fmpq_mat_clear(QxM[i]);
        }
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << now() << endl;
        return CMatF;
    }
    
    // together with boundary constant
    matrix NDEH::Series(matrix C, const ex & x0, const int xN, const lst & _las) { 
    
        if(!is_a<numeric>(x0)) throw Error("x0 is not a number.");
        
        if(!fuchsified) Fuchsify();
        if(!sheared) Shear();
        if(!is_sheared_form(a0)) throw Error("a0 is NOT sheared form.");

        if(mx.s<0) {    
            if(!In_GiNaC_Parallel && Verbose>10) cout << "     \\--LCM " << flush; 
            mx.lcm();
            if(!In_GiNaC_Parallel && Verbose>10) cout << "\r                 \r     \\--LCM s = " << mx.s << endl;
        }
        int matN = mx.n;
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series W.C. @ " << NN(x0,2) << endl;
                
        int s = mx.s;
        if(s>xN) s = xN;
        if(s<1) s = 1;

        auto _fp_ = fp+10000;
        auto _dp = dp+5000;
        acb_poly_t Qx, M[matN][matN];
        acb_poly_init(Qx);
        acb_poly_set_fmpz_poly(Qx,mx.Qx,_fp_);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
            acb_poly_init(M[r][c]);
            auto num = fmpz_poly_q_numref(mx.M[r][c]);
            acb_poly_set_fmpz_poly(M[r][c],num,_fp_); 
        }
        if(mx_clear) mx.clear();
        
        acb_t z0; acb_init(z0);
        acb_t z0n; acb_init(z0n);
        to_acb(z0,x0,_dp);
        
        acb_t Qs[s+1], q0;
        acb_mat_t QxM[s+1];
        acb_init(q0);
        acb_poly_get_coeff_acb(q0, Qx, 0);
        acb_one(z0n);
        for(int i=0; i<=s; i++) {
            acb_init(Qs[i]);
            acb_mat_init(QxM[i],matN,matN);
            acb_poly_get_coeff_acb(Qs[i], Qx, i);
            acb_mul(Qs[i],Qs[i],z0n,_fp_);
            acb_div(Qs[i],Qs[i],q0,_fp_);
            for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                auto item = acb_mat_entry(QxM[i],r,c);
                acb_poly_get_coeff_acb(item, M[r][c], i);
                acb_mul(item,item,z0n,_fp_);
                acb_div(item,item,q0,_fp_);
            }
            acb_mul(z0n,z0n,z0,_fp_);
        }
        
        // clear Qx & M
        acb_clear(z0);
        acb_clear(z0n);
        acb_clear(q0);
        acb_poly_clear(Qx);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) acb_poly_clear(M[r][c]);
        
        // now cycle w.r.t lambda set
        int nc = C.cols();
        acb_mat_t MatF, _MatF;
        acb_mat_init(MatF,matN,nc);
        acb_mat_init(_MatF,matN,nc);
        acb_mat_zero(MatF);
        
        int las_size = ST.las.size();
        if(_las.nops()>0) las_size = _las.nops();
        int _idx = 0;
        for(int idx=0; idx<ST.las.size(); idx++) { // cycle w.r.t lambda set
            auto la = ST.las[idx];
            if(_las.nops()>0) { // pick up the seleted la
                bool ok = false;
                for(auto _la : _las) {
                    ex diff = la-_la;
                    if(diff.info(info_flags::integer)) {
                        ok = true;
                        break;
                    }
                }
                if(!ok) continue;
            }
            _idx++;
            int kla = ST.K[la];
                        
            vector<acb_struct*> z_cls;
            vector<acb_mat_struct*> m_cls;
            
            acb_t z, la_;
            acb_init(z);
            z_cls.push_back(z);
            acb_init(la_); 
            z_cls.push_back(la_);
            to_acb(la_,la,_dp);
            
            // setup B0 & its inverse
            acb_mat_t B0;
            acb_mat_init(B0,matN,matN); 
            m_cls.push_back(B0);
            acb_mat_t iB0[kla+1]; 
            for(int i=0; i<=kla; i++) {
                acb_mat_init(iB0[i],matN,matN);
                m_cls.push_back(iB0[i]);
            }
            
            // setup CMats & init CMats[0][k]
            acb_t xla_; 
            acb_init(xla_);
            to_acb(xla_,pow(x0,la),_dp);
            acb_mat_set(_MatF,MatF); // get MatF from last la
            acb_mat_t CMats[s][kla+1], _CMats[kla+1];
            for(int k=0; k<=kla; k++) {
                for(int i=0; i<s; i++) {
                    acb_mat_init(CMats[i][k],matN,nc);
                    m_cls.push_back(CMats[i][k]);
                }
                acb_mat_init(_CMats[k],matN,nc);
                m_cls.push_back(_CMats[k]);
                to_acb(_CMats[k],ST.C0[la][k].mul(C),_dp);
                acb_mat_scalar_mul_acb(_CMats[k],_CMats[k],xla_,_fp_);
                if(k==0) acb_mat_add(_MatF,_MatF,_CMats[k],_fp_);
            }
            acb_clear(xla_);
            
            auto _fp = fp;
            if(_fp<500) _fp = 500;
            restart: ; // precision restart here
            for(int k=0; k<=kla; k++) acb_mat_set(CMats[0][k],_CMats[k]);
            acb_mat_set(MatF,_MatF);
            
            for(int cn=1; cn<=xN; cn++) {
                if(!In_GiNaC_Parallel && Verbose>10) {
                    cout << "\r                                  \r" << flush;
                    cout << "     \\--C Matrix [" << _idx << "/" << las_size << "][" << cn << "/" << xN << "]" << flush;
                }
                                
                // initialize iB0, depend only on la+n
                acb_add_si(z,la_,cn,_fp); // z=la+n, alpha in B
                acb_mat_one(B0);
                acb_mat_scalar_mul_acb(B0,B0,z,_fp); // q0 is set to 1
                acb_mat_sub(B0,QxM[0],B0,_fp);
                if(!acb_mat_inv(iB0[0],B0,_fp)) throw Error("acb_mat_inv gets wrong.");
                for(int i=1; i<=kla; i++) acb_mat_mul_threaded(iB0[i], iB0[i-1], iB0[0], _fp);
                
                for(int k=0; k<=kla; k++) { // k-th row
                    int tot = (cn<s ? cn : s);
                    acb_mat_t mat_vec[tot];
                    for(int i=0; i<tot; i++) acb_mat_init(mat_vec[i],matN,nc);
                    #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
                    for(int cm=1; cm<=tot; cm++) { // sum m from 1 to s in DESS
                        acb_mat_zero(mat_vec[cm-1]);
                        acb_mat_t Bm,T,V;
                        acb_t z;
                        acb_mat_init(Bm,matN,matN);
                        acb_mat_init(T,matN,matN);
                        acb_mat_init(V,matN,nc);
                        acb_init(z);
                        acb_add_si(z,la_,cn-cm,_fp); // z = la+n-m
                        acb_mul(z,z,Qs[cm],_fp);
                        acb_mat_one(Bm);
                        acb_mat_scalar_mul_acb(Bm,Bm,z,_fp);
                        acb_mat_sub(Bm,QxM[cm],Bm,_fp);
                        for(int j=0; j<=kla; j++) { // T is T[k,j] in DESS
                            // T[k,j] = -iB0[j-k].Bm + iB0[j-k-1].Qm, q0 set 1
                            if(j<k) continue; // zero
                            acb_mat_mul_threaded(T,iB0[j-k],Bm,_fp);
                            acb_mat_neg(T,T);
                            if(j>k) acb_mat_scalar_addmul_acb(T,iB0[j-k-1],Qs[cm],_fp);
                            // get Tkj.Cj
                            acb_mat_mul_threaded(V,T,CMats[(cn-cm)%s][j],_fp);
                            acb_mat_add(mat_vec[cm-1],mat_vec[cm-1],V,_fp);
                        }
                        acb_mat_clear(Bm);
                        acb_mat_clear(T);
                        acb_mat_clear(V);
                        acb_clear(z);
                        flint_cleanup();
                    } 
                    
                    acb_mat_t CMat;
                    acb_mat_init(CMat,matN,nc);
                    acb_mat_zero(CMat);
                    for(int i=0; i<tot; i++) {
                        acb_mat_add(CMat,CMat,mat_vec[i],_fp);
                        acb_mat_clear(mat_vec[i]);
                    }
                    
                    acb_mat_set(CMats[cn%s][k],CMat); // before here, CMats should not be modified
                    if(k==0) acb_mat_add(MatF,MatF,CMat,_fp);
                    acb_mat_clear(CMat);
                    
                    if(k==0) { // precision check
                        mag_t mag;
                        mag_init(mag);
                        bool ok = true;
                        for(int r=0; r<matN; r++) for(int c=0; c<nc; c++) {
                            auto item = acb_mat_entry(MatF,r,c);
                            auto ri = acb_realref(item);
                            if(arb_rel_error_bits(ri)<-rel_fp) continue;
                            arb_get_mag(mag,ri);
                            if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                //cout << endl; arb_printd(ri,5); cout << endl; 
                                ok = false;
                                goto done; 
                            }
                            if(arb_contains_si(ri,0)) arb_zero(ri);
                            ri = acb_imagref(item);
                            if(arb_rel_error_bits(ri)<-rel_fp) continue;
                            arb_get_mag(mag,ri);
                            if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                                //cout << endl; arb_printd(ri,5); cout << endl; 
                                ok = false;
                                goto done; 
                            }
                            if(arb_contains_si(ri,0)) arb_zero(ri);
                        }
                        done: ;
                        mag_clear(mag); 
                        if(!ok) {
                            _fp += 500;
                            flint_cleanup();
                            if(_fp>fp+10000) throw Error("still too low precision with high fp!");
                            goto restart;
                        }
                    }
                }
            }
            
            for(auto z : z_cls) acb_clear(z);
            for(auto m : m_cls) acb_mat_clear(m);
            if(!In_GiNaC_Parallel && Verbose>10 && xN>0) cout << endl;
            fp = _fp;
        }
        
        for(int i=0; i<=s; i++) {
            acb_clear(Qs[i]);
            acb_mat_clear(QxM[i]);
        }
        
        auto mat = acb_to_mat(MatF,_dp);
        acb_mat_clear(MatF);
        acb_mat_clear(_MatF);
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series W.C. @ " << now() << endl;
        return mat;
    }
    
    // together with boundary constant
    matrix NDEH::Taylor(matrix C, const ex & x0, const ex & dx, const int xN) {  
        if(mx.n<1) mx.init(Mat);
        if(mx.s<0) {    
            if(!In_GiNaC_Parallel && Verbose>10) cout << "     \\--LCM " << flush; 
            mx.lcm();
            if(!In_GiNaC_Parallel && Verbose>10) cout << "\r                 \r     \\--LCM s = " << mx.s << endl;
        }
                
        auto _fp = fp+10000;
        auto _dp = dp+5000;
        int matN = mx.n;        
        acb_t z,z0,dz,dzn;
        acb_init(z); acb_init(z0); acb_init(dz); acb_init(dzn);
        to_acb(dz,dx,_dp);

        acb_poly_t Qx, M[matN][matN];
        acb_poly_init(Qx);
        acb_poly_set_fmpz_poly(Qx,mx.Qx,_fp);
        acb_poly_shift_left(Qx,Qx,1); // note Q is the denominator of x*M after lcm, Mo = M'/(x*Q)
        to_acb(z,x0,_dp);
        acb_poly_taylor_shift(Qx,Qx,z,_fp);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
            acb_poly_init(M[r][c]);
            auto num = fmpz_poly_q_numref(mx.M[r][c]);
            acb_poly_set_fmpz_poly(M[r][c],num,_fp); 
            acb_poly_taylor_shift(M[r][c],M[r][c],z,_fp);
        }
        
        vector<acb_mat_struct*> mat_to_clear; // acb_mat_t to be cleared
        
        int s = mx.s+1; // note that we need s+1 now
        if(s>xN) s = xN;
        if(s<1) s = 1;
        
        acb_t Qs[s+1];
        acb_mat_t QxM[s+1];
        acb_one(dzn);
        for(int i=0; i<=s; i++) {
            acb_init(Qs[i]);
            acb_mat_init(QxM[i],matN,matN);
            mat_to_clear.push_back(QxM[i]);
            acb_poly_get_coeff_acb(Qs[i], Qx, i);
            acb_mul(Qs[i],Qs[i],dzn,_fp); // x -> dx*y
            if(i>0) { // due to x*M, QxM[0] is zero 
                for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                    acb_poly_get_coeff_acb(acb_mat_entry(QxM[i],r,c), M[r][c], i-1); // note i-1 here, because of x*M
                }
                acb_mat_scalar_mul_acb(QxM[i],QxM[i],dzn,_fp); // x -> dx*y
            }
            acb_mul(dzn,dzn,dz,_fp);
        }
        
        // clear Qx & M
        acb_poly_clear(Qx);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) acb_poly_clear(M[r][c]);
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Taylor W.C. @ " << NN(x0,2) << endl;
                
        int nc = C.cols();
        acb_mat_t MatF, _MatF;
        acb_mat_init(MatF,matN,nc); // column vector
        acb_mat_init(_MatF,matN,nc); // column vector
        
        // setup CMats & init CMats[0]
        acb_mat_t CMats[s];
        for(int i=0; i<s; i++) {
            acb_mat_init(CMats[i],matN,nc);
            mat_to_clear.push_back(CMats[i]);
        }
        to_acb(_MatF,C,_dp);

        _fp = fp;
        if(_fp<500) _fp = 500;
        restart: ;
        acb_mat_set(CMats[0],_MatF);
        acb_mat_set(MatF,_MatF);
        for(int cn=1; cn<=xN; cn++) {
            if(!In_GiNaC_Parallel && Verbose>10) {
                cout << "\r                                  \r" << flush;
                cout << "     \\--C Matrix [" << cn << "/" << xN << "]" << flush;
            }
            
            acb_set(z0,Qs[0]);
            acb_mul_si(z0,z0,cn,_fp);
            acb_inv(z0,z0,_fp); // z0 = -iB0 = 1/q0*n
            
            auto tot = (cn<s ? cn : s);
            acb_mat_t mat_vec[tot];
            for(int i=0; i<tot; i++) acb_mat_init(mat_vec[i],matN,nc);
            #pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
            for(int cm=1; cm<=tot; cm++) { // sum m from 1 to s in DESS
                acb_mat_t Bm;
                acb_mat_init(Bm,matN,matN);
                acb_t z;
                acb_init(z);
                acb_set(z,Qs[cm]);
                acb_mul_si(z,z,cn-cm,_fp);
                acb_mat_one(Bm);
                acb_mat_scalar_mul_acb(Bm,Bm,z,_fp);
                acb_mat_sub(Bm,QxM[cm],Bm,_fp); // Bm = QxM[m]-(n-m)*qm
                acb_mat_scalar_mul_acb(Bm,Bm,z0,_fp); // Bm->T = -iB0.Bm
                // get Tkj.Cj
                acb_mat_mul_threaded(mat_vec[cm-1],Bm,CMats[(cn-cm)%s],_fp); // Bm -> T
                acb_mat_clear(Bm);
                acb_clear(z);
                flint_cleanup();
            } 
            
            acb_mat_t CMat;
            acb_mat_init(CMat,matN,nc);
            acb_mat_zero(CMat);
            for(int i=0; i<tot; i++) {
                acb_mat_add(CMat,CMat,mat_vec[i],_fp);
                acb_mat_clear(mat_vec[i]);
            }
            
            acb_mat_set(CMats[cn%s],CMat); // before here, CMats should not be modified
            acb_mat_add(MatF,MatF,CMat,_fp);
            acb_mat_clear(CMat);

            mag_t mag;
            mag_init(mag);
            bool ok = true;
            for(int r=0; r<matN; r++) for(int c=0; c<nc; c++) {
                auto item = acb_mat_entry(MatF,r,c);
                auto ri = acb_realref(item);
                if(arb_rel_error_bits(ri)<-rel_fp) continue;
                arb_get_mag(mag,ri);
                if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                    //cout << endl; arb_printd(ri,5); cout << endl; 
                    ok = false;
                    goto done; 
                }
                if(arb_contains_si(ri,0)) arb_zero(ri);
                ri = acb_imagref(item);
                if(arb_rel_error_bits(ri)<-rel_fp) continue;
                arb_get_mag(mag,ri);
                if(mag_cmp_2exp_si(mag,-abs_fp)>0) { 
                    //cout << endl; arb_printd(ri,5); cout << endl; 
                    ok = false;
                    goto done; 
                }
                if(arb_contains_si(ri,0)) arb_zero(ri);
            }
            done: ;
            mag_clear(mag); 
            if(!ok) {
                _fp += 500;
                flint_cleanup();
                if(_fp>fp+10000) throw Error("still too low precision with high fp!");
                goto restart;
            }
        }
        
        acb_clear(z); acb_clear(z0); acb_clear(dz); acb_clear(dzn);
        for(int i=0; i<=s; i++) acb_clear(Qs[i]);
        for(auto m : mat_to_clear) acb_mat_clear(m);
        mat_to_clear.clear();
        if(!In_GiNaC_Parallel && Verbose>10 && xN>0) cout << endl;
                
        auto mat = acb_to_mat(MatF,_dp);
        acb_mat_clear(MatF);
        acb_mat_clear(_MatF);
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Taylor W.C. @ " << now() << endl;
        fp = _fp;
        return mat;
    }
    
    void NDEH::info() {
        if(Mat.nops()<1) return;
        
        int pr = prank(Mat,x);
        matrix a0 = a0_matrix(Mat, x, pr);
        auto ev_map = eigenvalues(a0);
        vector<exvector> ev_groups;
        for(auto &kv : ev_map) {
            auto ev = kv.first;
            bool in_g = false;
            for(auto &gi : ev_groups) {
                if(in_g) break;
                for(auto &ev_i : gi) {
                    auto diff_ev = normal(ev_i-ev);
                    if(is_a<numeric>(diff_ev) && ex_to<numeric>(diff_ev).is_integer()) {
                        in_g = true;
                        gi.push_back(ev);
                        break;
                    }
                }
            }
            if(!in_g) {
                exvector gi_new;
                gi_new.push_back(ev);
                ev_groups.push_back(gi_new);
            }
        }
        
        ostringstream ostr;
        for(auto &gi : ev_groups) {
            sort(gi.begin(),gi.end(),[&](const auto &a, const auto &b){
                return normal(a-b).info(info_flags::positive);
            });
            for(auto &ev_i : gi) {
                ostr << ev_i << "[" << ev_map[ev_i] << "], ";
            }
        }
        string eg_str = ostr.str();
        eg_str = eg_str.substr(0,eg_str.size()-2);
        
        ex den = matrix_den_lcm(Mat);
        exvector fvec;
        if(is_a<mul>(den)) for (const auto &f : den) fvec.push_back(f);
        else fvec.push_back(den);
        exset roots;
        for (const auto &f : fvec) {
            ex b = f; 
            int n = 1;
            if (is_a<power>(f)) {
                b = f.op(0);
                n = ex_to<numeric>(f.op(1)).to_int();
            } 
            int deg = b.degree(x);
            if (deg == 0) { }
            else if (deg == 1) {
                ex c0 = b.coeff(x,0);
                ex c1 = b.coeff(x,1);
                roots.insert(normal(-c0/c1));
            } else throw Error("Roots: higher powers found.");
        }
        lst rs;
        for(auto r : roots) rs.append(r);
        sort_lst(rs);
        ostr.clear();
        ostr.str("");
        for(auto r : rs) ostr << r << ", ";
        string rs_str = ostr.str();
        rs_str = rs_str.substr(0,rs_str.size()-2);
        
        cout << "---------------------------------------" << endl;
        cout << "Info Summary: " << endl;
        cout << "---------------------------------------" << endl;
        cout << "  Poincare Rank: " << pr << endl;
        cout << "  Roots of Mat: " << rs_str << endl;
        cout << "  EigenValues of Ao: " << eg_str << endl;
        cout << "---------------------------------------" << endl;
    }
    
    void NDEH::Reset() {
        ST.Reset();
        TT.T.clear();
    }
    
}

