/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {

    using namespace EoD;
    
    void SeriesT::Resize() {
        auto nla = las.size();
        T.clear();
        T.resize(nla);
        for(int idx=0; idx<nla; idx++) {
            int kla = K[las[idx]];
            vector<vector<vector<matrix>>> vcm(s+1);
            for(int cm=0; cm<=s; cm++) {
                vector<vector<matrix>> vi(kla+1);
                for(int i=0; i<=kla; i++) {
                    vector<matrix> vj(kla+1);
                    vi[i] = vj;
                }
                vcm[cm] = vi;
            }
            T[idx] = vcm;
        }
    }
    void SeriesT::Reset() {
        las.clear();
        K.clear();
        C0.clear();
    }
    
    DE::DE(const symbol & _x) : x(_x), scn("cn"), a("a") { }
    DE::DE(const matrix & _mat, const symbol & _x) : Mat(_mat), x(_x), scn("cn"), a("a") { }
    DE::DE(const symbol & _x, const matrix & _mat) : Mat(_mat), x(_x), scn("cn"), a("a") { }
    DE::DE(const DE & b) : Mat(b.Mat), x(b.x), Ts(b.Ts), scn("cn"), a("a") { }
        
    // Dx J = M.J ---> J = T.J' & Dx J' = M'.J' with M' = Ti.M.T - Ti.Dx T
    void DE::Apply(const matrix & t, bool st) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        if(st) Ts.push_back(t);
        Digits = oDigits;
    }
    
    void DE::Apply(const lst & diag, bool st) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        auto matN = diag.nops();
        matrix t(matN, matN);
        for(int i=0; i<matN; i++) t(i,i) = diag.op(i);
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        if(st) Ts.push_back(t);
        Digits = oDigits;
    }
    
    void DE::x2y(const ex & y) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
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
        Digits = oDigits;
    }
    
    void DE::xpow() {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        int n = Mat.nops();
        for(int i=0; i<n; i++) Mat.let_op(i) = HepLib::xpow(Mat.op(i),x);
        Digits = oDigits;
    }
    
    matrix DE::MatT() {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        int matN = Mat.rows();
        matrix t = ex_to<matrix>(unit_matrix(matN));
        for(auto ti : Ts) t = t.mul(ti);
        Digits = oDigits;
        return t;
    }
        
    void DE::subs(const ex & sub, unsigned opt) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
        Digits = oDigits;
    }
    
    void DE::subs(const exmap & sub, unsigned opt) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
        Digits = oDigits;
    }
    
    void DE::subs(const lst & l, const lst & r, unsigned opt) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(l, r, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(l, r, opt);
        }
        Digits = oDigits;
    }
    
    // fuchsify at x=0
    void DE::Fuchsify() {
        matrix m = Mat;
        auto a01 = a01_matrix(m, x);
        matrix a0 = a01.first;
        matrix a1 = a01.second;
        symbol lambda("L");
        matrix t= imatrix(m.rows());
        
        int pr = prank(m,x);
        if(pr < 1) {
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsify: Poincare rank is alreay " << pr << endl;
            Ts.push_back(t);
            return;
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsify: start @ " << now() << endl;
        
        while(true) {
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--Poincare rank: " << pr << ",  A0 rank: " << a0.rank() << endl;
            
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
            
            m = normal(with_balance_t(m, p, x));
            t = t.mul(imatrix(m.rows()).sub(p).sub(p.mul_scalar(1/x)));
            
            pr = prank(m,x);
            if (pr < 1) break;
            a01 = a01_matrix(m, x);
            a0 = a01.first;
            a1 = a01.second;
        }

        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsify: finished @ " << now() << endl;
        
        Mat = m;
        Ts.push_back(t);
    }
    
    // normalization at x=0 // TODO: still under developement
    void DE::Normalize() {
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Normalize: start @ " << now() << endl;

        matrix m = Mat;
        matrix t= imatrix(m.rows());
        if(prank(m,x)!=0) throw Error("Normalize: prank is NOT 0.");
        
        while(true) {
        
            matrix a0 = a0_matrix(m, x, 0);
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
            
            if(!In_GiNaC_Parallel && Verbose>0) cout << "     \\--EigenValues of A0 summary:" << endl;
            ostringstream ostr;
            ostr << "        ";
            for(auto &gi : ev_groups) {
                sort(gi.begin(),gi.end(),[&](const auto &a, const auto &b){
                    return normal(a-b).info(info_flags::positive);
                });
                for(auto &ev_i : gi) {
                    ostr << ev_i << "[" << ev_map[ev_i] << "], ";
                }
            }
            if(!In_GiNaC_Parallel && Verbose>1) {
                string str = ostr.str();
                cout << str.substr(0,str.size()-2) << endl;
            }
            
            bool is_normalized = true;
            for(auto &gi : ev_groups) {
                if(gi.size()<2) continue;
                is_normalized = false;
                auto last_ev = gi[gi.size()-1];
                for_each(gi.begin(), gi.end()-1, [&](const auto &ev) {
                    if((last_ev-ev).info(info_flags::positive)) {
                        if(!In_GiNaC_Parallel && Verbose>0) cout << "     \\--Raising ev = " << ev << " ..." << endl;
                        matrix vec_u = eigenvectors(a0, ev)[0];
                        matrix vec_v  = dual_basis(vec_u)[0];
                        matrix p = normal(vec_u.mul(vec_v));
                        m = with_balance_t(m, p, x);
                        t = t.mul(imatrix(m.rows()).sub(p).sub(p.mul_scalar(1/x)));
                    } else {
                        if(!In_GiNaC_Parallel && Verbose>0) cout << "     \\--Lowering ev = " << ev << " ..." << endl;
                        matrix vec_v = eigenvectors_left(a0, ev)[0];
                        matrix vec_u  = dual_basis_left(vec_v)[0];
                        matrix p = normal(vec_u.mul(vec_v));
                        m = with_balance_ti(m, p, x);
                        t = t.mul(imatrix(m.rows()).sub(p).sub(p.mul_scalar(x)));
                    }
                    m = normal(m);
                    t = normal(t);
                    a0 = a0_matrix(m, x, 0);
                    return 0;
                });
            }
            
            if(is_normalized) break;        
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Normalize: fininshed @ " << now() << endl;

        Mat = m;
        Ts.push_back(t);
    }
    
    // shearing transformation at x=0
    void DE::Shear() {
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: start @ " << now() << endl;

        matrix m = Mat;
        matrix t = imatrix(m.rows());
        matrix ti = imatrix(m.rows());
        if(prank(m,x)>0) throw Error("Shear: prank > 0.");
        
        matrix a0 = a0_matrix(m, x, 0);
        while(true) {
        
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--EigenValues of A0 summary:" << endl;
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
            
            bool is_normalized = true;
            for(auto &gi : ev_groups) {
                if(gi.size()>1) {
                    is_normalized = false;
                    break;
                }
            }
            
            int max_dep = -1;
            ostringstream ostr;
            ostr << "        ";
            for(auto &gi : ev_groups) {
                sort(gi.begin(),gi.end(),[&](const auto &a, const auto &b){
                    return normal(a-b).info(info_flags::positive);
                });
                if(!In_GiNaC_Parallel && Verbose>1) {
                    for(auto &ev_i : gi) {
                        ostr << ev_i << "[" << ev_map[ev_i] << "], ";
                    }
                }
                int tmp = ex_to<numeric>(gi[0]-gi[gi.size()-1]).to_int();
                if(max_dep < tmp) max_dep = tmp;
            }
            if(!In_GiNaC_Parallel && Verbose>1) {
                string str = ostr.str();
                cout << str.substr(0,str.size()-2) << endl;
            }
            
            if(is_normalized) break;
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--jordan decomposition ..." << endl;
            auto qj = jordan(a0);
            matrix smat = imatrix(m.rows());
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--shearing transformation ..." << endl;
            int cpos = 0;
            
            bool s_only_first = false;
            if(s_only_first) {
                for(auto kv : qj.second) {
                    bool is_first = false;
                    for(auto &gi : ev_groups) {
                        if(normal(gi[0] - kv.first).is_zero()) {
                            is_first = true;
                            break;
                        }
                    }
                    if(is_first) for(int j=0; j<kv.second; j++) smat(cpos+j, cpos+j) = x;
                    cpos += kv.second;
                }
            } else {
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
            }
            
            auto tm = qj.first.mul(smat);
            auto tmi = tm.inverse();
            m = transform(m, tmi, tm, x);
            t = t.mul(tm);
            ti = tmi.mul(ti);
            
            matrix_map_inplace(m, [this, &max_dep](const ex & e) { return series_ex(e, x, max_dep); });
            
            a0 = matrix_map(m, [this](const ex & e) { return normal(e.collect(x).coeff(x,-1)); });
            if(a0.has(x)) throw Error("a0 still has x.");
                    
            matrix_map_inplace(m, [](const ex & e) { return normal(e); });
        }
        
        auto qj = jordan(a0);
        matrix qi = normal(qj.first.inverse());
        a0 = normal(qi.mul(a0).mul(qj.first));
        t = t.mul(qj.first);
        ti = qi.mul(ti);
        matrix_map_inplace(t, [](const ex & e) { return normal(e); });
        matrix_map_inplace(ti, [](const ex & e) { return normal(e); });
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--DE Transformation  ..." << endl;
        m = transform(Mat, ti, t, x);
        matrix_map_inplace(m, [](const ex & e) { return normal(e); });
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: fininshed @ " << now() << endl;
        Mat = m;
        Ts.push_back(t);
    }
    
    BJF::BJF(matrix _A, ex _b, int K): A(_A), b(_b) { 
        B = ex_to<matrix>(unit_matrix(A.rows())).mul_scalar(b); 
    }
    
    matrix BJF::operator()(int i, int j) {
        if(i==j) return A;
        else if(i+1==j) return B;
        else throw Error("BJF: 0 matrix block.");
    }
       
    BJFinv::BJFinv(matrix _A, ex _b, int K): A(_A), b(_b) { 
        U = ex_to<matrix>(unit_matrix(A.rows()));
        matrix a1 = fermat_inv(A);
        Ain[1] = a1;
        matrix al = a1;
        for(int i=2; i<=K+1; i++) al = Ain[i] = al.mul(a1);
    }
    
    matrix BJFinv::operator()(int i, int j) {
        if(i>j) throw Error("BJFinv: 0 matrix block");
        pair<int, int> key = make_pair(i,j);
        auto f = cache.find(key);
        if(f!=cache.end()) return f->second;
        if(i==j) return cache[key] = Ain[1];
        else {
            int ij = j-i;
            return cache[key] = U.mul_scalar(pow(-b,ij)).mul(Ain[ij+1]);
        }
    }
    
    CMatrix DE::Series(const unsigned int xN) { 
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series" << endl;
        
        int matN = Mat.rows();
        
        if(!ST.inited) { // for cache
            if(prank(Mat,x)>0) Fuchsify();
            matrix a0 = normal(a0_matrix(Mat, x, 0));
            if(!is_jordan_form(a0)) {
                Shear();
                a0 = normal(a0_matrix(Mat, x, 0));
            }
            auto m = Mat;
            
            lst _las;
            vector<pair<ex,int>> js;            
            int lastN = 0;
            for(int r=0; r<a0.rows(); r++) {
                lastN++;
                if(r==a0.rows()-1 || a0(r,r+1)==0) { // the last row
                    auto la = a0(r,r);
                    _las.append(la);
                    js.push_back(make_pair(la, lastN));
                    if(lastN-1 > ST.K[la]) ST.K[la] = lastN-1;
                    lastN = 0;
                }
            }
            _las.sort().unique();
            sort_lst(_las);
            for(auto la : _las) {
                las.push_back(la);
                ST.las.push_back(la);
            }
            
            ST.C0.clear();
            for(int idx=0; idx<ST.las.size(); idx++) {
                auto la = ST.las[idx];
                int kla = ST.K[la];
                vector<matrix> vm;
                for(int li=0; li<=kla; li++) vm.push_back(matrix(matN, matN));
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
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--Q polynormial ..." << endl;
            
            matrix_map_inplace(m, [&](const ex & e){ return normal(e); });
            ex qlcm = matrix_den_lcm(m);
            if(normal(qlcm.subs(x==0,nopat)).is_zero()) qlcm = normal(qlcm/x); //make sure check
            matrix_map_inplace(m, [&](const ex & e){ return collect_ex(normal(qlcm*x*e),x); });
            
            // make sure m is polynomial in x
            for(int i=0; i<m.nops(); i++) if(!m.op(i).is_polynomial(x)) {
                cout << m << endl;
                throw Error("DE::Series, matrix is NOT a polynomial in x.");
            }
            
            qlcm = collect_ex(qlcm, x);
            int qmax = qlcm.degree(x);
            exvector qs;
            for(int i=0; i<=qmax; i++) {
                qs.push_back(qlcm.coeff(x,i));
            }
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--B Matrix ..." << endl;
            
            ST.s = qmax; // s in DESS
            for(int i=0; i<m.nops(); i++) {
                auto tmp = m.op(i).degree(x);
                if(ST.s<tmp) ST.s = tmp;
            }
            for(int i=0; i<matN; i++) m(i,i) -= a*qlcm;

            for(int i=qmax+1; i<=ST.s; i++) qs.push_back(0);
            vector<matrix> BMat(ST.s+1);
            for(int i=0; i<=ST.s; i++) {
                BMat[i] = matrix_map(m, [&](const ex & e){ return e.coeff(x,i); });
            }
            
            ST.Resize(); // resize vector
            
            // init ST.T
            for(int idx=0; idx<ST.las.size(); idx++) {
                auto la = ST.las[idx];
                int kla = ST.K[la];
                matrix zmat(matN, matN); // zero matrix
                for(int cm=1; cm<=ST.s; cm++) 
                for(int i=0; i<=kla; i++) 
                for(int j=0; j<=kla; j++) ST.T[idx][cm][i][j] = zmat;
                
                matrix B0(matN, matN);
                for(int i=0; i<matN*matN; i++) B0.let_op(i) = BMat[0].op(i).subs(a==la+scn,nopat);
                BJFinv iBJF(B0, -qs[0], kla);
                                
                for(int cm=1; cm<=ST.s; cm++) {
                    matrix Bm(matN, matN);
                    for(int i=0; i<matN*matN; i++) Bm.let_op(i) = BMat[cm].op(i).subs(a==la+scn-cm,nopat);
                    BJF mBJF(Bm, -qs[cm], kla);
                    
                    for(int k=0; k<=kla; k++) 
                    for(int j=0; j<=kla; j++) 
                    for(int i=0; i<=kla; i++) {
                        if(k>i) continue; // iBJF(k,i)
                        if(i!=j && i+1!=j) continue; // mBJF(i,j)
                        ST.T[idx][cm][k][j] = ST.T[idx][cm][k][j].sub(iBJF(k,i).mul(mBJF(i,j))); 
                        // T(n,m) = -iBJF(B0,-q0).mBJF(Bm,-qm)
                    }
                }
            }
            ST.inited = true;
        } 
        
        bool CMat_Parallel = (GiNaC_Parallel_NP.find("CMat")!=GiNaC_Parallel_NP.end()) && (GiNaC_Parallel_NP["CMat"]>0);
        if(CMat_Parallel && !In_GiNaC_Parallel && Verbose>1) cout << "     \\--C Matrix ..." << endl;
        
        int s = ST.s;
        if(s>xN) s = xN;
        if(s<1) s = 1;
        CMatrix CMatF; // CMatF[la][k][n] : coefficient of x^la*log(x)^k/k!
        for(int idx=0; idx<ST.las.size(); idx++) {
            auto la = ST.las[idx];
            int kla = ST.K[la];
            CMatF[la].resize(kla+1);
            matrix CMats[s][kla+1];
            for(int k=0; k<=kla; k++) {
                auto t = ST.C0[la][k];
                CMats[0][k] = t;
                CMatF[la][k].resize(xN+1);
                CMatF[la][k][0] = t;
            }
                        
            int nidx = 0;
            for(int cn=1; cn<=xN; cn++) {
                if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1) {
                    cout << "\r                                  \r" << flush;
                    cout << "     \\--C Matrix [" << idx+1 << "/" << ST.las.size() << "][" << cn << "/" << xN << "]" << flush;
                }
                exmap smap;
                smap[scn] = cn;
                //if(WDigits>0) smap[scn] = smap[scn].evalf();
                nidx = (nidx+1)%s;
                for(int k=0; k<=kla; k++) {
                    matrix cmat(matN, matN);
                    for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                        for(int j=0; j<=kla; j++) {
                            matrix Tkj = ST.T[idx][cm][k][j];
                            for(int i=0; i<Tkj.nops(); i++) Tkj.let_op(i) = Tkj.op(i).subs(smap,nopat);
                            cmat = cmat.add(Tkj.mul(CMats[(nidx-cm+s)%s][j])); 
                        }
                    }
                    if(WDigits<0) matrix_map_inplace(cmat, [&](const ex & e) { return normal(e); });
                    CMats[nidx][k] = cmat;
                    CMatF[la][k][cn] = cmat;
                }
            }
            if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1 && xN>0) cout << endl;
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << now() << endl;
        Digits = oDigits;
        return CMatF;
    }
    
    matrix DE::Series(const ex & x0, const unsigned int xN, const lst & _las) { 
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << NN(x0,2) << endl;
        
        int matN = Mat.rows();
        
        if(!ST.inited) { // for cache
            if(prank(Mat,x)>0) Fuchsify();
            matrix a0 = normal(a0_matrix(Mat, x, 0));
            if(!is_jordan_form(a0)) {
                Shear();
                a0 = normal(a0_matrix(Mat, x, 0));
            } else {
                for(int i=0; i<matN; i++) for(int j=i+1; j<matN; j++) {
                    ex diff = normal(a0(i,i)-a0(j,j));
                    if(!is_zero(diff) && diff.info(info_flags::integer)) {
                        Shear();
                        a0 = normal(a0_matrix(Mat, x, 0));
                    }
                }
            }
            
            auto m = Mat;
            
            lst _las;
            vector<pair<ex,int>> js;            
            int lastN = 0;
            for(int r=0; r<a0.rows(); r++) {
                lastN++;
                if(r==a0.rows()-1 || a0(r,r+1)==0) { // the last row
                    auto la = a0(r,r);
                    _las.append(la);
                    js.push_back(make_pair(la, lastN));
                    if(lastN-1 > ST.K[la]) ST.K[la] = lastN-1;
                    lastN = 0;
                }
            }
            _las.sort().unique();
            sort_lst(_las);
            for(auto la : _las) {
                las.push_back(la);
                ST.las.push_back(la);
            }
            
            ST.C0.clear();
            for(int idx=0; idx<ST.las.size(); idx++) {
                auto la = ST.las[idx];
                int kla = ST.K[la];
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
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--Q polynormial ..." << endl;
            
            matrix_map_inplace(m, [&](const ex & e){ return normal(e); });
            ex qlcm = matrix_den_lcm(m);
            if(normal(qlcm.subs(x==0,nopat)).is_zero()) qlcm = normal(qlcm/x); //make sure check
            matrix_map_inplace(m, [&](const ex & e){ return collect_ex(normal(qlcm*x*e),x); });
            
            // make sure m is polynomial in x
            for(int i=0; i<m.nops(); i++) if(!m.op(i).is_polynomial(x)) {
                cout << m << endl;
                throw Error("DE::Series, matrix is NOT a polynomial in x.");
            }
            
            qlcm = collect_ex(qlcm, x);
            int qmax = qlcm.degree(x);
            exvector qs;
            for(int i=0; i<=qmax; i++) {
                qs.push_back(qlcm.coeff(x,i));
            }
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--B Matrix ..." << endl;
            
            ST.s = qmax; // s in DESS
            for(int i=0; i<m.nops(); i++) {
                auto tmp = m.op(i).degree(x);
                if(ST.s<tmp) ST.s = tmp;
            }
            for(int i=0; i<matN; i++) m(i,i) -= a*qlcm;

            for(int i=qmax+1; i<=ST.s; i++) qs.push_back(0);
            vector<matrix> BMat(ST.s+1);
            for(int i=0; i<=ST.s; i++) {
                BMat[i] = matrix_map(m, [&](const ex & e){ return e.coeff(x,i); });
            }
            
            ST.Resize(); // resize vector
            
            // init ST.T
            for(int idx=0; idx<ST.las.size(); idx++) {
                auto la = ST.las[idx];
cout << la << endl;
                int kla = ST.K[la];
                matrix zmat(matN, matN); // zero matrix
                for(int cm=1; cm<=ST.s; cm++) 
                for(int i=0; i<=kla; i++) 
                for(int j=0; j<=kla; j++) ST.T[idx][cm][i][j] = zmat;
                
                matrix B0(matN, matN);
                for(int i=0; i<matN*matN; i++) B0.let_op(i) = BMat[0].op(i).subs(a==la+scn,nopat);
                BJFinv iBJF(B0, -qs[0], kla);
cout << B0 << endl;
                                
                for(int cm=1; cm<=ST.s; cm++) {
                    matrix Bm(matN, matN);
                    for(int i=0; i<matN*matN; i++) Bm.let_op(i) = BMat[cm].op(i).subs(a==la+scn-cm,nopat);
                    BJF mBJF(Bm, -qs[cm], kla);
                    
                    for(int k=0; k<=kla; k++) 
                    for(int j=0; j<=kla; j++) 
                    for(int i=0; i<=kla; i++) {
                        if(k>i) continue; // iBJF(k,i)
                        if(i!=j && i+1!=j) continue; // mBJF(i,j)
cout << iBJF(k,i) << endl << endl;
                        ST.T[idx][cm][k][j] = ST.T[idx][cm][k][j].sub(iBJF(k,i).mul(mBJF(i,j))); 
                        // T(n,m) = -iBJF(B0,-q0).mBJF(Bm,-qm)
                    }
                }
            }
            ST.inited = true;
        } 
        
        bool CMat_Parallel = (GiNaC_Parallel_NP.find("CMat")!=GiNaC_Parallel_NP.end()) && (GiNaC_Parallel_NP["CMat"]>0);
        if(CMat_Parallel && !In_GiNaC_Parallel && Verbose>1) cout << "     \\--C Matrix ..." << endl;
        
        int s = ST.s;
        if(s>xN) s = xN;
        if(s<1) s = 1;
        matrix MatF(matN, matN);
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
            int kla = ST.K[la];
            matrix CMats[s][kla+1];
            for(int k=0; k<=kla; k++) CMats[0][k] = ST.C0[la][k];
            
            for(int k=0; k<=kla; k++) {
                auto xterm = pow(x0,la)*pow(log(x0),k)/factorial(k);
                if(WDigits>0) xterm = xterm.evalf(); // OK
                MatF = MatF.add(CMats[0][k].mul_scalar(xterm));
            }
            
            int nidx = 0;
            for(int cn=1; cn<=xN; cn++) {
                if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1) {
                    cout << "\r                                  \r" << flush;
                    cout << "     \\--C Matrix [" << idx+1 << "/" << ST.las.size() << "][" << cn << "/" << xN << "]" << flush;
                }
                exmap smap;
                smap[scn] = cn;
                if(WDigits>0) smap[scn] = smap[scn].evalf();
                nidx = (nidx+1)%s;
                for(int k=0; k<=kla; k++) {
                    matrix cmat(matN, matN);
                    for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                        for(int j=0; j<=kla; j++) {
                            matrix Tkj = ST.T[idx][cm][k][j];
                            for(int i=0; i<Tkj.nops(); i++) Tkj.let_op(i) = Tkj.op(i).subs(smap,nopat);
                            cmat = cmat.add(Tkj.mul(CMats[(nidx-cm+s)%s][j])); 
                        }
                    }
                    if(WDigits<0) matrix_map_inplace(cmat, [&](const ex & e) { return normal(e); });
                    CMats[nidx][k] = cmat;
                    auto xterm = pow(x0,la)*pow(x0,cn)*pow(log(x0),k)/factorial(k);
                    if(WDigits>0) xterm = xterm.evalf(); // OK
                    MatF = MatF.add(cmat.mul_scalar(xterm));
                }

            }
            if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1 && xN>0) cout << endl;
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << now() << endl;
        Digits = oDigits;
        return MatF;
    }
    
    matrix DE::Taylor(const ex & x0, const ex & dx, const unsigned int xN) {        
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Taylor @ " << NN(x0,2) << endl;
        
        auto m = Mat;
        int matN = m.rows();
        
        if(!TT.inited) {
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--Q polynormial ..." << endl;
            
            matrix_map_inplace(m, [&](const ex & e){ return normal(e); });
            ex qlcm = matrix_den_lcm(m);
            if(normal(qlcm.subs(x==x0,nopat)).is_zero()) //make sure check
                throw Error("DE::Taylor only works for non-singualr point."); 
                
            matrix_map_inplace(m, [&](const ex & e) {
                ex res = normal(qlcm*e);
                res = res.subs(x==x+x0);
                res = collect_ex(x*res,x);
                return res; 
            });
            
            // make sure m is polynomial in x
            for(int i=0; i<m.nops(); i++) if(!m.op(i).is_polynomial(x)) {
                cout << m.op(i) << endl;
                throw Error("DE::Taylor, matrix is NOT a polynomial in x.");
            }
            qlcm = qlcm.subs(x==x+x0,nopat);
            qlcm = collect_ex(qlcm, x);
            if(!qlcm.is_polynomial(x)) throw Error("DE::Taylor, qlcm is NOT a polynomial in x.");
            int qmax = qlcm.degree(x);
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--B Matrix ..." << endl;
            
            TT.s = qmax; // s in DESS
            for(int i=0; i<m.nops(); i++) {
                auto tmp = m.op(i).degree(x);
                if(TT.s<tmp) TT.s = tmp;
            }
            for(int i=0; i<matN; i++) m(i,i) -= a*qlcm;

            vector<matrix> BMat(TT.s+1);
            for(int i=0; i<=TT.s; i++) {
                BMat[i] = matrix_map(m, [&](const ex & e){ return e.coeff(x,i); });
            }
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--T Matrix ..." << endl;
            
            TT.T.clear();
            TT.T.resize(TT.s+1);
            matrix zmat(matN, matN); // zero matrix
            for(int cm=1; cm<=TT.s; cm++) TT.T[cm] = zmat;
        
            matrix B0(matN, matN);
            for(int i=0; i<matN*matN; i++) B0.let_op(i) = BMat[0].op(i).subs(a==scn,nopat);
            B0 = B0.inverse();
                            
            for(int cm=1; cm<=TT.s; cm++) {
                matrix Bm(matN, matN);
                for(int i=0; i<matN*matN; i++) Bm.let_op(i) = BMat[cm].op(i).subs(a==scn-cm,nopat);
                TT.T[cm] = TT.T[cm].sub(B0.mul(Bm)); // T(n,m) = -(B0^-1).Bm
            }
            TT.inited = true;
        } 
        
        int s = TT.s;
        if(s>xN) s = xN;
        if(s<1) s = 1;
        matrix CMats[s]; // only keep last s CMats
        CMats[0] = ex_to<matrix>(unit_matrix(matN));
        matrix MatF = CMats[0];
            
        int nidx = 0; 
        for(int cn=1; cn<=xN; cn++) {
            if(!In_GiNaC_Parallel && Verbose>1) {
                cout << "\r                                  \r" << flush;
                cout << "     \\--C Matrix [" << cn << "/" << xN << "]" << flush;
            }
            exmap smap;
            smap[scn] = cn;
            if(WDigits>0) smap[scn] = smap[scn].evalf();
            matrix cmat(matN, matN);
            nidx = (nidx+1)%s;
            for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                matrix Tm = TT.T[cm];
                for(int i=0; i<Tm.nops(); i++) Tm.let_op(i) = Tm.op(i).subs(smap,nopat);
                cmat = cmat.add(Tm.mul(CMats[(nidx-cm+s)%s])); 
            }
            if(WDigits<0) matrix_map_inplace(cmat, [&](const ex & e) { return normal(e); });
            CMats[nidx] = cmat;
            
            auto xterm = pow(dx,cn);
            if(WDigits>0) xterm = xterm.evalf();
            MatF = MatF.add(cmat.mul_scalar(xterm));
        }
        if(!In_GiNaC_Parallel && Verbose>1) cout << endl;
            
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Taylor @ " << now() << endl;
        Digits = oDigits;
        return MatF;
    }
    
    void DE::info() {
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
    
    void DE::Reset() {
        ST.Reset();
        TT.T.clear();
    }

    //==============================================================================
    
    ex matrix_norm(const matrix & mat, unsigned opt) {
        if(opt==0) {
            ex res;
            for(int i=0; i<mat.nops(); i++) {
                ex t = abs(mat.op(i));
                res += t*t;
            }
            return sqrt(res);
        }
        throw Error("matrix_norm: option NOT supported yet.");
        return 0;
    }

    ex matrix_den_lcm(const matrix & m) {
        exmap pn_map;
        for(int i=0; i<m.nops(); i++) {
            auto den = m.op(i).denom();
            den = exfactor(den);
            if(!is_a<mul>(den)) den = lst{den};
            for(auto item : den) {
                ex p = item;
                ex n = 1;
                if(item.match(pow(w1,w2)) && item.op(1).info(info_flags::integer)) {
                    p = item.op(0);
                    n = item.op(1);
                }
                auto kv = pn_map.find(p);
                if(kv==pn_map.end() || kv->second<n) pn_map[p] = n;
            }
        }
        ex res = 1;
        for(auto kv : pn_map) res *= pow(kv.first, kv.second);
        return res;
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
    
    matrix PolynomialFit(const exvector & xs, const exvector & ys, unsigned int k, int k0) {
        unsigned int n = xs.size();
        if(ys.size() != n) throw Error("PolynomialFit: the size of xs is not the same as ys.");
        matrix X(k+1,n);
        for(int c=0; c<n; c++) {
            ex xp = 1;
            for(int r=0; r<=k; r++) {
                X(r,c) = xp;
                xp *= xs[c];
            }
        }
        matrix Y(n,1);
        for(int r=0; r<n; r++) {
            if(k0==0) Y(r,0) = ys[r];
            else Y(r,0) = ys[r]/pow(xs[r], k0);
        }
        auto mat = X.mul(X.transpose()).inverse().mul(X).mul(Y);
        return mat;
    }
    
    matrix C2Mat(const CMatrix & cmat, const ex & x0) {
        // Cmat[la][k][n] : coefficient of x^la*log(x)^k/k!;
        matrix mat;
        bool first = true;
        for(const auto kv : cmat) {
            ex la = kv.first;
            const auto & vvm = kv.second;
            for(int k=0; k<vvm.size(); k++) {
                const auto & vm = vvm[k];
                for(int n=0; n<vm.size(); n++) {
                    auto xterm = pow(x0,la)*pow(x0,n)*pow(log(x0),k)/factorial(k);
                    if(first) {
                        mat = vm[n].mul_scalar(xterm);
                        first = false;
                    } else mat = mat.add(vm[n].mul_scalar(xterm));
                }
            }
        }
        return mat;
    }
    
}

