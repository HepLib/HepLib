/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {

    using namespace EoD;
    
    DE::DE(const symbol & _x) : x(_x), a("a") { }
    DE::DE(const matrix & _mat, const symbol & _x) : Mat(_mat), x(_x), a("a") { }
    DE::DE(const symbol & _x, const matrix & _mat) : Mat(_mat), x(_x), a("a") { }
    DE::DE(const DE & b) : Mat(b.Mat), x(b.x), Ts(b.Ts), a("a") { }
        
    // Dx J = M.J ---> J = T.J' & Dx J' = M'.J' with M' = Ti.M.T - Ti.Dx T
    void DE::Apply(const matrix & t, bool st) {
        if(WDigits>0) set_precision(WDigits);
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        if(st) Ts.push_back(t);
        if(WDigits>0) reset_precision();
    }
    
    void DE::Apply(const lst & diag, bool st) {
        if(WDigits>0) set_precision(WDigits);
        auto matN = diag.nops();
        matrix t(matN, matN);
        for(int i=0; i<matN; i++) t(i,i) = diag.op(i);
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        if(st) Ts.push_back(t);
        if(WDigits>0) reset_precision();
    }
    
    void DE::x2y(const ex & y) {
        if(WDigits>0) set_precision(WDigits);
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
        if(WDigits>0) reset_precision();
    }
    
    void DE::xpow() {
        if(WDigits>0) set_precision(WDigits);
        int n = Mat.nops();
        for(int i=0; i<n; i++) Mat.let_op(i) = HepLib::xpow(Mat.op(i),x);
        if(WDigits>0) reset_precision();
    }
    
    matrix DE::MatT() {
        if(WDigits>0) set_precision(WDigits);
        int matN = Mat.rows();
        matrix t = ex_to<matrix>(unit_matrix(matN));
        for(auto ti : Ts) t = t.mul(ti);
        if(WDigits>0) reset_precision();
        return t;
    }
        
    void DE::subs(const ex & sub, unsigned opt) {
        if(WDigits>0) set_precision(WDigits);
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
        if(WDigits>0) reset_precision();
    }
    
    void DE::subs(const exmap & sub, unsigned opt) {
        if(WDigits>0) set_precision(WDigits);
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
        if(WDigits>0) reset_precision();
    }
    
    void DE::subs(const lst & l, const lst & r, unsigned opt) {
        if(WDigits>0) set_precision(WDigits);
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(l, r, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(l, r, opt);
        }
        if(WDigits>0) reset_precision();
    }
    
    // fuchsify at x=0
    void DE::Fuchsify() {
        matrix m = Mat;
        auto a01 = a01_matrix(m, x);
        a0 = a01.first;
        matrix a1 = a01.second;
        symbol lambda("L");
        matrix t= imatrix(m.rows());
        
        int pr = prank(m,x);
        if(pr < 1) {
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsified: Poincare rank is alreay " << pr << endl;
            Ts.push_back(t);
            return;
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsifying @ " << now() << endl;
        
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

        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Fuchsified @ " << now() << endl;
        
        Mat = m;
        a0 = a0_matrix(m, x, 0);
        Ts.push_back(t);
        fuchsified = true;
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
        fuchsified = true;
    }
    
    // shearing transformation at x=0
    void DE::Shear() {
        if(!fuchsified) Fuchsify();
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: start @ " << now() << endl;

        matrix m = Mat;
        matrix t = imatrix(m.rows());
        matrix ti = imatrix(m.rows());
        
        while(true) {
        
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--eigen values of Ao ..." << endl;
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
            auto qj = jordan(a0,ev_map);
            matrix smat = imatrix(m.rows());
            
            if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--shearing transformation ..." << endl;
            
            int cpos = 0;
            int stype = 1; 
            if(stype==0) {
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
            } else if(stype==1) {
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
            auto tmi = smat.inverse().mul(qj.first.inverse());
            m = transform(m, tmi, tm, x);
            t = t.mul(tm);
            ti = tmi.mul(ti);
            
            GiNaC_Parallel_Verb["Shear"] = 0;
            auto mat_a0_vec = GiNaC_Parallel(m.nops(), [&](int idx)->ex {
                ex mi = m.op(idx);
                mi = series_ex(mi,x,max_dep);
                mi = normal(mi);
                ex a0i = series_ex(mi,x,-1).coeff(x,-1);
                return lst{mi,a0i};
            }, "Shear");
            
            for(int i=0; i<mat_a0_vec.size(); i++) {
                m.let_op(i) = mat_a0_vec[i].op(0);
                a0.let_op(i) = mat_a0_vec[i].op(1);
            }
            if(a0.has(x)) throw Error("a0 still has x.");
        }
        
        auto qj = jordan(a0);
        matrix qi = (qj.first.inverse());
        t = t.mul(qj.first);
        ti = qi.mul(ti);
        matrix_map_inplace(t, [](const ex & e) { return normal(e); });
        matrix_map_inplace(ti, [](const ex & e) { return normal(e); });
        matrix_map_inplace(Mat, [](const ex & e) { return normal(e); });
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "     \\--DE Transformation  ..." << endl;
        m = transform(Mat, ti, t, x);
        a0 = normal(a0_matrix(m, x, 0));
        matrix_map_inplace(m, [](const ex & e) { return normal(e); });
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Shear: fininshed @ " << now() << endl;
        Mat = m;
        Ts.push_back(t);
        sheared = true;
    }
    
    void DE::STInit() { // for cache
        if(!fuchsified) Fuchsify();
        if(!sheared) Shear();
        if(!is_sheared_form(a0)) throw Error("a0 is NOT sheared form.");
        auto m = Mat;
        int matN = Mat.rows();
        
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
        
        ST.kmax = -1;
        ST.C0.clear();
        for(int idx=0; idx<ST.las.size(); idx++) {
            auto la = ST.las[idx];
            int kla = ST.K[la];
            if(ST.kmax<kla) ST.kmax = kla;
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
        int kmax = ST.kmax;
        matrix zmat(matN, matN); // zero matrix
        for(int cm=1; cm<=ST.s; cm++) 
        for(int i=0; i<=kmax; i++) 
        for(int j=0; j<=kmax; j++) ST.T[cm][i][j] = zmat;
        
        matrix B0(matN, matN);
        for(int i=0; i<matN*matN; i++) B0.let_op(i) = BMat[0].op(i);
        invBJF iBJF(B0, -qs[0], kmax);
                        
        for(int cm=1; cm<=ST.s; cm++) {
            matrix Bm(matN, matN);
            for(int i=0; i<matN*matN; i++) {
                static symbol b("b");
                Bm.let_op(i) = BMat[cm].op(i).subs(a==b-cm,nopat).subs(b==a,nopat);
            }
            BJF mBJF(Bm, -qs[cm], kmax);
            
            for(int k=0; k<=kmax; k++) 
            for(int j=0; j<=kmax; j++) 
            for(int i=0; i<=kmax; i++) {
                if(k>i) continue; // iBJF(k,i)
                if(i!=j && i+1!=j) continue; // mBJF(i,j)
                ST.T[cm][k][j] = ST.T[cm][k][j].sub(iBJF(k,i).mul(mBJF(i,j))); 
                // T(n,m) = -iBJF(B0,-q0).mBJF(Bm,-qm)
            }
        }
        ST.inited = true;
    }
    
    CMatrix DE::Series(const int xN) { 
        if(WDigits>0) set_precision(WDigits);
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series" << endl;
        
        int matN = Mat.rows();
        
        if(!ST.inited) STInit(); 
        
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
                        
            for(int cn=1; cn<=xN; cn++) {
                if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1) {
                    cout << "\r                                  \r" << flush;
                    cout << "     \\--C Matrix [" << idx+1 << "/" << ST.las.size() << "][" << cn << "/" << xN << "]" << flush;
                }
                exmap smap;
                smap[a] = la+cn;
                //if(WDigits>0) smap[a] = smap[a].evalf();
                for(int k=0; k<=kla; k++) {
                    matrix cmat(matN, matN);
                    for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                        for(int j=0; j<=kla; j++) {
                            matrix Tkj = ST.T[cm][k][j];
                            for(int i=0; i<Tkj.nops(); i++) Tkj.let_op(i) = Tkj.op(i).subs(smap,nopat);
                            cmat = cmat.add(Tkj.mul(CMats[(cn-cm)%s][j])); 
                        }
                    }
                    if(WDigits<0) matrix_map_inplace(cmat, [&](const ex & e) { return normal(e); });
                    CMats[cn%s][k] = cmat;
                    CMatF[la][k][cn] = cmat;
                }
            }
            if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1 && xN>0) cout << endl;
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << now() << endl;
        if(WDigits>0) reset_precision();
        return CMatF;
    }
    
    matrix DE::Series(const ex & x0, const int xN, const lst & _las) { 
        if(WDigits>0) set_precision(WDigits);
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << NN(x0,2) << endl;
        
        int matN = Mat.rows();
        if(!ST.inited) STInit();  
        
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
            
            for(int cn=1; cn<=xN; cn++) {
                if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1) {
                    cout << "\r                                  \r" << flush;
                    cout << "     \\--C Matrix [" << idx+1 << "/" << ST.las.size() << "][" << cn << "/" << xN << "]" << flush;
                }
                exmap smap;
                smap[a] = la+cn;
                if(WDigits>0) smap[a] = smap[a].evalf();
                for(int k=0; k<=kla; k++) {
                    matrix cmat(matN, matN);
                    for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                        for(int j=0; j<=kla; j++) {
                            matrix Tkj = ST.T[cm][k][j];
                            for(int i=0; i<Tkj.nops(); i++) Tkj.let_op(i) = Tkj.op(i).subs(smap,nopat);
                            cmat = cmat.add(Tkj.mul(CMats[(cn-cm)%s][j])); 
                        }
                    }
                    if(WDigits<0) matrix_map_inplace(cmat, [&](const ex & e) { return normal(e); });
                    CMats[cn%s][k] = cmat;
                    auto xterm = pow(x0,la)*pow(x0,cn)*pow(log(x0),k)/factorial(k);
                    if(WDigits>0) xterm = xterm.evalf(); // OK
                    MatF = MatF.add(cmat.mul_scalar(xterm));
                }

            }
            if(!CMat_Parallel && !In_GiNaC_Parallel && Verbose>1 && xN>0) cout << endl;
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Series @ " << now() << endl;
        if(WDigits>0) reset_precision();
        return MatF;
    }
    
    matrix DE::Taylor(const ex & x0, const ex & dx, const int xN) {        
        if(WDigits>0) set_precision(WDigits);
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
            for(int i=0; i<matN*matN; i++) B0.let_op(i) = BMat[0].op(i);
            B0 = B0.inverse();
                            
            for(int cm=1; cm<=TT.s; cm++) {
                matrix Bm(matN, matN);
                for(int i=0; i<matN*matN; i++) {
                    static symbol b("b");
                    Bm.let_op(i) = BMat[cm].op(i).subs(a==b-cm,nopat).subs(b==a,nopat);
                }
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
            
        for(int cn=1; cn<=xN; cn++) {
            if(!In_GiNaC_Parallel && Verbose>1) {
                cout << "\r                                  \r" << flush;
                cout << "     \\--C Matrix [" << cn << "/" << xN << "]" << flush;
            }
            exmap smap;
            smap[a] = cn;
            if(WDigits>0) smap[a] = smap[a].evalf();
            matrix cmat(matN, matN);
            for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                matrix Tm = TT.T[cm];
                for(int i=0; i<Tm.nops(); i++) Tm.let_op(i) = Tm.op(i).subs(smap,nopat);
                cmat = cmat.add(Tm.mul(CMats[(cn-cm)%s])); 
            }
            if(WDigits<0) matrix_map_inplace(cmat, [&](const ex & e) { return normal(e); });
            CMats[cn%s] = cmat;
            
            auto xterm = pow(dx,cn);
            if(WDigits>0) xterm = xterm.evalf();
            MatF = MatF.add(cmat.mul_scalar(xterm));
        }
        if(!In_GiNaC_Parallel && Verbose>1) cout << endl;
            
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Taylor @ " << now() << endl;
        if(WDigits>0) reset_precision();
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
    
}

