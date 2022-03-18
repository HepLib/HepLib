/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {

    using namespace EoD;
    
    DE::DE(const symbol & _x) : x(_x) { }
    DE::DE(const matrix & _mat, const symbol & _x) : Mat(_mat), x(_x) { }
    DE::DE(const symbol & _x, const matrix & _mat) : Mat(_mat), x(_x) { }
    DE::DE(const DE & b) : Mat(b.Mat), x(b.x), Ts(b.Ts) { }
    
    // Dx J = M.J ---> J = T.J' & Dx J' = M'.J' with M' = Ti.M.T - Ti.Dx T
    void DE::Apply(const matrix & t) {
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        Ts.push_back(t);
    }
    
    void DE::Apply(const lst & diag) {
        auto matN = diag.nops();
        matrix t(matN, matN);
        for(int i=0; i<matN; i++) t(i,i) = diag.op(i);
        Mat = t.inverse().mul(Mat.mul(t).sub(matrix_diff(t,x)));
        Ts.push_back(t);
    }
    
    void DE::x2y(const ex & y) {
        ex det = diff_ex(y,x);
        int matN = Mat.rows();
        for(int r=0; r<matN; r++) {
            for(int c=0; c<matN; c++) {
                Mat(r,c) = det * Mat(r,c).subs(x==y, nopat);
            }
        }
        auto n = Ts.size();
        for(int i=0; i<n; i++) {
            for(int r=0; r<matN; r++) {
                for(int c=0; c<matN; c++) {
                    Ts[i](r,c) = Ts[i](r,c).subs(x==y, nopat);
                }
            }
        }
    }
    
    void DE::xpow() {
        int n = Mat.nops();
        for(int i=0; i<n; i++) Mat.let_op(i) = xpow(Mat.op(i));
    }
    
    ex DE::xpow(const ex & e) {
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
    
    matrix DE::MatT() {
        int matN = Mat.rows();
        matrix t = ex_to<matrix>(unit_matrix(matN));
        for(auto ti : Ts) t = t.mul(ti);
        return t;
    }
    
    void DE::subs(const ex & sub, unsigned opt) {
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
    }
    
    void DE::subs(const exmap & sub, unsigned opt) {
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(sub, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(sub, opt);
        }
    }
    
    void DE::subs(const lst & l, const lst & r, unsigned opt) {
        auto mN = Mat.nops();
        for(int i=0; i<mN; i++) Mat.let_op(i) = Mat.op(i).subs(l, r, opt);
        for(auto & ti : Ts) {
            for(int i=0; i<mN; i++) ti.let_op(i) = ti.op(i).subs(l, r, opt);
        }
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
            if(Verbose>1) cout << "  \\--Fuchsify: Poincare rank is alreay " << pr << endl;
            Ts.push_back(t);
            return;
        }
        
        if(Verbose>1) cout << "  \\--Fuchsify: start @ " << now() << endl;
        
        while(true) {
            if(Verbose>1) cout << "     \\--Poincare rank: " << pr << ",  A0 rank: " << a0.rank() << " @ " << now() << endl;
            
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

        if(Verbose>1) cout << "  \\--Fuchsify: finished @ " << now() << endl;
        
        Mat = m;
        Ts.push_back(t);
    }
    
    // normalization at x=0
    void DE::Normalize() {
        if(Verbose>1) cout << "  \\--Normalize: start @ " << now() << endl;

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
            
            if(Verbose>0) cout << "     \\--EigenValues of A0 summary:" << endl;
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
            if(Verbose>1) {
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
                        if(Verbose>0) cout << "     \\--Raising ev = " << ev << " ..." << endl;
                        matrix vec_u = eigenvectors(a0, ev)[0];
                        matrix vec_v  = dual_basis(vec_u)[0];
                        matrix p = normal(vec_u.mul(vec_v));
                        m = with_balance_t(m, p, x);
                        t = t.mul(imatrix(m.rows()).sub(p).sub(p.mul_scalar(1/x)));
                    } else {
                        if(Verbose>0) cout << "     \\--Lowering ev = " << ev << " ..." << endl;
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
        
        if(Verbose>1) cout << "  \\--Normalize: fininshed @ " << now() << endl;

        Mat = m;
        Ts.push_back(t);
    }
    
    // shearing at x=0
    void DE::Shear() {
        if(Verbose>1) cout << "  \\--Shear: start @ " << now() << endl;

        matrix m = Mat;
        matrix t = imatrix(m.rows());
        matrix ti = imatrix(m.rows());
        if(prank(m,x)>0) throw Error("Shear: prank > 0.");
        
        matrix a0 = a0_matrix(m, x, 0);
        while(true) {
        
            if(Verbose>1) cout << "     \\--EigenValues of A0 summary:" << endl;
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
                if(Verbose>1) {
                    for(auto &ev_i : gi) {
                        ostr << ev_i << "[" << ev_map[ev_i] << "], ";
                    }
                }
                int tmp = ex_to<numeric>(gi[0]-gi[gi.size()-1]).to_int();
                if(max_dep < tmp) max_dep = tmp;
            }
            if(Verbose>1) {
                string str = ostr.str();
                cout << str.substr(0,str.size()-2) << endl;
            }
            
            if(is_normalized) break;
            
            if(Verbose>1) cout << "     \\--jordan decomposition ..." << endl;
            auto qj = jordan(a0);
            matrix smat = imatrix(m.rows());
            
            if(Verbose>1) cout << "     \\--shearing transformation ..." << endl;
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
        
        if(Verbose>1) cout << "     \\--DE Transformation  ..." << endl;
        m = transform(Mat, ti, t, x);
        matrix_map_inplace(m, [](const ex & e) { return normal(e); });
        
        if(Verbose>1) cout << "  \\--Shear: fininshed @ " << now() << endl;
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
        
    pair<matrix,matrix> DE::Series(const ex & x0, const int xN) {
        if(xN<0 && !is_a<numeric>(d0)) throw Error("DE::Series, xN<0 only works with numeric d0.");
    
        auto oDigits = Digits;
        if(Precision>0) Digits = Precision+ExDigits;
        if(Verbose>1) cout << "  \\--Series: start @ " << now() << endl;
        
        int matN = Mat.rows();
        
        if(prank(Mat,x)>0) Fuchsify();
        matrix a0 = normal(a0_matrix(Mat, x, 0));
        if(!is_jordan_form(a0)) {
            Shear();
            a0 = normal(a0_matrix(Mat, x, 0));
        }
        auto m = Mat;
        
        lst las;
        vector<pair<ex,int>> js;
        map<ex, unsigned, ex_is_less> laK;
        map<ex, vector<matrix>, ex_is_less> laC0;
        
        int lastN = 0;
        for(int r=0; r<a0.rows(); r++) {
            lastN++;
            if(r==a0.rows()-1 || a0(r,r+1)==0) { // the last row
                auto la = a0(r,r);
                las.append(la);
                js.push_back(make_pair(la, lastN));
                if(lastN-1 > laK[la]) laK[la] = lastN-1;
                lastN = 0;
            }
        }
        las.sort().unique();
        
        for(auto la : las) {
            int kla = laK[la];
            vector<matrix> vm;
            for(int li=0; li<=kla; li++) vm.push_back(matrix(matN, matN));
            laC0[la] = vm;
        }
        
        //init Cla0[la]
        int cur_pos = 0;
        for(int ji=0; ji<js.size(); ji++) {
            auto la = js[ji].first;
            int size = js[ji].second;
            for(int ir=0;ir<size;ir++) 
                for(int ic=0;ic<size;ic++)
                    if(ic-ir>=0) laC0[la][ic-ir](cur_pos+ir, cur_pos+ic) = 1; 
            cur_pos += size;
        }
        
        // save CoMat
        matrix CoMat(matN, matN);
        for(auto la : las) {
            int kla = laK[la];
            for(int k=0; k<=kla; k++) {
                auto mat = laC0[la][k];
                auto xterm = pow(x,la)*pow(log(x),k)/factorial(k);
                for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) CoMat(r,c) += mat(r,c)*xterm;
            }
        }
        
        if(Verbose>1) cout << "     \\--Q polynormial ..." << endl;
        
        matrix_map_inplace(m, [&](const ex & e){ return normal(e); });
        ex qlcm = matrix_den_lcm(m);
        if(normal(qlcm.subs(x==0,nopat)).is_zero()) qlcm = normal(qlcm/x); //make sure check
        matrix_map_inplace(m, [&](const ex & e){ return collect_ex(normal(qlcm*x*e),x); });
        
        qlcm = collect_ex(qlcm, x);
        int qmax = qlcm.degree(x);
        exvector qs;
        for(int i=0; i<=qmax; i++) {
            if(Precision>0 && i>0) qs.push_back(qlcm.coeff(x,i).evalf());
            else qs.push_back(qlcm.coeff(x,i));
        }
        
        if(Verbose>1) cout << "     \\--B Matrix ..." << endl;
        
        int s = qmax; // s in DESS
        for(int i=0; i<m.nops(); i++) {
            auto tmp = m.op(i).degree(x);
            if(s<tmp) s = tmp;
        }
        symbol a("a"); // alpha
        for(int i=0; i<matN; i++) m(i,i) -= a*qlcm;

        for(int i=qmax+1; i<=s; i++) qs.push_back(0);
        vector<matrix> BMat(s+1);
        for(int i=0; i<=s; i++) {
            if(Precision>0 && i>0)
                BMat[i] = matrix_map(m, [&](const ex & e){ return e.coeff(x, i).evalf(); });
            else 
                BMat[i] = matrix_map(m, [&](const ex & e){ return e.coeff(x, i); });
        }
        
        if(Verbose>1) cout << "     \\--C Matrix ..." << endl;
        
        auto res = GiNaC_Parallel(las.nops(), [&](int idx)->ex {
            auto la = las.op(idx);
            int kla = laK[la];
            vector<vector<matrix>> CMats(xN>=0 ? xN+1 : 1);
            CMats[0] = laC0[la];
            
            matrix MatF(matN, matN);
            for(int k=0; k<=kla; k++) {
                auto xterm = pow(x0,la)*pow(log(x0),k)/factorial(k);
                for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                    MatF(r,c) += CMats[0][k](r,c)*xterm;
                }
            }
            
            for(int cn=1; (xN<0 || cn<=xN); cn++) {
                if(xN<0) CMats.push_back(vector<matrix>());
                for(int i=0; i<=kla; i++) CMats[cn].push_back(matrix(matN, matN));
                matrix B0(matN, matN);
                for(int i=0; i<matN*matN; i++) B0.let_op(i) = BMat[0].op(i).subs(a==la+cn,nopat);
                BJFinv iBJF(B0, -qs[0], kla);
                for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                    matrix Bm(matN, matN);
                    for(int i=0; i<matN*matN; i++) Bm.let_op(i) = BMat[cm].op(i).subs(a==la+cn-cm,nopat);
                    BJF mBJF(Bm, -qs[cm], kla);
                    // T(n,m) . C(n-m, 0..kla) = -iBJF(B0,-q0).mBJF(Bm,-qm).CMats
                    for(int k=0; k<=kla; k++)
                    for(int i=0; i<=kla; i++)
                    for(int j=0; j<=kla; j++) {
                        if(k>i) continue; // iBJF(k,i)
                        if(i!=j && i+1!=j) continue; // mBJF(i,j)
                        CMats[cn][k] = CMats[cn][k].sub(iBJF(k,i).mul(mBJF(i,j)).mul(CMats[cn-cm][j])); 
                    }
                }
                
                if(xN>=0) {
                    for(int k=0; k<=kla; k++)  // normalize
                        matrix_map_inplace(CMats[cn][k], [&](const ex & e) { return normal(e); });
                }
                
                bool chk = true;
                for(int k=0; k<=kla; k++) {
                    auto xterm = pow(x0,la)*pow(x0,cn)*pow(log(x0),k)/factorial(k);
                    for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                        auto item = CMats[cn][k](r,c)*xterm;
                        MatF(r,c) += item;
                        if(xN<0 && chk) {
                            if(abs(MatF(r,c))>pow(10,Digits+3)) chk = abs(item/MatF(r,c)) < pow(10,-Precision);
                            else chk = abs(item) < pow(10,-Precision);
                        }
                    }
                }
                if(xN<0 && chk) {
                    if(Debug && Verbose>1) cout << "  \\--Taylor: Precision Reached after " << cn << " times." << endl;
                    break;
                }
            }
            
            return MatF;
        }, "CMat");
        matrix CMat(matN, matN);
        for(auto item : res) CMat = CMat.add(ex_to<matrix>(item));
        
        if(Verbose>1) cout << "  \\--Series: fininshed @ " << now() << endl;
        Digits = oDigits;
        return make_pair(CMat, CoMat);
    }
    
    matrix DE::Taylor(const ex & x0, const int xN) {
        if(xN<0 && !is_a<numeric>(d0)) throw Error("DE::Series, xN<0 only works with numeric d0.");
        
        auto oDigits = Digits;
        if(Precision>0) Digits = Precision+ExDigits;
        if(Verbose>1) cout << "  \\--Taylor: start @ " << now() << endl;
        
        auto m = Mat;
        int matN = m.rows();
        if(prank(m,x)>-1) throw Error("DE::Taylor only works for non-singualr point.");
        
        if(Verbose>1) cout << "     \\--Q polynormial ..." << endl;
        
        matrix_map_inplace(m, [&](const ex & e){ return normal(e); });
        ex qlcm = matrix_den_lcm(m);
        if(normal(qlcm.subs(x==0,nopat)).is_zero()) //make sure check
            throw Error("DE::Taylor only works for non-singualr point."); 
        matrix_map_inplace(m, [&](const ex & e){ return collect_ex(normal(qlcm*x*e),x); });
        
        qlcm = collect_ex(qlcm, x);
        int qmax = qlcm.degree(x);
        exvector qs;
        for(int i=0; i<=qmax; i++) {
            if(Precision>0) qs.push_back(qlcm.coeff(x,i).evalf());
            else qs.push_back(qlcm.coeff(x,i));
        }
        
        if(Verbose>1) cout << "     \\--B Matrix ..." << endl;
        
        int s = qmax; // s in DESS
        for(int i=0; i<m.nops(); i++) {
            auto tmp = m.op(i).degree(x);
            if(s<tmp) s = tmp;
        }
        symbol a("a"); // alpha
        for(int i=0; i<matN; i++) m(i,i) -= a*qlcm;

        for(int i=qmax+1; i<=s; i++) qs.push_back(0);
        vector<matrix> BMat(s+1);
        for(int i=0; i<=s; i++) {
            if(Precision>0)
                BMat[i] = matrix_map(m, [&](const ex & e){ return e.coeff(x, i).evalf(); });
            else 
                BMat[i] = matrix_map(m, [&](const ex & e){ return e.coeff(x, i); });
        }
        
        if(Verbose>1) cout << "     \\--C Matrix ..." << endl;
        
        vector<matrix> CMats(xN>=0 ? xN+1 : 1);
        CMats[0] = ex_to<matrix>(unit_matrix(matN));
        
        matrix CMat(matN, matN);
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
            CMat(r,c) += CMats[0](r,c);
        }
        for(int cn=1; (xN<0 || cn<=xN); cn++) {
            if(xN<0) CMats.push_back(matrix());
            CMats[cn] = matrix(matN, matN);
            matrix B0(matN, matN);
            for(int i=0; i<matN*matN; i++) B0.let_op(i) = BMat[0].op(i).subs(a==cn,nopat);
            auto iB0 = B0.inverse();
            for(int cm=1; (cm<=cn) && (cm<=s); cm++) {
                matrix Bm(matN, matN);
                for(int i=0; i<matN*matN; i++) Bm.let_op(i) = BMat[cm].op(i).subs(a==cn-cm,nopat);
                CMats[cn] = CMats[cn].sub(iB0.mul(Bm).mul(CMats[cn-cm])); 
            }
            
            if(xN>=0) matrix_map_inplace(CMats[cn], [&](const ex & e) { return normal(e); });
            
            bool chk = true;
            auto xterm = pow(x0,cn);
            for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) {
                ex item = CMats[cn](r,c)*xterm;
                CMat(r,c) += item;
                if(xN<0 && chk) {
                    if(abs(CMat(r,c))>pow(10,Digits+3)) chk = abs(item/CMat(r,c)) < pow(10,-Precision);
                    else chk = abs(item) < pow(10,-Precision);
                }
            }
            if(xN<0 && chk) {
                if(Verbose>1) cout << "  \\--Taylor: Precision Reached after " << cn << " times." << endl;
                break;
            }
        }
            
        if(Verbose>1) cout << "  \\--Taylor: fininshed @ " << now() << endl;
        Digits = oDigits;
        return CMat;
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
        if (is_a<mul>(den)) for (const auto &f : den) fvec.push_back(f);
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
                ex c0 = b.coeff(x, 0);
                ex c1 = b.coeff(x, 1);
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

    //==============================================================================

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
    
    
}

