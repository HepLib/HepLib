#include <fuchsia.cpp>

/*-----------------------------------------------------*/
/*               Functions by Feng                       */
/*-----------------------------------------------------*/

//----------------------------------------
// poincare rank at x = 0
//----------------------------------------
int prank(const matrix &mat, const symbol &x) {
    int p = 5;
    symbol y;
    bool first = true;
    while(true) {
        matrix m = matrix_map(mat, [&](auto &&e){ return normal(series_to_poly((e+y*pow(x,-p+1)+y*pow(x,-p+2)).series(x==0, -p))); });
        if(!m.is_zero_matrix()) {
            if(!first) return p;
            p++;
        } else {
            p--;
        }
        first = false;
    }
    return p;
}

//----------------------------------------
// matrix a0 and a1 at x = 0
//----------------------------------------
pair<matrix,matrix> a01_matrix(const matrix &mat, const symbol &x, int pr=19790923) {
    int p = (pr == 19790923 ? prank(mat, x) : pr);
    matrix a1 = matrix_map(mat, [&](auto &&e){ return collect(series_to_poly(e.series(x==0, -p+1)),x); });
    matrix a0 = matrix_map(a1, [&](auto &&e){ return e.coeff(x, -p-1); });
    matrix_map_inplace(a1, [&](auto &&e){ return e.coeff(x, -p); });
    assert(!a0.has(x) && !a1.has(x));
    return make_pair(a0, a1);
}

//----------------------------------------
// matrix a0 at x = 0
//----------------------------------------
matrix a0_matrix(const matrix &mat, const symbol &x, int pr=19790923) {
    int p = (pr == 19790923 ? prank(mat, x) : pr);
    matrix a0 = matrix_map(mat, [&](auto &&e){ return collect(series_to_poly(e.series(x==0, -p)),x); });
    matrix_map_inplace(a0, [&](auto &&e){ return e.coeff(x, -p-1); });
    assert(!a0.has(x));
    return a0;
}

//----------------------------------------
// balance transformation at x = 0 and x=oo
//----------------------------------------
matrix with_balance_t(const matrix &m, const matrix &P, const symbol &x) {
    matrix coP = imatrix(m.rows()).sub(P);
    matrix res = coP.sub(P.mul_scalar(x)).mul(m).mul(coP.sub(P.mul_scalar(1/x))).add(P.mul_scalar(1/x));
    return res;
}

//----------------------------------------
// balance transformation at x = oo and x=0, the inverse of 0 and oo
//----------------------------------------
matrix with_balance_ti(const matrix &m, const matrix &P, const symbol &x) {
    matrix coP = imatrix(m.rows()).sub(P);
    matrix res = coP.sub(P.mul_scalar(1/x)).mul(m).mul(coP.sub(P.mul_scalar(x))).sub(P.mul_scalar(1/x));
    return res;
}

//----------------------------------------
// dual basis matrix for u with x*u=1
//----------------------------------------
vector<matrix> dual_basis(const matrix &u) {
    matrix identity = imatrix(u.cols());
    vector<matrix> results;
    
    matrix tmp(u.cols(), u.rows());
    exmap tmpz;
    for (unsigned i = 0; i < tmp.nops(); i++) {
        symbol t;
        tmp.let_op(i) = t;
        tmpz[t] = 0;
    }
    
    try {
        // Solve x*u=1.
        matrix x = matrix_solve_left(normal(u), tmp, identity);
        matrix_map_inplace(x, [&](auto &&e) { return e.subs(tmpz); });
        if (PARANOID) {
            matrix z = normal(x.mul(u)).sub(imatrix(u.cols()));
            assert(z.is_zero_matrix());
        }
        results.push_back(x);
    } catch (const std::runtime_error &e) {
        if (e.what() != std::string("matrix::solve(): inconsistent linear system")) {
            throw;
        }
    }
    return results;
}

//----------------------------------------
// dual basis for left matrix u with u*x=1
//----------------------------------------
vector<matrix> dual_basis_left(const matrix &v) {
    vector<matrix> results;
    for(auto item : dual_basis(v.transpose())) {
        results.push_back(item.transpose());
    }
    return results;
}

//----------------------------------------
// date time string
//----------------------------------------
string now(bool use_date = true) {
    time_t timep;
    time (&timep);
    char tmp[64];
    if(use_date) strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    else strftime(tmp, sizeof(tmp), "%H:%M:%S",localtime(&timep) );
    return tmp;
}

//----------------------------------------
// matrix exponential
//----------------------------------------
matrix matrix_exp(const matrix &mat, const ex &lnx = 1) {
    matrix m(mat.rows(), mat.cols());
    auto qj = jordan(mat);
    int cur_pos = 0;
    for(auto kv : qj.second) {
        for(int ir=0;ir<kv.second;ir++) {
            for(int ic=0;ic<kv.second;ic++) {
                if(ir == ic) {
                    m(cur_pos+ir, cur_pos+ic) = exp(lnx*kv.first);
                } else if(ic - ir > 0) {
                    m(cur_pos+ir, cur_pos+ic) = exp(lnx*kv.first)*pow(lnx, ic-ir)/factorial(ic-ir);
                }
            }
        }
        cur_pos += kv.second;
    }
    m = qj.first.mul(m).mul(qj.first.inverse());
    return m;
}

//----------------------------------------
// collect function
//----------------------------------------
template <typename F>
ex collect_map(const ex &expr, const ex &vars, F f) {
    ex ret = 0;
    auto tmp = collect(expand(expr), vars);
    if(is_a<add>(tmp)) {
        for(auto const &ti : tmp) ret += f(ti);
    } else {
        ret = f(tmp);
    }
    return ret;
}

ex collect_normal(const ex &expr, const ex &vars) {
    ex ret = 0;
    auto tmp = collect(expand(expr), vars);
    if(is_a<add>(tmp)) {
        for(auto const &ti : tmp) ret += normal(ti);
    } else {
        ret = normal(tmp);
    }
    return ret;
}

//----------------------------------------
// fuchsify at x=0
//----------------------------------------
pair<matrix,matrix> fuchsify(const matrix &mat, const symbol &x) {
    matrix m = mat;
    auto a01 = a01_matrix(m, x);
    matrix a0 = a01.first;
    matrix a1 = a01.second;
    symbol lambda("L");
    matrix t= imatrix(m.rows());
    
    int pr = prank(m,x);
    if(pr < 1) {
        cout << "---------------------------------------" << endl;
        cout << "fichsify skipped: Poincare rank is " << pr << endl;
        return make_pair(m,t);
    }
    
    cout << "---------------------------------------" << endl;
    cout << "fuchsify start @ " << now() << endl;
    cout << "---------------------------------------" << endl;
    
    while(true) {
    
        cout << "  " << now(false) << " - Poincare rank: " << pr << "\tA0 rank: " << a0.rank() << endl;
        
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
                    s(deg, j) = factor(normal(e.coeff(lambda, deg)));
                }
            }
            ws.add_rows(s);
        }
        
        assert(!ws.basis_rows().has(lambda));
        
        ws.normalize();
        matrix wscols = ws.basis_cols();
        if (ws.dim() == 0) {
            throw fuchsia_error("matrix is Moser-irreducible");
        }
        if (PARANOID) {
            // ws is a subset of nullspace(a0)
            matrix z = normal(a0.mul(wscols));
            assert(z.is_zero_matrix());
        }
        
        auto dualbs = dual_basis(wscols);
        
        assert(dualbs.size() > 0);
        
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

    cout << "---------------------------------------" << endl;
    cout << "fuchsify finished @ " << now() << endl;
    
    return make_pair(m, t);
}

//----------------------------------------
// normalization at x=0
//----------------------------------------
pair<matrix,matrix> normalize(const matrix &mat, const symbol &x) {
    cout << "---------------------------------------" << endl;
    cout << "normalize start @ " << now() << endl;
    cout << "---------------------------------------" << endl;

    matrix m = mat;
    matrix t= imatrix(m.rows());
    assert(prank(m,x)==0);
    
    while(true) {
    
        matrix a0 = a0_matrix(m, x, 0);
        auto ev_map = eigenvalues(a0);
        vector<vector<ex>> ev_groups;
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
                vector<ex> gi_new;
                gi_new.push_back(ev);
                ev_groups.push_back(gi_new);
            }
        }
        
        cout << "  eigen values of A0 summary:" << endl;
        for(auto &gi : ev_groups) {
            sort(gi.begin(),gi.end(),[&](const auto &a, const auto &b){
                return normal(a-b).info(info_flags::positive);
            });
            ostringstream ostr;
            cout << "    ";
            for(auto &ev_i : gi) {
                ostr << ev_i << "[" << ev_map[ev_i] << "],  ";
            }
            string str = ostr.str();
            cout << str.substr(0,str.size()-3) << endl;
        }
        cout << "---------------------------------------" << endl;
        
        bool is_normalized = true;
        for(auto &gi : ev_groups) {
            if(gi.size()<2) continue;
            is_normalized = false;
            auto last_ev = gi[gi.size()-1];
            for_each(gi.begin(), gi.end()-1, [&](const auto &ev) {
                if((last_ev-ev).info(info_flags::positive)) {
                    cout << "  " << now(false) << " - raising ev = " << ev << " ..." << endl;
                    matrix vec_u = eigenvectors(a0, ev)[0];
                    matrix vec_v  = dual_basis(vec_u)[0];
                    matrix p = normal(vec_u.mul(vec_v));
                    m = with_balance_t(m, p, x);
                    t = t.mul(imatrix(m.rows()).sub(p).sub(p.mul_scalar(1/x)));
                } else {
                    cout << "  " << now(false) << " - lowering ev = " << ev << " ..." << endl;
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
        cout << "---------------------------------------" << endl;
    
    }
    
    cout << "normalize fininshed @ " << now() << endl;

    return make_pair(m, t);
}

//----------------------------------------
// shearing at x=0
//----------------------------------------
pair<matrix, matrix> shearing(const matrix &mat, const symbol &x, const symbol ep, const int epN = -19790923) {
    cout << "---------------------------------------" << endl;
    cout << "shearing start @ " << now() << endl;
    cout << "---------------------------------------" << endl;

    matrix m = mat;
    matrix t = imatrix(m.rows());
    matrix ti = imatrix(m.rows());
    assert(prank(m,x)==0);
    
    matrix a0 = normal(a0_matrix(m, x, 0));
    while(true) {
    
        cout << "  eigen values of A0 summary:" << endl;
        auto ev_map = eigenvalues(a0);
        vector<vector<ex>> ev_groups;
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
                vector<ex> gi_new;
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
        for(auto &gi : ev_groups) {
            sort(gi.begin(),gi.end(),[&](const auto &a, const auto &b){
                return normal(a-b).info(info_flags::positive);
            });
            ostringstream ostr;
            ostr << "    ";
            for(auto &ev_i : gi) {
                ostr << ev_i << "[" << ev_map[ev_i] << "],  ";
            }
            string str = ostr.str();
            cout << str.substr(0,str.size()-3) << endl;
            int tmp = ex_to<numeric>(gi[0]-gi[gi.size()-1]).to_int();
            if(max_dep < tmp) max_dep = tmp;
        }
        cout << "---------------------------------------" << endl;
        
        if(is_normalized) break;
        
        cout << "  " << now(false) << " - jordan decomposition ..." << endl;
        auto qj = jordan(a0);
        matrix smat = imatrix(m.rows());
        
        cout << "  " << now(false) << " - shearing ..." << endl;
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
        
        auto tm = normal(qj.first.mul(smat));
        auto tmi = normal(tm.inverse());
        m = transform(m, tmi, tm, x);
        t = normal(t.mul(tm));
        ti = normal(tmi.mul(ti));
        
        int current = 0;
        int total = m.rows() * m.cols();
        matrix_map_inplace(m, [&](auto &&e){
            cout << "\r  [ series ] : " << (++current)*100/total << "%" << flush;
            return series_to_poly(normal(e).series(x==0, max_dep-1));
        });
        
        a0 = matrix_map(m, [&](auto &&e){
            return normal(e.collect(x).coeff(x,-1));
        });
        assert(!a0.has(x));
                
        current = 0;
        cout << "\r                       \r" << flush;
        matrix_map_inplace(m, [&](auto &&e){
            cout << "\r  [ normal ] : " << (++current)*100/total << "%" << flush;
            return normal(e);
        });
        cout << "\r" << flush;
        
        cout << "---------------------------------------" << endl;
        
    }
    
    cout << "  " << now(false) << " - normal t/ti  ..." << endl;
    auto qj = jordan(a0);
    matrix qi = normal(qj.first.inverse());
    a0 = normal(qi.mul(a0).mul(qj.first));
    t = t.mul(qj.first);
    ti = qi.mul(ti);
    int current = 0;
    int total = t.rows() * t.rows();
    matrix_map_inplace(t, [&](auto &&e){
        cout << "\r  [ normal 1/2 ] : " << (++current)*100/total << "%" << flush;
        return normal(e);
    });
    
    current = 0;
    cout << "\r                       \r" << flush;
    matrix_map_inplace(ti, [&](auto &&e){
        cout << "\r  [ normal 2/2 ] : " << (++current)*100/total << "%" << flush;
        return normal(e);
    });
    
    cout << "\r" << flush;
    cout << "  " << now(false) << " - final transformation  ..." << endl;
    m = transform(mat, ti, t, x);
    current = 0;
    matrix_map_inplace(m, [&](auto &&e){
        cout << "\r  [ normal ] : " << (++current)*100/total << "%" << flush;
        return epN == -19790923 ? normal(e) : collect_normal(series_to_poly(e.expand().series(ep==0, epN)), lst{ep});
    });
    cout << "\r" << flush;
    
    
    cout << "---------------------------------------" << endl;
    cout << "shearing fininshed @ " << now() << endl;

    
    return make_pair(m,t);
}

//----------------------------------------
// dess at x = 0
//----------------------------------------
bool is_jordan_form(matrix mat) {
    for(int r=0; r<mat.rows(); r++) {
        for(int c=0; c<mat.cols(); c++) {
            if(r==c) continue;
            if(r>c) {
                if(mat(r,c)!=0) return false;
            } else if(r<c-1) {
                if(mat(r,c)!=0) return false;
            } else {
                if(mat(r,c)==1) {
                    if(mat(r,c-1)!=mat(r+1,c)) return false;
                } else if(mat(r,c)!=0)
                    return false;
            }
        }
    }
    return true;
}

matrix BJF(matrix const &a, ex const &b, int Klambda) {
    int sa = a.rows();
    int sm = (Klambda+1) * sa;
    matrix ret(sm, sm);
    for(int k=0; k<=Klambda; k++) {
        for(int r=0; r<sa; r++) {
            if(k<Klambda) ret(k*sa+r, (k+1)*sa+r) = b;
            for(int c=0; c<sa; c++)
                ret(k*sa+r, k*sa+c) = a(r,c);
        }
    }
    return ret;
}

matrix dess(const matrix mat, const symbol &x, const symbol &ep,
    const ex x0, const int xN = 10, const int epN = -19790923
) {
    cout << "---------------------------------------" << endl;
    cout << "dess start @ " << now() << endl;
    cout << "---------------------------------------" << endl;
    
    int matN = mat.rows();
    bool ep_series = !(epN == -19790923);
    
    int current = 0;
    int total = matN * matN;
    matrix mm = mat;
    matrix a0 = normal(a0_matrix(mm, x, 0));
    symbol a("a");
    
    assert(is_jordan_form(a0));
    
    lst lambda_list;
    vector<pair<ex,int>> js;
    map<ex, unsigned, ex_is_less> dessK;
    map<ex, vector<matrix>, ex_is_less> CMatrix;
    
    int lastN = 0;
    for(int r=0; r<a0.rows()-1; r++) {
        lastN++;
        if(a0(r,r+1)==0) {
            auto la = a0(r,r);
            lambda_list.append(la);
            js.push_back(make_pair(la, lastN));
            if(lastN-1 > dessK[la]) dessK[la] = lastN-1;
            lastN = 0;
        }
    }
    lambda_list.sort().unique();
    
    for(auto lambda : lambda_list) {
        matrix cmat(matN*(dessK[lambda]+1), matN);
        vector<matrix> tmp;
        tmp.push_back(cmat);
        CMatrix[lambda] = tmp;
    }
    
    //init CMatrix[la][0]
    int cur_pos = 0;
    for(int ji=0; ji<js.size(); ji++) {
        auto la = js[ji].first;
        int size = js[ji].second;
        for(int ir=0;ir<size;ir++) {
            for(int ic=0;ic<size;ic++) {
                if(ir == ic) {
                    CMatrix[la][0](cur_pos+ir, cur_pos+ic) = 1;
                } else if(ic - ir > 0) {
                    CMatrix[la][0]( (ic-ir)*matN + (cur_pos+ir), cur_pos+ic) = 1/factorial(ic-ir);
                }
            }
        }
        cur_pos += size;
    }
    
    if(ep_series) {
        matrix_map_inplace(mm, [&](auto &&e){
            cout << "\r  [ ep series ] : " << (++current)*100/total << "%" << flush;
            return collect_normal(series_to_poly(e.expand().series(ep==0, epN)), lst{ep});
        });
    }
    
    cout << "\r" << flush;
    cout << "  " << now(false) << " - Q-polynormial ..." << endl;
    ex qlcm = 1;
    current = 0;
    matrix_map(mm, [&](auto &&e) {
        cout << "\r  [ lcm ] : " << (++current)*100/total << "%" << flush;
        return qlcm = lcm(qlcm, e.denom());
    });
    
    if(normal(qlcm.subs(x==0)).is_zero()) qlcm = normal(qlcm/x); //make sure check
    qlcm = qlcm.expand().collect(x);
    int qmax = qlcm.degree(x);
    vector<ex> dessq;
    for(int ii=0; ii<=qmax; ii++) {
        dessq.push_back(qlcm.coeff(x, ii));
    }
    
    cout << "\r" << flush;
    cout << "  " << now(false) << " - B-Matrix ..." << endl;
    current = 0;
    matrix_map_inplace(mm, [&](auto &&e){
        cout << "\r  [ normal ] : " << (++current)*100/total << "%" << flush;
        return normal(qlcm * x * e);
    });
    
    for(int ii=0; ii<matN; ii++) mm(ii,ii) -= a*qlcm;
    
    int bmax = 0;
    matrix_map_inplace(mm, [&](auto &&e){
        auto tmp = e.expand().collect(x);
        if(bmax < tmp.degree(x)) bmax = tmp.degree(x);
        return tmp;
    });
    
    for(int ii=qmax+1; ii<=bmax; ii++) dessq.push_back(0);
    
    vector<matrix> BMatrix;
    for(int ii=0; ii<=bmax; ii++) {
        BMatrix.push_back(matrix_map(mm, [&](auto &&e){
            return e.coeff(x, ii);
        }));
    }
    
    cout << "\r" << flush;
    cout << "  " << now(false) << " - C-Matrix ..." << endl;
    
    matrix dessMat(matN, matN);
    for(auto lambda : lambda_list) {
        matrix lamMat(matN, matN);
        for(int cn=0; cn<=xN; cn++) {
            cout << "\r                       \r" << flush;
            cout << "\r  [ \u03BB = " << lambda << " ] : " << cn << "/" << xN << flush;
            matrix csum(matN*(dessK[lambda]+1), matN);
            for(int cm=1; (cm<=cn) && (cm<=bmax); cm++) {
                matrix Tnm = BJF(
                    matrix_map(BMatrix[0],[&](auto &&e){ return e.subs(a==lambda+cn); }),
                    -dessq[0], dessK[lambda]
                ).inverse().mul_scalar(-1).mul(BJF(
                    matrix_map(BMatrix[cm],[&](auto &&e){ return e.subs(a==lambda+cn-cm); }),
                    -dessq[cm], dessK[lambda]
                ));
                auto ctmp = Tnm.mul(CMatrix[lambda][cn-cm]);
                
                matrix_map_inplace(ctmp, [&](auto &&e) {
                    return ep_series ? collect_normal(series_to_poly(e.expand().series(ep==0, epN)), lst{ep}) : normal(e);
                });
                csum = csum.add(ctmp);
            }

            matrix_map_inplace(csum, [&](auto &&e) {
                return ep_series ? collect_normal(e, lst{ep}) : normal(e);
            });
            
            if(cn>0) CMatrix[lambda].push_back(csum);
            else csum = CMatrix[lambda][0];
            
            matrix kMat(matN, matN);
            for(int kk=0; kk<=dessK[lambda]; kk++) {
                kMat = kMat.add(ex_to<matrix>(sub_matrix(csum, kk*matN,matN, 0,matN)).mul_scalar(pow(log(x0),kk)/factorial(kk)));
            }
            lamMat = lamMat.add(kMat.mul_scalar(pow(x0, cn)));
        }
        dessMat = dessMat.add(lamMat.mul_scalar(pow(x0,lambda)));
    }
    
    cout << "\r" << flush;
    cout << "  " << now(false) << " - finalizing ..." << endl;

    current = 0;
    matrix_map_inplace(dessMat, [&](auto &&e) {
        cout << "\r  [ Evaluated ] : " << (++current)*100/total << "%" << flush;
        return ep_series ? collect_normal(series_to_poly(e.expand().series(ep==0, epN)), lst{ep}) : normal(e);
    });
    
    cout << "\r" << flush;
    cout << "---------------------------------------" << endl;
    cout << "dess fininshed @ " << now() << endl;
    return dessMat;
}
