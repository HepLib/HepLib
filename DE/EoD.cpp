/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib::EoD {

    // DE transformation
    matrix matrix_diff(const matrix &mat, const symbol &x, const int n) {
        return matrix_map(mat, [&](auto &&e){ return e.diff(x,n); });
    }

    // Dx J = M.J ---> J = T.J' & Dx J' = M'.J' with M' = Ti.M.T - Ti.Dx T
    matrix transform(const matrix &m, const matrix &t, const symbol &x) {
        return t.inverse().mul(m.mul(t).sub(matrix_diff(t,x)));
    }
    
    matrix transform(const matrix &m, const matrix &tinverse, const matrix &t, const symbol &x) {
        return tinverse.mul(m.mul(t).sub(matrix_diff(t,x)));
    }

    matrix normal(const matrix &m) {
        return matrix_map(m, [&](auto &&e) { return normal(e); });
    }

    matrix matrix_cut(const matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc) {
        matrix res(nr, nc);
        for (unsigned i = 0; i < nr; i++) {
            for (unsigned j = 0; j < nc; j++) {
                res(i, j) = m(r + i, c + j);
            }
        }
        return res;
    }
    
    void rescale_submatrix(matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc) {
        exvector n(nr*nc), d(nr*nc);
        ex mul = 0;
        ex div = 0;
        auto it_n = n.begin();
        auto it_d = d.begin();
        for (unsigned i = r; i < r + nr; i++) {
            for (unsigned j = c; j < c + nc; j++) {
                ex nd = m(i, j).normal().numer_denom();
                ex numer = nd.op(0);
                ex denom = nd.op(1);
                *it_n++ = numer;
                *it_d++ = denom;
                div = div.is_zero() ? numer : gcd(div, numer);
                mul = mul.is_zero() ? denom : lcm(mul, denom);
            }
        }
        if (div.is_zero()) return;
        if (mul.is_zero()) return;
        // It would be tempting to exit here if div=mul=1, but we'd
        // then discard the normal() call results.
        it_n = n.begin();
        it_d = d.begin();
        for (unsigned i = r; i < r + nr; i++) {
            for (unsigned j = c; j < c + nc; j++) {
                ex nn, dd;
                bool ok1 = divide(*it_n++, div, nn);
                bool ok2 = divide(mul, *it_d++, dd);
                if(!ok1 || !ok2) throw Error("Error");
                m(i, j) = nn*dd;
            }
        }
    }
    
    void echelon_form_gauss(matrix &m) {
        unsigned nr = m.rows();
        unsigned nc = m.cols();
        exvector &mv = ((matrix_hack*)&m)->mvec();
        for (unsigned i = 0; i < nr*nc; i++) {
            mv[i] = mv[i].normal();
        }
        unsigned r0 = 0;
        unsigned c0 = 0;
        for (; (c0 < nc) && (r0 < nr - 1); c0++) {
            for (unsigned r = r0; r < nr; r++) {
                mv[r*nc + c0] = mv[r*nc + c0].normal();
            }
            // No normalization before is_zero() here, because
            // we maintain the matrix normalized throughout the
            // algorithm.
            unsigned pivot = r0;
            while ((pivot < nr) && mv[pivot*nc + c0].is_zero()) {
                pivot++;
            }
            if (pivot == nr) {
                // The whole column below r0:c0 is zero, let's skip
                // it.
                continue;
            }
            if (pivot > r0) {
                // Found a non-zero row somewhere below r0; let's
                // swap it in.
                for (unsigned c = c0; c < nc; c++) {
                    mv[pivot*nc + c].swap(mv[r0*nc + c]);
                }
            }
            ex a = mv[r0*nc + c0];
            for (unsigned r = r0 + 1; r < nr; r++) {
                ex b = mv[r*nc + c0];
                if (!b.is_zero()) {
                    ex k = b/a;
                    mv[r*nc + c0] = 0;
                    for (unsigned c = c0 + 1; c < nc; c++) {
                        mv[r*nc + c] = normal(mv[r*nc + c] - k*mv[r0*nc + c]);
                    }
                }
            }
            r0++;
        }
        // Zero out the remaining rows (just in case).
        for (unsigned r = r0 + 1; r < nr; r++) {
            for (unsigned c = 0; c < nc; c++) {
                mv[r*nc + c] = 0;
            }
        }
    }

    void matrix_hack::append_rows(const matrix &src) {
        if(col != src.cols()) throw Error("error");
        row += src.rows();
        m.insert(m.end(), ((const matrix_hack*)&src)->m.begin(), ((const matrix_hack*)&src)->m.end());
    }

    exvector & matrix_hack::mvec() {
        ensure_if_modifiable();
        return m;
    }

    void matrix_hack::resize(unsigned nrows) {
        row = nrows;
        m.resize(row*col);
    }

    namespace {
        ex ratcan(const ex &e) {
            ex nd = e.normal().numer_denom();
            return nd.op(0).expand()/nd.op(1).expand();
        }
    }

    matrix imatrix(unsigned n) {
        return ex_to<matrix>(unit_matrix(n));
    }
    
    map<ex, unsigned, ex_is_less> eigenvalues(const matrix &m) {
        map<ex, unsigned, ex_is_less> eigenvalues;
        symbol lambda("L");
        ex charpoly = fermat_Det(imatrix(m.rows()).mul_scalar(lambda).sub(m));
        ex num = exfactor(normal(charpoly).numer());
        exvector fvec;
        if (is_a<mul>(num)) for (const auto &f : num) fvec.push_back(f);
        else fvec.push_back(num);
        for (const auto &f : fvec) {
            ex b = f; 
            int n = 1;
            if (is_a<power>(f)) {
                b = f.op(0);
                n = ex_to<numeric>(f.op(1)).to_int();
            } 
            int deg = b.degree(lambda);
            if (deg == 0) { }
            else if (deg == 1) {
                ex c0 = b.coeff(lambda, 0);
                ex c1 = b.coeff(lambda, 1);
                eigenvalues[ratcan(-c0/c1)] += n;
            } else {
                cout << charpoly << endl;
                throw Error("eigenvalues: higher powers found.");
            }
        }
        return eigenvalues;
    }
    
    int prank(const matrix &mat, const symbol &x) {
        int p = 5;
        while(true) {
            int pr = -100000;
            for (int i=0; i<mat.nops(); i++) {
                auto ir = series_ex(mat.op(i), x, -p);
                if(!is_zero(ir)) {
                    int cpr = -ir.ldegree(x)-1;
                    if(pr<cpr) pr = cpr;
                }
            }
            if(pr!=-100000) return pr;
            p--;
        }
        throw Error("prank");
    }
    
    pair<matrix,matrix> a01_matrix(const matrix &mat, const symbol &x, int pr) {
        int p = (pr == 19790923 ? prank(mat, x) : pr);
        int nr=mat.rows(), nc=mat.cols();
        matrix a0(nr,nc), a1(nr,nc);
        auto res = GiNaC_Parallel(mat.nops(), [&mat,&p,&x](int idx) {
            auto item = mat.op(idx);
            item = series_ex(item, x, -p);
            return lst{ item.coeff(x, -p-1), item.coeff(x, -p) };
        }, "a01");
        
        for(int i=0; i<mat.nops(); i++) {
            a0.let_op(i) = res[i].op(0);
            a1.let_op(i) = res[i].op(1);
        }
        if(a0.has(x) || a1.has(x)) throw Error("a0/a1 staill has x.");
        return make_pair(a0, a1);
    }


    matrix a0_matrix(const matrix &mat, const symbol &x, int pr) {
        int p = (pr == 19790923 ? prank(mat, x) : pr);
        int nr=mat.rows(), nc=mat.cols();
        GiNaC_Parallel_Verb["a0"] = 0;
        auto res = GiNaC_Parallel(mat.nops(), [&mat,&p,&x](int idx) {
            auto item = mat.op(idx);
            item = series_ex(item, x, -p-1);
            return item.coeff(x,-p-1);
        }, "a0");
        matrix a0(nr,nc);
        for(int i=0; i<mat.nops(); i++) a0.let_op(i) = res[i];
        if(a0.has(x)) throw Error("a0 still has x.");
        return a0;
    }
    
    // balance transformation at x = 0 and x=oo
    matrix with_balance_t(const matrix &m, const matrix &P, const symbol &x) {
        matrix coP = imatrix(m.rows()).sub(P);
        matrix res = coP.sub(P.mul_scalar(x)).mul(m).mul(coP.sub(P.mul_scalar(1/x))).add(P.mul_scalar(1/x));
        return res;
    }
    
    // balance transformation at x = oo and x=0, the inverse of 0 and oo
    matrix with_balance_ti(const matrix &m, const matrix &P, const symbol &x) {
        matrix coP = imatrix(m.rows()).sub(P);
        matrix res = coP.sub(P.mul_scalar(1/x)).mul(m).mul(coP.sub(P.mul_scalar(x))).sub(P.mul_scalar(1/x));
        return res;
    }
    
    // dual basis matrix for u with x*u=1
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
            results.push_back(x);
        } catch (const std::runtime_error &e) {
            if (e.what() != std::string("matrix::solve(): inconsistent linear system")) {
                throw;
            }
        }
        return results;
    }
    
    vector<matrix> dual_basis_left(const matrix &v) {
        vector<matrix> results;
        for(auto item : dual_basis(v.transpose())) {
            results.push_back(item.transpose());
        }
        return results;
    }
    
    // matrix exponential
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
    
    bool is_jordan_form(const matrix & mat) {
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
    
    matrix matrix_solve_left(const matrix &m, const matrix &vars, const matrix &rhs) {
        return m.transpose().solve(vars.transpose(), rhs.transpose()).transpose();
    }

    void matrix_map_inplace(matrix &m, std::function<ex(const ex &)> f) {
        GiNaC_Parallel_Verb["MMI"] = 0;
        auto res = GiNaC_Parallel(m.nops(), [&m, &f](int idx) {
            return f(m.op(idx));
        }, "MMI");
        for(unsigned i=0; i<m.nops(); i++) m.let_op(i) = res[i];
    }

    matrix matrix_map(const matrix &m, std::function<ex(const ex &)> f) {
        matrix r(m.rows(), m.cols());
        for(unsigned i=0; i<m.nops(); i++) r.let_op(i) =f(m.op(i));
        return r;
    }
        
}

