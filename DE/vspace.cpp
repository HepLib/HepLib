
#include "vspace.h"

// https://github.com/magv/fuchsia.cpp

namespace HepLib::ED {

    namespace {
        ex ratcan(const ex &e) {
            ex nd = HepLib::exnormal(e).numer_denom();
            return nd.op(0).expand()/nd.op(1).expand();
        }
    }

    vspace::vspace(unsigned n) : basis(0, n) { }

    vspace::vspace(const matrix &b) : basis(b) {
        for (unsigned i = 0; i < basis.rows(); i++) {
            rescale_submatrix(basis, i, 1, 0, basis.cols());
        }
        normalize();
    }

    void vspace::add_rows(const matrix &v) {
        ((matrix_hack*)&basis)->append_rows(v);
    }

    void vspace::normalize() {
        echelon_form_gauss(basis);
        unsigned nrows = basis.rows();
        for (; nrows > 0; nrows--) {
            for (unsigned c = 0; c < basis.cols(); c++) {
                if (!HepLib::exnormal(basis(nrows - 1, c)).is_zero()) goto done;
            }
        }
    done:;
        ((matrix_hack*)&basis)->resize(nrows);
        for (unsigned i = 0; i < basis.rows(); i++) {
            rescale_submatrix(basis, i, 1, 0, basis.cols());
        }
    }

    unsigned vspace::dim() const {
        return basis.rows();
    }

    matrix vspace::basis_col(unsigned i) const {
        if(i >= basis.rows()) throw Error("Error");
        matrix v(basis.cols(), 1);
        for (unsigned j = 0; j < basis.cols(); j++) {
            v[j] = basis(i, j);
        }
        return v;
    }

    matrix vspace::basis_row(unsigned i) const {
        if(i >= basis.rows()) throw Error("Error");
        matrix v(1, basis.cols());
        for (unsigned j = 0; j < basis.cols(); j++) {
            v[j] = basis(i, j);
        }
        return v;
    }

    const matrix & vspace::basis_rows() const {
        return basis;
    }

    const matrix vspace::basis_cols() const {
        return basis.transpose();
    }

    bool vspace::contains(const matrix &v) const {
        if(v.nops() != basis.cols()) throw Error("Error");
        matrix vv = v;
        rescale_submatrix(vv, 0, v.rows(), 0, v.cols());
        unsigned p = 0;
        // Division-free subtraction of basis vectors from v.
        for (unsigned i = 0; i < basis.rows(); i++, p++) {
            // Advance p to the first non-zero column of basis[i].
            for (;;) {
                // This assertion should only fail if the normalize()
                // was not called between add_rows() and contains().
                if(p >= basis.cols()) throw Error("Error");
                if (!basis(i, p).is_zero()) break;
                vv[p] = HepLib::exnormal(vv.op(p));
                // If vv has non-zero columns before p, it's not in
                // the basis.
                if (!vv.op(p).is_zero())
                    return false;
                p++;
            }
            // Subtract basis[i] from vv, if vv[p] != 0.
            const ex &vv_p = vv.op(p);
            if (!vv_p.is_zero()) {
                const ex b_ip = basis(i, p);
                vv[p] = 0;
                for (unsigned j = p + 1; j < basis.cols(); j++) {
                    vv[j] = HepLib::exnormal(vv.op(j)*b_ip - basis(i, j)*vv_p);
                }
            }
        }
        for (unsigned i = p; i < basis.cols(); i++) {
            vv[i] = HepLib::exnormal(vv.op(i));
            if (!vv.op(i).is_zero())
                return false;
        }
        return true;
    }
    
    vspace nullspace(const matrix &m) {
        unsigned nrows = m.rows();
        unsigned ncols = m.cols();
        matrix v(ncols, 1);
        std::vector<symbol> tmp;
        for (unsigned i = 0; i < ncols; i++) {
            symbol t;
            v(i, 0) = t;
            tmp.push_back(t);
        }
        // Solve M*V = 0
        matrix zero(nrows, 1);
        matrix s = exnormal(m).solve(v, zero);
        matrix coeff(ncols, ncols);
        for (unsigned k = 0; k < ncols; k++) {
            for (unsigned i = 0; i < ncols; i++) {
                coeff(k, i) = s(i, 0).coeff(tmp[k]);
            }
        }
        return vspace(coeff);
    }

    vspace eigenspace(const matrix &m, const ex &eval) {
        matrix mm = m;
        for (unsigned i = 0; i < m.rows(); i++) {
            mm(i, i) -= eval;
        }
        return nullspace(mm);
    }

    vspace eigenspace_left(const matrix &m, const ex &eval) {
        return eigenspace(m.transpose(), eval);
    }

    vector<matrix> eigenvectors(const matrix &m, const ex &eval) {
        vspace es = eigenspace(m, eval);
        vector<matrix> evectors(es.dim());
        for (unsigned i = 0; i < es.dim(); i++) {
            evectors[i] = es.basis_col(i);
        }
        return evectors;
    }

    vector<matrix> eigenvectors_left(const matrix &m, const ex &eval) {
        vspace es = eigenspace_left(m, eval);
        vector<matrix> evectors(es.dim());
        for (unsigned i = 0; i < es.dim(); i++) {
            evectors[i] = es.basis_row(i);
        }
        return evectors;
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
        // added by Feng
        for (unsigned i = r; i < r + nr; i++) {
            for (unsigned j = c; j < c + nc; j++) {
                if(!m(i, j).info(info_flags::numeric)) goto done;
            }
        }
        return; 
        done: ;
        // added END
        
        exvector n(nr*nc), d(nr*nc);
        ex mul = 0;
        ex div = 0;
        auto it_n = n.begin();
        auto it_d = d.begin();
        for (unsigned i = r; i < r + nr; i++) {
            for (unsigned j = c; j < c + nc; j++) {
                ex nd = HepLib::exnormal(m(i, j)).numer_denom();
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
            mv[i] = HepLib::exnormal(mv[i]);
        }
        unsigned r0 = 0;
        unsigned c0 = 0;
        for (; (c0 < nc) && (r0 < nr - 1); c0++) {
            for (unsigned r = r0; r < nr; r++) {
                mv[r*nc + c0] = HepLib::exnormal(mv[r*nc + c0]);
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
                        mv[r*nc + c] = HepLib::exnormal(mv[r*nc + c] - k*mv[r0*nc + c]);
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
    
    // dual basis matrix for u with x*u=1
    vector<matrix> dual_basis(const matrix &u) {
        matrix identity = ex_to<matrix>(unit_matrix(u.cols()));
        vector<matrix> results;
        
        matrix tmp(u.cols(), u.rows());
        exmap tmpz;
        for (unsigned i = 0; i < tmp.nops(); i++) {
            symbol t;
            tmp[i] = t;
            tmpz[t] = 0;
        }
        
        try {
            // Solve x*u=1.
            matrix x = matrix_solve_left(exnormal(u), tmp, identity);
            for(int i=0; i<x.nops(); i++) x[i] = x.op(i).subs(tmpz);
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
    
    matrix matrix_solve_left(const matrix &m, const matrix &vars, const matrix &rhs) {
        return m.transpose().solve(vars.transpose(), rhs.transpose()).transpose();
    }

    matrix exnormal(const matrix & m0) {
        matrix m = m0;
        for(unsigned i=0; i<m.nops(); i++) m[i] = HepLib::exnormal(m.op(i));
        return m;
    }    
    
    //=*********************************************************************=
    
    ex fmpq_mat_charpoly(const matrix & mat, const symbol & lambda) {
        auto nr = mat.rows();
        auto nc = mat.cols();
        fmpq_mat_t qmat;
        fmpq_mat_init(qmat, nr, nc);
        _to_(qmat, mat);
        fmpq_poly_t pol;
        fmpq_poly_init(pol);
        fmpq_mat_charpoly(pol, qmat);
        fmpq_mat_clear(qmat);
        ex charpoly = _to_(lambda, pol);
        fmpq_poly_clear(pol);
        return charpoly;
    }
    
    ev_am_t ev_am(const matrix &mat) { // eigenvalue_algebric_multiplicity
        ev_am_t ev2am;
        symbol lambda("x");
        matrix m = ex_to<matrix>(unit_matrix(mat.rows()));
        m = m.mul_scalar(lambda).sub(mat); // lambda I - M
        //ex charpoly = fermat_Det(m);
        ex charpoly = fmpq_mat_charpoly(mat, lambda);
        ex num = factor_flint(HepLib::exnormal(charpoly).numer());
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
                ev2am[ratcan(-c0/c1)] += n;
            } else {
                cout << endl << "charpoly: " << charpoly << endl;
                cout << "item: " << f << endl;
                throw Error("ev_am: higher powers found.");
            }
        }
        return ev2am;
    }

    pair<matrix, jcf_t> jordan(const matrix &m) {
        if(m.cols() != m.rows()) throw Error("jordan");
        auto ev2am = ev_am(m);
        return jordan(m, ev2am);
    }
        
    pair<matrix, jcf_t> jordan(const matrix &m, ev_am_t & ev2am) {
        unsigned n = m.rows();
        matrix q(n, n);
        vector<pair<ex,int>> jcs;
        int idxq = 0;
        for (const auto &kv : ev2am) {
            const auto &eval = kv.first;
            const auto &almul = kv.second;
            matrix mm = m;
            for (unsigned i = 0; i < n; i++) {
                mm(i, i) -= eval;
            }
            vector<vspace> nspace;
            vspace xspace(n);
            matrix mmk = mm;
            for (;;) {
                auto ns = nullspace(mmk);
                nspace.push_back(ns);
                if (ns.dim() == almul) break;
                mmk = mmk.mul(mm);
            }
            int kk = 0;
            for (int i = nspace.size() - 1; i >= 0; i--) {
                /* We're expecting kt - kk chains to start at rank i.
                 */
                int kt = (i == 0) ?
                    nspace[i].dim() : nspace[i].dim() - nspace[i-1].dim();
                for (int k = nspace[i].dim() - 1; kk < kt; k--) {
                    /* Find a vector in nspace[i] that is not in
                     * nspace[i-1] (and thus, a chain head), and not
                     * in xspace (not excluded).
                     */
                    if(k < 0) throw Error("jordan");
                    auto v = nspace[i].basis_col(k);
                    if ((i != 0) && nspace[i-1].contains(v)) continue;
                    if (xspace.contains(v)) continue;
                    /* Exclude v from the future results (so no
                     * other vector is linearly dependent on it),
                     * and build a Jordan chain starting from v.
                     */
                    for (unsigned r = 0; r < n; r++) {
                        q(r, idxq + i) = v.op(r);
                    }
                    for (int j = i - 1; j >= 0; j--) {
                        v = mm.mul(v);
                        for (unsigned r = 0; r < n; r++) {
                            q(r, idxq + j) = v.op(r);
                        }
                    }
                    if (xspace.contains(v)) {
                        // The last element of a chain is in xspace.
                        // We'll skip this chain.
                        continue;
                    }
                    /* Only rescale chains of length > 1; single
                     * eigenvectors are already well scaled by
                     * vspace().
                     */
                    if (i != 0) {
                        rescale_submatrix(q, 0, n, idxq, i+1);
                    }
                    xspace.add_rows(matrix_cut(q, 0, n, idxq, i+1).transpose());
                    xspace.normalize();
                    idxq += i + 1;
                    jcs.push_back(make_pair(eval, i + 1));
                    kk++;
                }
                if(kk != kt) throw Error("jordan");
            }
        }
        return make_pair(q, jcs);
    }
    
    matrix proj_mat(const matrix & a0, const matrix & a1) {
        symbol lambda("L");
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
                    s(deg, j) = exfactor(HepLib::exnormal(e.coeff(lambda, deg)));
                }
            }
            ws.add_rows(s);
        }
        
        if(ws.basis_rows().has(lambda)) throw Error("Error: ws.basis_rows().has(lambda)");
        ws.normalize();
        matrix wscols = ws.basis_cols();
        if (ws.dim() == 0) {
            cout << endl;
            cout << "a0: " << a0 << endl;
            cout << "a1: " << a1 << endl;
            throw Error("Error: matrix is Moser-irreducible");
        }
        
        auto dualbs = dual_basis(wscols);
        
        if(dualbs.size()<=0) throw Error("Error: dualbs.size()<=0");
        
        auto dualb = dualbs[0];
        matrix p = exnormal(wscols.mul(dualb));
        return p;
    }

}






