
#include "DE.h"

// https://github.com/magv/fuchsia.cpp 

namespace HepLib::D_E {

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
                if (!basis(nrows - 1, c).normal().is_zero()) goto done;
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
            v.let_op(j) = basis(i, j);
        }
        return v;
    }

    matrix vspace::basis_row(unsigned i) const {
        if(i >= basis.rows()) throw Error("Error");
        matrix v(1, basis.cols());
        for (unsigned j = 0; j < basis.cols(); j++) {
            v.let_op(j) = basis(i, j);
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
                vv.let_op(p) = normal(vv.op(p));
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
                vv.let_op(p) = 0;
                for (unsigned j = p + 1; j < basis.cols(); j++) {
                    vv.let_op(j) = normal(vv.op(j)*b_ip - basis(i, j)*vv_p);
                }
            }
        }
        for (unsigned i = p; i < basis.cols(); i++) {
            vv.let_op(i) = normal(vv.op(i));
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
        matrix s = normal(m).solve(v, zero);
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

    pair<matrix, vector<pair<ex,int>>> jordan(const matrix &m) {
        if(m.cols() != m.rows()) throw Error("jordan");
        unsigned n = m.rows();
        map<ex, unsigned, ex_is_less> eval2almul = eigenvalues(m);
        matrix q(n, n);
        vector<pair<ex,int>> jcs;
        int idxq = 0;
        for (const auto &kv : eval2almul) {
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

}






