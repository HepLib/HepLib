#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <chrono>
#include <fstream>
#include <tuple>
#include <list>

using namespace GiNaC;
using namespace std;

/*-----------------------------------------------------*/
/*          Functions from fuchsia.cpp                   */
/*          https://github.com/magv/fuchsia.cpp          */
/*-----------------------------------------------------*/

static bool PARANOID = false;

//----------------------------------------
// fuchsia_error
//----------------------------------------
class fuchsia_error : public exception {
    const char *message;
    public:
    fuchsia_error(const char *message) : message(message) {};
    virtual const char* what() const throw() { return message; }
};

//----------------------------------------
// matrix load & save
//----------------------------------------
matrix load_matrix(istream &i, parser &reader) {
    for (int depth = 0;;) {
        int c = i.get();
        if (c == EOF)
            throw parse_error("got EOF before any data");
        if (isspace(c)) continue;
        if (c == '(') {
            int c = i.get();
            if (c == EOF)
                throw parse_error("got EOF before any data");
            if (c == '*') { depth++; continue; }
            throw parse_error("matrix file should start with a '{' or a (* comment *)");
        }
        if (depth == 0) {
            if (c == '{') { i.putback(c); break; }
            throw parse_error("matrix file should start with a '{' or a (* comment *)");
        } else {
            if (c == '*') {
                int c = i.get();
                if (c == EOF)
                    throw parse_error("EOF before any data");
                if (c == ')') { depth--; continue; }
            }
        }
    }
    ex x = reader(i);
    // TODO: check that the input is indeed a matrix of rational
    // expressions, with no imaginary or irrational numbers; signal
    // errors if this is not the case.
    matrix m = ex_to<matrix>(lst_to_matrix(ex_to<lst>(x)));
    return m;
}

matrix load_matrix(const char *filename, parser &reader) {
    ifstream i(filename);
    if (!i) throw parse_error("the matrix file was not found");
    return load_matrix(i, reader);
}

void save_matrix(ostream &f, const matrix &m) {
    f << "{";
    for (unsigned i = 0; i < m.rows(); i++) {
        if (i != 0) f << ",\n";
        f << "{";
        for (unsigned j = 0; j < m.cols(); j++) {
            if (j != 0) f << ",";
            f << m(i, j);
        }
        f << "}";
    }
    f << "}";
}

void save_matrix(const char *filename, const matrix &m) {
    // Write to a temporary file, rename it to the correct location
    // at the end. This is done so that if the filesystem has
    // atomic renames, the resulting file would never contain
    // partial output, even if Fuchsia was Ctrl-C'ed.
    string tmpfilename = string(filename) + string(".tmp~");
    ofstream f(tmpfilename);
    save_matrix(f, m);
    f.close();
    int r = rename(tmpfilename.c_str(), filename);
    if (r != 0) {
        throw system_error(r, system_category());
    }
}

//----------------------------------------
// matrix functions
//----------------------------------------
matrix imatrix(unsigned n) {
    matrix id(n, n);
    for (unsigned i = 0; i < n; i++) {
        id(i, i) = 1;
    }
    return id;
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

matrix matrix_solve_left(const matrix &m, const matrix &vars, const matrix &rhs) {
    return m.transpose().solve(vars.transpose(), rhs.transpose()).transpose();
}

template <typename F>
void matrix_map_inplace(matrix &m, F f) {
    for (unsigned i = 0; i < m.nops(); i++) {
        m.let_op(i) = f(m.op(i));
    }
}

template <typename F>
matrix matrix_map(const matrix &m, F f) {
    matrix r(m.rows(), m.cols());
    for (unsigned i = 0; i < m.nops(); i++) {
        r.let_op(i) = f(m.op(i));
    }
    return r;
}

matrix normal(const matrix &m) {
    return matrix_map(m, [&](auto &&e) { return normal(e); });
}

ex ratcan(const ex &e) {
    ex nd = e.normal().numer_denom();
    return nd.op(0).expand()/nd.op(1).expand();
}

matrix ratcan(const matrix &m) {
    return matrix_map(m, [&](auto &&e) { return ratcan(e); });
}

//----------------------------------------
// factor functions
//----------------------------------------
template <typename F>
void factor_iter(const ex &e, F yield) {
    if (is_a<mul>(e)) {
        for (const auto &f : e) {
            if (is_a<power>(f)) {
                yield(f.op(0), ex_to<numeric>(f.op(1)).to_int());
            } else {
                yield(f, 1);
            }
        }
    } else {
        if (is_a<power>(e)) {
            yield(e.op(0), ex_to<numeric>(e.op(1)).to_int());
        } else {
            yield(e, 1);
        }
    }
}

template<typename F>
void factor_and_iter(const ex &e, F yield) {
    factor_iter(e,
        [&](const ex &f1, int k1) {
            factor_iter(factor(f1),
                [&](const ex &f2, int k2) {
                    yield(f2, k1*k2);
                });
        });
}

ex factor_fixed(const ex &e) {
    ex res = 1;
    factor_and_iter(e,
        [&](const ex &f, int k) {
            res *= pow(f, k);
        });
    return res;
}

template <typename F>
void term_iter(const ex &e, F yield) {
    if (is_a<add>(e)) {
        for (const auto &t : e) {
            yield(t);
        }
    } else {
        yield(e);
    }
}

//----------------------------------------
// eigenvalues
//----------------------------------------
map<ex, unsigned, ex_is_less> eigenvalues(const matrix &m, bool skip_roots=false) {
    map<ex, unsigned, ex_is_less> eigenvalues;
    symbol lambda("L");
    
    //ex charpoly = m.charpoly(lambda);
    // determinant_algo :: automatic, gauss, divfree, laplace, bareiss
    ex charpoly = imatrix(m.rows()).mul_scalar(lambda).sub(m).determinant(1);
    
    factor_iter(factor_fixed(normal(charpoly).numer()),
        [&](const ex &f, int k) {
            int deg = f.degree(lambda);
            if (deg == 0) {
                // No roots in λ here.
            }
            else if (deg == 1) {
                // f == c0 + c1*λ
                ex c0 = f.coeff(lambda, 0);
                ex c1 = f.coeff(lambda, 1);
                eigenvalues[ratcan(-c0/c1)] += k;
            }
            else {
                if (skip_roots) {
                    cout << "skipping eigenvalues with roots: RootOf" << endl;
                } else {
                    cout << "could not factor this part of charpoly" << endl;
                    throw fuchsia_error("eigenvalues(): eigenvalues contain roots of 2nd degree or higher");
                }
            }
        });
    return eigenvalues;
}

//----------------------------------------
// matrix_hack
//----------------------------------------
class matrix_hack : public matrix {
    public:
    void append_rows(const matrix &src);
    exvector &mvec();
    void resize(unsigned nrows);
};

void matrix_hack::append_rows(const matrix &src) {
    assert(col == src.cols());
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

//----------------------------------------
// rescale_submatrix
//----------------------------------------
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
            assert(ok1 && ok2);
            m(i, j) = nn*dd;
        }
    }
}

//----------------------------------------
// ex_to_matrix
//----------------------------------------
matrix ex_to_matrix(const ex &m) {
    return ex_to<matrix>(m.evalm());
}

//----------------------------------------
// echelon_form_gauss
//----------------------------------------
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

//----------------------------------------
// vspace
//----------------------------------------
struct vspace {
    matrix basis;
    vspace(unsigned n);
    vspace(const matrix &b);
    unsigned dim() const;
    matrix basis_col(unsigned i) const;
    matrix basis_row(unsigned i) const;
    const matrix basis_cols() const;
    const matrix &basis_rows() const;
    bool contains(const matrix &v) const;
    void add_rows(const matrix &v);
    void normalize();
};

vspace::vspace(unsigned n) : basis(0, n) {
}

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
    assert(i < basis.rows());
    matrix v(basis.cols(), 1);
    for (unsigned j = 0; j < basis.cols(); j++) {
        v.let_op(j) = basis(i, j);
    }
    return v;
}

matrix vspace::basis_row(unsigned i) const {
    assert(i < basis.rows());
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
    assert(v.nops() == basis.cols());
    matrix vv = v;
    rescale_submatrix(vv, 0, v.rows(), 0, v.cols());
    unsigned p = 0;
    // Division-free subtraction of basis vectors from v.
    for (unsigned i = 0; i < basis.rows(); i++, p++) {
        // Advance p to the first non-zero column of basis[i].
        for (;;) {
            // This assertion should only fail if the normalize()
            // was not called between add_rows() and contains().
            assert(p < basis.cols());
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

//----------------------------------------
// nullspace
//----------------------------------------
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

//----------------------------------------
// eigenspace
//----------------------------------------
vspace eigenspace(const matrix &m, const ex &eval) {
    matrix mm = m;
    for (unsigned i = 0; i < m.rows(); i++) {
        mm(i, i) -= eval;
    }
    return nullspace(mm);
}

//----------------------------------------
// eigenspace_left
//----------------------------------------
vspace eigenspace_left(const matrix &m, const ex &eval) {
    return eigenspace(m.transpose(), eval);
}

//----------------------------------------
// eigenvectors_right
//----------------------------------------
vector<matrix> eigenvectors(const matrix &m, const ex &eval) {
    vspace es = eigenspace(m, eval);
    vector<matrix> evectors(es.dim());
    for (unsigned i = 0; i < es.dim(); i++) {
        evectors[i] = es.basis_col(i);
    }
    return evectors;
}

//----------------------------------------
// eigenvectors_left
//----------------------------------------
vector<matrix> eigenvectors_left(const matrix &m, const ex &eval) {
    vspace es = eigenspace_left(m, eval);
    vector<matrix> evectors(es.dim());
    for (unsigned i = 0; i < es.dim(); i++) {
        evectors[i] = es.basis_row(i);
    }
    return evectors;
}

//----------------------------------------
// jordan
//----------------------------------------
pair<matrix, vector<pair<ex,int>>> jordan(const matrix &m) {
    assert(m.cols() == m.rows());
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
                assert(k >= 0);
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
            assert(kk == kt);
        }
    }
    if (PARANOID) {
        assert(q.rank() == n);
    }
    return make_pair(q, jcs);
}

//----------------------------------------
// dot
//----------------------------------------
ex dot(const matrix &v1, const matrix &v2) {
    assert(v1.nops() == v2.nops());
    ex res = 0;
    for (unsigned i = 0; i < v1.nops(); i++) {
        res += v1.op(i)*v2.op(i);
    }
    return res;
}

//----------------------------------------
// cross
//----------------------------------------
matrix cross(const matrix &v1, const matrix &v2) {
    matrix res(v1.nops(), v2.nops());
    for (unsigned i = 0; i < v1.nops(); i++) {
        const ex &a = v1.op(i);
        for (unsigned j = 0; j < v2.nops(); j++) {
            res(i, j) = a*v2.op(j);
        }
    }
    return res;
}

//----------------------------------------
// DE transformation
//----------------------------------------
matrix matrix_diff(const matrix &mat, const symbol &x, const int n=1) {
    return matrix_map(mat, [&](auto &&e){ return e.diff(x,n); });
}

matrix transform(const matrix &m, const matrix &t, const symbol &x) {
    return t.inverse().mul(m.mul(t).sub(matrix_diff(t,x)));
}

matrix transform(const matrix &m, const matrix &tinverse, const matrix &t, const symbol &x) {
    return tinverse.mul(m.mul(t).sub(matrix_diff(t,x)));
}

