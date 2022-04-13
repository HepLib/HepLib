/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {

    using namespace EoD;
    
    BJF::BJF(matrix _A, ex _b, int K): A(_A), b(_b) { 
        B = ex_to<matrix>(unit_matrix(A.rows())).mul_scalar(b); 
    }
    
    matrix BJF::operator()(int i, int j) {
        if(i==j) return A;
        else if(i+1==j) return B;
        else throw Error("BJF: 0 matrix block.");
    }
       
    invBJF::invBJF(matrix _A, ex _b, int K): A(_A), b(_b) { 
        U = ex_to<matrix>(unit_matrix(A.rows()));
        matrix a1 = fermat_inv(A);
        Ain[1] = a1;
        matrix al = a1;
        for(int i=2; i<=K+1; i++) al = Ain[i] = al.mul(a1);
    }
    
    matrix invBJF::operator()(int i, int j) {
        if(i>j) throw Error("invBJF: 0 matrix block");
        pair<int, int> key = make_pair(i,j);
        auto f = cache.find(key);
        if(f!=cache.end()) return f->second;
        if(i==j) return cache[key] = Ain[1];
        else {
            int ij = j-i;
            return cache[key] = U.mul_scalar(pow(-b,ij)).mul(Ain[ij+1]);
        }
    }
    
    void SeriesT::Resize() {  
        T.clear();
        T.resize(s+1);
        for(int cm=0; cm<=s; cm++) {
            vector<vector<matrix>> vi(kmax+1);
            for(int i=0; i<=kmax; i++) {
                vector<matrix> vj(kmax+1);
                vi[i] = vj;
            }
            T[cm] = vi;
        }
    }
    void SeriesT::Reset() {
        las.clear();
        K.clear();
        C0.clear();
    }
        
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
        auto den_vec = GiNaC_Parallel(m.nops(), [&m](int idx) {
            auto den = m.op(idx).denom();
            return factor_flint(den);
        });
        
        exmap pn_map;
        for(int i=0; i<den_vec.size(); i++) {
            auto den = den_vec[i];
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

