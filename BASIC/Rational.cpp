// Rational Polynomial ReConstruction
#include "BASIC.h"
#include <vector>
#include <execution>

namespace HepLib {

    // FROM FIRE6.m
    numeric RationalReconstruct(numeric a, numeric p) {
        numeric g = a, s = 1, t = 0, g1 = p, s1 = 0, t1 = 1;
        while(g*g>p) {
            numeric o[6];
            int i = 0;
            o[i++] = g1;
            o[i++] = s1;
            o[i++] = t1; 
            o[i++] = g - iquo(g,g1)*g1;
            o[i++] = s - iquo(g,g1)*s1;
            o[i++] = t - iquo(g,g1)*t1;
            i = 0;
            g = o[i++]; 
            s = o[i++]; 
            t = o[i++]; 
            g1 = o[i++]; 
            s1 = o[i++]; 
            t1 = o[i++];
        }
        return g/s;
    }
    
    // FROM FIRE6.m
    numeric mulInv(numeric a, numeric b) {
        numeric b0 = b;
        numeric x0 = 0;
        numeric x1 = 1;

        if (b == 1) return 1;

        while (a > 1) {
            numeric q = iquo(a, b);
            numeric amb = mod(a, b);   
            a = b;
            b = amb;

            numeric xqx = x1 - q * x0;
            x1 = x0;
            x0 = xqx;
        } 
        if (x1 < 0) x1 += b0;

        return x1;
    }

    numeric ChineseRemainder(std::vector<numeric> a, std::vector<numeric> n) {
        numeric prod = 1;
        for (int i = 0; i < n.size(); i++) prod *= n[i];
     
        numeric sm = 0;
        for (int i = 0; i < n.size(); i++) {
            numeric p = iquo(prod, n[i]);
            sm += a[i] * mulInv(p, n[i]) * p;
        }
     
        return mod(sm, prod);
    }
    
    numeric RationalReconstruct(vector<numeric> aa, vector<numeric> pp) { 
        numeric prod = 1, res = 0, resQ = 1, i, temp, steps = 0, resQprev = 0;
        for(int i=0; i<aa.size(); i++) {
            res = ChineseRemainder({res, aa[i]}, {prod, pp[i]});
            prod = prod*pp[i];
            resQ = RationalReconstruct(res, prod);
            if(resQ == resQprev) return resQ;
            resQprev = resQ;
            steps++;
        }
        throw Error("RationalReconstruct unstable!");
        return resQ;
    }

    // FROM FIRE6.m
    ex Thiele(const exvector & keys, const exvector & values, const ex & d) {
        int n = keys.size();
        exvector coeffs(n);
        for(int i=0; i<n; i++) coeffs[i] = 0;
        for(int i=0; i<n; i++) {
            auto temp = values[i];
            for(int j=0; j<i; j++) {
                if(temp.is_equal(coeffs[j])) { 
                    n = i-1;
                    goto out;
                }
                temp = (keys[i] - keys[j])/(temp - coeffs[j]);
            }
            coeffs[i] = temp;
        }
        throw Error("Thiele unstable!");
        out: ;
        ex result = coeffs[n];
        for(int i=n-1; i>=0; i--) result = coeffs[i] + (d - keys[i])/result;
        result = normal(result);
        //result = exnormal(result);
        return result;
    }
    
    ex Newton(const exvector & keys, const exvector & values, const ex & d, const ex factor) {
        int n = keys.size();
        exvector coeffs(n);
        for(int i=0; i<n; i++) coeffs[i] = 0;
        for(int i=0; i<n; i++) {
            auto temp = values[i] * factor.subs(d==keys[i]);
            temp = normal(temp); 
            //temp = exnormal(temp); 
            for(int j=0; j<i; j++) {
                if(temp.is_equal(coeffs[j])) { 
                    n = i-1;
                    goto out;
                }
                temp = (temp - coeffs[j])/(keys[i] - keys[j]);
                temp = normal(temp);
                //temp = exnormal(temp);
            }
            coeffs[i] = temp;
        }
        throw Error("Newton unstable!");
        out: ;
        ex result = coeffs[n];
        for(int i=n-1; i>=0; i--) result = coeffs[i] + (d - keys[i]) * result;
        result = normal(result/factor);
        //result = exnormal(result/factor);
        return result;
    }
    
}
