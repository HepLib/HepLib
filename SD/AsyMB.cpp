/**
 * @file
 * @brief Functions to do Asy
 */
 
#include "SD.h"
#include <math.h>
#include <cmath>

namespace HepLib::SD {

    /**
     * @brief PRank from FIESTA
     * @param m input matrix
     * @return PRank for matrix
     */
    int SecDec::PRank(matrix m) {
        int nr = m.rows();
        int nc = m.cols();
        int pr = 0;
        if(nr>1) {
            matrix pr_mat(nr, nc);
            for(int ic=0; ic<nc; ic++) {
                for(int ir=0; ir<nr; ir++) {
                    pr_mat(ir,ic) = m(ir,ic)-m(0,ic);
                }
            }
            pr = pr_mat.rank();
        }
        return pr;
    }

    /**
     * @brief PExpand from asy2.1.1.m
     * @param xpol polynomial
     * @param delta has delta or not
     * @return PExpand result
     */
    ex SecDec::PExpand(ex xpol, bool delta) {
        lst nlst;
        ex pol = collect_common_factors(xpol.expand());
        if(is_a<mul>(pol)) {
            ex tmp = 1;
            for(int i=0; i<pol.nops(); i++) {
                if(is_a<numeric>(pol.op(i))) continue;
                if(pol.op(i).match(x(w))) continue;
                if(pol.op(i).match(pow(x(w1),w2))) continue;
                if(pol.op(i).match(pow(y(w1),w2))) continue;
                if(pol.op(i).match(pow(vs,w))) continue;
                tmp *= pol.op(i);
            }
            pol = tmp;
        }
        pol = pol.expand();
        auto xs = get_xy_from(pol);
        int nx = xs.size();
        int id = delta ? nx-1 : nx;
        
        lst lxs;
        for(auto item : xs) lxs.append(item);
        lxs.append(vs);
        pol = pol.expand().collect(lxs, true);
        int np = is_a<add>(pol) ? pol.nops() : 1;
        matrix rs_mat(np, nx+1);
            
        for(int n=0; n<np; n++) {
            ex tmp = pol.op(n);
            rs_mat(n,0) = tmp.degree(vs);
            for(int ix=0; ix<nx; ix++) {
                rs_mat(n, ix+1) = collect_ex(tmp, xs[ix]).degree(xs[ix]);
            }
        }

        int pr = PRank(rs_mat);
        if(pr-1 != id) return nlst;
        
        matrix rp_mat(np, id+1);
        for(int n=0; n<np; n++) {
            for(int ix=0; ix<id+1; ix++) {
                rp_mat(n,ix) = rs_mat(n,ix);
            }
        }
        int rp = PRank(rp_mat);
        
        if(pr != rp) {
            cout << "pr=" << pr << ", rp=" << rp << endl;
            throw Error("PExpand: projection method does not work.");
        }
        
        auto fs = SecDecG::RunQHull(rp_mat);
        lst vs, ret;
        for(int i=0; i<id+1; i++) {
            symbol v;
            vs.append(v);
        }
        auto vv = vs.op(id);
        for(auto fi : fs) {
            lst eqns;
            for(auto pi : fi) {
                ex eqn = rs_mat(pi,0);
                for(int i=0; i<id; i++) {
                    eqn += vs.op(i) * rs_mat(pi,i+1);
                }
                eqns.append(eqn==vv);
            }
            auto lss = lsolve(eqns, vs);
            if(lss.nops()>1) {
                bool bf = true;
                auto vs2 = subs(vs,lss);
                for(int r=0; r<np; r++) {
                    ex eqn = rs_mat(r,0);
                    for(int i=0; i<id; i++) {
                        eqn += vs2.op(i) * rs_mat(r,i+1);
                    }
                    ex chk = eqn-vs2.op(id);
                    if(!is_a<numeric>(chk)) {
                        cerr << ErrColor << "chk is NOT a number: " << chk << RESET << endl;
                        throw Error("PExpand: chk is NOT a number.");
                    }
                    if(chk<0) {
                        bf = false;
                        break;
                    }
                }
                if(bf) {
                    vs2[id] = 0;
                    ex min = vs2.op(0);
                    for(auto vsi : vs2) {
                        if(vsi<min) min = vsi;
                    }
                    for(int i=0; i<vs2.nops(); i++) vs2[i] = vs2.op(i)-min;
                    ret.append(vs2);
                }
            }
        }
        
        lst lxs2;
        for(auto item : xs) lxs2.append(item);
        ret.prepend(lxs2);
        return ret;
    }

    /**
     * @brief DoAsy 
     */
    void SecDec::DoAsy() {
        
        while(true) {
            exvector funexp;
            for(auto fe : FunExp) {
                funexp.push_back(fe);
            }
            FunExp.clear();
            FunExp.shrink_to_fit();
            
            bool ret = false;
            for(auto fe : funexp) {
                if(fe.nops()<=2) {
                    cerr << ErrColor << "needs 3 elements: " << fe << RESET << endl;
                    throw Error("DoAsy: we needs 3 elements.");
                }
                ex ft = fe.op(0).op(1).subs(vs==0);
                ex eqn;
                bool ok2 = true;
                auto xs = get_x_from(ft);
                for(int i=0; i<xs.size(); i++) { // keep only 2 xi's
                    for(int j=i+1; j<xs.size(); j++) { // keep only 2 xi's
                        bool delta_ij = false;
                        for(int di=0; di<fe.op(2).nops(); di++) {
                            if(fe.op(2).op(di).has(xs[i]) && fe.op(2).op(di).has(xs[j])) {
                                delta_ij = true;
                                break;
                            }
                        }
                        if(!delta_ij) continue;
                        
                        symbol xi("xi"), xj("xj");
                        auto ftt = ft.subs(lst{xs[i]==xi, xs[j]==xj});
                        ftt = Factor(ftt);
                        lst fts2;
                        if(is_a<mul>(ftt)) {
                            for(auto item : ftt) fts2.append(item);
                        } else {
                            fts2.append(ftt);
                        }
                        lst tmp_lst;
                        for(auto item : fts2) {
                            auto tmp = Factor(item.subs(x(w)==0));
                            if(is_a<mul>(tmp)) {
                                for(auto it : tmp) tmp_lst.append(it);
                            } else {
                                tmp_lst.append(tmp);
                            }
                        }
                        fts2 = tmp_lst;
                                            
                        for(auto item : fts2) {
                            auto tmp = item;
                            if(item.match(pow(w1,w2))) tmp = item.op(0);
                            if(tmp.degree(xi)==1 && tmp.degree(xj)==1) {
                                ex ci = tmp.coeff(xi);
                                ex cj = tmp.coeff(xj);
                                if((ci*xi+cj*xj-tmp).is_zero() && is_a<numeric>(ci*cj) && (ci*cj)<0) {
                                    eqn = tmp.subs(lst{xi==xs[i], xj==xs[j]});
                                    ok2 = false;
                                    goto OK2;
                                }
                            }
                        }
                    }
                }
                
                OK2:
                if(!ok2) {
                    auto xij = get_x_from(eqn);
                    ex xi = xij[0];
                    ex xj = xij[1];
                    
                    ex ci = eqn.coeff(xi);
                    ex cj = eqn.coeff(xj);
                    
                    // handle eqn==ci xi - cj xj
                    ret = true;
                    ci = abs(ci);
                    cj = abs(cj);
                    symbol yi,yj;
                    // Part I: ci xi-cj xj>0, i.e., xi>cj/ci xj
                    auto f1 = ex_to<lst>(fe.op(0));
                    auto e1 = ex_to<lst>(fe.op(1));
                    ex c1 = cj/ci;
                    for(int i=0; i<f1.nops(); i++) {
                        f1[i] = f1.op(i).subs(lst{xi==c1*yj/(1+c1)+yi,xj==yj/(1+c1)}).subs(lst{yi==xi,yj==xj});
                    }
                    if(e1.op(0)==1) f1[0] = f1.op(0)/(1+c1); // Jaccobi
                    else if(e1.op(1)==1) f1[1] = f1.op(1)/(1+c1);
                    else {
                        f1.append(1/(1+c1));
                        e1.append(1);
                    }
                    
                    FunExp.push_back(lst{f1,e1,fe.op(2)});
                    // Part II: ci xi-cj xj<0, i.e., i.e., xj>ci/cj xi
                    auto f2 = ex_to<lst>(fe.op(0));
                    auto e2 = ex_to<lst>(fe.op(1));
                    ex c2 = ci/cj;
                    for(int i=0; i<f2.nops(); i++) {
                        f2[i] = f2.op(i).subs(lst{xj==c2*yi/(1+c2)+yj,xi==yi/(1+c2)}).subs(lst{yi==xi,yj==xj});
                    }
                    if(e2.op(0)==1) f2[0] = f2.op(0)/(1+c2); // Jaccobi
                    else if(e2.op(1)==1) f2[1] = f2.op(1)/(1+c2);
                    else {
                        f2.append(1/(1+c2));
                        e2.append(1);
                    }
                    FunExp.push_back(lst{f2,e2,fe.op(2)});
                } else {
                    FunExp.push_back(fe);
                }
            }
            if(!ret) break;
        }
        
        
        for(auto fe : FunExp) {
            ex expn = 0;
            symbol s;
            for(int i=0; i<fe.op(0).nops(); i++) {
                auto item = fe.op(0).op(i).subs(x(w)==s*y(w)).subs(y(w)==x(w));
                if(!item.has(s)) continue;
                item = collect_ex(item, s);
                if(item.ldegree(s)!=item.degree(s)) {
                    cerr << ErrColor << "Not Homogeneous: " << s << RESET << endl;
                    throw Error("DoAsy: Not Homogeneous");
                }
                expn += item.degree(s) * fe.op(1).op(i);
            }
            auto xsize = get_x_from(fe.op(0)).size();
            if(!normal(expn+xsize).is_zero()) {
                cout << ErrColor << "expn=" << expn << ", xsize=" << xsize << RESET << endl;
                throw Error("DoAsy: expn + xsize != 0.");
            }
        }
        
        
        auto fes = FunExp;
        FunExp.clear();
        FunExp.shrink_to_fit();
        
        for(auto fe : fes) {
            ex xpol = 1;
            lst uf;
            for(int i=0; i<fe.op(0).nops(); i++) {
                if(fe.op(1).op(i).info(info_flags::nonnegint)) continue;
                uf.append(fe.op(0).op(i));
            }
            uf.sort();
            uf.unique();
            //sort_lst(uf);
            for(auto item : uf) xpol *= item;
            if(!xpol.has(vs)) {
                FunExp.push_back(fe);
                continue;
            }
            
            bool has_delta = true;
            auto rs = PExpand(xpol, has_delta);
            if(rs.nops()<1) {
                throw Error("PExpand returned with nothing, even without hard region!");
            }
            if(Verbose>10) {
                cout << PRE << "\\--Asy Regions:" << (rs.nops()-1) << endl;
                if(rs.nops()>1) {
                    for(auto ri : rs) cout << "     " << ri << endl;
                }
            }
            
            auto r0 = rs.op(0);
            auto r0y = subs(r0,x(w)==y(w));
            for(int i=1; i<rs.nops(); i++) {
                lst srepl;
                auto ri = rs.op(i);
                ex vs_pow = 0;
                for(int j=0; j<r0.nops(); j++) {
                    srepl.append(r0.op(j)==r0y.op(j) * pow(vs, ri.op(j)));
                    vs_pow += ri.op(j);
                }

                auto fs = subs(fe.op(0), srepl);
                fs = subs(fs, y(w)==x(w));
                auto es = fe.op(1);
                ex fpre = fs.op(0);
                if(!(es.op(0)-1).is_zero()) {
                    cerr << ErrColor << "op(0) is Not 1: " << es << RESET << endl;
                    throw Error("DoAsy: op(0) is Not 1");
                }
                
                lst fs2, es2;
                for(int j=1; j<fs.nops(); j++) {
                    auto fj = fs.op(j);
                    auto tmp = fj.expand();
                    auto vsp = 0;
                    try {
                        vsp = tmp.ldegree(vs);
                    } catch(exception &e) {
                        cout << e.what() << endl;
                        throw Error("DoAsy: non-integer exponent.");
                    }
                    vs_pow += vsp * es.op(j);
                    tmp = collect_common_factors(tmp)/pow(vs,vsp);
                    fs2.append(tmp);
                    es2.append(es.op(j));
                }
                exmap eps_map;
                for(auto epi : eps_lst) eps_map[epi.op(0)] = 0;

                auto vsn0 = vsRank(fpre); // maybe need to expand ep/eps first
                auto vsn = vsn0 + vs_pow.subs(eps_map);
                int di=0;
                lst fss, ess;
                fss.append(fs2);
                ess.append(es2);
                while(di<=vsN-vsn) { // fss, ess will get updated
                    lst fss2, ess2;
                    for(int ife=0; ife<fss.nops(); ife++) {
                        lst fs3 = ex_to<lst>(fss.op(ife));
                        lst es3 = ex_to<lst>(ess.op(ife));
                        if(di<vsN-vsn) {
                            for(int ii=0; ii<fs3.nops(); ii++) {
                                lst fs4 = fs3;
                                lst es4 = es3;
                                auto dit = diff_ex(fs4.op(ii),vs);
                                if(!dit.is_zero()) {
                                    if((es4.op(ii)-1).is_zero()) {
                                        fs4[ii] = dit;
                                    } else {
                                        fs4.append(es4.op(ii)*dit);
                                        es4[ii] = es4.op(ii)-1;
                                        es4.append(1);
                                    }
                                    fss2.append(fs4);
                                    ess2.append(es4);
                                }
                            }
                        }
                        fs3 = ex_to<lst>(subs(fs3,vs==0));
                        fs3.prepend(fpre/factorial(di) * pow(vs,di+vs_pow));
                        es3.prepend(1);
                        FunExp.push_back(lst{fs3,es3,fe.op(2)});
                    }
                    fss = fss2;
                    ess = ess2;
                    di++;
                }
            }
        }
    }
    
    /**
     * @brief MB
     */
    void SecDec::MB() {
        for(auto &fe : FunExp) {
            if(fe.has(vz) || !fe.has(vs)) continue; // 2nd entrance
            
            // check epz
            if(fe.has(epz)) {
                cout << ErrColor << "MB: epz found at fe = " << fe << RESET << endl;
                throw Error("MB: epz found at fe.");
            }
            
            exmap eps_map;
            ex epn = ex(1)/111;
            for(auto epi : eps_lst) {
                eps_map[epi.op(0)] = epn;
                epn = epn / 13;
            }
            
            // check variables besides x or PL
            // CV should only appear at kv.op(0).op(0), i.e., the prefactor
            for(int i=2; i<fe.op(0).nops(); i++) {
                // make sure only Constant/F terms can contain small variable: vs
                if(i!=1 && fe.op(0).op(i).has(vz)) {
                    cout << "vz Found @ " << i << " of " << fe.op(0) << endl;
                    throw Error("MB: vz Found.");
                }
                
                auto tmp = fe.op(0).op(i).subs(lst{WRA(w)==0,x(w)==1,PL(w)==1}).subs(eps_map);
                if(!is_a<numeric>(NN(tmp,20))) {
                    cout << ErrColor << tmp << RESET << endl;
                    cout << fe.op(0) << endl;
                    throw Error("MB: Extra Variable(^[ep,eps,PL,x]) Found.");
                }
            }
        
            ex ft = fe.op(0).op(1);
            if(ft.has(vs)) {
                ft = collect_ex(ft, vs);
                if(!ft.is_polynomial(vs) || (ft.degree(vs)-1)!=0) {
                    cout << ErrColor << "Not supported F-term with s: " << ft << RESET << endl;
                    throw Error("MB: Not supported F-term with s.");
                }
                ex expn = ex(0)-fe.op(1).op(1);
                // (2*Pi*I) dropped out, since we will take residue later.
                fe[0][0] = fe.op(0).op(0) * tgamma(expn+vz)*tgamma(ex(0)-vz)/tgamma(expn)*pow(vs,vz);
                ex w1 = ft.coeff(vs);
                ex w2 = ft.subs(vs==0);
                if(!w2.is_zero()) {
                    if(xPositive(w1)) {
                        let_op_append(fe, 0, w1);
                        let_op_append(fe, 1, vz);
                        fe[0][1] = w2;
                        fe[1][1] = fe.op(1).op(1)-vz-epz;
                    } else if(xPositive(w2)) {
                        let_op_append(fe, 0, w2);
                        let_op_append(fe, 1, fe.op(1).op(1)-vz-epz);
                        fe[0][1] = w1;
                        fe[1][1] = vz;
                    } else if(xPositive(ex(0)-w1)) {
                        let_op_append(fe, 0, ex(0)-w1);
                        let_op_append(fe, 1, vz);
                        let_op_append(fe, 0, exp(-I*Pi*vz));
                        let_op_append(fe, 1, 1);
                        fe[0][1] = w2;
                        fe[1][1] = fe.op(1).op(1)-vz-epz;
                    } else if(xPositive(ex(0)-w2)) {
                        let_op_append(fe, 0, ex(0)-w2);
                        let_op_append(fe, 1, fe.op(1).op(1)-vz-epz);
                        let_op_append(fe, 0, exp(-I*Pi*(fe.op(1).op(1)-vz-epz)));
                        let_op_append(fe, 1, 1);
                        fe[0][1] = w1;
                        fe[1][1] = vz;
                    } else {
                        cout << "w1=" << w1 << endl;
                        cout << "w2=" << w2 << endl;
                        throw Error("MB: Neither w1 nor w2 is xPositive.");
                    }
                }
            }
        }
    }
    
    
}
