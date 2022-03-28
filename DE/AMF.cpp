/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"
#include <cmath>

namespace HepLib {
    
    using namespace EoD;
        
    AMF::AMF(IBP & _ibp) : ibp(_ibp), x(iet) { } 
    
    void AMF::InitDE() {
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Generating DE @ " << now() << endl;
        
        if(ibp.MIntegrals.nops()<1) {
            ibp.Reduce();
            ibp.RM(false); // since Propagators will be changed
        }
        if(ibp.MIntegrals.nops()<1) throw Error("No MI Found, maybe No need DE!");
        ibp.FindRules(true);
        Rules = ibp.Rules;
        for(int i=0; i<ibp.Propagators.nops(); i++) {
            ibp.Propagators.let_op(i) = ibp.Propagators.op(i) + x;
        }
        while(true) {
            ibp.Rules.remove_all();
            MIntegrals = ibp.MIntegrals;
            sort_lst(MIntegrals);
            ibp.MIntegrals.remove_all();
            ibp.Integrals.remove_all();
            lst dmis;
            for(auto mi : MIntegrals) {
                ex dmi = 0;
                auto ns0 = mi.op(1);
                for(int i=0; i<ns0.nops(); i++) {
                    if(is_zero(ns0.op(i))) continue;
                    auto ns = ns0;
                    ns.let_op(i) = ns0.op(i)+1;
                    dmi -= ns0.op(i)*F(mi.op(0),ns);
                    ibp.Integrals.append(ns);
                }
                dmis.append(dmi);
            }
            ibp.Reduce();
            ibp.RM(true); // keep .start & .config
            ibp.FindRules(true);
            sort_lst(ibp.MIntegrals);

            if(ibp.MIntegrals==MIntegrals) {
                int matN = MIntegrals.nops();
                matrix mat(matN,matN);
                for(int r=0; r<matN; r++) {
                    auto dmi = dmis.op(r).subs(ibp.Rules,nopat);
                    dmi = collect_ex(dmi, F(w1,w2));
                    ex chk = 0;
                    for(int c=0; c<matN; c++) {
                        mat(r,c) = dmi.coeff(MIntegrals.op(c));
                        chk += mat(r,c) * MIntegrals.op(c);
                    }
                    if(!is_zero(normal(chk-dmi))) throw Error("AMFlow::InitDE, Check failed");
                }
                Mat = mat;
                matrix_map_inplace(Mat, [](const ex & e) { return normal(e); });
                break;
            }
        }
        system(("rm -rf "+ibp.WorkingDir).c_str());
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Generated DE @ " << now() << endl;
        
        if(true) { // BC 
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Generating BC @ " << now() << endl;
            exmap lmap;
            static symbol y("y");
            for(auto lp : ibp.Internal) lmap[lp] = lp*y;
            lmap[x] = y*y;
            lst props;
            for(int i=0; i<ibp.Propagators.nops(); i++) {
                props.append(ibp.Propagators.op(i).subs(lmap).expand().coeff(y,2));
            }
            
            FIRE bc;
            bc.ProblemNumber = ibp.ProblemNumber;
            bc.Propagators = props;
            bc.Propagators.sort().unique();
            for(auto mi : MIntegrals) {
                lst intg;
                exmap p2n;
                for(int i=0; i<props.nops(); i++) p2n[props.op(i)] += mi.op(1).op(i);
                for(int i=0; i<bc.Propagators.nops(); i++) intg.append(p2n[bc.Propagators.op(i)]);
                bc.Integrals.append(intg);
                _MIntegrals.append(F(bc.ProblemNumber,intg));
            }
            bc.Internal = ibp.Internal;
            bc.Reduce();
            bc.FindRules(true);
            for(int i=0; i<_MIntegrals.nops(); i++) _MIntegrals.let_op(i) = _MIntegrals.op(i).subs(bc.Rules);
            system(("rm -rf "+bc.WorkingDir).c_str());
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Generated BC @ " << now() << endl;
        }
        
        if(!In_GiNaC_Parallel && Verbose>10) cout << "  \\--DE Poles starting @ " << now() << endl;
        if(true) { // get poles of Mat
            ex den = matrix_den_lcm(Mat);
            exvector fvec;
            if (is_a<mul>(den)) for(const auto &f : den) fvec.push_back(f);
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
                if (deg==0) { }
                else if (deg==1) {
                    ex c0 = b.coeff(x, 0);
                    ex c1 = b.coeff(x, 1);
                    roots.insert(normal(-c0/c1));
                } else if(deg==2) {
                    ex c0 = b.coeff(x, 0);
                    ex c1 = b.coeff(x, 1);
                    ex c2 = b.coeff(x, 2);
                    roots.insert(Rationalize(evalf((-c1+sqrt(c1*c1-4*c2*c0))/(2*c2)),10));
                    roots.insert(Rationalize(evalf((-c1-sqrt(c1*c1-4*c2*c0))/(2*c2)),10));
                } else {
                    cout << "current factor: " << f << endl;
                    throw Error("AMF::InitDE, higher powers found.");
                }
            }
            lst rs;
            for(auto ri : roots) rs.append(ri);
            sort_lst(rs);
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--DE Poles: " << NN(rs,2) << endl;
            pts.remove_all();
            if(roots.size()>0) {
                ex max = -1, min = -1;
                for(auto r : roots) {
                    if(is_zero(r)) continue;
                    ex ar = abs(r);
                    if(max<0 || ar>max) max = ar;
                    if(min<0 || ar<min) min = ar;
                }
                pts.remove_all();
                min /= 2;
                max *= 2;
                ex x0 = I*max; // the last point
                x0 = Rationalize(x0, 20);
                pts.append(x0);
                
                while(true) {
                    ex mm = -1;
                    for(auto r : roots) {
                        ex ar = abs(r-x0);
                        if(mm<0 || ar<mm) mm = ar;
                    }
                    x0 -= I*mm/2;
                    x0 = Rationalize(x0, 20);
                    pts.prepend(x0);
                    if(x0/I<0) throw Error("AMF::InitDE, et<0 FOUND.");
                    if(abs(x0)<min) break;
                }
            } else throw Error("AMF::InitDE, NO root found.");
        }
        if(!In_GiNaC_Parallel && Verbose>10) {
            lst npts;
            for(auto pi : pts) npts.append(NN(pi/I,2));
            cout << "  \\--AMF Points: I*" << npts << endl;
            cout << "  \\--DE Poles finished @ " << now() << endl;
        }
    }
    
    matrix AMF::RU(const ex & x1, const ex & x2) {
        ex dis = x1-x2;
        DE de(Mat,x);
        if(d!=d0) de.subs(d==d0, nopat);
        de.WDigits = WDigits;
        auto mat = de.Taylor(x2,dis,xN);
        return mat;
    }
    
    lst AMF::Evaluate() {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        int matN = Mat.rows();
        int nloop = ibp.Internal.nops();
        
        // DE at infinity
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ infinity ..." << RESET << endl;
        DE oo(Mat,x);
        if(d!=d0) oo.subs(d==d0, nopat);
        oo.WDigits = WDigits;
        oo.x2y(1/x); // infinity to origin
        oo.xpow();
        int npts = pts.nops();
        ex xoo = 1/pts.op(npts-1);
        auto ooU = oo.Series(xoo,xN);
        auto ooT = oo.MatT();
              
        // ooC 
        matrix ooC(matN,1); // F in ooC should be treated as Bubble type
        if(true) {
            for(int i=0; i<matN; i++) ooC(i,0) = iWF(i);
            
            ex ila; // select the lambda of type n-L*d/2 
            bool first = true;
            for(const auto & la : oo.las) {
                if(!normal(la+nloop*d0/2).info(info_flags::integer)) continue;
                if(!first) throw Error("something may be wrong here.");
                ila = la;
                first = false;
            }
            lst dlst;
            int dim = ibp.Propagators.nops();
            for(int i=0; i<matN; i++) {
                ex v = -nloop*d0/2-ila;
                ex ns = MIntegrals.op(i).op(1);
                for(int j=0; j<dim; j++) v += ns.op(j);
                dlst.append(pow(x,v));
            }
            auto mat = ex_to<matrix>(diag_matrix(dlst)).inverse().mul(ooT);
            auto pr = prank(mat,x);
            pr++; if(pr<0) pr = 0;
            
            auto ooCMat = oo.Series(pr);
            if(ooCMat[ila].size()!=1) throw Error("ooCMat[ila].size()!=1"); // log term exist, not consider yet
            
            for(const auto & kv : ooCMat) {
                ex la = kv.first;
                if(!is_zero(la-ila)) {
                    matrix C00 = kv.second[0][0]; // k=0, n=0;
                    for(int i=0; i<matN; i++) {
                        if(!is_zero(C00(i,i))) ooC(i,0) = 0;
                    }
                } 
            }
            
            matrix cmat = ooCMat[ila][0][0];
            for(int n=1; n<=pr; n++) cmat = cmat.add(ooCMat[ila][0][n].mul_scalar(pow(x,n)));
            mat = mat.mul(cmat);
            for(int i=0; i<mat.nops(); i++) mat.let_op(i) = series_ex(mat.op(i),x,0);
            matrix fmat(matN, matN+1);
            for(int r=0; r<matN; r++) {
                for(int c=0; c<matN; c++) fmat(r,c) = mat(r,c);
                fmat(r,matN) = _MIntegrals.op(r).subs(d==d0);
            }
            fmat = fermat_Redrowech(fmat);

            lst feqns;
            for(int r=0; r<matN; r++) {
                int idx = -1;
                for(int c=0; c<matN; c++) {
                    if(!is_zero(fmat(r,c)) && !is_zero(fmat(r,c)-1)) {
                        cout << matN << endl;
                        throw Error("something may be wrong here.");
                    }
                    if(is_zero(fmat(r,c)-1)) {
                        if(idx!=-1) throw Error("something may be wrong here.");
                        idx = c;
                    }
                }
                if(idx!=-1) ooC(idx,0) = fmat(r,matN);
                else {
                    feqns.append(fmat(r,matN));
                    if(!is_zero(fmat(r,matN))) {
                        cout << feqns << endl;
                        throw Error("fmat(r,matN) is NOT zero.");
                    }
                }
            }

            if(false && feqns.nops()>0) { // false to disable
                exset fset;
                find(feqns, F(w1,w2), fset);
                lst fs;
                for(auto fi : fset) fs.append(fi);
                sort_lst(fs);
                int neqs = feqns.nops();
                int nfs = fs.nops();
                matrix fmat(neqs,nfs+1);
                for(int r=0; r<neqs; r++) {
                    for(int c=0; c<nfs; c++) fmat(r,c) = feqns.op(r).coeff(fs.op(c));
                }
                fmat = fermat_Redrowech(fmat);
                exmap fsol;
                for(int r=0; r<neqs; r++) {
                    int first_index = -1;
                    for(int c=0; c<nfs; c++) {
                        if(!is_zero(fmat(r,c))) { first_index = c; break;}
                    }
                    if(first_index!=-1) {
                        ex res = 0;
                        for(int c=first_index+1; c<nfs; c++) res += fmat(r,c) * fs.op(c);
                        fsol[fs.op(first_index)] = res/fmat(r,first_index);
                    }
                }
                for(int r=0; r<matN; r++) ooC(r,0) = ooC(r,0).subs(fsol);
            }
            
            if(ooC.has(iWF(w))) {
                cout << fmat << endl;
                cout << ooC << endl;
                throw Error("ooC is not determined completely.");
            }
        }
            
        // J(xoo)
        matrix ooTUC = ooU.mul(ooC);
        if(WDigits>0) ooTUC = ex_to<matrix>(subs(ooT, x==xoo.evalf())).mul(ooTUC);
        else ooTUC = ex_to<matrix>(subs(ooT, x==xoo)).mul(ooTUC);
        
        // Middle U matrix
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ regular points ..." << RESET << endl;
        matrix MatU = ex_to<matrix>(unit_matrix(matN));
        for(int i=0; i<npts-1; i++) MatU = MatU.mul(RU(pts.op(i),pts.op(i+1)));
                
        // DE at origin
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ origin ..." << RESET << endl;
        DE o(Mat,x);
        if(d!=d0) o.subs(d==d0, nopat);
        o.WDigits = WDigits;
        ex xo = pts.op(0);
        auto oU = o.Series(xo,xN);
        matrix oT = o.MatT();
        matrix ioT = fermat_inv(oT); // oT.inverse();
        if(WDigits>0) ioT = ex_to<matrix>(subs(ioT,x==xo.evalf()));
        else ioT = ex_to<matrix>(subs(ioT,x==xo));
        matrix ioTU;
        if(WDigits>0) ioTU = oU.inverse().mul(ioT);
        else ioTU = fermat_inv(oU).mul(ioT);
        
        // Final C at origin
        matrix oC = ioTU.mul(MatU.mul(ooTUC));
        
        // take x->0 limit at origin
        auto pr = prank(oT,x);
        pr++;
        if(pr<0) pr = 0;
        for(auto ev : o.las) {
            if(!ev.info(info_flags::integer)) continue;
            pr -= ex_to<numeric>(ev).to_int();
            break;
        }
        matrix oU0 = o.Series(x,pr,lst{0}); // only pick up x^integer
        auto oTU = oT.mul(oU0);
        xpow(oTU,x);
        if(WDigits>0) oTU = ex_to<matrix>(oTU.evalf());
        HepLib::subs(oTU,x==0,nopat);
        
        lst res; // result for master integrals
        oC = oTU.mul(oC);
        for(int i=0; i<matN; i++) res.append(oC(i,0));
        Digits = oDigits;
        return res;
    }
    
    lst AMF::FitEps(const lst & eps, int lp, bool parallel) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        exvector eps_vec(eps.begin(), eps.end());
        int nmi = MIntegrals.nops();
        exvector mis_vec[nmi];
        if(parallel && eps.nops()>1) {
            auto od0 = d0;
            auto res_vec = GiNaC_Parallel(eps.nops(), [&](int idx)->ex {
                d0 = 4-2*eps.op(idx);
                return Evaluate();
            }, "AMF");
            d0 = od0;
            for(auto mis : res_vec) {
                for(int i=0; i<nmi; i++) mis_vec[i].push_back(mis.op(i));
            }
        } else {
            auto od0 = d0;
            for(auto epi : eps) {
                d0 = 4-2*epi;
                if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "ep = " << epi << RESET << endl;
                auto mis = Evaluate();
                for(int i=0; i<nmi; i++) mis_vec[i].push_back(mis.op(i));
            }
            d0 = od0;
        }
        
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--Final PolynomialFit ..." << endl;
        if(WDigits>0) Digits = 2*WDigits;
        lst mis_lst;
        for(int i=0; i<nmi; i++) {
            int tn = eps.nops()-1;
            auto cs = PolynomialFit(eps_vec, mis_vec[i], tn, lp);
            ex mi = 0;
            for(int i=0; i<tn; i++) mi += pow(ep,lp+i) * cs.op(i);
            mis_lst.append(mi);
        }
        Digits = oDigits;
        return mis_lst;
    }
    
    lst AMF::FitEps(int goal, int order, bool parallel) { // form AMFlow
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        int nloop = ibp.Internal.nops();
        numeric nn = 5*order/numeric(2)+2*nloop;
        int n = ceil(nn.to_double());
        if(n>100) throw Error("FitEps: too large order.");
        ex eps0 = pow(10, -nloop/numeric(2)-goal/numeric(order+1));
        eps0 = Rationalize(evalf(eps0));
        int lp = -2*nloop;
        lst eps;
        for(int i=0; i<n; i++) eps.append(eps0 + eps0*ex(i+1)/100);
        numeric nsp = (n+2*nloop)*(nloop/numeric(2)+goal/numeric(order+1));
        int sp = ceil(nsp.to_double());
        if(sp<30) sp = 30;
        WDigits = 2*sp;
        xN = 4*sp;
        auto res = FitEps(eps,lp,parallel);
        Digits = oDigits;
        return res;
    }  
        
}

