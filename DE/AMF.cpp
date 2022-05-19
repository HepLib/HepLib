/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {

    namespace {
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
    }
    
    //=*********************************************************************=

    AMF::AMF(IBP & _ibp) : ibp(_ibp), x(iet) { } 
    
    void AMF::ExportDE(const string fn) {
        lst res = { Rules, MIntegrals, pts, Mat, _MIntegrals };
        garWrite(res, fn);
    }
    
    void AMF::ImportDE(const string fn) {
        auto res = garRead(fn);
        Rules = ex_to<lst>(res.op(0));
        MIntegrals = ex_to<lst>(res.op(1)); 
        pts = ex_to<lst>(res.op(2)); 
        Mat = ex_to<matrix>(res.op(3));
        _MIntegrals = ex_to<lst>(res.op(4));
    }

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
                int N = MIntegrals.nops();
                matrix mat(N,N);
                for(int r=0; r<N; r++) {
                    auto dmi = dmis.op(r).subs(ibp.Rules,nopat);
                    dmi = collect_ex(dmi, F(w1,w2));
                    ex chk = 0;
                    for(int c=0; c<N; c++) {
                        mat(r,c) = dmi.coeff(MIntegrals.op(c));
                        chk += mat(r,c) * MIntegrals.op(c);
                    }
                    if(!is_zero(normal(chk-dmi))) throw Error("AMFlow::InitDE, Check failed");
                }
                Mat = mat;
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
                ex base = f; 
                int n = 1;
                if (is_a<power>(f)) {
                    base = f.op(0);
                    n = ex_to<numeric>(f.op(1)).to_int();
                } 
                int deg = base.degree(x);
                if (deg==0) { }
                else if (deg==1) {
                    ex c0 = base.coeff(x, 0);
                    ex c1 = base.coeff(x, 1);
                    ex x1 = normal(-c0/c1);
                    x1 = Rationalize(x1.evalf(), 10);
                    roots.insert(x1);
                } else if(deg==2) {
                    ex c = base.coeff(x, 0);
                    ex b = base.coeff(x, 1);
                    ex a = base.coeff(x, 2);
                    ex x1 = (-b+sqrt(b*b-4*a*c))/(2*a);
                    ex x2 = (-b-sqrt(b*b-4*a*c))/(2*a);
                    x1 = Rationalize(x1.evalf(), 10);
                    x2 = Rationalize(x2.evalf(), 10);
                    roots.insert(x1);
                    roots.insert(x2);
                } else if(deg==3) { // 卡丹公式
                    ex d = base.coeff(x, 0);
                    ex c = base.coeff(x, 1);
                    ex b = base.coeff(x, 2);
                    ex a = base.coeff(x, 3);
                    
                    ex u = (9*a*b*c-27*a*a*d-2*b*b*b)/(54*a*a*a);
                    ex v = sqrt(3*(4*a*c*c*c-b*b*c*c-18*a*b*c*d+27*a*a*d*d+4*b*b*b*d))/(18*a*a);
                    ex m;
                    if(abs(u+v)>=abs(u-v)) m = pow(u+v,1/ex(3));
                    else m = pow(u-v,1/ex(3));
                    ex n = 0;
                    if(!is_zero(m)) n = (b*b-3*a*c)/(9*a*a*m);
                    ex w1 = -1/ex(2)+sqrt(ex(3))/2*I;
                    ex w2 = -1/ex(2)-sqrt(ex(3))/2*I;
                    ex x1 = m+n-b/(3*a);
                    ex x2 = w1*m+w2*n-b/(3*a);
                    ex x3 = w2*m+w1*n-b/(3*a);
                    
                    x1 = Rationalize(x1.evalf(), 10);
                    x2 = Rationalize(x2.evalf(), 10);
                    x3 = Rationalize(x3.evalf(), 10);
    
                    roots.insert(x1);
                    roots.insert(x2);
                    roots.insert(x3);
                } else if(deg==4) { 
                    ex e = base.coeff(x, 0);
                    ex d = base.coeff(x, 1);
                    ex c = base.coeff(x, 2);
                    ex b = base.coeff(x, 3);
                    ex a = base.coeff(x, 4);
                    
                    ex P = (c*c+12*a*e-3*b*d)/9;
                    ex Q = (27*a*d*d+2*c*c*c+27*b*b*e-72*a*c*e-9*b*c*d)/54;
                    ex D = sqrt(Q*Q-P*P*P);
                    ex u;
                    if(abs(Q+D)>=abs(Q-D)) u = pow(Q+D, 1/ex(3));
                    else u = pow(Q-D, 1/ex(3));
                    ex v = 0;
                    if(!is_zero(u)) v = P/u;
                    ex w = -1/ex(2)+sqrt(ex(3))/2*I;
                    
                    ex m=0, S=0;
                    for(int k=1; k<4; k++) {
                        ex wk = pow(w,k-1);
                        ex w4k = pow(w,4-k);
                        ex mk = sqrt(b*b-8*a*c/3+4*a*(wk*u+w4k*v));
                        ex Sk = 2*b*b-16*a*c/3-4*a*(wk*u+w4k*v);
                        if(abs(mk).evalf()>abs(m).evalf()) { m = mk; S = Sk; }
                    }
                    
                    ex T;
                    if(!is_zero(m)) T = (8*a*b*c-16*a*a*d-2*b*b*b)/m;
                    else {
                        m = 0;
                        S = b*b-8*a*c/3;
                        T = 0;
                    }
                    
                    ex x1 = (-b-m+sqrt(S-T))/(4*a);
                    ex x2 = (-b-m-sqrt(S-T))/(4*a);
                    ex x3 = (-b+m+sqrt(S+T))/(4*a);
                    ex x4 = (-b+m-sqrt(S+T))/(4*a);
                    
                    x1 = Rationalize(x1.evalf(), 10);
                    x2 = Rationalize(x2.evalf(), 10);
                    x3 = Rationalize(x3.evalf(), 10);
                    x4 = Rationalize(x4.evalf(), 10);
                    
                    roots.insert(x1);
                    roots.insert(x2);
                    roots.insert(x3);
                    roots.insert(x4);
                } else {
                    cout << "current factor: " << f << endl;
                    throw Error("AMF::InitDE, higher powers found.");
                }
            }
            lst rs;
            for(auto ri : roots) rs.append(ri);
            sort_lst(rs);
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Total DE Poles: " << rs.nops() << endl;
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
                    if(x0/I<0) x0 = min/2;;
                    if(abs(x0)<min) break;
                }
            } else throw Error("AMF::InitDE, NO root found.");
        }
        if(!In_GiNaC_Parallel && Verbose>10) {
            cout << "  \\--Total AMF Points: " << pts.nops() << endl;
            cout << "  \\--DE Poles finished @ " << now() << endl;
        }
    }
    
    lst AMF::Evaluate(const ex & d0, int xn, int dp) {
        set_precision(dp);
        int N = Mat.rows();
        auto nmat = ex_to<matrix>(subs(Mat,lst{d==d0}));
        int nloop = ibp.Internal.nops();
        int npts = pts.nops();
        
        //--------------------------------------------------------------------------------------
        // DE at origin
        //--------------------------------------------------------------------------------------
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ origin ..." << RESET << endl;
        matrix oTU, ioTU;
        if(true) {
            DEX de(x);
            de.init(nmat);
            de.fuchsify(); // fuchsify & shear
            matrix T = de.T();
            
            // take x->0 limit at origin
            auto pr = prank(T,x);
            pr++;
            if(pr<0) pr = 0;
            for(auto ev : de.las) {
                if(!ev.info(info_flags::integer)) continue;
                pr -= ex_to<numeric>(ev).to_int();
                break;
            }
            auto CMat = de.series(pr,lst{0});
            auto C00 = CMat[0][0]; // also drop ln^k x
            matrix U(N,N);
            for(int n=0; n<C00.size(); n++) U = U.add(C00[n].mul_scalar(pow(x,n)));
            oTU = T.mul(U);
            xpow(oTU,x);
//            subs(oTU,x==0,nopat);

            ex x0 = pts.op(0);
            matrix iT = ex_to<matrix>(subs(T,x==x0)).inverse();
            U = de.series(xn,dp,x0);
            ioTU = U.inverse().mul(iT);
        }
        
        //--------------------------------------------------------------------------------------
        // DE at infinity
        //--------------------------------------------------------------------------------------
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ infinity ..." << RESET << endl;
        matrix fBC;
        matrix iBC;
        if(true) {
            DEX de(x);
            auto mat = nmat;
            mat = x2y(mat,1/x,x); // infinity to origin
            //xpow(mat,x);
            de.init(mat);
            
            ex x0 = 1/pts.op(npts-1);
            matrix T;              
            matrix C(N,1); // F in ooC should be treated as Bubble type
            if(true) {
                auto CMat = de.series(0);
                T = de.T();
                for(int i=0; i<N; i++) C(i,0) = iWF(i);
                
                ex ila; // select the lambda of type n-L*d/2 
                bool first = true;
                for(const auto & la : de.las) {
                    if(!normal(la+nloop*d0/2).info(info_flags::integer)) continue;
                    if(!first) throw Error("!first: something may be wrong here.");
                    ila = la;
                    first = false;
                }
                if(first) throw Error("first: something may be wrong here.");
                
                for(const auto & kv : CMat) {
                    ex la = kv.first;
                    if(!is_zero(la-ila)) {
                        matrix C00 = kv.second[0][0]; // k=0, n=0;
                        for(int i=0; i<N; i++) {
                            if(!is_zero(C00(i,i))) C(i,0) = 0;
                        }
                    } 
                }
                            
                lst dlst;
                int dim = MIntegrals.op(0).op(1).nops();
                for(int i=0; i<N; i++) {
                    ex v = -nloop*d0/2-ila;
                    ex ns = MIntegrals.op(i).op(1);
                    for(int j=0; j<dim; j++) v += ns.op(j);
                    dlst.append(pow(x,v));
                }
                auto mat = ex_to<matrix>(diag_matrix(dlst)).inverse().mul(T);
                auto pr = prank(mat,x);
                pr++; 
                if(pr<0) pr = 0;
                if(pr>xn) throw Error("pr>xn");
                CMat = de.series(pr, lst{ila});
                vector<matrix> fmat_vec;
                auto ks = CMat[ila].size();
                for(int k=0; k<ks; k++) {
                    matrix cmat = CMat[ila][k][0];
                    for(int n=1; n<=pr; n++) cmat = cmat.add(CMat[ila][k][n].mul_scalar(pow(x,n)));
                    mat = mat.mul(cmat);
                    int ldeg = 1;
                    for(int i=0; i<mat.nops(); i++) {
                        mat.let_op(i) = series_ex(mat.op(i),x,0);
                        if(ldeg>mat.op(i).ldegree(x)) ldeg = mat.op(i).ldegree(x);
                    }
                    for(int l=ldeg; l<=0; l++) {
                        matrix fmat(N, N+1);
                        for(int r=0; r<N; r++) {
                            for(int c=0; c<N; c++) fmat(r,c) = mat(r,c).coeff(x,l);
                            if(k==0 && l==0) fmat(r,N) = _MIntegrals.op(r).subs(d==d0);
                            else fmat(r,N) = 0;
                        }
                        fmat_vec.push_back(fmat);
                    }
                }
                    
                matrix fmat(N*fmat_vec.size(),N+1); 
                for(int n=0; n<fmat_vec.size(); n++) {
                    for(int r=0; r<N; r++) for(int c=0; c<=N; c++) fmat(n*N+r,c) = fmat_vec[n](r,c);
                }
                fmat = fermat_Redrowech(fmat);
cout << endl << fmat << endl << endl;
                lst feqns;
                for(int r=0; r<N; r++) {
                    int idx = -1;
                    for(int c=0; c<N; c++) {
                        if(!is_zero(fmat(r,c)) && !is_zero(fmat(r,c)-1)) {
                            cout << "fmat(r,c)=" << fmat(r,c) << endl;
                            throw Error("fmat(r,c) not 0/1: something may be wrong here.");
                        }
                        if(is_zero(fmat(r,c)-1)) {
                            if(idx!=-1) throw Error("idx!=-1: something may be wrong here.");
                            idx = c;
                        }
                    }
                    if(idx!=-1) C(idx,0) = fmat(r,N);
                    else {
                        feqns.append(fmat(r,N));
                        if(false && !is_zero(fmat(r,N))) { // TODO: check again
                            cout << feqns << endl;
                            throw Error("fmat(r,N) is NOT zero."); 
                        }
                    }
                }

                if(feqns.nops()>0) { // false to disable
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
                    for(int r=0; r<N; r++) C(r,0) = C(r,0).subs(fsol);
                }
                
                if(C.has(iWF(w))) {
                    cout << fmat << endl;
                    cout << C << endl;
                    throw Error("ooC is not determined completely.");
                }
            }
            
            exset fset;
            find(C,F(w1,w2),fset);
            if(fset.size()<1) throw Error("only zero solution?!");
            exvector fs(fset.begin(),fset.end());
            int nc = fs.size();
            fBC = matrix(nc,1);
            iBC = matrix(N,nc);
            for(int c=0; c<nc; c++) {
                fBC(c,0) = fs[c];
                for(int r=0; r<N; r++) iBC(r,c) = C(r,0).coeff(fs[c]);                
            }
            iBC = de.series(xn,iBC,dp,x0); // TODO: add x0
            T = ex_to<matrix>(subs(T,x==x0));
            iBC = T.mul(iBC);
        }
        
        //--------------------------------------------------------------------------------------
        // from infinity to origin
        //--------------------------------------------------------------------------------------
        if(!In_GiNaC_Parallel && Verbose>0) {
            cout << "  \\--" << WHITE << "AMF @ middle " << npts-1 << " points ... " << RESET << endl;
        }
        if(true) {
            DEX de(x);
            de.init(nmat);
            auto nc = iBC.cols();
            for(int i=npts-2; i>=0; i--) {
                auto x1 = pts.op(i);
                auto x2 = pts.op(i+1);
                auto dx = x1-x2;
                iBC = de.taylor(xn,iBC,x2,dp,dx);
            }
        }

        // Final C at origin
        matrix oC = ioTU.mul(iBC).mul(fBC);
cout << oC << endl;
        lst res; // result for master integrals
        oC = oTU.mul(oC);
        for(int i=0; i<N; i++) res.append(oC(i,0));
        reset_precision();
cout << NN(res,20) << endl;
        return res;
    }
    
//    lst AMF::FitEps(const lst & eps, int lp, bool parallel) {
//        if(WDigits>0) set_precision(WDigits);
//        exvector eps_vec(eps.begin(), eps.end());
//        int nmi = MIntegrals.nops();
//        exvector mis_vec[nmi];
//        if(parallel && eps.nops()>1) {
//            auto od0 = d0;
//            auto res_vec = GiNaC_Parallel(eps.nops(), [&](int idx)->ex {
//                d0 = 4-2*eps.op(idx);
//                return Evaluate();
//            }, "AMF");
//            for(auto mis : res_vec) {
//                for(int i=0; i<nmi; i++) mis_vec[i].push_back(mis.op(i));
//            }
//            d0 = od0;
//        } else {
//            auto od0 = d0;
//            for(auto epi : eps) {
//                d0 = 4-2*epi;
//                if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "ep = " << epi << RESET << endl;
//                auto mis = Evaluate();
//                for(int i=0; i<nmi; i++) mis_vec[i].push_back(mis.op(i));
//            }
//            d0 = od0;
//        }
//        if(WDigits>0) reset_precision();
//        set_precision(100*WDigits);
//        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--Final PolynomialFit ..." << endl;
//        lst mis_lst;
//        for(int i=0; i<nmi; i++) {
//            int tn = eps.nops()-3;
//            if(tn<1) tn = 1;
//            auto cs = PolynomialFit(eps_vec, mis_vec[i], tn, lp);
//            ex mi = 0;
//            for(int i=0; i<tn; i++) mi += pow(ep,lp+i) * cs.op(i);
//            mis_lst.append(mi);
//        }
//        reset_precision();
//        return mis_lst;
//    }
    
//    lst AMF::FitEps(int goal, int order, bool parallel) { // form AMFlow
//        if(WDigits>0) set_precision(WDigits);
//        int nloop = ibp.Internal.nops();
//        ex nn = 5*order/numeric(2)+2*nloop;
//        int n = cln_ceiling(nn);
//        if(n>100) throw Error("FitEps: too large order.");
//        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ " << n << " eps points ..." << RESET << endl;
//        
//        nn = nloop/numeric(2)+goal/numeric(order+1);
//        auto tn = cln_ceiling(nn);        
//        ex eps0 = GiNaC::pow(10,-tn);
//        
//        int lp = -2*nloop;
//        lst eps;
//        for(int i=0; i<n; i++) eps.append(eps0 + eps0*i/ex(100));
//        auto sp = (n+2*nloop)*tn;
//        if(sp<30) sp = 30;
//        if(WDigits<2*sp) WDigits = 2*sp;
//        if(dp<WDigits) dp = WDigits;
//        xN = 4*sp;
//        auto res = FitEps(eps,lp,parallel);
//        if(WDigits>0) reset_precision();
//        return res;
//    }  
        
    //=*********************************************************************=
    
    ex AMF::Vacuum(int nl, int np) {
        set_precision(100);
        ex res;
        static ex J = tgamma(-1+ep);
        if(nl==1 && np==1) res = 1;
        else if(nl==2 && np==2) res = 1;
        else if(nl==2 && np==3) {
            res = str2ex("-1.500000000000000000000000000000000000000000000000 -1.500000000000000000000000000000000000000000000000*ep +0.515860858034188335902343433308415603643104514453*ep^2 -8.540503339614544671799894997792116772367413851777*ep^3 +1.039200541451345629937997428402437565814452745044*ep^4 -34.02412109418437876206777042875976448646874597234*ep^5");
        } else if(nl==3 && np==4) {
            res = str2ex("-2.000000000000000000000000000000000000000000000000 -1.666666666666666666666666666666666666666666666666*ep -0.499999999999999999999999999999999999999999999999*ep^2 +8.583333333333333333333333333333333333333333333333*ep^3 +2.664875615375146678409775303572533678107178419288*ep^4 +196.7353782591730433858563053732030434664159925326*ep^5");
        } else if(nl==3 && np==5) {
            res = str2ex("1.000000000000000000000000000000000000000000000000 +2.666666666666666666666666666666666666666666666666*ep +1.301611617264956661528646466716502126047124304425*ep^2 +16.17027687753648029273000749323203572619203425393*ep^3 +50.36368751002766464712459970887058692375220280352*ep^4 +72.00897461295336034290529765044698186018294035915*ep^5");
        } else if(nl==3 && np==6) {
            res = str2ex("-2.404113806319188570799476323022899981529972584680*ep^2 +17.24761989872635488431312965422760018324025125004*ep^3 -73.26296589040362104788617737106101541072605775907*ep^4 +259.4946671222559246930353806588203939311375233114*ep^5 -855.0640324263683182684972461824631640925159683337*ep^6 +2715.946776452544387893443991909756653155929639372*ep^7");
        } else throw Error("Not Supported Yet.");
        
        
        res *= pow(J,nl) * exp(-I*Pi*(2-ep)*nl);
        res = series_ex(res, ep, 5-nl).evalf();
        reset_precision();
        return res;
    }
    
    //=*********************************************************************=
}

