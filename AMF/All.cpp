/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "HEP.h"
#include "AMF.h"
#include "flint/fmpz_mod.h"
#include "flint/fmpz_mod_mat.h"
#include <cmath>

namespace HepLib {
    
    void All::InitDE() {
        MIntegral.remove_all();
        xMIntegral.remove_all();
        
        // from now on, we use iPropagator, the internal propagator
        iPropagator = Propagator;
        for(int i=0; i<iPropagator.nops(); i++) {
            auto pi = iPropagator.op(i).expand();
            bool sq = false;
            for(auto lp : Internal) {
                if(pi.has(lp)) {
                    sq = !pi.coeff(lp,2).is_zero();
                    break;
                }
            }
            if(!sq) {
                // make sure the i-th index as numerators
                for(auto fi : Integral) {
                    if(fi.op(i)>0) {
                        cout << "Integral: " << Integral << endl;
                        throw Error("the index is positive for Non-squared propagator!");
                    }
                }
                if(pi.match(w1*w2) && pi.nops()==2) {
                    iPropagator.let_op(i) = pow(pi.op(0), 2);
                    if(!isComplete(iPropagator, Internal, External)) {
                        iPropagator.let_op(i) = pow(pi.op(1), 2);
                        if(!isComplete(iPropagator, Internal, External)) {
                            iPropagator.let_op(i) = pow(pi.op(0)+pi.op(1), 2);
                            if(!isComplete(iPropagator, Internal, External)) throw Error("isComplete failed, somthing may be wrong here.");
                        }
                    }
                } else {
                    cout << "pi: " << pi << endl;
                    throw Error("Not yet supported Non-squared propagator.");
                }
            }
        }
        
        map<int,ex> pi_by_prop;
        for(int i=0; i<iPropagator.nops(); i++) {
            if(!iPropagator.op(i).is_equal(Propagator.op(i))) {
                auto ret = express_by(Propagator.op(i), iPropagator, Internal, External);
                ret = ex_to<lst>(subs(ret,Replacement));
                ex res = 0;
                for(int j=0; j<iPropagator.nops(); j++) res += ret.op(j) * iWF(j);
                res += ret.op(iPropagator.nops());
                pi_by_prop[i] = res;
            }
        }

        lst ns0, _Integral_;
        for(int i=0; i<iPropagator.nops(); i++) ns0.append(0);
        for(int i=0; i<Integral.nops(); i++) {
            ex ei = F(0, ns0);
            ex ni = Integral.op(i);
            for(int j=0; j<iPropagator.nops(); j++) {
                ei = MapFunction([&](const ex & e, MapFunction &self)->ex {
                    if(!e.has(F(w1,w2))) return e;
                    else if(e.match(F(w1,w2))) {
                        auto fns = e.op(1);
                        auto itr = pi_by_prop.find(j);
                        if(itr == pi_by_prop.end()) {
                            fns.let_op(j) = fns.op(j) + ni.op(j);
                            return F(e.op(0), fns);
                        } else {
                            if(ni.op(j)>0) throw Error("ni.op(j)>0.");
                            ex res = pow(itr->second, -ni.op(j));
                            auto cvs = collect_lst(res, iWF(w));
                            res = 0;
                            for(auto cv : cvs) {
                                lst fns2 = ex_to<lst>(fns);
                                int deg, pos;
                                if(!cv.op(1).has(iWF(w))) {
                                    deg = 0;
                                    pos = -1;
                                } else if(cv.op(1).match(iWF(w))) {
                                    deg = 1;
                                    pos = ex2int(cv.op(1).op(0));
                                } else {
                                    deg = ex2int(cv.op(1).op(1));
                                    pos = ex2int(cv.op(1).op(0).op(0));
                                }
                                if(pos>=0) fns2.let_op(pos) = fns2.op(pos) - deg;
                                res += cv.op(0) * F(e.op(0), fns2);
                            }
                            return res;
                        }
                    } else return e.map(self);
                })(ei);
            }
            iIntegral.append(ei);
            exset fs;
            find(ei, F(w1,w2), fs);
            for(auto fi : fs) _Integral_.append(fi.op(1));
        }
        _Integral_.sort();
        _Integral_.unique();
        
        // TODO: reorder, put numerator to the right
            
        auto pre_o = PRE;
        PRE = pre;
        FIRE fire;
        fire.T1 = T1;
        fire.T2 = T2;
        fire.LT1 = LT1;
        fire.LT2 = LT2;
        fire.PosPref = PosPref;
        fire.Internal = Internal;
        fire.External = External;
        fire.NVariables[d] = d0;
        fire.Integral = _Integral_;
        fire.Propagator = iPropagator;
        
        auto rr = Replacement;
        for(auto r : Replacement) rr.append(w*r.op(0) == w*r.op(1));
        rr.sort();
        rr.unique();
        fire.Replacement.remove_all();
        for(auto ri : rr) {
            auto ri1 = ri.op(1).subs(rr);
            while(true) {
                auto ri2 = ri1.subs(rr);
                if(ri2==ri1) break;
                ri1 = ri2;
            }
            fire.Replacement.append(ri.op(0)==ri1);
        }
        Replacement = fire.Replacement;
    
        if(!In_GiNaC_Parallel && Verbose>1) cout << pre << "\\--Generating DE @ " << now() << endl;
        
        fire.RM(false); // make sure delete all files
        fire.Reduce();
        fire.RM(false); // since Propagator will change
        system(("rm -rf "+fire.WorkingDir).c_str());
        
        if(fire.MIntegral.nops()<1) throw Error("No MI Found, maybe No need DE!");
        fire.FindRules(true);
        Rules = fire.Rules;
        
        MIntegral = fire.MIntegral;
        auto nvar_o = fire.NVariables;
        auto exe_o = fire.Execute;
        fire.PIntegral.remove_all();
        for(auto mi : MIntegral) fire.PIntegral.append(mi.op(1));
        
        for(int i=0; i<iPropagator.nops(); i++) {
            bool any = false;
            for(auto mi : MIntegral) {
                if(mi.op(1).op(i)>0) {
                    any = true;
                    break;
                }
            }
            if(any) TopTopo.append(1);
            else TopTopo.append(0);
        }
                
        fire.NVariables[x] = 1/ex(19);
        fire.Execute = InstallPrefix + "/FIRE/P/FIRE";
        
        // add ieta to all propagators
        for(int i=0; i<iPropagator.nops(); i++) {
            fire.Propagator.let_op(i) = iPropagator.op(i) + x;
        }
        auto cMIntegral = MIntegral;
        while(true) {
            fire.Rules.remove_all();
            fire.MIntegral.remove_all();
            fire.Integral.remove_all();
            for(auto mi : cMIntegral) {
                auto ns0 = mi.op(1);
                for(int j=0; j<ns0.nops(); j++) {
                    if(is_zero(ns0.op(j))) continue;
                    auto ns = ns0;
                    ns.let_op(j) = ns.op(j)+1;
                    fire.Integral.append(ns);
                }
            }
            for(auto mi : cMIntegral) fire.Integral.append(mi.op(1));
            for(auto mi : MIntegral) fire.Integral.append(mi.op(1));
            fire.Integral.sort();
            fire.Integral.unique();
            fire.Reduce();
            fire.RM(true); // keep .start & .config
            if(using_FR) fire.FindRules(true);
            sort_lst(fire.MIntegral);

            if(fire.MIntegral==cMIntegral) {
                fire.RM(false);
                xMIntegral = cMIntegral;
                break;
            }
            cMIntegral = fire.MIntegral;
        }
        
        // check TopSector
        for(int i=0; i<iPropagator.nops(); i++) {
            bool any = false;
            for(auto mi : xMIntegral) {
                if(mi.op(1).op(i)>0) {
                    any = true;
                    break;
                }
            }
            if(any && TopTopo.op(i)==0) cout << RED << "Warn: TopSector Changed!" << RESET << endl;
        }
                
        if(true) { // Final Reducation -> DE
            fire.NVariables = nvar_o;
            fire.Execute = exe_o;
            fire.Rules.remove_all();
            fire.MIntegral.remove_all();
            fire.Integral.remove_all();
            lst dmis;
            for(auto mi : xMIntegral) {
                ex dmi = 0;
                auto ns0 = mi.op(1);
                for(int i=0; i<ns0.nops(); i++) {
                    if(is_zero(ns0.op(i))) continue;
                    auto ns = ns0;
                    ns.let_op(i) = ns.op(i)+1;
                    fire.Integral.append(ns);
                    dmi -= ns0.op(i)*F(mi.op(0),ns);
                }
                dmis.append(dmi);
            }
            for(auto mi : xMIntegral) fire.Integral.append(mi.op(1));
            for(auto mi : MIntegral) fire.Integral.append(mi.op(1));
            fire.Integral.sort();
            fire.Integral.unique();
            fire.Reduce();
            if(using_FR) fire.FindRules(true);
            sort_lst(fire.MIntegral);
            if(xMIntegral != fire.MIntegral) {
                cout << "xMIntegral size: " << xMIntegral.nops() << endl;
                cout << "ibp.MIntegral size: " << fire.MIntegral.nops() << endl;
                throw Error("xMIntegral != ibp.MIntegral");
            }
            MIxMI = ex_to<lst>(subs(MIntegral,fire.Rules));
            if(MIxMI.has(x)) {
                cout << endl << MIxMI << endl;
                throw Error("MIxMI still has iet.");
            }
            int n = xMIntegral.nops();
            matrix mat(n,n);
            for(int r=0; r<n; r++) {
                auto dmi = dmis.op(r).subs(fire.Rules,nopat);
                dmi = collect_ex(dmi, F(w1,w2));
                ex chk = 0;
                for(int c=0; c<n; c++) {
                    mat(r,c) = dmi.coeff(xMIntegral.op(c));
                    chk += mat(r,c) * xMIntegral.op(c);
                }
                if(!is_zero(normal(chk-dmi))) throw Error("AMFlow::InitDE, Check failed");
            }
            Mat = mat;
            system(("rm -rf "+fire.WorkingDir).c_str());
        }
        
        if(!In_GiNaC_Parallel && Verbose>1) cout << pre << "\\--Generated DE @ " << now() << endl;
        PRE = pre_o;
        
        if(true) { // BC
            if(!In_GiNaC_Parallel && Verbose>1) cout << pre << "\\--Generating BC @ " << now() << endl;
            exmap lmap;
            static symbol y("y");
            for(auto lp : Internal) lmap[lp] = lp*y;
            lst props;
            for(int i=0; i<iPropagator.nops(); i++) {
                auto e = iPropagator.op(i).subs(lmap).expand().coeff(y,2)+1; // +1 refers to +x
                props.append(e);
            }
            
            FIRE bc;
            bc.T1 = T1;
            bc.T2 = T2;
            bc.LT1 = LT1;
            bc.LT2 = LT2;
            bc.NVariables = nvar_o;
            bc.ProblemNumber = fire.ProblemNumber;
            bc.Propagator = props;
            bc.Propagator.sort().unique();
            for(auto mi : xMIntegral) {
                lst intg;
                exmap p2n;
                for(int i=0; i<props.nops(); i++) p2n[props.op(i)] += mi.op(1).op(i);
                for(int i=0; i<bc.Propagator.nops(); i++) intg.append(p2n[bc.Propagator.op(i)]);
                bc.Integral.append(intg);
                ooMIntegral.append(F(bc.ProblemNumber, intg));
            }
            bc.Internal = Internal;
            bc.Reduce();
            bc.FindRules(true);
            for(int i=0; i<ooMIntegral.nops(); i++) ooMIntegral.let_op(i) = ooMIntegral.op(i).subs(bc.Rules).subs(F(w1,w2)==F(bc.Propagator, w2));
            system(("rm -rf "+bc.WorkingDir).c_str());
            if(!In_GiNaC_Parallel && Verbose>1) cout << pre << "\\--Generated BC @ " << now() << endl;
        }
    }
    
    void All::Poles() {
        auto pre_o = PRE;
        PRE = pre;
        if(!In_GiNaC_Parallel && Verbose>10) cout << pre << "\\--Poles(1/" << WHITE << rr << RESET << ") starting @ " << now() << endl;

        int ndig = 30;
        lst irs;
        ex den = matrix_den_lcm(Mat);
        ex pex = 1;
        if(!is_a<mul>(den)) den = lst{den};
        for(auto di : den) {
            start:
            if(!di.has(x)) continue;
            else if(di.match(pow(x,w)) || di.is_equal(x)) continue; // skip zero root
            else if(di.match(pow(w1,w2))) {
                di = di.op(0);
                goto start;
            } else if(di.degree(x)==2 && di.coeff(x).is_zero()) {
                ex cc = di.coeff(x,0) / di.coeff(x,2);
                if(cc>0) {
                    irs.append(NN(sqrt(cc), ndig));
                    continue;
                }
            } else if(di.degree(x)==4 && di.coeff(x,3).is_zero() && di.coeff(x,3).is_zero() && di.coeff(x).is_zero()) {
                ex cc = di.coeff(x,0) / di.coeff(x,4);
                if(cc<0) {
                    irs.append(NN(pow(-cc,1/ex(4)), ndig));
                    continue;
                }
            }
            pex *= di;
        }
        lst rs = poly_roots(pex, 200); // 2^100 ~ 10^30
        sort_lst(rs);
        sort_lst(irs, false); // large first
        lst hd_irs;
        for(auto item : irs) {
            ex min = item;
            for(auto r : rs) {
                ex ar = abs(r-I*item); // note I here
                if(min<0 || ar<min) min = ar;
            }
            for(auto r : irs) {
                if(r==item) continue;
                ex ar = abs(r-item); // no I here
                if(ar<min) min = ar;
            }
            hd_irs.append(min/2);
        }
        if(irs.nops()>0 && !In_GiNaC_Parallel) cout << pre << RED << "Note: imaginary poles found, using contour: ie->ie-e'." << RESET << endl;
        
        pts.remove_all();
        if(rs.nops()>0) {
            ex max = -1, min = -1;
            for(auto r : rs) {
                ex ar = abs(r);
                if(max<0 || ar>max) max = ar;
                if(min<0 || ar<min) min = ar;
            }
            if(irs.nops()>0 && min>irs.op(irs.nops()-1)) min = irs.op(irs.nops()-1);
            min /= rr;
            max *= rr;
            
            double beta = 2*std::asin(1/(2.0L*rr));
            int nn = std::floor(3.14159265358979323846L/beta);
            int sign = -1;
            
            ex x0;
            int cur_pos = 0;
            // the 1st point near infinity
            if(irs.nops()>cur_pos && irs.op(cur_pos)+hd_irs.op(cur_pos)>max) {
                for(int i=0; i<=nn; i++) {
                    x0 = I*irs.op(cur_pos)+hd_irs.op(cur_pos)*(I*cos(i*beta)+sign*sin(i*beta));
                    x0 = Rationalize(x0, ndig);
                    pts.prepend(x0);
                }
                x0 = I*(irs.op(cur_pos)-hd_irs.op(cur_pos));
                x0 = Rationalize(x0, ndig);
                pts.prepend(x0);
                cur_pos++;
            } else {
                x0 = I*max; // the last point
                x0 = Rationalize(x0, ndig);
                pts.prepend(x0);
            }
            
            // other points from infinity to origin
            while(true) {
                ex mm = abs(x0); // distance between x0 and origin
                for(auto r : rs) {
                    ex ar = abs(r-x0);
                    if(mm<0 || ar<mm) mm = ar;
                }
                for(auto r : irs) {
                    ex ar = abs(I*r-x0); // note I here
                    if(mm<0 || ar<mm) mm = ar;
                }
                if(pts.nops()>1000) {
                    cout << "pex=" << pex << endl;
                    throw Error("maybe, pole is pure imaginary!");
                }
                x0 -= I*mm/rr;

                if(irs.nops()>cur_pos && irs.op(cur_pos)+hd_irs.op(cur_pos)>(x0/I)) {
                    for(ex i=0; i<=nn; i++) {
                        x0 = I*irs.op(cur_pos)+hd_irs.op(cur_pos)*(I*cos(i*beta)+sign*sin(i*beta));
                        x0 = Rationalize(x0, ndig);
                        pts.prepend(x0);
                    }
                    x0 = I*(irs.op(cur_pos)-hd_irs.op(cur_pos));
                    x0 = Rationalize(x0, ndig);
                    pts.prepend(x0);
                    cur_pos++;
                } else {
                    if(x0/I<0) x0 = I*9*min/10; // x0 -> 0.9min
                    x0 = Rationalize(x0, ndig);
                    pts.prepend(x0);
                }
                if(abs(x0)<min) break;
            }
        } else throw Error("All::InitDE, NO root found.");
    
        if(!In_GiNaC_Parallel && Verbose>10) {
            cout << pre << "\\--Total AMF Points: " << WHITE << pts.nops() << RESET << endl;
            cout << pre << "\\--Poles(1/" << WHITE << rr << RESET << ") finished @ " << now() << endl;
        }
        PRE = pre_o;
    }
    
    matrix All::oo() {
        if(!In_GiNaC_Parallel && Verbose>0) cout << pre << "\\--" << WHITE << "All Mode @ infinity ..." << RESET << endl;

        int N = Mat.rows();
        int nloop = Internal.nops();
        ex x0 = 1/pts.op(pts.nops()-1);
    
        DEX de(x, pre);
        auto mat = Mat;
        mat = x2y(mat,1/x,x); // infinity to origin
        //xpow(mat,x); // seems no need, at least for DEX
        
        de.init(mat);
        auto CMat = de.series(0);
        auto T = de.T();
        
        ex ila; // select the lambda of type n-L*d/2
        bool first = true;
        for(const auto & la : de.las) {
            if(!normal(la+nloop*d0/2).info(info_flags::integer)) continue;
            if(!first) {
                cout << endl << "las=" << de.las << endl;
                throw Error("!first: something may be wrong here.");
            }
            ila = la;
            first = false;
        }
        if(first) throw Error("first: something may be wrong here.");
        
        list<matrix> fmat_vec; // boundar equations
        list<lst> fcol_vec; // the last column
        
        lst zero_cl;
        for(int i=0; i<N; i++) zero_cl.append(0);
        
        for(const auto & kv : CMat) {
            ex la = kv.first;
            if(!is_zero(la-ila)) {
                int kmax = kv.second.size();
                for(int k=0; k<kmax; k++) {
                    auto cmat = kv.second[k][0];
                    matrix fmat(N, N);
                    for(int r=0; r<N; r++) for(int c=0; c<N; c++) fmat(r,c) = cmat(r,c);
                    fmat_vec.push_back(fmat);
                    fcol_vec.push_back(lst{zero_cl , 0}); // 0 will not used
                }
            }
        }
                    
        lst dlst;
        int dim = iPropagator.nops();
        for(int i=0; i<N; i++) {
            ex v = -nloop*d0/2-ila;
            ex ns = xMIntegral.op(i).op(1);
            for(int j=0; j<dim; j++) v += ns.op(j);
            dlst.append(pow(x,v));
        }
        auto t_mat = ex_to<matrix>(diag_matrix(dlst)).inverse(solve_algo::gauss).mul(T);
        auto pr = prank(t_mat,x);
        pr++;
        if(pr<0) pr = 0;
        CMat = de.series(pr, lst{ila});
        
        auto ks = CMat[ila].size();
        for(int k=0; k<ks; k++) {
            matrix cmat = CMat[ila][k][0];
            for(int n=1; n<=pr; n++) cmat = cmat.add(CMat[ila][k][n].mul_scalar(pow(x,n)));
            auto mat = t_mat.mul(cmat);
            int ldeg = 0;
            for(int i=0; i<mat.nops(); i++) {
                mat.let_op(i) = series_ex(mat.op(i),x,0);
                if(ldeg>mat.op(i).ldegree(x)) ldeg = mat.op(i).ldegree(x);
            }
            for(int l=ldeg; l<=0; l++) {
                matrix fmat(N, N);
                for(int r=0; r<N; r++) {
                    for(int c=0; c<N; c++) fmat(r,c) = mat(r,c).coeff(x,l);
                }
                if(l>=0 && k==0) {
                    fmat_vec.push_back(fmat);
                    fcol_vec.push_back(lst{ooMIntegral, 0}); // 0 will not used
                } else {
                    fmat_vec.push_front(fmat);
                    fcol_vec.push_front(lst{zero_cl, 0}); // 0 will not used
                }
            }
        }
        
        matrix fmat(N,N*fmat_vec.size()); // note transpose
        auto itr = fmat_vec.begin();
        for(int n=0; n<fmat_vec.size(); n++,itr++) {
            for(int r=0; r<N; r++) for(int c=0; c<N; c++) fmat(c,n*N+r) = (*itr)(r,c);
        }

        //auto rref_fmat = fermat_Redrowech_Sparse(fmat);
        fmpq_mat_t qmat;
        fmpq_mat_init(qmat, fmat.rows(), fmat.cols());
        _to_(qmat, fmat);
        fmpz_mat_t zmat;
        fmpz_t z;
        fmpz_mat_init(zmat, fmat.rows(), fmat.cols());
        fmpz_init(z);
        fmpq_mat_get_fmpz_mat_matwise(zmat, z, qmat);
        fmpq_mat_clear(qmat);
        
        //auto rank = fmpq_mat_rref(qmat, qmat);
        fmpz_set_ui(z, 99999999u);
        fmpz_nextprime(z, z, 1);
        
        fmpz_mod_ctx_t mod_ctx;
        fmpz_mod_ctx_init(mod_ctx, z);
        fmpz_mod_mat_t zmat_mod;
        fmpz_mod_mat_init(zmat_mod, fmat.rows(), fmat.cols(), mod_ctx);
        fmpz_mod_mat_set_fmpz_mat(zmat_mod, zmat, mod_ctx);
        auto rank = fmpz_mod_mat_rref(zmat_mod, zmat_mod, mod_ctx);
        fmpz_mod_mat_clear(zmat_mod, mod_ctx);
        fmpz_mod_ctx_clear(mod_ctx);
        
        if(rank!=N) {
            cout << "rank=" << rank << " != "<< "N=" << N << endl;
            throw Error("oo: BC is NOT enough!");
        }
        auto rref_fmat = _to_(zmat);
        fmpz_clear(z);
        fmpz_mat_clear(zmat);
        
        //  ========== Solve Boundary Conditions ==========
        
        matrix col_bc(N,1);
        matrix inv_mat(N, N);
        int c = 0;
        for(int r=0; r<N; r++) {
            for( ; c<N*fmat_vec.size(); c++) {
                if(!rref_fmat(r,c).is_zero()) {
                    if(rref_fmat(r,c)!=1) throw Error("pivot is NOT 1 in RREF.");
                    for(int j=0; j<N; j++) inv_mat(r, j) = fmat(j, c);
                    auto itr = fcol_vec.begin();
                    std::advance(itr, c/N);
                    col_bc(r,0) = itr->op(0).op(c%N);
                    break;
                }
            }
        }
        inv_mat = inv_mat.inverse(solve_algo::gauss); // TODO: try Fermat or Flint
        
        exset fs;
        find(col_bc, F(w1,w2), fs);
        int nc = fs.size();
        if(!using_fBC) nc = 1;
        matrix iBC(N,nc);
        if(using_fBC) {
            fBC = matrix(nc,1);
            int c = 0;
            for(auto fi : fs) {
                fBC(c,0) = fi;
                for(int r=0; r<N; r++) iBC(r,c) = col_bc(r,0).coeff(fi);
                c++;
            }
        } else {
            exmap f2f;
            for(auto fi : fs) {
                for(auto it : fi.op(1)) {
                    if(!it.is_zero() && !it.is_equal(1)) throw Error("oo: exponent is Not either 0 or 1.");
                }
                int nL = Internal.nops();
                lst mom_lst;
                for(int i=0; i<fi.op(0).nops(); i++) {
                    if(fi.op(1).op(i).is_zero()) continue;
                    auto pi = fi.op(0).op(i)-1; // p^+1-1 -> p^2
                    pi = exfactor(pi);
                    if(pi.match(pow(w,2))) {
                        mom_lst.append(pi.op(0));
                    } else {
                        cout << "pi: " << pi << endl;
                        throw Error("pi not match w^2!");
                    }
                }
                mom_lst.sort();
                mom_lst.unique();
                
                lst loops;
                exmap loops_back;
                for(int j=0; j<nL; j++) {
                    auto item = symbol("_q"+to_string(j+1));
                    loops.append(item);
                    loops_back[item] = Internal.op(j);
                }
                int g_max = 0;
                lst g_last;
                Combinations(mom_lst.nops(), nL, [&](const int* is) {
                    lst eqs;
                    for(int j=0; j<nL; j++) {
                        eqs.append(mom_lst.op(is[j])==loops.op(j));
                    }
                    auto sols = lsolve(eqs, Internal);
                    if(sols.nops()>0) {
                        auto vlp = subs(mom_lst, sols).subs(loops_back);
                        lst grp;
                        for(auto vi : vlp) grp.append(lst{vi});
                        for(auto li : Internal) {
                            lst grp2, gi;
                            for(auto item : grp) {
                                if(item.has(li)) for(auto it : item) gi.append(it);
                                else grp2.append(item);
                            }
                            grp2.prepend(gi);
                            grp = grp2;
                        }
                        if(g_max<grp.nops()) {
                            g_max = grp.nops();
                            g_last = grp;
                        }
                    }
                });
                
                ex res = 1;
                for(auto gi : g_last) {
                    int nl=0, np = gi.nops();
                    for(auto li : Internal) if(gi.has(li)) nl++;
                    res *= Vacuum(nl, np, d0);
                }
            
                f2f[fi] = res;
            }
        
            iBC = ex_to<matrix>(MapFunction([&](const ex & e, MapFunction &self)->ex{
                if(!e.has(F(w1,w2))) return e;
                else if(e.match(F(w1,w2))) {
                    auto itr = f2f.find(e);
                    if(itr==f2f.end()) throw Error("A::oo, F not found.");
                    return itr->second;
                } else return e.map(self);
            })(col_bc));
        }
        iBC = inv_mat.mul(iBC);
        
        auto fp = dp2fp(dp);
        gr_ctx_t ctx;
        gr_ctx_init_complex_float_acf(ctx, fp);
        gr_ptr z0 = gr_heap_init(ctx);
        _to_(z0, x0, ctx);
        iBC = de.series(xn, iBC, z0, ctx);
        gr_heap_clear(z0, ctx);
        gr_ctx_clear(ctx);
        
        T = ex_to<matrix>(subs(T,x==x0));
        iBC = T.mul(iBC);
        
        return iBC;
    }
    
    void All::oo2o(matrix & iBC) {
        int N = Mat.rows();
        int npts = pts.nops();
        if(!In_GiNaC_Parallel && Verbose>0) cout << pre << "\\--" << WHITE << "All Mode @ middle " << npts-1 << " points ... " << RESET << endl;

        DEX de(x, pre);
        de.init(Mat); // note T matrix - permutation to triangular block
        auto T = de.T();
        iBC = T.inverse(solve_algo::gauss).mul(iBC);

        auto fp = dp2fp(dp);
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_ctx_init_complex_float_acf(ctx, fp);
        gr_ptr z0 = gr_heap_init(ctx);
        gr_ptr dz = gr_heap_init(ctx);
        gr_mat_t imat;
        gr_mat_init(imat, iBC.rows(), iBC.cols(), ctx);
        _to_(imat, iBC, ctx);
        
        for(int i=npts-2; i>=0; i--) {
            if(!In_GiNaC_Parallel && Verbose>0 && Verbose<=5) {
                cout << "\r                                          \r" << flush;
                cout << pre << "\\--" << "All Mode @ middle points: " << npts-1 << "|" << npts-1-i << flush;
            }
            auto x1 = pts.op(i);
            auto x2 = pts.op(i+1);
            auto dx = x1-x2;
            _to_(z0, x2, ctx);
            _to_(dz, dx, ctx);
            string es = "(" + to_string(npts-2) + "/" + to_string(npts-2-i) + ")";
            de.taylor(xn+xxn, imat, z0, dz, ctx, es);
        }

        gr_heap_clear(z0, ctx);
        gr_heap_clear(dz, ctx);

        iBC = _to_(imat, ctx);
        gr_mat_clear(imat, ctx);
        gr_ctx_clear(ctx);
        if(!In_GiNaC_Parallel && Verbose>0 && Verbose<=5) cout << " @ " << now(false) << endl;
        iBC = T.mul(iBC);
    }
    
    matrix All::o() {
        if(!In_GiNaC_Parallel && Verbose>0) cout << pre << "\\--" << WHITE << "All Mode @ origin ..." << RESET << endl;
        int N = Mat.rows();
        matrix oTUi;
    
        DEX de(x, pre);
        de.init(Mat);
        de.fuchsify(); // fuchsify & shear
        matrix T = de.T();
        
        // take x->0 limit at origin
        auto pr = prank(T,x);
        pr++;
        if(pr<0) pr = 0;
        ex key = 0;
        for(auto ev : de.las) {
            if(!ev.info(info_flags::integer)) continue;
            pr -= ex2int(ev);
            key = ev;
            break;
        }
        auto CMat = de.series(pr,lst{0});
        auto itr = CMat.find(key);
        if(itr==CMat.end()) throw Error("All::o, integer key is not found.");
        auto C00 = itr->second[0]; // also drop ln^k x
        matrix U(N,N);
        for(int n=0; n<C00.size(); n++) U = U.add(C00[n].mul_scalar(pow(x,n)));
        auto TU = T.mul(U);
        //xpow(TU,x); // TODO: maybe no need for DEX ?
        for(int i=0; i<TU.nops(); i++) TU.let_op(i) = series_ex(TU.op(i),x,0);

        ex x0 = pts.op(0);
        matrix iT = ex_to<matrix>(subs(T,x==x0)).inverse(solve_algo::gauss);
        
        //----------------------------------------------
        auto fp = dp2fp(dp);
        gr_ctx_t ctx;
        gr_ctx_init_complex_float_acf(ctx, fp);
        gr_ptr z0 = gr_heap_init(ctx);
        _to_(z0, x0, ctx);
        U = de.series(xn, z0, ctx);
        gr_heap_clear(z0, ctx);
        gr_ctx_clear(ctx);
        //----------------------------------------------
        
        auto iTU = U.inverse(solve_algo::gauss).mul(iT);
        oTUi = TU.mul(iTU);
        return oTUi;
    }

    void All::operator()() {
        if(!In_GiNaC_Parallel && Verbose>0) cout << WHITE << pre << "\\--All Mode Start @ " << now() << RESET << endl;
        if(!init_de) {
            InitDE();
            Poles();
            init_de = true;
        }
        
        auto oTUi = o();
        auto iBC = oo();
        oo2o(iBC);
        
        // Final C at origin
        matrix oC = oTUi.mul(iBC);
        if(using_fBC) oC = oC.mul(fBC);
        int n = Mat.rows();
        NxMIntegral.remove_all();
        exmap x_nmi_map;
        for(int i=0; i<n; i++) {
            x_nmi_map[xMIntegral.op(i)] = oC(i,0);
            NxMIntegral.append(oC(i,0));
        }
        
        NMIntegral.remove_all();
        exmap nmi_map;
        for(int i=0; i<MIxMI.nops(); i++) {
            auto nmi = MIxMI.op(i).subs(x_nmi_map);
            NMIntegral.append(nmi);
            nmi_map[MIntegral.op(i)] = nmi;
        }
        
        NIntegral.remove_all();
        n = Integral.nops();
        for(int i=0; i<n; i++) {
            ex e = iIntegral.op(i);
            e = e.subs(Rules).subs(nmi_map);
            if(!using_fBC && e.has(F(w1,w2))) {
                cout << xMIntegral << endl;
                cout << MIntegral << endl;
                cout << MIxMI << endl;
                cout << "iPropagator: " << iPropagator << endl;
                cout << "Rules: " << Rules << endl;
                cout << "item has F: " << e << endl;
                throw Error("All Mode: F still exists!");
            }
            NIntegral.append(e);
        }
        
        if(!In_GiNaC_Parallel && Verbose>0) cout << WHITE << pre << "\\--All Mode Done @ " << now() << RESET << endl;
    }
    
    ex All::Vacuum(int nl, int np, const ex & nd) {
        ex J = tgamma(-1+ep);
        ex res = pow(J,nl) * exp(-I*Pi*(2-ep)*nl);
        if(nl==np)  return res.subs(ep==(4-nd)/2);
        
        static exmap cache;
        ex key = lst{nl,np,nd};
        auto itr = cache.find(key);
        if(itr!=cache.end()) return itr->second;
        if(nl==2 && np==3) {
            Single amf(nd);
            amf.xn = xn;
            amf.dp = dp;
            amf.rr = rr;
            amf.pre = pre + "   ";
            amf.Internal = str2lst("{ q1, q2 }");
            amf.Propagator = str2lst("{ q1^2+1, q2^2+1, (q1+q2)^2+1 }");
            amf.Integral = str2lst("{ {1, 1, 1} }");
            amf();
            return cache[key] = amf.NIntegral.op(0);
        } else if(nl==3 && np==4) {
            Single amf(nd);
            amf.xn = xn;
            amf.dp = dp;
            amf.rr = rr;
            amf.pre = pre + "   ";
            amf.Internal = str2lst("{ q1, q2, q3 }");
            amf.Propagator = str2lst("{ q1^2+1, q2^2+1, q3^2+1, (q1+q2+q3)^2+1, (q1+q2)^2, (q2+q3)^2 }");
            amf.Integral = str2lst("{ {1, 1, 1, 1, 0, 0} }");
            amf();
            return cache[key] = amf.NIntegral.op(0);
        } else if(nl==3 && np==5) {
            Single amf(nd);
            amf.xn = xn;
            amf.dp = dp;
            amf.rr = rr;
            amf.pre = pre + "   ";
            amf.Internal = str2lst("{ q1, q2, q3 }");
            amf.Propagator = str2lst("{ q1^2+1, q2^2+1, q3^2+1, (q1+q2+q3)^2+1, (q1+q2)^2+1, (q2+q3)^2 }");
            amf.Integral = str2lst("{ {1, 1, 1, 1, 1, 0} }");
            amf();
            return cache[key] = amf.NIntegral.op(0);
        } else if(nl==3 && np==6) {
            Single amf(nd);
            amf.xn = xn;
            amf.dp = dp;
            amf.rr = rr;
            amf.pre = pre + "   ";
            amf.Internal = str2lst("{ q1, q2, q3 }");
            amf.Propagator = str2lst("{ q1^2+1, q2^2+1, q3^2+1, (q1+q2)^2+1, (q2+q3)^2+1, (q1+q2+q3)^2+1 }");
            amf.Integral = str2lst("{ {1, 1, 1, 1, 1, 1} }");
            amf();
            return cache[key] = amf.NIntegral.op(0);
        } else if(nl==4 && np==5) {
            Single amf(nd);
            amf.xn = xn;
            amf.dp = dp;
            amf.rr = rr;
            amf.pre = pre + "   ";
            amf.Internal = str2lst("{ q1, q2, q3, q4 }");
            amf.Propagator = str2lst("{ q1^2+1, q2^2+1, q3^2+1, q4^2+1, (q1+q2)^2+1, (q1+q3)^2, (q1+q4)^2, (q2+q3)^2, (q2+q4)^2, (q3+q4)^2 }");
            amf.Integral = str2lst("{ {1, 1, 1, 1, 1, 0, 0, 0, 0, 0} }");
            amf();
            return cache[key] = amf.NIntegral.op(0);
        }
        
        cout << "nl=" << nl << ", np=" << np << endl;
        throw Error("All::Vacuum, Not Supported Yet.");
        
    }
    
    ex All::Vacuum(int nl, int np) {
        ex res;
        static ex J = tgamma(-1+ep);
        if(nl==1 && np==1) res = 1;
        else if(nl==2 && np==2) res = 1;
        else if(nl==2 && np==3) {
            res = str2ex("-1.500000000000000000000000000000000000000000000000 -1.500000000000000000000000000000000000000000000000*ep +0.515860858034188335902343433308415603643104514453*ep^2 -8.540503339614544671799894997792116772367413851777*ep^3 +1.039200541451345629937997428402437565814452745044*ep^4 -34.02412109418437876206777042875976448646874597234*ep^5");
        } else if(nl==3 && np==3) res = 1;
        else if(nl==3 && np==4) {
            res = str2ex("-2.000000000000000000000000000000000000000000000000 -1.666666666666666666666666666666666666666666666666*ep -0.499999999999999999999999999999999999999999999999*ep^2 +8.583333333333333333333333333333333333333333333333*ep^3 +2.664875615375146678409775303572533678107178419288*ep^4 +196.7353782591730433858563053732030434664159925326*ep^5");
        } else if(nl==3 && np==5) {
            res = str2ex("1.000000000000000000000000000000000000000000000000 +2.666666666666666666666666666666666666666666666666*ep +1.301611617264956661528646466716502126047124304425*ep^2 +16.17027687753648029273000749323203572619203425393*ep^3 +50.36368751002766464712459970887058692375220280352*ep^4 +72.00897461295336034290529765044698186018294035915*ep^5");
        } else if(nl==3 && np==6) {
            res = str2ex("-2.404113806319188570799476323022899981529972584680*ep^2 +17.24761989872635488431312965422760018324025125004*ep^3 -73.26296589040362104788617737106101541072605775907*ep^4 +259.4946671222559246930353806588203939311375233114*ep^5 -855.0640324263683182684972461824631640925159683337*ep^6 +2715.946776452544387893443991909756653155929639372*ep^7");
        } else throw Error("Not Supported Yet.");
        
        res *= pow(J,nl) * exp(-I*Pi*(2-ep)*nl);
        res = series_ex(res, ep, 10-nl).evalf();
        return res;
    }
    bool All::isVacuum(int nl, int np) { // single mass
        if(nl==1 && np==1) return true;
        else if(nl==2 && np==3) return true;
        else if(nl==3 && np==4) return true;
        else if(nl==3 && np==5) return true;
        else if(nl==3 && np==6) return true;
        else if(nl==4 && np==5) return true;
        return false;
    }
    
}

