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
    
    void Single::InitDE() {
    
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
                
        fire.NVariables[x] = 1/ex(19);
        fire.Execute = InstallPrefix + "/FIRE/P/FIRE";
        
        bool fix = (prop_idx>=0);
        if(prop_idx+1>iPropagator.nops()) {
            cout << "prop_idx = " << prop_idx << endl;
            cout << " #(iPropagator) = " << iPropagator.nops() << endl;
            throw Error("iPropagator Position should be < #(iPropagator)-1.");
        }
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
        for(int i=0; i<iPropagator.nops(); i++) {
            if(fix && i!=prop_idx) continue;
            if(TopTopo.op(i)==0) continue;
            if(!In_GiNaC_Parallel && Verbose>10) cout << pre << "\\--Try Position @ " << i << endl;
            fire.Propagator = iPropagator;
            fire.Propagator.let_op(i) = iPropagator.op(i) + x;
            auto cMIntegral = MIntegral;
            sort_lst(cMIntegral);
            while(true) {
                fire.Rules.remove_all();
                fire.MIntegral.remove_all();
                fire.Integral.remove_all();
                for(auto mi : cMIntegral) {
                    auto ns = mi.op(1);
                    if(is_zero(ns.op(i))) continue;
                    ns.let_op(i) = ns.op(i)+1;
                    fire.Integral.append(ns);
                }
                if(fire.Integral.nops()<1) break;
                for(auto mi : cMIntegral) fire.Integral.append(mi.op(1));
                for(auto mi : MIntegral) fire.Integral.append(mi.op(1));
                fire.Integral.sort();
                fire.Integral.unique();
                fire.Reduce();
                fire.RM(true); // keep .start & .config
                if(using_FR) fire.FindRules(true);
                sort_lst(fire.MIntegral);

                if(fire.MIntegral==cMIntegral) {
                    if(cMIntegral.nops()<xMIntegral.nops() || xMIntegral.nops()<1) {
                        xMIntegral = cMIntegral;
                        prop_idx = i;
                    }
                    fire.RM(false);
                    break;
                }
                
                cMIntegral = fire.MIntegral;
            }
        }
        if(prop_idx<0) throw Error("prop_idx Not found!");
        
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
        
        if(!In_GiNaC_Parallel && Verbose>5) cout << pre << WHITE << "\\--Use Position @ " << prop_idx << RESET << endl;
        if(true) { // Final Reducation -> DE
            fire.RM(false);
            fire.NVariables = nvar_o;
            fire.Execute = exe_o;
            fire.Propagator = iPropagator;
            fire.Propagator.let_op(prop_idx) = iPropagator.op(prop_idx) + x;
            fire.Rules.remove_all();
            fire.MIntegral.remove_all();
            fire.Integral.remove_all();
            lst dmis;
            for(auto mi : xMIntegral) {
                auto ns = mi.op(1);
                auto n0 = ns.op(prop_idx);
                if(!n0.is_zero()) {
                    ns.let_op(prop_idx) = n0+1;
                    fire.Integral.append(ns);
                    dmis.append(-n0*F(mi.op(0),ns));
                } else dmis.append(0);
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
                throw Error("InitDE: MIxMI still has iet.");
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
    }
    
    void Single::Poles() {
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
        } else throw Error("Single::InitDE, NO root found.");

        if(!In_GiNaC_Parallel && Verbose>10) {
            cout << pre << "\\--Total AMF Points: " << WHITE << pts.nops() << RESET << endl;
            cout << pre << "\\--Poles(1/" << WHITE << rr << RESET << ") finished @ " << now() << endl;
        }
        PRE = pre_o;

    }
    
    lst Single::apart(const ex & lps, const ex & eps, const lst & expr) {
        int nl = lps.nops();
        int ne = eps.nops();
        int nle = nl*(nl+1)/2 + nl*ne;
        
        lst sps, sqs; // sqs: squared list
        matrix mat_all(nle+2, nle);
        for(int i=0; i<nl; i++) {
            sps.append(lps.op(i)*lps.op(i));
            sqs.append(lps.op(i)*lps.op(i));
            for(int j=i+1; j<nl; j++) {
                sps.append(lps.op(i)*lps.op(j));
                sqs.append(pow(lps.op(i)+lps.op(j),2));
            }
            for(int j=0; j<ne; j++) {
                sps.append(lps.op(i)*eps.op(j));
                sqs.append(pow(lps.op(i)+eps.op(j),2));
            }
        }
        
        matrix mat_sq(nle+2, nle);
        for(int c=0; c<nle; c++) {
            int r = 0;
            auto cc = sqs.op(c).expand();
            for(int i=0; i<nl; i++) {
                mat_sq(r, c) = cc.coeff(lps.op(i), 2);
                r++;
                auto ci = cc.coeff(lps.op(i));
                for(int j=i+1; j<nl; j++) {
                    mat_sq(r, c) = ci.coeff(lps.op(j));
                    r++;
                }
                for(int j=0; j<ne; j++) {
                    mat_sq(r, c) = ci.coeff(eps.op(j));
                    r++;
                }
                cc = cc.coeff(lps.op(i), 0);
            }
            mat_sq(r, c) = cc.subs(Replacement);
            r++;
            mat_sq(r, c) = 0; // note 0 power
        }
        
        auto res = MapFunction([&](const ex & e, MapFunction & self)->ex {
            if(!e.has(F(w1,w2))) return e;
            else if(e.match(F(w1,w2))) {
                exmap ps_ns;
                int n = e.op(0).nops();
                for(int i=0; i<n; i++) ps_ns[e.op(0).op(i)] += e.op(1).op(i);
                n = ps_ns.size();
                matrix mat(nle+2, n);
                
                auto itr = ps_ns.begin();
                for(int c=0; c<n; c++, itr++) {
                    if(itr->first.is_zero() && itr->second.is_zero()) throw Error("apart: 0^0 found!");
                    int r = 0;
                    auto cc = itr->first.expand();
                    for(int i=0; i<nl; i++) {
                        mat(r, c) = cc.coeff(lps.op(i), 2);
                        r++;
                        auto ci = cc.coeff(lps.op(i));
                        for(int j=i+1; j<nl; j++) {
                            mat(r, c) = ci.coeff(lps.op(j));
                            r++;
                        }
                        for(int j=0; j<ne; j++) {
                            mat(r, c) = ci.coeff(eps.op(j));
                            r++;
                        }
                        cc = cc.coeff(lps.op(i), 0);
                    }
                    mat(r, c) = cc.subs(Replacement);
                    r++;
                    mat(r, c) = -itr->second; // F convention
                }
                auto air = Apart(mat);
                
                // Apart Complete
                air = MapFunction([&](const ex & e, MapFunction &self)->ex {
                    if(!e.has(ApartIR(w))) return e;
                    else if(e.match(ApartIR(w))) {
                        if((e.op(0)).is_equal(1)) return ApartIR(mat_sq);
                        if(!is_a<matrix>(e.op(0))) throw Error("ApartIRC: Not matrix : " + ex2str(e.op(0)));
                        auto mat = ex_to<matrix>(e.op(0));
                        int cc = mat.cols();
                        if(cc==nle) return ApartIR(mat);
                        else {
                            matrix cmat(nle+2,nle);
                            for(int r=0; r<nle+2; r++) {
                                for(int c=0; c<cc; c++) cmat(r,c) = mat(r,c);
                            }
                            for(int i=0; i<nle; i++) {
                                for(int j=0; j<nle+2; j++) cmat(j, cc) = mat_sq(j, i);
                                auto r = ex_to<matrix>(sub_matrix(cmat,0,nle,0,cc+1)).rank();
                                if(r==nle) break;
                                else if(r==cc+1) cc++;
                            }
                            if(cmat.rank()!=nle) throw Error("ApartIRC failed, NOT full rank.");
                            return ApartIR(cmat);
                        }
                    } else return e.map(self);
                })(air);
                
                // ApartIR to F
                air = MapFunction([sps](const ex & e, MapFunction &self)->ex{
                    if(e.match(ApartIR(1))) return 1;
                    else if(!e.has(ApartIR(w))) return e;
                    else if(e.match(ApartIR(w))) {
                        matrix mat = ex_to<matrix>(e.op(0));
                        lst ps, ns;
                        auto nr = mat.rows();
                        auto nc = mat.cols();
                        for(int c=0; c<nc; c++) {
                            ex sum=0;
                            for(int r=0; r<nr-2; r++) sum += mat(r,c) * sps.op(r);
                            sum += mat(nr-2,c);
                            ps.append(sum);
                            ns.append( -mat(nr-1,c));
                        }
                        return F(ps, ns);
                    } else return e.map(self);
                })(air);
                return air;
            } else return e.map(self);
        })(expr);
        return ex_to<lst>(res);
    }
    
    // ns can be {n1,n2,...} or { {n1,n2,..}, ...}
    ex Single::region(const lst & ls, const lst & es, const lst & ps, const lst & ns) {
        if(ns.nops()<1) return 0;
        auto pre_o = PRE;
        PRE = pre;
        lst fs;
        bool is_lst = is_a<lst>(ns.op(0));
        if(is_lst) for(auto mi : ns) fs.append(F(ps, mi));
        else fs.append(F(ps, ns));
        fs = apart(ls, es, fs);
        map<ex, exset, ex_is_less> ps_int;
        map<ex, int, ex_is_less> ps_pn;
        map<int, ex> pn_ps;
        int pn = 0;
        ex res = MapFunction([&](const ex & e, MapFunction &self)->ex{
            if(!e.has(F(w1,w2))) return e;
            else if(e.match(F(w1,w2))) {
                ex ps = e.op(0);
                ex ns = e.op(1);
                auto itr = ps_int.find(ps);
                if(itr==ps_int.end()) {
                    ps_int[ps].insert(ns);
                    ps_pn[ps] = pn;
                    pn_ps[pn] = ps;
                    pn++;
                } else itr->second.insert(ns);
                return F(ps_pn[ps], ns);
            } else return e.map(self);
        })(fs);
        
        int verb = Verbose;
        Verbose = 0;
        for(auto kv : ps_int) {
            FIRE bc;
            bc.T1 = T1;
            bc.T2 = T2;
            bc.LT1 = LT1;
            bc.LT2 = LT2;
            bc.PosPref = PosPref;
            bc.Internal = ls;
            bc.External = es;
            bc.NVariables[d] = d0;
            bc.Replacement = Replacement;
            bc.Propagator = ex_to<lst>(kv.first);
            bc.ProblemNumber = ps_pn[kv.first];
            for(auto mi : kv.second) bc.Integral.append(mi);
            bc.Reduce();
            bc.FindRules(true);
            system(("rm -rf "+bc.WorkingDir).c_str());
            res = res.subs(bc.Rules);
        }
        Verbose = verb;
        
        res = MapFunction([&](const ex & e, MapFunction &self)->ex{
            if(!e.has(F(w1,w2))) return e;
            else if(e.match(F(w1,w2))) return F(pn_ps[ex2int(e.op(0))], e.op(1));
            else return e.map(self);
        })(res);
        
        PRE = pre_o;
        return is_lst ? res : res.op(0);
    }
    
    matrix Single::oo() {
        if(!In_GiNaC_Parallel && Verbose>0) cout << pre << "\\--" << WHITE << "Single Mode @ infinity ..." << RESET << endl;
        
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
        
        //  ========== to Find all Regions ==========
        
        auto xPropagator = iPropagator;
        xPropagator.let_op(prop_idx) = iPropagator.op(prop_idx) + x;
        
        static symbol y("y");
        int nL = Internal.nops();
        exmap int_map;
        for(auto ii : Internal) int_map[ii] = ii * y;
        lst vLoopPat;
        for(int i=0; i<xPropagator.nops(); i++) {
            if(TopTopo.op(i).is_zero()) continue;
            auto pi = xPropagator.op(i);
            auto item = pi.subs(int_map).expand().coeff(y,2);
            item = exfactor(item);
            if(item.match(pow(w,2))) vLoopPat.append(item.op(0));
            else {
                cout << "item: " << item << endl;
                cout << "pi: " << pi << endl;
                throw Error("item not match w^2!");
            }
        }
        vLoopPat.sort();
        vLoopPat.unique();
        
        lst loops;
        exmap loops_back;
        for(int j=0; j<nL; j++) {
            auto item = symbol("_q"+to_string(j+1));
            loops.append(item);
            loops_back[item] = Internal.op(j);
        }
        
        exmap region_lsp;
        Combinations(vLoopPat.nops(), nL, [&](const int* is) {
            lst eqs;
            for(int j=0; j<nL; j++) {
                eqs.append(vLoopPat.op(is[j])==loops.op(j));
            }
            auto sols = lsolve(eqs, Internal);
            if(sols.nops()>0) {
                auto vlp = subs(vLoopPat, sols).subs(loops_back);
                auto xprop = subs(xPropagator, sols).subs(loops_back).subs(x==y*y);
                PermutationsR(2, nL, [&](const int* pis) {
                    exmap int_map;
                    lst setL, setS;
                    for(int j=0; j<nL; j++) {
                        if(pis[j]) {
                            setL.append(Internal.op(j));
                            int_map[Internal.op(j)] = Internal.op(j) * y;
                        } else setS.append(Internal.op(j));
                    }
                    lst vlp_scale;
                    for(auto item : vlp) {
                        if(item.subs(int_map).expand().coeff(y).normal().is_zero()) vlp_scale.append(0);
                        else vlp_scale.append(1);
                    }
                    if(region_lsp.find(vlp_scale) == region_lsp.end()) {
                        lst pp = ex_to<lst>(xprop.subs(int_map));
                        int nP = xprop.nops();
                        lst powers;
                        for(int i=0; i<nP; i++) {
                            pp.let_op(i) = pp.op(i).normal();
                            if(pp.op(i).has(y)) {
                                powers.append(1);
                                pp.let_op(i) = pp.op(i).expand().coeff(y,2);
                                if(is_zero(pp.op(i))) throw Error("oo: linear propagator with non-zero exponent!");
                            } else powers.append(0);
                        }
                        region_lsp[vlp_scale] = lst{ setL, setS, pp, powers };
                    }
                });
            }
        });
                           
        //  ========== to Select Boundary Conditions ==========
        
        list<matrix> fmat_vec; // boundar equations
        list<lst> fcol_vec; // the last column
        
        lst zero_cl;
        for(int i=0; i<N; i++) zero_cl.append(0);

        for(const auto & kv : CMat) { // not a pattern: l*d0/2=n
            ex la = kv.first;
            bool any = false;
            for(int ll=nloop; ll>=0; ll--) {
                if(normal(la+ll*d0/2).info(info_flags::integer)) {
                    any = true;
                    break;
                }
            }
            if(!any) {
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

        for(int ll=nloop; ll>=0; ll--) {

            ex ila; // select the lambda of type n-L*d/2
            bool first = true;
            for(const auto & la : de.las) {
                if(!normal(la+ll*d0/2).info(info_flags::integer)) continue;
                if(!first) {
                    cout << endl << "las=" << de.las << endl;
                    throw Error("!first: something may be wrong here.");
                }
                ila = la;
                first = false;
            }
            if(first) continue;

            lst dlst;
            int dim = iPropagator.nops();
            lst reg_lst;
            for(int i=0; i<N; i++) {
                ex lp = 0;
                bool first = true;
                ex ns = xMIntegral.op(i).op(1);
                lst reg;
                for(auto kv : region_lsp) {
                    auto lsp = kv.second;
                    if(lsp.op(0).nops() != ll) continue;
                    ex lp2 = 0;
                    auto powers = lsp.op(3);
                    for(int j=0; j<powers.nops(); j++) lp2 += powers.op(j) * ns.op(j);
                    if(first) { first=false; lp = lp2; reg = lst{ lsp }; }
                    else if(lp>lp2) { lp = lp2;  reg = lst{ lsp }; }
                    else if(lp==lp2) reg.append(lsp);
                }
                if(reg.nops()<1) throw Error("reg: region is null.");
                reg_lst.append(reg);
                ex v = lp-ll*d0/2-ila;
                dlst.append(pow(x,v));
            }

            auto t_mat = ex_to<matrix>(diag_matrix(dlst)).inverse(solve_algo::gauss).mul(T);
            auto pr = prank(t_mat,x);
            pr++;
            if(pr<0) pr = 0;
            int order = 0;
            if(ll==0) order = 3;
            pr += order;
            CMat = de.series(pr, lst{ila});

            auto ks = CMat[ila].size();
            for(int k=0; k<ks; k++) {
                matrix cmat = CMat[ila][k][0];
                for(int n=1; n<=pr; n++) cmat = cmat.add(CMat[ila][k][n].mul_scalar(pow(x,n)));
                auto mat = t_mat.mul(cmat);
                int ldeg = 0;
                for(int r=0; r<N; r++) for(int c=0; c<N; c++) {
                    mat(r,c) = series_ex(mat(r,c),x,order); // up to order
                    if(ldeg>mat(r,c).ldegree(x)) ldeg = mat(r,c).ldegree(x);
                }
                for(int l=ldeg; l<=order; l++) {
                    matrix fmat(N, N);
                    for(int r=0; r<N; r++) {
                        for(int c=0; c<N; c++) fmat(r,c) = mat(r,c).coeff(x,l);
                    }
                    if(l>=0 && k==0) {
                        fmat_vec.push_back(fmat);
                        fcol_vec.push_back(lst{reg_lst, l}); // note the order: l
                    } else {
                        fmat_vec.push_front(fmat);
                        fcol_vec.push_front(lst{zero_cl, 0}); // 0 will not used
                    }
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
        
        matrix iBC(N,1);
        matrix inv_mat(N, N);
        int c = 0;
        lst sub_intg;
        for(int r=0; r<N; r++) {
            for( ; c<N*fmat_vec.size(); c++) {
                if(!rref_fmat(r,c).is_zero()) {
                    if(rref_fmat(r,c)!=1) throw Error("pivot is NOT 1 in RREF.");
                    for(int j=0; j<N; j++) inv_mat(r, j) = fmat(j, c);
                    auto itr = fcol_vec.begin();
                    std::advance(itr, c/N);
                    auto lcs = itr->op(0).op(c%N);
                    int xon = ex2int(itr->op(1)); // the expanded order-th
                    auto mi = xMIntegral.op(c%N);
                    if(lcs.is_zero()) iBC(r,0) = 0;
                    else {
                        for(auto lc : lcs) {
                            ex ibc;
                            if(lc.op(0).nops()==0) { // all SMALL region
                                lst ns = ex_to<lst>(mi.op(1));
                                ns.let_op(prop_idx) = -xon;
                                ex cc = 1;
                                if(xon>0) {
                                    symbol p("p");
                                    cc = series_ex(pow(p+1, -mi.op(1).op(prop_idx)), p, xon);
                                    cc = cc.coeff(p, xon);
                                }
                                ibc = cc * F(mi.op(0), ns);
                                sub_intg.append(ns);
                            } else {
                                if(xon!=0) throw Error("Single::oo, xon=!0 NOT supported Yet!");
                                lst lps, sps, lns, sns;
                                auto props = lc.op(2);
                                auto pows = lc.op(3);
                                auto ns = mi.op(1);
                                for(int i=0; i<pows.nops(); i++) {
                                    if(pows.op(i).is_zero()) { sps.append(props.op(i)); sns.append(ns.op(i)); }
                                    else { lps.append(props.op(i)); lns.append(ns.op(i)); }
                                }
                                
                                // Large region
                                lst ls = ex_to<lst>(lc.op(0));
                                lst es = {};
                                auto res = region(ls, es, lps, lns);
                                res = MapFunction([&](const ex & e, MapFunction &self)->ex{
                                    if(!e.has(F(w1,w2))) return e;
                                    else if(e.match(F(w1,w2))) {
                                        int sum = 0;
                                        for(auto item : e.op(1)) {
                                            if(item!=1 && !item.is_zero()) {
                                                cout << "Vacuum: " << e << endl;
                                                throw Error("Single::oo, Unknow Vacuum Integral Found.");
                                            }
                                            sum += ex2int(item);
                                        }
                                        if(isVacuum(ls.nops(), sum)) return Vacuum(ls.nops(), sum).subs(ep==(4-d0)/2);
                                        cout << "Vacuum: " << e << endl;
                                        throw Error("Single::oo, Unknow Vacuum Integral Found.");
                                    } else return e.map(self);
                                })(res);
                                ibc = res;
                                
                                // Small region
                                ls = ex_to<lst>(lc.op(1));
                                es = External;
                                if(ls.nops()>0 && !is_zero(res)) {
                                    res = region(ls, es, sps, sns);
                                    exset fs;
                                    find(res, F(w1,w2), fs);
                                    map<ex,lst,ex_is_less> prop_intg;
                                    for(auto fi : fs) prop_intg[fi.op(0)].append(fi.op(1));
                                    exmap ni_map;
                                    for(auto kv : prop_intg) {
                                        Single amf(*this);
                                        amf.using_FR = true;
                                        amf.Propagator = iPropagator; // using iPropagator
                                        amf.pre = pre + "   ";
                                        amf.Internal = ls;
                                        amf.External = es;
                                        amf.xn += xxn;
                                        amf.Propagator = ex_to<lst>(kv.first);
                                        amf.Integral = ex_to<lst>(kv.second);
                                        amf.Replacement = Replacement;
                                        amf();
                                        int n = kv.second.nops();
                                        for(int i=0; i<n; i++) {
                                            ni_map[F(kv.first, kv.second.op(i))] = amf.NIntegral.op(i);
                                        }
                                    }
                                    res = MapFunction([&](const ex & e, MapFunction &self)->ex{
                                        if(!e.has(F(w1,w2))) return e;
                                        else if(e.match(F(w1,w2))) {
                                            auto itr = ni_map.find(e);
                                            if(itr==ni_map.end()) throw Error("oo: ni_map not found.");
                                            return itr->second;
                                        } else return e.map(self);
                                    })(res);
                                    ibc *= res;
                                }
                            }
                            iBC(r,0) += ibc;
                        }
                    }
                    break;
                }
            }
        }

        if(sub_intg.nops()>0) {
            Single amf(*this);
            amf.using_FR = true;
            amf.Propagator = iPropagator; // using iPropagator
            amf.xn += xxn;
            amf.pre = pre + "   ";
            amf.Integral = sub_intg;
            amf();
            int n = sub_intg.nops();
            exmap ni_map;
            for(int i=0; i<n; i++) {
                ni_map[F(0,sub_intg.op(i))] = amf.NIntegral.op(i);
            }
            n = iBC.rows();
            for(int r=0; r<n; r++) iBC(r,0) = iBC(r,0).subs(ni_map);
        }
        
        inv_mat = inv_mat.inverse(solve_algo::gauss); // TODO: try Fermat or Flint
        iBC = inv_mat.mul(iBC);
        if(iBC.has(F(w1,w2))) {
            cout << iBC << endl;
            throw Error("oo: iBC still has iet or F!");
        }
        if(iBC.has(x)) {
            for(int r=0; r<N; r++) {
                auto item = collect_ex(iBC(r,0), x);
                auto ldeg = item.ldegree(x);
                for(int k=ldeg; k<0; k++) {
                    auto nn = item.coeff(x,k);
                    if(!In_GiNaC_Parallel) {
                        if(abs(nn)>1E-25) {
                            cout << pre << RED << "Note: c/x^" << (-k) << " found with c=" << NN(nn,3) << ", assuming c=0." << RESET << endl;
                        } else {
                            cout << pre << WHITE << "Note: c/x^" << (-k) << " found, but |c| < 10^(-25), assuming c=0." << RESET << endl;
                        }
                    }
                }
                if(ldeg!=0) iBC(r,0) = item.coeff(x,0);
            }
        }
        
        //  ========== to use iBC ==========
        
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
    
    void Single::oo2o(matrix & iBC) {
        int N = Mat.rows();
        int npts = pts.nops();
        if(!In_GiNaC_Parallel && Verbose>0) cout << pre << "\\--" << WHITE << "Single Mode @ middle " << npts-1 << " points ... " << RESET << endl;
        
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
                cout << pre << "\\--" << "Single Mode @ middle points: " << npts-1 << "|" << npts-1-i << flush;
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
    
    matrix Single::o() {
        if(!In_GiNaC_Parallel && Verbose>0) cout << pre << "\\--" << WHITE << "Single Mode @ origin ..." << RESET << endl;
        
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
        if(itr==CMat.end()) throw Error("Single::o, integer key is not found.");
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

    void Single::operator()() {
        if(!In_GiNaC_Parallel && Verbose>0) cout << WHITE << pre << "\\--Single Mode Start @ " << now() << RESET << endl;
        
        if(!init_de) {
            InitDE();
            Poles();
            init_de = true;
        }
        
        auto oTUi = o();
        auto iBC = oo();
        oo2o(iBC);

        matrix oC = oTUi.mul(iBC); // Final result at origin
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
            if(e.has(F(w1,w2))) {
                cout << xMIntegral << endl;
                cout << MIntegral << endl;
                cout << MIxMI << endl;
                cout << "iPropagator: " << iPropagator << endl;
                cout << "Rules: " << Rules << endl;
                cout << "item has F: " << e << endl;
                throw Error("Single Mode: F still exists!");
            }
            NIntegral.append(e);
        }
                        
        if(!In_GiNaC_Parallel && Verbose>0) cout << WHITE << pre << "\\--Single Mode Done @ " << now() << RESET << endl;
    }
        
    ex Single::Vacuum(int nl, int np) { // single mass
        ex res = exp(-I*Pi*(2-ep)*nl)*std::pow(-1,np);
        if(nl==1 && np==1) res *= -tgamma(-1+ep);
        else if(nl==2 && np==3) res *= -pow(tgamma(1-ep),2)*tgamma(ep)*tgamma(-1+2*ep)/tgamma(2-ep);
        else if(nl==3 && np==4) res *= pow(tgamma(1-ep),3)*tgamma(-1+2*ep)*tgamma(-2+3*ep)/tgamma(2-ep);
        else if(nl==3 && np==5) res *= -tgamma(2-3*ep)*pow(tgamma(1-ep),4)*pow(tgamma(ep),2)*tgamma(-1+3*ep)/(pow(tgamma(2-2*ep),2)*tgamma(2-ep));
        else if(nl==4 && np==5) res *= -pow(tgamma(1-ep),4)*tgamma(-2+3*ep)*tgamma(-3+4*ep)/tgamma(2-ep);
        else {
            cout << "nl=" << nl << ", np=" << np << endl;
            throw Error("Vacuum: undefined single-mass vacuum integral");
        }
        return res;
    }
    bool Single::isVacuum(int nl, int np) { // single mass
        if(nl==1 && np==1) return true;
        else if(nl==2 && np==3) return true;
        else if(nl==3 && np==4) return true;
        else if(nl==3 && np==5) return true;
        else if(nl==4 && np==5) return true;
        return false;
    }
    
}

