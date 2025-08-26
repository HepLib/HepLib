/**
 * @file
 * @brief Others
 */
 
#include "HEP.h"
#include <cmath>

namespace HepLib {
    
    namespace {
        int ep_ldegree(const ex & e) {
            auto ee = exnormal(e);
            auto nd = ee.numer_denom();
            if(nd.op(0).is_polynomial(ep) && nd.op(1).is_polynomial(ep)) {
                return nd.op(0).ldegree(ep) - nd.op(1).ldegree(ep);
            }
            cout << "numer or denom is NOT a polynormial of ep." << endl;
            abort();
        }
        
        
        class Family {
        public:
            int pn;
            lst intg;
            Family(int pn_) : pn(pn_) { }
        };
            
        exmap _sp2x_(const ex & ips, const ex & eps) {
            int ni = 0;
            exmap sp2x;
            for(auto ip : ips) {
                for(auto ipx : ips) {
                    if(sp2x.find(ip*ipx)==sp2x.end()) {
                        auto xi = x(ni);
                        sp2x[ip*ipx] = xi;
                        sp2x[w*ip*ipx] = w*xi;
                        ni++;
                    }
                }
                for(auto ep : eps) {
                    auto xi = x(ni);
                    sp2x[ep*ip] = xi;
                    sp2x[w*ep*ip] = w*xi;
                    ni++;
                }
            }
            return sp2x;
        }

        bool is_complete(const ex & props, const ex & ips, const ex & eps) {
            auto sp2x = _sp2x_(ips, eps);
            auto ps = subs(props, sp2x);
            int n = ps.nops();
            matrix mat(n,n);
            for(int r=0; r<n; r++) {
                auto xi = x(r);
                for(int c=0; c<n; c++) {
                    mat(r,c) = ps.op(c).coeff(xi);
                }
            }
            return mat.rank()==n;
        }
        bool is_complete(const ex & props, const exmap & sp2x) {
            auto ps = subs(props, sp2x);
            int n = ps.nops();
            matrix mat(n,n);
            for(int r=0; r<n; r++) {
                auto xi = x(r);
                for(int c=0; c<n; c++) {
                    mat(r,c) = ps.op(c).coeff(xi);
                }
            }
            return mat.rank()==n;
        }
    }
    
    bool Exp2AMF::is_loop(const ex & e) {
        for(auto ei : Internal) if(ei.is_equal(e)) return true;
        return false;
    }
    
    void Exp2AMF::Export(const ex & expr, const string & dir) {
        static string wlo =
R"EOF(
If[Length[$ScriptCommandLine]<2, Print["Usage: wolframescript -f "<>ToString[$ScriptCommandLine[[1]]]<>" <n>"];Quit[]];
current = DirectoryName@If[$FrontEnd===Null,$InputFileName,NotebookFileName[]];

<<AMFlow`AMFlow`

SetReductionOptions["IBPReducer"->"Blade"];
ID = ToExpression[$ScriptCommandLine[[2]]];
PropsInts=Get[FileNameJoin[{current,ToString[ID]<>".m"}]];

AMFlowInfo["Family"] = Symbol["I"<>ToString[ID]];
AMFlowInfo["Loop"] = <<Loop>>;
AMFlowInfo["Leg"] = <<Leg>>;
AMFlowInfo["Conservation"] = { };
AMFlowInfo["Replacement"] = <<Replacement>>;
AMFlowInfo["Propagator"] = PropsInts[[1]]/.<<M2M2>>;
AMFlowInfo["Numeric"] = <<Numeric>>;
AMFlowInfo["NThread"] = <<NThread>>;

integrals = Map[(j[Symbol["I"<>ToString[ID]], Sequence@@#])&, PropsInts[[2]]];
precision = <<Precision>>;
epsorder = 2*Length@AMFlowInfo["Loop"]-PropsInts[[3]] + <<Order>>;
sol = SolveIntegrals[integrals, precision, epsorder];
Put[sol, FileNameJoin[{current, "sol"<>ToString[ID]<>".m"}]];

Quit[];
)EOF";

        string sLoop = ex2str(Internal);
        string_replace_all(sLoop, "==", "->");
        string sLeg = ex2str(External);
        string_replace_all(sLeg, "==", "->");
        string sReplacement = ex2str(Replacement);
        string_replace_all(sReplacement, "==", "->");
        string sNumeric = ex2str(Numeric);
        string_replace_all(sNumeric, "==", "->");
        string sM2M2 = ex2str(M2M2);
        string_replace_all(sM2M2, "==", "->");
        string wl = wlo;
        string_replace_all(wl, "<<Loop>>", sLoop);
        string_replace_all(wl, "<<Leg>>", sLeg);
        string_replace_all(wl, "<<Replacement>>", sReplacement);
        string_replace_all(wl, "<<Numeric>>", sNumeric);
        string_replace_all(wl, "<<NThread>>", to_string(NThread));
        string_replace_all(wl, "<<M2M2>>", sM2M2);
        string_replace_all(wl, "<<Order>>", ex2str(Order));
        string_replace_all(wl, "<<Precision>>", ex2str(Precision));

        system(("mkdir -p "+dir).c_str());
        str2file(wl, dir+"/run.wl");
        
        auto res = expr;
        
        // handle (qi*pi)^-n
        cout << "Replacing linear propagators ..." << endl;
        sp2x = _sp2x_(Internal, External);
        exset pq_fs;
        find(res, F(w1,w2), pq_fs);
        exmap pq_f2f;
        for(auto f : pq_fs) {
            lst ps = ex_to<lst>(f.op(0));
            lst ns = ex_to<lst>(f.op(1));
            lst ns_p; // problematic ns
            for(int i=0; i<ps.nops(); i++) {
                bool ok = false;
                for(auto ip : Internal) {
                    if(ps.op(i).has(ip*ip)) {
                        ok = true;
                        break;
                    }
                }
                if(!ok) {
                    ex pi = ps.op(i);
                    if(!is_a<mul>(pi) || !pi.nops()==2) {
                        cout << "prop is NOT the form of a*b" << endl;
                        abort();
                    }
                
                    if(ns.op(i)==-1) {
                        ns_p.append(i);
                    } else if(ns.op(i)!=0) {
                        cout << f << endl;
                        cout << "n of qi.pi is NOT -1" << endl;
                        abort();
                    } else { // ns.op(i)==0
                        auto pa = pi.op(0), pb = pi.op(1);
                        ex pab = expand(pow(pa+pb,2)).subs(Replacement);
                        ps[i] = pab;
                        bool ic = is_complete(ps, sp2x);
                        if(!ic && is_loop(pa)) {
                            ps[i] = expand(pow(pa,2)).subs(Replacement);
                            ic = is_complete(ps, sp2x);
                        }
                        if(!ic && is_loop(pb)) {
                            ps[i] = expand(pow(pb,2)).subs(Replacement);
                            ic = is_complete(ps, sp2x);
                        }
                        if(!ic) {
                            cout << "1st: ic is NOT true!" << endl;
                            abort();
                        }
                    }
                }
            }
            if(ns_p.nops()>1) {
                cout << f << endl;
                abort();
            } else if(ns_p.nops()>0) {
                auto idx = ex2int(ns_p.op(0));
                ns[idx] = 0;
                auto pidx = ps.op(idx);
                
                ex pa = pidx.op(0);
                ex pb = pidx.op(1);
                
                ex pab = expand(pow(pa+pb,2)).subs(Replacement);
                ps[idx] = pab;
                bool ic = is_complete(ps, sp2x);
                if(!ic && is_loop(pa)) {
                    ps[idx] = expand(pow(pa,2)).subs(Replacement);
                    ic = is_complete(ps, sp2x);
                }
                if(!ic && is_loop(pb)) {
                    ps[idx] = expand(pow(pb,2)).subs(Replacement);
                    ic = is_complete(ps, sp2x);
                }
                if(!ic) {
                    cout << "2nd: ic is NOT true!" << endl;
                    abort();
                }
                
                auto xps = subs(ps,sp2x);
                auto xpidx = pidx.subs(sp2x);
                for(auto ip : Internal) {
                    if(xps.has(ip)) {
                        cout << "xps still has ip" << endl;
                        cout << "xps: " << xps << endl;
                        abort();
                    }
                    if(xpidx.has(ip)) {
                        cout << "xpidx still has ip" << endl;
                        cout << "xpidx: " << xpidx << endl;
                        abort();
                    }
                }
                
                ex xeq = 0, eq = 0, fres = 0;
                lst ys;
                for(int ii=0; ii<ps.nops(); ii++) {
                    symbol yi;
                    xeq += yi*xps.op(ii);
                    eq += yi*ps.op(ii);
                    lst nn = ns;
                    nn[ii] = ns.op(ii)-1;
                    fres += yi*F(ps, nn);
                    ys.append(yi);
                }
                lst eqs;
                for(int ii=0; ii<ps.nops(); ii++) {
                    auto xi = x(ii);
                    eqs.append(xeq.coeff(xi)==xpidx.coeff(xi));
                }
                auto sol = lsolve(eqs,ys);
                
                if(sol.nops()<1 || sol.has(x(w))) {
                    cout << "wrong sol: " << sol << endl;
                }
                
                auto rem = pidx-eq.subs(sol);
                fres = fres.subs(sol);
                if(!is_zero(rem)) {
                    fres += rem * F(ps, ns);
                }
                
                pq_f2f[f] = fres;
                ex chk = exnormal(F2ex(f-pq_f2f[f]));
                if(chk!=0) {
                    cout << "Failed: chk=" << chk << endl;
                }
            } else if(!f.op(0).is_equal(ps)) {
                pq_f2f[f] = F(ps, ns);
                ex chk = exnormal(F2ex(f-pq_f2f[f]));
                if(chk!=0) {
                    cout << "Failed: chk=" << chk << endl;
                }
            }
        }
        
        cout << "Collecting ..." << endl;
        
        res = MapFunction([&](const ex & e, MapFunction &self)->ex{
            if(e.match(F(w1, w2))) {
                auto itr = pq_f2f.find(e);
                if(itr==pq_f2f.end()) {
                    cout << "F Not Found" << endl;
                    abort();
                }
                return itr->second;
            } else return e.map(self);
        })(res);
        
        res = collect_ex(res, F(w1,w2));

        exset fs;
        find(res, F(w1,w2), fs);
        map<ex,Family,ex_is_less> p2f;
        exmap f2f;
        int pn = 0;
        for(auto f : fs) {
            auto itr = p2f.find(f.op(0));
            if(itr == p2f.end()) {
                itr = p2f.insert({f.op(0), Family(pn)}).first;
                pn++;
            }
            itr->second.intg.append(f.op(1));
            f2f[f] = F(itr->second.pn, f.op(1));
        }
        
        res = MapFunction([&](const ex & e, MapFunction &self)->ex{
            if(e.match(F(w1, w2))) {
                auto itr = f2f.find(e);
                if(itr==f2f.end()) {
                    cout << "F Not Found" << endl;
                    abort();
                }
                return itr->second;
            } else return e.map(self);
        })(res);
        
        cout << "Exporting res.txt ..." << endl;
        ex2file(res.subs(Numeric), dir+"/res.txt");
        
        string sres =
R"EOF(
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;
<<HepLib`
Module[{sols, res},
sols = Flatten@Table[Import[FileNameJoin[{current,"sol"<>ToString[i]<>".m"}]], {i,0,<<NN>>}];
res = C2M@Import[FileNameJoin[{current, "res.txt"}]];
res = res/.{F[p_, n_] :> j[Symbol["I"<>ToString[p]], Sequence@@n]} /. Dispatch@sols;
res
]
)EOF";
        
        string_replace_all(sres, "<<NN>>", to_string(pn-1));
        str2file(sres, dir+"/res.m");
        
        int tot = p2f.size(), cur = 0;
        for(auto kv : p2f) {
            cur++;
            cout << "[ " <<  cur << " / " << tot << "] pn = " << kv.second.pn << endl;
            cout << kv.first << endl;
            cout << kv.second.intg << endl;
            cout << endl;
            
            lst prop = ex_to<lst>(kv.first);
            auto & f = kv.second;
            for(int i=0; i<prop.nops(); i++) {
                bool ok = false;
                for(auto ip : Internal) {
                    if(prop.op(i).has(ip*ip)) {
                        ok = true;
                        break;
                    }
                }
                if(!ok) {
                    cout << "Propagators: " << prop.op(i) << " @ pos: " << i << endl;
                    cout << "Integrals: " << f.intg << endl;
                    abort();
                }
            }
            
            GiNaC_Parallel_Verb["Exp2AMF"] = 0;
            auto order_vec = GiNaC_Parallel(f.intg.nops(), [&](int idx)->ex{
                auto ns = f.intg.op(idx);
                ex cc = res.coeff(F(f.pn, ns)).subs(Numeric);
                cc = cc.subs(d==4-2*ep);
                auto ldeg = ep_ldegree(cc);
                return ldeg;
            }, "Exp2AMF");
            
            int order = 0;
            for(auto ldeg : order_vec) {
                if(order>ldeg) order = ex2int(ldeg);
            }
            
            ex2file(lst{ prop, f.intg, order }, dir+"/"+to_string(f.pn)+".m");
        }
                
    }

}
