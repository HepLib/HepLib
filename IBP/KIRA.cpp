#include "Process.h"
#include "IBP.h"

namespace HepLib::IBP {

    namespace {
        string Fout(const ex expr) {
            ex f = expr;
            if(is_a<lst>(f)) f = F(f);
            string fstr = ex2str(f);
            string_replace_all(fstr,"{","");
            string_replace_all(fstr,"}","");
            string_replace_all(fstr,"(","[");
            string_replace_all(fstr,")","]");
            return fstr;
        }
        
        ex Fin(const string & expr) {
            string fstr = expr;
            string_replace_all(fstr,"[","({");
            string_replace_all(fstr,"]","})");
            return str2ex(fstr);
        }
    }

    void KIRA::Export() {

        if(Integrals.nops()<1) return;
        int pn = 0; // to avoid unsigned short overflow in FIRE
        int pdim = Propagators.nops();
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);
        
        if(Pairs.nops()<1) {
            for(auto it : Internal) {
                for(auto ii : InExternal) Pairs.append(it*ii);
            }
            Pairs.sort();
            Pairs.unique();
            
            if(Pairs.nops()<pdim) {
                lst sps_ext;
                for(auto it : External) {
                    for(auto ii : External) sps_ext.append(it*ii);
                }
                sps_ext.sort();
                sps_ext.unique();
                for(auto item : sps_ext) {
                    auto item2 = subs_all(item,Replacements);
                    if(is_zero(item-item2)) Pairs.append(item);
                }
                Pairs.sort();
                Pairs.unique();
            }
        }
        
        if(Pairs.nops() != pdim) {
            cout << "Pairs = " << Pairs << endl;
            cout << "Propagators = " << Propagators << endl;
            throw Error("KIRA::Export: Pairs failed.");
        }
        
        lst sp2s, s2sp, ss;
        for(auto item : Pairs) {
            symbol si;
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst ibp_eqns;
        for(int i=0; i<Propagators.nops(); i++) {
            auto eq = Propagators.op(i).expand();
            eq = eq.subs(sp2s, subs_options::algebraic);
            eq = eq.subs(Replacements, subs_options::algebraic);
            ibp_eqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(ibp_eqns, ss);
        if(s2p.nops() != pdim) throw Error("KIRA::Export: lsolve failed.");

        lst ibps;
        lst ns;
        for(int i=0; i<Propagators.nops(); i++) ns.append(a(i));
        for(auto l : Internal) {
            ex ibp = 0;
            symbol sl;
            for(int i=0; i<Propagators.nops(); i++) {
                auto ns_tmp = ns;
                ns_tmp.let_op(i) = ns.op(i) + 1;
                auto dp = Propagators.op(i).subs(l==sl).diff(sl).subs(sl==l);
                ibp -= a(i) * F(ns_tmp) * dp;
            }
            
            for(auto ii : External) {
                auto ibp_tmp = ibp * ii;
                ibp_tmp = ibp_tmp.expand();
                ibp_tmp = ibp_tmp.subs(sp2s, subs_options::algebraic);
                ibp_tmp = ibp_tmp.subs(Replacements, subs_options::algebraic);
                ibp_tmp = ibp_tmp.subs(s2p, subs_options::algebraic);
                ex res = 0;
                for(int i=0; i<Propagators.nops(); i++) {
                    auto ci = ibp_tmp.coeff(iWF(i), 1);
                    ci = MapFunction([i](const ex & e, MapFunction &self)->ex{
                        if(e.match(F(w))) {
                            auto tmp = e.op(0);
                            tmp.let_op(i) = tmp.op(i)-1;
                            return F(tmp);
                        } else if(!e.has(F(w))) return e;
                        else return e.map(self);
                    })(ci);
                    res += ci;
                }
                res += ibp_tmp.subs(lst{iWF(w)==0});
                ibps.append(res);
            }
            
            for(auto ii : Internal) {
                auto ibp_tmp = ibp * ii;
                ibp_tmp = ibp_tmp.expand();
                ibp_tmp = ibp_tmp.subs(sp2s, subs_options::algebraic);
                ibp_tmp = ibp_tmp.subs(Replacements, subs_options::algebraic);
                ibp_tmp = ibp_tmp.subs(s2p, subs_options::algebraic);
                ex res = 0;
                for(int i=0; i<Propagators.nops(); i++) {
                    auto ci = ibp_tmp.coeff(iWF(i), 1);
                    ci = MapFunction([i](const ex &e, MapFunction &self)->ex {
                        if(e.match(F(w))) {
                            auto tmp = e.op(0);
                            tmp.let_op(i) = tmp.op(i)-1;
                            return F(tmp);
                        } else if(!e.has(F(w))) return e;
                        else return e.map(self);
                    })(ci);
                    res += ci;
                }
                res += ibp_tmp.subs(lst{iWF(w)==0});
                if(ii==l) res += d*F(ns);
                ibps.append(res);
            }
        }
        
        // seeds generation
        exvector eqns;
        exset seeds;
        for(auto seed : Integrals) {
            for(auto const & item : ibps) {
                exset fs;
                item.find(F(w), fs);
                for(auto fi : fs) {
                    lst sol, as;
                    for(int i=0; i<fi.op(0).nops(); i++) {
                        auto expn = fi.op(0).op(i);
                        auto ca = expn.coeff(a(i),1);
                        auto c0 = expn.coeff(a(i),0);
                        auto ai = (seed.op(i)-c0)/ca;
                        sol.append(a(i)==ai);
                        as.append(ai);
                    }

                    seeds.insert(as);
                    auto ii = item.subs(sol);
                    if(ii.is_zero()) continue;
                    eqns.push_back(ii);
                }
            }
        }
        
        for(int r=0; r<SeedRuns; r++) {
            auto seeds_tmp = seeds;
            for(auto const & item : ibps) {
                exset fs;
                item.find(F(w), fs);
                for(auto fi : fs) {
                    int tot = fi.op(0).nops();
                    ex cc[tot][2];
                    for(int i=0; i<tot; i++) {
                        auto expn = fi.op(0).op(i);
                        auto ca = expn.coeff(a(i),1);
                        auto c0 = expn.coeff(a(i),0);
                        cc[i][0] = ca;
                        cc[i][1] = c0;
                    }
                    
                    for(auto seed : seeds_tmp) {
                        lst sol, as;
                        for(int i=0; i<tot; i++) {
                            auto ai = (seed.op(i)-cc[i][1])/cc[i][0];
                            sol.append(a(i)==ai);
                            as.append(ai);
                        }
                        seeds.insert(as);
                        auto ii = item.subs(sol);
                        if(ii.is_zero()) continue;
                        eqns.push_back(ii);
                    }
                }
            }
        }
        
        if(Cuts.nops()>1) {
            int total = eqns.size();
            for(int i=0; i<total; i++) {
                auto tmp = eqns[i];
                exset fs;
                find(tmp, F(w), fs);
                exmap repl;
                for(auto fi : fs){
                    for(auto ic : Cuts) {
                        int j = ex_to<numeric>(ic).to_int();
                        if(fi.op(0).op(j)<=0) {
                            repl[fi]=0;
                            break;
                        }
                    }
                }
                eqns[i] = tmp.subs(repl);
            }
        }
        
        if(WorkingDir.length()<1) WorkingDir = to_string(getpid());
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        if(!dir_exists(job_dir)) system(("mkdir -p "+job_dir).c_str());
        
        ostringstream oss;
        for(auto eqn : eqns) {
            if(is_zero(eqn)) continue;
            auto cvs = mma_collect_lst(eqn,F(w));
            for(auto cv : cvs) {
                oss << Fout(cv.op(1)) << " * (" << cv.op(0) << ")" << endl;
            }
            oss << endl;
        }
        
        ofstream eqn_out(job_dir+"/equations");
        string ostr = oss.str();
        string_replace_all(ostr, "{", "");
        string_replace_all(ostr, "}", "");
        eqn_out << ostr << endl;
        eqn_out.close();
        
        oss.str("");
        oss.clear();
        oss << "jobs:" << endl;
        oss << "    - reduce_user_defined_system:" << endl;
        oss << "        input_system: equations" << endl;
        oss << "    - kira2file:" << endl;
        oss << "        target:" << endl;
        oss << "            - [F,integrals]" << endl;
        ofstream job_out(job_dir+"/jobs");
        job_out << oss.str() << endl;
        job_out.close();
        
        oss.str("");
        oss.clear();
        for(auto integral : Integrals) oss << Fout(integral) << endl;
        ofstream intg_out(job_dir+"/integrals");
        intg_out << oss.str() << endl;
        intg_out.close();
        
    }
    
    void KIRA::Run() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream cmd;
        cmd << "cd " << job_dir << " && kira jobs > /dev/null 2>/dev/null";
        system(cmd.str().c_str());
    }

    void KIRA::Import() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream fn;
        fn << job_dir << "/results/F/kira_integrals.kira";
        auto strvec = file2vec(fn.str());
        
        ex exL=0, exR=0;
        for(auto line : strvec) {
            if(line.size()==0) {
                if(!is_zero(exL)) {
                    Rules.append(exL==exR);
                }
                exL = exR = 0;
            } else if(is_zero(exL)) {
                exL -= Fin(line);
                if(!is_zero(exL.subs(F(w)==1)-1)) {
                    cout << line << endl;
                    throw Error("KIRA::Import error found.");
                }
            } else {
                exR += Fin(line);
            }
        }
    }

}
