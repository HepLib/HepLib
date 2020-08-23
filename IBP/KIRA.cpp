#include "Process.h"
#include "IBP.h"

namespace HepLib::IBP {

    namespace {
        string Fout(const ex expr) {
            ex f = expr;
            if(is_a<lst>(f)) f = F(f);
            else if(expr.match(F(w1,w2))) f = F(f.op(1));
            
            string fstr = ex2str(f);
            string_replace_all(fstr,"{","");
            string_replace_all(fstr,"}","");
            string_replace_all(fstr,"(","[");
            string_replace_all(fstr,")","]");
            return fstr;
        }
        
        ex Fin(const string & expr, int pn=0) {
            string fstr = expr;
            string_replace_all(fstr,"[","("+to_string(pn)+",{");
            string_replace_all(fstr,"]","})");
            return str2ex(fstr);
        }
    }

    void KIRA::Export() {

        if(Integrals.nops()<1) return;
        
        if(Round==0 || ibps.nops()<1) {
            _Integrals = Integrals;
            for(auto intg : Integrals) RIntegrals.append(F(ProblemNumber, intg));
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
            
            if(Pairs.nops() > pdim) {
                cout << "Pairs = " << Pairs << endl;
                cout << "Propagators = " << Propagators << endl;
                throw Error("KIRA::Export: Pairs more than Propagators.");
            }
            pdim = Pairs.nops();
            
            lst sp2s, s2sp, ss;
            for(auto item : Pairs) {
                symbol si;
                ss.append(si);
                sp2s.append(item==si);
                s2sp.append(si==item);
            }
            
            lst ibp_eqns;
            for(int i=0; i<pdim; i++) {
                auto eq = Propagators.op(i).expand();
                eq = eq.subs(sp2s, subs_options::algebraic);
                eq = eq.subs(Replacements, subs_options::algebraic);
                ibp_eqns.append(eq == iWF(i));
            }
            auto s2p = lsolve(ibp_eqns, ss);
            if(s2p.nops() != pdim) {
                cout << ibp_eqns << endl;
                cout << s2p << endl;
                throw Error("KIRA::Export: lsolve failed.");
            }

            ibps.remove_all();
            for(auto l : Internal) {
                lst ns;
                for(int i=0; i<Propagators.nops(); i++) ns.append(a(i));
                ex ibp = 0;
                symbol sl;
                for(int i=0; i<Propagators.nops(); i++) {
                    auto ns_tmp = ns;
                    ns_tmp.let_op(i) = ns.op(i) + 1;
                    auto dp = Propagators.op(i).subs(l==sl).diff(sl).subs(sl==l);
                    ibp -= (a(i)+Shift[i]) * F(ns_tmp) * dp;
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
        }
        
        // seeds generation
        exvector eqns;
        exset seeds;
        for(auto seed : _Integrals) {
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
        
        if(true) {
            exset fs;
            find(exvec2lst(eqns), F(w), fs);
            for(auto fi : fs) seeds.insert(fi.op(0));
        }
        
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
                
                for(auto seed : seeds) {
                    lst sol, as;
                    for(int i=0; i<tot; i++) {
                        auto ai = (seed.op(i)-cc[i][1])/cc[i][0];
                        sol.append(a(i)==ai);
                        as.append(ai);
                    }
                    auto ii = item.subs(sol);
                    if(ii.is_zero()) continue;
                    eqns.push_back(ii);
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
                for(auto fi : fs) {
                    lst ns = ex_to<lst>(fi.op(0));
                    for(auto ic : Cuts) {
                        int j = ex_to<numeric>(ic).to_int()-1;
                        if(ns.op(j)<=0) {
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
        system(("rm -rf "+job_dir).c_str());
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
        if(mi_pref.nops()>0) {
            ostringstream oss2;
            int nn = mi_pref.nops();
            for(int i=0; i<nn; i++) oss2 << Fout(mi_pref.op(i)) << endl;
            ofstream pref_out(job_dir+"/preferred");
            pref_out << oss2.str() << endl;
            pref_out.close();
            oss << "        preferred_masters: preferred" << endl;
        }
        oss << "    - kira2file:" << endl;
        oss << "        target:" << endl;
        oss << "            - [F,integrals]" << endl;
        ofstream job_out(job_dir+"/jobs");
        job_out << oss.str() << endl;
        job_out.close();
        
        oss.str("");
        oss.clear();
        for(auto integral : _Integrals) oss << Fout(integral) << endl;
        ofstream intg_out(job_dir+"/integrals");
        intg_out << oss.str() << endl;
        intg_out.close();
        
    }
    
    void KIRA::Run() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream cmd;
        cmd << "cd " << job_dir << " && kira " << cmd_args << " --silent jobs >/dev/null 2>&1";
        system(cmd.str().c_str());
    }

    void KIRA::Import() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream fn;
        fn << job_dir << "/results/F/kira_integrals.kira";
        auto strvec = file2vec(fn.str());
        
        ex exL=0, exR=0;
        map<ex,int,ex_is_less> flags;
        lst exRs;
        for(auto intg : _Integrals) flags[F(ProblemNumber,intg)] = 1;
        _Rules.remove_all();
        for(auto line : strvec) {
            if(line.size()==0) {
                if(!is_zero(exL)) {
                    _Rules.append(exL==exR);
                    flags[exL] = 0;
                    exRs.append(exR);
                }
                exL = exR = 0;
            } else if(is_zero(exL)) {
                exL -= Fin(line,ProblemNumber);
                if(!exL.match(F(w1,w2))) {
                    cout << line << endl;
                    throw Error("KIRA::Import error found.");
                }
            } else {
                exR += Fin(line,ProblemNumber).subs(d==D);
            }
        }
        if(!is_zero(exL)) {
            _Rules.append(exL==exR);
            flags[exL] = 0;
            exRs.append(exR);
        }
        MasterIntegrals.remove_all();
        for(auto kv : flags) if(kv.second!=0) MasterIntegrals.append(kv.first);
        exset miset;
        find(exRs,F(w1,w2),miset);
        for(auto mi : miset) MasterIntegrals.append(mi);
        MasterIntegrals.sort();
        MasterIntegrals.unique();
        
        auto _RIntegrals = RIntegrals;
        for(int i=0; i< RIntegrals.nops(); i++) {
            RIntegrals.let_op(i) = RIntegrals.op(i).subs(_Rules);
        }
        
        auto integrals = _Integrals;
        _Integrals.remove_all();
        for(auto mi : MasterIntegrals) _Integrals.append(mi.op(1));
        
        bool red = _Rules.nops()>0 && _Integrals.nops()>0;
        red = red && !is_zero(_Integrals-integrals);
        red = red && !is_zero(_RIntegrals-RIntegrals);
        
        Round++;
        if(Round<Rounds && red) Reduce();
        else {
            Rules.remove_all();
            for(int i=0; i<Integrals.nops(); i++) {
                auto ii = F(ProblemNumber,Integrals.op(i));
                auto ri = RIntegrals.op(i);
                if(!is_zero(ii-ri)) Rules.append(ii==ri);
            }
        }
        
    }

}
