/**
 * @file
 * @brief IBP with KIRA
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib::IBP {

    /**
     * @brief export integral to KIRA form
     */
    string UKIRA::Fout(const ex & expr) {
        if(!using_uw) {
            ex f = expr;
            if(is_a<lst>(f)) f = F(f);
            else if(expr.match(F(w1,w2))) f = F(f.op(1));
            string fstr = ex2str(f);
            string_replace_all(fstr,"{","");
            string_replace_all(fstr,"}","");
            string_replace_all(fstr,"(","[");
            string_replace_all(fstr,")","]");
            return fstr;
        } else {
            ex idx;
            if(expr.match(F(w1,w2))) idx = expr.op(1);
            else if(expr.match(F(w))) idx = expr.op(0);
            else idx = expr;
            return to_string(i2w[idx]);
        }
    }
    
    /**
     * @brief import integral from KIRA
     */
    ex UKIRA::Fin(const string & expr) {
        if(!using_uw) {
            string fstr = expr;
            string_replace_all(fstr,"[","("+to_string(ProblemNumber)+",{");
            string_replace_all(fstr,"]","})");
            return str2ex(fstr);
        } else {
            auto cpos = expr.find("*");
            if(cpos==string::npos) {
                if(expr=="0") return 0;
                throw Error("UKIRA::Fin with 0 or * NOT Found.");
            }
            auto wstr = expr.substr(0,cpos);
            unsigned long long weight = stoull(wstr,NULL,0);
            auto oex = str2ex(expr.substr(cpos+1,string::npos));
            return F(ProblemNumber, w2i[weight]) * oex;
        }
    }
    
    /**
     * @brief Export input data for KIRA
     */
    void UKIRA::Export() {

        if(Integrals.nops()<1) return;
        
        int pdim = Propagators.nops();
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);
        
        if(ISP.nops()<1) {
            for(auto it : Internal) {
                for(auto ii : InExternal) ISP.append(it*ii);
            }
            ISP.sort();
            ISP.unique();
        }
        
        if(ISP.nops() > pdim) {
            cout << "ISP = " << ISP << endl;
            cout << "Propagators = " << Propagators << endl;
            throw Error("UKIRA::Export: #(ISP) > #(Propagators).");
        }
        
        lst sp2s, s2sp, ss;
        int _pic=0;
        for(auto item : ISP) {
            _pic++;
            Symbol si("P"+to_string(_pic));
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst leqns;
        for(int i=0; i<ISP.nops(); i++) { // note NOT pdim
            auto eq = Propagators.op(i).expand().subs(iEpsilon==0); // drop iEpsilon
            eq = eq.subs(sp2s, algbr);
            eq = eq.subs(Replacements, algbr);
            if(eq.has(iWF(w))) throw Error("UKIRA::Export, iWF used in eq.");
            leqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(leqns, ss);
        if(s2p.nops() != ISP.nops()) throw Error("UKIRA::Export: lsolve failed.");
        
        if(DSP.nops()<1) {
            for(auto p1 : Internal)
            for(auto p2 : InExternal)
            DSP.append(lst{p1,p2});
        }

        ibps.remove_all(); // no need
        lst nsa;
        for(int i=0; i<pdim; i++) nsa.append(a(i));
        for(auto sp : DSP) {
            auto ilp = sp.op(0);
            auto iep = sp.op(1);
            
            ex ibp = 0;
            symbol ss;
            for(int i=0; i<pdim; i++) {
                auto ns = nsa;
                ns.let_op(i) = nsa.op(i) + 1;
                auto dp = Propagators.op(i).subs(ilp==ss).diff(ss).subs(ss==ilp);
                ibp -= (a(i)+Shift[i+1]) * F(ns) * dp;
            }
            
            ibp = ibp * iep;
            ibp = ibp.expand();
            ibp = ibp.subs(sp2s, algbr);
            ibp = ibp.subs(Replacements, algbr);
            ibp = ibp.subs(s2p, algbr);
            
            ex res = 0;
            for(int i=0; i<pdim; i++) {
                auto ci = ibp.coeff(iWF(i), 1);
                ci = MapFunction([i](const ex &e, MapFunction &self)->ex {
                    if(!e.has(F(w))) return e;
                    else if(e.match(F(w))) {
                        lst tmp = ex_to<lst>(e.op(0));
                        tmp.let_op(i) = tmp.op(i)-1;
                        return F(tmp);
                    } else return e.map(self);
                })(ci);
                res += ci;
            }
            res += ibp.subs(lst{iWF(w)==0});
            auto cilp = iep.coeff(ilp);
            if(!is_zero(cilp)) res += d*cilp*F(nsa);
            ibps.append(res);
        }
        
        // seeds generation
        int rmax = -1, smax = -1;
        int rrmax[pdim], ssmax[pdim];
        for(int i=0; i<pdim; i++) rrmax[i] = ssmax[i] = -1;
        for(auto integral : Integrals) {
            if(integral.nops()!=pdim) throw Error("UKIRA::Export, integral dimension not match propagators.");
            int rr = 0;
            int ss = 0;
            for(int i=0; i<pdim; i++) {
                auto item = ex2int(integral.op(i));
                if(item>0) {
                    if(rrmax[i]<item) rrmax[i] = item;
                    rr += item;
                } else {
                    item = 0-item;
                    if(ssmax[i]<item) ssmax[i] = item;
                    ss += item;
                }
            }
            if(rmax<rr) rmax = rr;
            if(smax<ss) smax = ss;
        }
        
        if(seed_option == 0) {
            int _rrmax = -1, _ssmax = -1;
            for(int i=0; i<pdim; i++) {
                if(_rrmax<rrmax[i]) _rrmax = rrmax[i];
                if(_ssmax<ssmax[i]) _ssmax = ssmax[i];
            }
            for(int i=0; i<pdim; i++) {
                rrmax[i] = _rrmax + rap;
                ssmax[i] = _ssmax + sap;
            }
            rmax += ra;
            smax += sa;
        } else {
            for(int i=0; i<pdim; i++) {
                rrmax[i] += rap;
                ssmax[i] += sap;
            }
            rmax += ra;
            smax += sa;
        }
        
        // generate IBP equations
        int as[pdim];
        for(int i=0; i<pdim; i++) as[i] = -ssmax[i];
        vector<vector<int>> asvec;
        while(true) {
            for(int i=0; i<pdim; i++) {
                if(as[i]+1>rrmax[i]) {
                    if(i+1 == pdim) goto done;
                    as[i] = -ssmax[i];
                } else {
                    as[i] = as[i] + 1;
                    break;
                }
            }
            int _rmax = 0, _smax = 0;
            for(int i=0; i<pdim; i++) {
                if(as[i]>0) _rmax += as[i];
                else _smax -= as[i];
            }
            if(_rmax>rmax || _smax>smax) continue;
            vector<int> asv;
            for(int i=0; i<pdim; i++) asv.push_back(as[i]);
            asvec.push_back(asv);
        }
        done: ;
        
        int nCut = Cuts.nops();
        bool hasCut = (nCut>1);
        int iCuts[nCut+1];
        for(int i=0; i<nCut; i++) iCuts[i] = ex_to<numeric>(Cuts.op(i)).to_int();
        auto eqns_result =
        GiNaC_Parallel(asvec.size(), [this,&asvec,pdim,hasCut,nCut,&iCuts](int idx)->ex {
            auto as = asvec[idx];
            exmap sol;
            for(int i=0; i<pdim; i++) sol[a(i)]=as[i];
            lst eqns;
            for(auto item : ibps) {
                auto ii = item.subs(sol);
                if(ii.is_zero()) continue;
                exset fs;
                find(ii, F(w), fs);
                if(hasCut) {
                    exmap repl;
                    for(auto fi : fs) {
                        lst ns = ex_to<lst>(fi.op(0));
                        for(int nc=0; nc<nCut; nc++) {
                            int j = iCuts[nc]-1;
                            if(ns.op(j)<=0) {
                                repl[fi]=0;
                                break;
                            }
                        }
                    }
                    ii = ii.subs(repl);
                }
                if(ii.is_zero()) continue;
                
                exmap repl;
                for(auto li : Internal) {
                    for(auto fi : fs) {
                        lst ns = ex_to<lst>(fi.op(0));
                        bool has = false;
                        for(int i=0; i<pdim; i++) {
                            if(is_zero(Shift[i+1]) && ns.op(i)<=0) continue;
                            if(Propagators.op(i).has(li)) {
                                has = true;
                                break;
                            }
                        }
                        if(!has) repl[fi]=0;
                    }
                }
                ii = ii.subs(repl);
                if(ii.is_zero()) continue;
                
                eqns.append(ii);
            }
            return eqns;
        }, "SEED");
        
        lst eqns;
        for(auto ilst : eqns_result) {
            if(ilst.nops()<1) continue;
            for(auto eqn : ilst) eqns.append(eqn);
        }

        if(using_uw) {
            exset fs;
            for(auto eqn : eqns) find(eqn, F(w), fs);
            for(auto intg : Integrals) fs.insert(F(intg));
            for(auto intg : PIntegrals) fs.insert(F(intg));
            exvector intg_vec;
            for(auto fi : fs) {
                lst rs,ss;
                int sid=0, rsum=0, ssum=0;
                auto idx_lst = fi.op(0);
                for(int i=0; i<idx_lst.nops(); i++) {
                    auto idx = ex_to<numeric>(idx_lst.op(i)).to_int();
                    if(idx!=0) sid += std::pow(2,idx_lst.nops()-i-1);
                    if(idx>0) {
                        rs.append(idx);
                        rsum += idx;
                    } else {
                        ss.append(idx);
                        ssum -= idx;
                    }
                }
                
                lst item;
                if(sort_option==0) { // {r+s,r,s}
                    item.append(rsum+ssum);
                    item.append(rsum);
                    item.append(ssum);
                    for(auto ii : ss) item.append(-ii);
                    for(auto ii : rs) item.append(ii);
                } else if(sort_option==1) { // {S,r,s,ss,rr}
                    item.append(rsum+ssum);
                    item.append(rsum);
                    item.append(ssum);
                    item.append(sid);
                    for(auto ii : ss) item.append(ii);
                    for(auto ii : rs) item.append(ii);
                } else if(sort_option==2) { // {S,s,r,rr,ss}
                    item.append(rsum+ssum);
                    item.append(ssum);
                    item.append(rsum);
                    item.append(sid);
                    for(auto ii : rs) item.append(ii);
                    for(auto ii : ss) item.append(ii);
                } else if(sort_option==-1) { // {S,r,s,-ss,rr}
                    item.append(rsum+ssum);
                    item.append(rsum);
                    item.append(ssum);
                    item.append(sid);
                    for(auto ii : ss) item.append(-ii);
                    for(auto ii : rs) item.append(ii);
                } else if(sort_option==-2) { // {S,s,r,rr,-ss}
                    item.append(rsum+ssum);
                    item.append(ssum);
                    item.append(rsum);
                    item.append(sid);
                    for(auto ii : rs) item.append(ii);
                    for(auto ii : ss) item.append(-ii);
                }
                item.append(fi);
                intg_vec.push_back(item);
            }

            sort_vec(intg_vec);
            unsigned long long int64 = 100000000000000;
            for(auto intg : intg_vec) {
                int64++;
                unsigned long long weight = int64;
                auto idx = intg.op(intg.nops()-1).op(0);
                for(int nc=0; nc<nCut; nc++) {
                    int j = iCuts[nc]-1;
                    if(idx.op(j)>1) {
                        weight += 100000000000000;
                        break;
                    }
                }
                i2w[idx] = weight;
                w2i[weight] = idx;
            }
        }
        
        if(WorkingDir.length()<1) WorkingDir = to_string(getpid());
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        system(("rm -rf "+job_dir).c_str());
        if(!dir_exists(job_dir)) system(("mkdir -p "+job_dir).c_str());
        
        ostringstream oss;
        for(auto eqn : eqns) {
            if(is_zero(eqn)) continue;
            auto cvs = collect_lst(eqn,F(w));
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
        oss << "  - reduce_user_defined_system:" << endl;
        oss << "      input_system: " << endl;
        oss << "        config: false" << endl;
        oss << "        files: [equations]" << endl;
        if(PIntegrals.nops()>0) {
            ostringstream oss2;
            int nn = PIntegrals.nops();
            for(int i=0; i<nn; i++) oss2 << Fout(PIntegrals.op(i)) << endl;
            ofstream pref_out(job_dir+"/preferred");
            pref_out << oss2.str() << endl;
            pref_out.close();
            oss << "      preferred_masters: preferred" << endl;
        }
        oss << "  - kira2file:" << endl;
        oss << "      target:" << endl;
        if(!using_uw) oss << "        - [F,integrals]" << endl;
        else oss << "        - [Tuserweight,integrals]" << endl;
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
    
    /**
     * @brief invoke kira program for reduction
     */
    void UKIRA::Run() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream cmd;
        cmd << "cd " << job_dir << " && kira " << KArgs << " --silent jobs >/dev/null 2>&1";
        system(cmd.str().c_str());
    }

    /**
     * @brief import kira result
     */
    void UKIRA::Import() {
        string job_dir = WorkingDir + "/" + to_string(ProblemNumber);
        ostringstream fn;
        if(!using_uw) fn << job_dir << "/results/F/kira_integrals.kira";
        else fn << job_dir << "/results/Tuserweight/kira_integrals.kira";
        auto strvec = file2strvec(fn.str());
        
        ex exL=0, exR=0;
        map<ex,int,ex_is_less> flags;
        lst exRs;
        for(auto intg : Integrals) flags[F(ProblemNumber,intg)] = 1;
        Rules.remove_all();
        for(auto line : strvec) {
            if(line.size()==0) {
                if(!is_zero(exL)) {
                    Rules.append(exL==exR);
                    flags[exL] = 0;
                    exRs.append(exR);
                }
                exL = exR = 0;
            } else if(is_zero(exL)) {
                exL -= Fin(line);
                if(!exL.match(F(w1,w2))) {
                    cout << line << endl;
                    throw Error("UKIRA::Import error found.");
                }
            } else {
                exR += Fin(line);
            }
        }
        if(!is_zero(exL)) {
            Rules.append(exL==exR);
            flags[exL] = 0;
            exRs.append(exR);
        }
        MIntegrals.remove_all();
        for(auto kv : flags) {
            if(kv.second!=0) MIntegrals.append(kv.first);
        }
        exset miset;
        find(exRs,F(w1,w2),miset);
        for(auto mi : miset) MIntegrals.append(mi);
        MIntegrals.sort();
        MIntegrals.unique();
        
    }

}
