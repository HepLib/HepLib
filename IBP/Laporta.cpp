/**
 * @file
 * @brief IBP with Leporta
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib::IBP {

    /**
     * @brief export integral to KIRA form
     */
    string Laporta::Fout(const ex & expr) {
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
    ex Laporta::Fin(const string & expr) {
        if(!using_uw) {
            string fstr = expr;
            string_replace_all(fstr,"[","("+to_string(ProblemNumber)+",{");
            string_replace_all(fstr,"]","})");
            return str2ex(fstr);
        } else {
            auto cpos = expr.find("*");
            if(cpos==string::npos) {
                if(expr=="0") return 0;
                throw Error("KIRA::Fin with 0 or * NOT Found.");
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
    void Laporta::Export() {

        if(Integrals.nops()<1) return;
        
        if(Round==0 || ibps.nops()<1) {
            _Integrals = Integrals;
            for(auto intg : Integrals) RIntegrals.append(F(ProblemNumber, intg));
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
                
                if(ISP.nops()<pdim) {
                    lst sps_ext;
                    for(auto it : External) {
                        for(auto ii : External) sps_ext.append(it*ii);
                    }
                    sps_ext.sort();
                    sps_ext.unique();
                    for(auto item : sps_ext) {
                        auto item2 = subs_all(item,Replacements);
                        if(is_zero(item-item2)) ISP.append(item);
                    }
                    ISP.sort();
                    ISP.unique();
                }
            }
            
            if(ISP.nops() > pdim) {
                cout << "ISP = " << ISP << endl;
                cout << "Propagators = " << Propagators << endl;
                throw Error("KIRA::Export: ISP more than Propagators.");
            }
            pdim = ISP.nops();
            
            lst sp2s, s2sp, ss;
            for(auto item : ISP) {
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
            if(s2p.nops() != pdim) throw Error("Laporta::Export: lsolve failed.");

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
        int pgDIM = Propagators.nops();
        int rmax = -1, smax = -1;
        int rrmax[pgDIM], ssmax[pgDIM];
        for(int i=0; i<pgDIM; i++) rrmax[i] = ssmax[i] = -1;
        for(auto integral : _Integrals) {
            if(integral.nops()!=pgDIM) throw Error("UKIRA::Export, integral dimension not match propagators.");
            int rr = 0;
            int ss = 0;
            for(int i=0; i<pgDIM; i++) {
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
            for(int i=0; i<pgDIM; i++) {
                if(_rrmax<rrmax[i]) _rrmax = rrmax[i];
                if(_ssmax<ssmax[i]) _ssmax = ssmax[i];
            }
            for(int i=0; i<pgDIM; i++) {
                rrmax[i] = _rrmax + rap;
                ssmax[i] = _ssmax + sap;
            }
            rmax += ra;
            smax += sa;
        } else {
            for(int i=0; i<pgDIM; i++) {
                rrmax[i] += rap;
                ssmax[i] += sap;
            }
            rmax += ra;
            smax += sa;
        }
        
        // generate IBP equations
        int as[pgDIM];
        for(int i=0; i<pgDIM; i++) as[i] = -ssmax[i];
        vector<vector<int>> asvec;
        while(true) {
            for(int i=0; i<pgDIM; i++) {
                if(as[i]+1>rrmax[i]) {
                    if(i+1 == pgDIM) goto done;
                    as[i] = -ssmax[i];
                } else {
                    as[i] = as[i] + 1;
                    break;
                }
            }
            int _rmax = 0, _smax = 0;
            for(int i=0; i<pgDIM; i++) {
                if(as[i]>0) _rmax += as[i];
                else _smax -= as[i];
            }
            if(_rmax>rmax || _smax>smax) continue;
            vector<int> asv;
            for(int i=0; i<pgDIM; i++) asv.push_back(as[i]);
            asvec.push_back(asv);
        }
        done: ;
        
        int nCut = Cuts.nops();
        bool hasCut = (nCut>1);
        int iCuts[nCut+1];
        for(int i=0; i<nCut; i++) iCuts[i] = ex_to<numeric>(Cuts.op(i)).to_int();
        auto verb = Verbose;
        Verbose = 0;
        auto eqns_result =
        GiNaC_Parallel(asvec.size(), [&](int idx)->ex {
            auto as = asvec[idx];
            exmap sol;
            for(int i=0; i<pgDIM; i++) sol[a(i)]=as[i];
            lst eqns;
            for(auto const & item : ibps) {
                auto ii = item.subs(sol);
                if(ii.is_zero()) continue;
                if(hasCut) {
                    exset fs;
                    find(ii, F(w), fs);
                    exmap repl;
                    for(auto fi : fs) {
                        lst ns = ex_to<lst>(fi.op(0));
                        for(auto ic : iCuts) {
                            int j = ic-1;
                            if(ns.op(j)<=0) {
                                repl[fi]=0;
                                break;
                            }
                        }
                    }
                    ii = ii.subs(repl);
                }
                if(ii.is_zero()) continue;
                eqns.append(ii);
            }
            return eqns;
        }, "Seeds");
        Verbose = verb;
        
        lst eqns;
        for(auto ilst : eqns_result) for(auto eqn : ilst) eqns.append(eqn);

        exset fs;
        for(auto eqn : eqns) find(eqn, F(w), fs);
        for(auto intg : _Integrals) fs.insert(F(intg));
        exvector intg_vec;
        for(auto fi : fs) {
            lst rs,ss;
            int sid=0, rsum=0, ssum=0, sn=0, rn=0;
            auto idx_lst = fi.op(0);
            for(int i=0; i<idx_lst.nops(); i++) {
                auto idx = ex_to<numeric>(idx_lst.op(i)).to_int();
                if(idx==0) continue;
                if(idx!=0) sid += std::pow(2,idx_lst.nops()-i-1);
                if(idx>0) {
                    rs.append(idx);
                    rsum += idx;
                    rn++;
                } else {
                    ss.append(idx);
                    ssum -= idx;
                    sn++;
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
        int rows = 0;
        int rows2 = 0;
        for(auto intg : intg_vec) {
            auto idx = intg.op(intg.nops()-1).op(0);
            bool isCut = false;
            for(auto ic : iCuts) {
                int j = ic-1;
                if(idx.op(j)>1) {
                    isCut = true;
                    break;
                }
            }
            if(isCut) {
                rows2--;
                i2w[idx] = rows2;
            } else {
                rows++;
                i2w[idx] = rows;
                w2i[rows] = idx;
            }
        }
        
        for(auto &kv : i2w) {
            if(kv.second<0) {
                kv.second = rows - kv.second;
                w2i[kv.second] = kv.first;
            }
        }
        int cols = rows - rows2;
        
        exvector eqns_cvs;
        for(auto eqn : eqns) {
            if(is_zero(eqn)) continue;
            auto cvs = collect_lst(eqn,F(w));
            lst cv_lst;
            for(auto cv : cvs) {
                cv_lst.append(lst{ (i2w[cv.op(1).op(0)]), cv.op(0) });
            }
            sort_lst(cv_lst,false);
            lst cv1, cv2;
            for(auto item : cv_lst) {
                cv1.append(item.op(0));
                cv2.append(item.op(1));
            }
            eqns_cvs.push_back(lst{cv1,cv2});
        }
        sort_vec(eqns_cvs,false);
        rows = eqns_cvs.size();
        
        //------
        cout << " - Fermating @ " << now(false) << " :> " << flush;
    
        Fermat fermat;
        fermat.Init();
        
        auto Variables = gather_symbols(ibps);
        ostringstream ss;
        for(auto j : Variables) ss << "&(J="<<j<<");" << endl;
        fermat.Execute(ss.str());
cout << ss.str() << endl;
        ss.clear();
        ss.str("");
        ss << "Array mat[" << rows << "," << cols << "] Sparse;" << endl;
        ostringstream oss;
        for(int r=0; r<rows; r++) {
            auto cv1 = eqns_cvs[r].op(0);
            auto cv2 = eqns_cvs[r].op(1);
            for(int j=0; j<cv1.nops(); j++) {
                ss << "mat["<<(rows-r)<<","<<cv1.op(j)<<"] := " << cv2.op(j) << ";" << endl;
cout << ss.str();
                fermat.Execute(ss.str());
                ss.clear();
                ss.str("");
            }
            oss << endl;
        }

        ss << "Redrowech([mat]);" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
           
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![mat" << endl;
        auto ostr = fermat.Execute(ss.str());
        fermat.Exit();
        cout << ostr << endl;
    }
    
    /**
     * @brief invoke kira program for reduction
     */
    void Laporta::Run() {
        
    }

    /**
     * @brief import kira result
     */
    void Laporta::Import() {
        
    }

}

