/**
 * @file
 * @brief Functions to perform partial fraction
 */

#include "HEP.h"

namespace HepLib {

    namespace {

        void A2FSave(string garfn, exvector air_vec, ex IntFs, vector<IBP*> &ibp_vec) {
            archive ar;
            
            ar.archive_ex(air_vec.size(), "air_vec_size");
            for(int i=0; i<air_vec.size(); i++) {
                ar.archive_ex(air_vec[i], ("air_"+to_string(i)).c_str());
            }
            
            ar.archive_ex(IntFs, "IntFs");
            
            ar.archive_ex(ibp_vec.size(), "ibp_vec_size");
            for(int i=0; i<ibp_vec.size(); i++) {
                ar.archive_ex(ibp_vec[i]->TO(), ("ibp_"+to_string(i)).c_str());
            }

            ofstream out(garfn);
            out << ar;
            out.close();
        }

        void A2FGet(string garfn, exvector &air_vec, ex &IntFs, vector<IBP*> &ibp_vec, int IBPmethod) {
            archive ar;
            ifstream in(garfn);
            in >> ar;
            in.close();
            
            map<string, ex> dict;
            for(int i=0; i<ar.num_expressions(); i++) {
                string name;
                ex res = ar.unarchive_ex(name, i);
                dict[name] = res;
            }
            
            auto size = ex_to<numeric>(dict["air_vec_size"]).to_int();
            if(size != air_vec.size()) throw Error("ApartIBP: ari_vec size is wrong!");
            air_vec.clear();
            air_vec.shrink_to_fit();
            air_vec.resize(size);
            for(int i=0; i<size; i++) air_vec[i] = dict["air_"+to_string(i)];
            
            size = ex_to<numeric>(dict["ibp_vec_size"]).to_int();
            ibp_vec.resize(size);
            for(int i=0; i<size; i++) {
                IBP* ibp;
                if(IBPmethod==0) ibp = new IBP();
                else if(IBPmethod==1) ibp = new FIRE();
                else if(IBPmethod==2) ibp = new KIRA();
                else if(IBPmethod==3) ibp = new UKIRA();
                else ibp = new IBP();
                ibp->FROM(dict["ibp_"+to_string(i)]);
                ibp_vec[i] = ibp;
            }
            
            IntFs = dict["IntFs"];
        }
        
    }
    
    /**
     * @brief convert F(ps, ns) to normal ex, ns is like FIRE convention
     * @param expr_in expression contains F
     * @return F(ps, ns) converted into normal expression
     */
     ex F2ex(const ex & expr_in) {
        ex ret = expr_in;
        ret = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(!e.has(F(w1,w2))) return e;
            else if(e.match(F(w1, w2))) {
                auto ps = e.op(0);
                auto ns = e.op(1);
                ex res = 1;
                for(int i=0; i<ps.nops(); i++) res *= pow(ps.op(i), ex(0)-ns.op(i));
                return res;
            } else return e.map(self);
        })(ret);
        return ret;
     }

    /**
     * @brief perform IBP reduction on the Aparted input
     * @param IBPmethod ibp method used, 0-No IBP, 1-FIRE, 2-KIRA
     * @param air_vec vector contains aparted input, ApartIRC will be call internally
     * @param aio AIOption for ApartIBP input
     * @return nothing returned, the input air_vec will be updated
     */
    void ApartIBP(exvector &air_vec, AIOption aio) {
        if(aio.smap.size()<1) aio.init_smap();
        int IBPmethod = aio.IBPmethod;
        
        lst lmom = ex_to<lst>(aio.Internal);
        lst emom = ex_to<lst>(aio.External);
        
        string wdir;
        if(aio.SaveDir != "") {
            if(IBPmethod==1) wdir = aio.SaveDir + "/FIRE";
            else if(IBPmethod==2) wdir = aio.SaveDir + "/KIRA";
            else if(IBPmethod==3) wdir = aio.SaveDir + "/UKIRA";
        } else {
            wdir = to_string(getpid());
            if(IBPmethod==1) wdir = wdir + "_FIRE";
            else if(IBPmethod==2) wdir = wdir + "_KIRA";
            else if(IBPmethod==3) wdir = wdir + "_UKIRA";
        }
        
        ex IntFs;
        vector<IBP*> ibp_vec;
        if(aio.SaveDir != "" && file_exists(aio.SaveDir+"/AIR2F.gar")) {
            if(Verbose > 1) cout << "  \\--Read AIR2F.gar" << flush; 
            A2FGet(aio.SaveDir+"/AIR2F.gar", air_vec, IntFs, ibp_vec, IBPmethod);
            for(auto ibp : ibp_vec) ibp->WorkingDir = wdir; // update working directory
            if(Verbose > 1) cout << " @ " << now(false) << endl; 
            goto AIR2F_Done;
        }
        
        if(aio.SaveDir != "") {
            if(file_exists(aio.SaveDir+"/AP.gar")) {
                if(Verbose > 1) cout << "  \\--Read AP.gar" << flush;
                garRead(air_vec, aio.SaveDir+"/AP.gar");
                if(Verbose > 1) cout << " @ " << now(false) << endl; 
                goto Apart_Done;
            } else system(("mkdir -p "+aio.SaveDir).c_str());
        } 
        
        if(true) {
            int av_size = air_vec.size();
            air_vec = GiNaC_Parallel(av_size, [air_vec,lmom] (int idx) {
                return collect_lst(air_vec[idx],lmom, o_normalF);
            }, "ApPre");
            
            exset vset;
            for(int i=0; i<av_size; i++) {
                for(auto cv : air_vec[i]) vset.insert(cv.op(1));
            }
            exvector vvec(vset.begin(), vset.end());
            vset.clear();
            
            if(true) {
                exmap v2ap;
                for(int i=0; i<vvec.size(); i++) v2ap[vvec[i]] = i;
                for(int i=0; i<av_size; i++) {
                    int n = air_vec[i].nops();
                    for(int j=0; j<n; j++) {
                        auto & v = air_vec[i][j][1]; // reference here
                        v = v2ap[v];
                    }
                }
                v2ap.clear();
            }
        
            auto ap_vec = GiNaC_Parallel(vvec.size(), [&vvec,lmom,emom,aio] (int idx) {
                auto air = vvec[idx];
                air = Apart(air,lmom,emom,aio.smap);
                air = collect_lst(air, ApartIR(w1,w2), o_normalF);
                return air;
            }, "Apart");
            vvec.clear();
            
            exset ap_set;
            for(auto cvs : ap_vec) for(auto cv : cvs) ap_set.insert(cv.op(1));
            exvector ap_ir_vec(ap_set.begin(), ap_set.end());
            ap_set.clear();
            auto ap_rules = ApartRules(ap_ir_vec); // including ApartIRC
            ap_ir_vec.clear();
            
            ap_vec = GiNaC_Parallel(ap_vec.size(), [&ap_vec,&ap_rules] (int idx) {
                auto cvs = ap_vec[idx];
                ex res = 0;
                for(auto cv : cvs) {
                    auto fi = ap_rules.find(cv.op(1));
                    if(fi==ap_rules.end()) res += cv.op(0) * cv.op(1);
                    else res += cv.op(0) * fi->second;
                }
                return res;
            }, "ApRule");
            ap_rules.clear();
            
            air_vec = GiNaC_Parallel(av_size, [&air_vec,&ap_vec] (int idx) {
                auto & cvs = air_vec[idx];
                ex res = 0;
                for(auto const & cv : cvs) {
                    int idx = ex_to<numeric>(cv.op(1)).to_int();
                    res += cv.op(0) * ap_vec[idx];
                }
                res = collect_ex(res, ApartIR(w1,w2), false, false, o_normalF);
                return res; // air_vec updated to ApartIR
            }, "ApPost");
            ap_vec.clear();
            
            if(aio.SaveDir != "") {
                garWrite(air_vec,aio.SaveDir+"/AP.gar");
                garRead(air_vec,aio.SaveDir+"/AP.gar");
            }
        }
        Apart_Done: ;
        
        if(true) {
        
            exvector AIR;
            if(true) {
                auto ret = GiNaC_Parallel(air_vec.size(), [&air_vec](int idx)->ex {
                    auto air = air_vec[idx];
                    exset airs;
                    find(air, ApartIR(w1,w2), airs);
                    lst ret;
                    for(auto item : airs) ret.append(item);
                    return ret;
                }, "ApIRC");
                exset intg;
                for(auto airs : ret) for(auto air : airs) intg.insert(air);
                AIR = exvector(intg.begin(), intg.end());
            }
            
            for(auto sp : aio.CSP) SP_map.erase(sp);
            // from here, Vector will be replaced by its name Symbol
            
            lst repls;
            auto sps = sp_map();
            for(auto kv : sps) repls.append(kv.first == kv.second);
                    
            lst loops, exts; // to match FIRE notation, not Vector, just Symbol
            for(auto li : lmom) {
                if(is_a<Vector>(li)) loops.append(ex_to<Vector>(li).name);
                else loops.append(li);
            }
            for(auto li : emom) {
                if(is_a<Vector>(li)) exts.append(ex_to<Vector>(li).name);
                else exts.append(li);
            }
        
            if(Verbose>0) cout <<  "  \\--Prepare " << WHITE << "IBP" << RESET << " reduction @ " << now(false) << flush;
            
            exmap AIR2F;
            std::map<ex, IBP*, ex_is_less> p2IBP;
            int pn=1;
            int ntot = AIR.size();
            for(int i=0; i<ntot; i++) {
                if(Verbose>0 && (((i+1)%1000)==0 || i+1==ntot)) cout << "\r                                 \r" << "  \\--Prepare " << WHITE << "IBP" << RESET << " reduction [" << (i+1) << "/" << ntot << "] @ " << now(false) << flush;
                auto const & ir = AIR[i];
                auto mat = ex_to<matrix>(ir.op(0));
                auto vars = ex_to<lst>(ir.op(1));
                lst pns;
                int nrow = mat.rows();
                for(int c=0; c<mat.cols(); c++) {
                    ex pc = 0;
                    for(int r=0; r<nrow-2; r++) pc += mat(r,c) * vars.op(r);
                    pc += mat(nrow-2,c);
                    pc = SP2sp(pc);
                    pns.append(lst{ pc, ex(0)-mat(nrow-1,c) }); // note the convension
                }
                sort_lst(pns); // use sort_lst
                
                int nCuts = aio.Cuts.nops();
                if(nCuts>0) {
                    ex cuts = aio.Cuts;
                    cuts = cuts.subs(SP_map,nopat);
                    if(aio.CutFirst) for(auto cut : cuts) pns.prepend(lst{ SP2sp(cut), 1 });
                    else for(auto cut : cuts) pns.append(lst{ SP2sp(cut), 1 });
                }
                
                lst props, ns;
                for(auto item : pns) {
                    props.append(item.op(0));
                    ns.append(item.op(1));
                }

                if(p2IBP.find(props)==p2IBP.end()) {
                    IBP* ibp;
                    if(IBPmethod==0) ibp = new IBP();
                    else if(IBPmethod==1) ibp = new FIRE();
                    else if(IBPmethod==2) ibp = new KIRA();
                    else if(IBPmethod==3) ibp = new UKIRA();
                    else {
                        ibp = new IBP();
                        IBPmethod = 0;
                    }
                    
                    p2IBP[props] = ibp;
                    ibp->Propagators = props;
                    ibp->Internal = loops;
                    ibp->External = exts;
                    ibp->Replacements = repls;
                    if(aio.ISP.nops()>0) for(auto item : aio.ISP) ibp->ISP.append(SP2sp(item));
                    if(aio.DSP.nops()>0) {
                        for(auto item : aio.DSP) {
                            lst sp = ex_to<lst>(item);
                            if(is_a<Vector>(sp.op(0))) sp.let_op(0) = (ex_to<Vector>(sp.op(0)).name);
                            if(is_a<Vector>(sp.op(1))) sp.let_op(1) = (ex_to<Vector>(sp.op(1)).name);
                            ibp->DSP.append(sp);
                        }
                    }
                    ibp->WorkingDir = wdir;
                    ibp->ProblemNumber = pn;
                    pn++;
                    if(nCuts>0) {
                        if(aio.CutFirst) for(int i=0; i<nCuts; i++) ibp->Cuts.append(i+1);
                        else for(int i=0; i<nCuts; i++) ibp->Cuts.append(nCuts-i);
                    }
                    ibp_vec.push_back(ibp);
                }
                IBP* ibp = p2IBP[props];
                ibp->Integrals.append(ns);

                AIR2F[AIR[i]] = F(ibp->ProblemNumber, ns);
            }
            if(Verbose>0) cout << endl;

            if(Verbose>0) cout << "  \\--Total Ints/Pros: " << WHITE << ntot << "/" << ibp_vec.size() << RESET << " @ " << now(false) << endl;
        
            if(true) {
                vector<IBP*> ibp_vec;
                for(auto ibp : ibp_vec) ibp_vec.push_back(ibp);
                auto int_fr = FindRules(ibp_vec, false, aio.UF);
                IntFs = int_fr.second;
                air_vec = GiNaC_Parallel(air_vec.size(), [&air_vec,&AIR2F,&int_fr] (int idx) {
                    auto air = air_vec[idx];
                    air = air.subs(AIR2F,nopat);
                    air = air.subs(int_fr.first,nopat);
                    air = collect_ex(air, F(w1,w2), false, false, o_normalF);
                    return air;
                }, "AIR2F");
                if(aio.SaveDir != "") A2FSave(aio.SaveDir+"/AIR2F.gar", air_vec, IntFs, ibp_vec);
            }
        }
        AIR2F_Done: ;
        
        MapFunction _F2ex([&ibp_vec,aio](const ex &e, MapFunction &self)->ex {
            if(!e.has(F(w1,w2))) return e;
            else if(e.match(F(w1,w2))) {
                int pn = ex_to<numeric>(e.op(0)).to_int();
                auto pso = ex_to<lst>(ibp_vec[pn-1]->Propagators);
                auto nso = ex_to<lst>(e.op(1));
                lst ps, ns;
                for(int i=0; i<pso.nops(); i++) {
                    if(!aio.keep0F && nso.op(i).is_zero()) continue;
                    ps.append(pso.op(i));
                    ns.append(nso.op(i));
                }
                return F(ps,ns);
            } else return e.map(self);
        });
        
        if(IBPmethod==0) { // no IBP reduction
            air_vec = GiNaC_Parallel(air_vec.size(), [&air_vec,&_F2ex](int idx)->ex {
                auto res = air_vec[idx];
                return _F2ex(res);
            }, "F2F");
            for(auto fp : ibp_vec) delete fp;
            if(aio.SaveDir == "") system(("rm -rf "+wdir).c_str());
            return;
        }

        exmap ibpRules; // IBP rules for problem pn
        if(aio.SaveDir != "" && file_exists(aio.SaveDir+"/Rules.gar")) {
            goto Rules_Done;
        }
        if(true) { 
            vector<IBP*> ibp_vec_re;
            if(true) {
                map<int,lst> pn_ints_map;
                for(auto item : IntFs) {
                    int pn = ex_to<numeric>(item.op(0)).to_int();
                    pn_ints_map[pn].append(item.op(1));
                }

                int nints = 0;
                for(auto pi : pn_ints_map) {
                    auto ibp = ibp_vec[pi.first-1];
                    ibp->Integrals = pi.second;
                    nints += ibp->Integrals.nops();
                    ibp_vec_re.push_back(ibp);
                }

                if(Verbose>0) cout << "  \\--Refined Ints/Pros: " << WHITE << nints << "/" << ibp_vec_re.size() << RESET << " @ " << now(false) << endl;
            } 
            
            if(IBPmethod==1) {
                if(GiNaC_Parallel_NB.find("Expo")==GiNaC_Parallel_NB.end() || GiNaC_Parallel_NB["Expo"]>100) GiNaC_Parallel_NB["Expo"] = 100;
                auto pRes = GiNaC_Parallel(ibp_vec_re.size(), [&ibp_vec_re](int idx)->ex {
                    ibp_vec_re[idx]->Export();
                    auto ret = lst{ ibp_vec_re[idx]->IsAlwaysZero ? 1 : 0, ibp_vec_re[idx]->Rules };
                    return ret;
                }, "Expo");
                for(int i=0; i<ibp_vec_re.size(); i++) {
                    ibp_vec_re[i]->IsAlwaysZero = (pRes[i].op(0)==1 ? true : false);
                    ibp_vec_re[i]->Rules = ex_to<lst>(pRes[i].op(1));
                }
                
                int nproc = aio.NIBP;
                if(nproc<1) nproc = CpuCores()/FIRE::Threads;
                int cproc = 0;
                if(nproc<2) nproc = 2;
                #pragma omp parallel for num_threads(nproc) schedule(dynamic, 1)
                for(int pi=0; pi<ibp_vec_re.size(); pi++) {
                    if(Verbose>1) {
                        #pragma omp critical
                        cout << "\r                                        \r" << "  \\--" << WHITE << "FIRE" << RESET << " Reduction [" << (++cproc) << "/" << ibp_vec_re.size() << "] " << flush;
                    }
                    ibp_vec_re[pi]->Run();
                }
                if(Verbose>1) cout << "@" << now(false) << endl;
                
                if(ibp_vec_re.size()>100) {
                    auto ret = GiNaC_Parallel(ibp_vec_re.size(), [&ibp_vec_re,wdir](int idx)->ex {
                        ibp_vec_re[idx]->Import();
                        return ibp_vec_re[idx]->TO();
                    }, "Impo");
                    for(int i=0; i<ibp_vec_re.size(); i++) ibp_vec_re[i]->FROM(ret[i]);
                } else {
                    cproc = 0;
                    for(auto item : ibp_vec_re) {
                        if(Verbose>1) cout << "\r                                        \r" << "  \\--" << WHITE << "FIRE" << RESET << " Import [" << (++cproc) << "/" << ibp_vec_re.size() << "] " << flush;
                        item->Import();
                    }
                    if(Verbose>1) cout << "@" << now(false) << endl;
                }
                //IBP::ReShare(ibp_vec_re);
                
                if(aio.SaveDir == "") system(("rm -rf "+wdir).c_str());
            } else if(IBPmethod==2 || IBPmethod==3) {
                for(auto ibp : ibp_vec_re) ibp->Reduce();
                if(aio.SaveDir == "") system(("rm -rf "+wdir).c_str());
            }
            
            // Find Rules in MIs
            exmap miRules = FindRules(ibp_vec_re, true, aio.UF).first;
            if(true) { // scope for ret
                if(aio.SaveDir != "") system(("mkdir -p "+aio.SaveDir+"/Rules").c_str());
                auto rules_vec = GiNaC_Parallel(ibp_vec_re.size(), [&ibp_vec_re,&miRules,&aio](int idx)->ex {
                    lst rules = ex_to<lst>(ibp_vec_re[idx]->Rules);
                    lst res;
                    for(auto ri : rules) res.append(lst { 
                        ri.op(0),  
                        ri.op(1).subs(miRules,nopat).subs(d==D,nopat)
                    });
                    for(auto mi : ibp_vec_re[idx]->MIntegrals) {
                        auto fi = miRules.find(mi);
                        if(fi!=miRules.end()) res.append(lst{ mi, fi->second });
                    }
                    auto pn = ibp_vec_re[idx]->ProblemNumber;
                    if(aio.SaveDir != "") {
                        garWrite(aio.SaveDir+"/Rules/"+to_string(pn)+".gar", res);
                        return 0;
                    } else return res;
                }, "FR2MI");
                if(aio.SaveDir != "") {
                    garWrite(aio.SaveDir+"/Rules.gar", 1);
                } else {
                    for(auto rs : rules_vec) {
                        for(auto ri : rs) if(ri.op(0)!=ri.op(1)) ibpRules[ri.op(0)] = ri.op(1);
                    }
                }
            }
        }
        Rules_Done: ;
        
        air_vec =
        GiNaC_Parallel(air_vec.size(), [&air_vec,&ibpRules,&_F2ex,&aio](int idx)->ex {
            ex res = air_vec[idx];
            exmap rules;
            if(aio.SaveDir != "") {
                exset fs;
                find(res, F(w1,w2), fs);
                for(auto fi : fs) {
                    auto pn = fi.op(0);
                    auto rs = garRead(aio.SaveDir+"/Rules/"+ex2str(pn)+".gar");
                    for(auto ri : rs) if(ri.op(0)!=ri.op(1)) rules[ri.op(0)] = ri.op(1);
                }
            } else rules = ibpRules;
            res = res.subs(rules,nopat);
            if(aio.pat.nops()>0) {
                auto cvs = collect_lst(res, aio.pat);
                res = 0;
                for(auto cv : cvs) {
                    auto c = cv.op(0);
                    auto v = cv.op(1);
                    if(aio.cv!=nullptr) {
                        auto _cv = aio.cv(c,v);
                        c = _cv.op(0);
                        v = _cv.op(1);
                    }
                    res += c * v;
                }
            }
            return _F2ex(res);
        }, "F2MI");
                                    
        for(auto fp : ibp_vec) delete fp;
    }
    
    /**
     * @brief perform IBP reduction on the Aparted input
     * @param IBPmethod ibp method used, 0-No IBP, 1-FIRE, 2-KIRA
     * @param air_vec vector contains aparted input, ApartIRC will be call internally
     * @param loops loop vectors
     * @param exts external vectors
     * @param cut_props cut propagators, default is { }
     * @param uf the function to compute UF polynomial
     * @return nothing returned, the input air_vec will be updated
     */
    void ApartIBP(exvector &air_vec, int IBPmethod, const lst & loops, const lst & exts, const lst & cut_props,
        std::function<lst(const IBP &, const ex &)> uf) {
                
        AIOption aio;
        aio.IBPmethod = IBPmethod;
        aio.Internal = loops;
        aio.External = exts;
        aio.Cuts = cut_props;
        if(cut_props.nops()>0) {
            for(auto p1 : loops) {
                for(auto p2 : loops) aio.CSP.append(SP(p1,p2));
                for(auto p2 : exts) aio.CSP.append(SP(p1,p2));
            }
            aio.CSP.sort();
            aio.CSP.unique();
        }
        for(auto li : loops) aio.smap[SP(li)] = 1;
        aio.UF = uf;
        ApartIBP(air_vec, aio);
    }

}
