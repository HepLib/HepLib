/**
 * @file
 * @brief IBP with FIRE
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib::IBP {

    /**
     * @brief Export start config intgral etc. files
     */
    void FIRE::Export() {
        ProblemDimension = Propagators.nops();
        if(Integrals.nops()<1) return;
        int pn = 0; // to avoid unsigned short overflow in FIRE
        int pdim = ProblemDimension;
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
            throw Error("FIRE::Export: Pairs failed.");
        }
        
        lst sp2s, s2sp, ss;
        for(auto item : Pairs) {
            symbol si;
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst eqns, iWFs;
        for(int i=0; i<pdim; i++) {
            auto eq = Propagators.op(i).expand().subs(iEpsilon==0); // drop iEpsilon
            eq = eq.subs(sp2s, subs_options::algebraic);
            eq = eq.subs(Replacements, subs_options::algebraic);
            if(eq.has(iWF(w))) throw Error("FIRE::Export, iWF used in eq.");
            eqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(eqns, ss);
        if(s2p.nops() != pdim) throw Error("FIRE::Export: lsolve failed.");

        vector<exmap> ibps;
        exvector IBPvec;
        lst ns0;
        for(int i=0; i<pdim; i++) ns0.append(0);
        for(auto loop : Internal) {
            lst dp_lst;
            for(int i=0; i<pdim; i++) { 
                auto s = ex_to<Symbol>(loop);
                dp_lst.append(Propagators.op(i).diff(s));
            } 
            
            for(auto iep : InExternal) { 
                exmap nc_map;
                for(int i=0; i<pdim; i++) { // diff on each propagator
                    auto ns = ns0;
                    ns.let_op(i) = ns.op(i)+1; // note the covention
                    auto tmp = dp_lst.op(i) * iep;
                    tmp = mma_collect(tmp, InExternal);
                    tmp = tmp.subs(Replacements, subs_options::algebraic);
                    tmp = tmp.subs(sp2s, subs_options::algebraic);
                    tmp = tmp.subs(s2p, subs_options::algebraic);
                    tmp = ex(0) - a(i+1)*tmp;

                    for(int j=0; j<pdim; j++) {
                        auto cj = tmp.coeff(iWF(j));
                        if(is_zero(cj)) continue;
                        auto cns = ns;
                        cns.let_op(j) = cns.op(j)-1; // note the covention
                        nc_map[cns] = nc_map[cns] + cj;
                    }
                    tmp = tmp.subs(iWF(w)==0); // constant term
                    if(!is_zero(tmp)) nc_map[ns] = nc_map[ns] + tmp;
                }
                
                if(is_zero(loop-iep)) nc_map[ns0] = nc_map[ns0] + VectorDimension;
                bool ok = false;
                for(auto nc : nc_map) {
                    if(!is_zero(nc.second)) {
                        ok = true;
                        IBPvec.push_back(nc.second.subs(D==d)); // replace D to d
                    }
                }
                if(ok) ibps.push_back(nc_map);
            }
        }

        auto Variables = gather_symbols(IBPvec);
        
        ostringstream start;

        start << "ExampleDimension[" << pn << "]=" << pdim << endl << endl;
        start << "ProblemNumber=" << pn << endl << endl;
        
        // .start - SBasisL
        if(Version==5) {
            PermutationsR(2, pdim, [pdim,pn,&start](const int *ns) {
                start << "SBasisL[" << pn << ",{";
                for(int i=0; i<pdim; i++) start << (ns[i]<1 ? -1 : 1) << (i<pdim-1 ? "," : "");
                start << "}]=0" << endl << endl;
            });
        }
        
        // .start - SBasis0L, SBasis0D && SBasis0C
        start << "SBasis0L[" << pn << "]=" << ibps.size() << endl << endl;
        ostringstream oss;
        for(int i=0; i<ibps.size(); i++) {
            start << "SBasis0D[" << pn << "," << (i+1) << "]=";
            lst items;
            for(auto kv : ibps[i]) {
                items.append(kv.first);
                if(Version==5) {
                    oss << "SBasis0C[" << pn << "," << (i+1) << "," << kv.first << "]=" << 
                    collect_common_factors(kv.second.normal()) << endl << endl;
                } else {
                    lst olst;
                    auto cv_lst = mma_collect_lst(kv.second, a(w));
                    for(auto item : cv_lst) {
                        auto cc = item.op(0);
                        auto cv = item.op(1);
                        cc = collect_common_factors(cc.normal());
                        if(is_zero(cc)) continue;
                        if(is_zero(cv-1)) cv=0;
                        else cv = cv.subs(a(w)==w);
                        olst.append(lst{cc, cv});
                    }
                    oss << "SBasis0C[" << pn << "," << (i+1) << "," << kv.first << "]=" << olst << endl << endl;
                }
            }
            start << items << endl << endl;
        }
        start << oss.str();
        
        // .start - SBasisS
        start << "SBasisS[" << pn << "]={{{";
        // 1,2,3,...
        for(int i=0; i<pdim; i++) start << (i+1) << (i<pdim-1 ? "," : "");
        start << "},{";
        // 1,1,1,...
        for(int i=0; i<pdim; i++) start << 1 << (i<pdim-1 ? "," : "");
        start << "},{";
        // 0,0,0,...
        for(int i=0; i<pdim; i++) start << 0 << (i<pdim-1 ? "," : "");
        start << "}}}" << endl << endl;
        
        // .start - SBasisR
        lst Rlst;
        Rlst.append(lst{});
        for(int i=0; i<pdim; i++) {
            let_op_append(Rlst, 0, -1);
        }
        for(auto lpi : Internal) {
            vector<int> ns_vec;
            lst ns0;
            for(int i=0; i<pdim; i++) ns0.append(1);
            for(int i=0; i<pdim; i++) {
                if(Propagators.op(i).has(lpi)) ns0.let_op(i) = -1;
                else ns_vec.push_back(i);
            }
            for(int n=0; n<std::pow(2,ns_vec.size()); n++) {
                int cn = n;
                lst ns1 = ns0;
                for(int j=0; j<ns_vec.size(); j++) {
                    if((cn%2)==1) ns1.let_op(ns_vec[j]) = -1;
                    cn /= 2;
                }
                Rlst.append(ns1);
            } 
        }
        
        // handle Cut Propagators
        if(Cuts.nops()>0) {
            //Rlst.remove_all();  // TODO: check this one
            for(auto cx : Cuts) {
                int ci = ex_to<numeric>(cx-1).to_int(); // start from 1 in Cuts
                lst ns0;
                for(int i=0; i<pdim; i++) ns0.append(1);
                ns0.let_op(ci) = -1;
                for(int n=0; n<std::pow(2,pdim-1); n++) {
                    int cn = n;
                    lst ns1 = ns0;
                    for(int j=0; j<pdim; j++) {
                        if(ci==j) continue;
                        if((cn%2)==1) ns1.let_op(j) = -1;
                        cn /= 2;
                    }
                    Rlst.append(ns1);
                } 
            }
        }
        
        Rlst.sort();
        Rlst.unique();
        
        for(auto iR : Rlst) {
            start << "SBasisR[" << pn << ",{";
            for(int i=0; i<pdim; i++) start << iR.op(i) << (i<pdim-1 ? "," : "");
            start << "}]=True" << endl << endl;
        }
        
        // .start - Others
        start << "SBasisRL[" << pn << "]=0" << endl << endl;
        start << "HPI[" << pn << "]={}" << endl << endl;
        
        string sss = start.str();
        string_replace_all(sss, "=", " = ");
        string_replace_all(sss, ",", ", ");
        
        string spn = to_string(ProblemNumber);
        if(!dir_exists(WorkingDir)) system(("mkdir -p " + WorkingDir).c_str());
        
        // .config
        ostringstream config;
        if(Version>5) config << "#compressor none" << endl;
        config << "#threads " << Threads << endl;
        //config << "#fthreads 4" << endl;
        //config << "#fermat fer64" << endl;
        config << "#variables ";
        bool first = true;
        for(auto v : Variables) { 
            if(!islower(ex_to<symbol>(v).get_name()[0])) 
                throw Error("Reduce: Fermat requires a name must begin with a lower case letter.");
            config << (first ? "" : ",") << v; 
            first=false; 
        }
        config << endl;
        config << "#database db" << ProblemNumber << endl;
        if(Version>5 && pos_pref!=1) config << "#pos_pref "<< pos_pref << endl;
        if(Version==5) config << "#bucket 20" << endl;
        config << "#start" << endl;
        config << "#problem " << pn << " " << ProblemNumber << ".start" << endl;
        if(mi_pref.nops()>0) {
            ostringstream oss;
            oss << "{";
            int nn = mi_pref.nops();
            for(int i=0; i<nn; i++) oss << "{" << pn << "," << mi_pref.op(i) << (i<nn-1 ? "}," : "}");
            oss << "}";
            ofstream pref_out(WorkingDir+"/"+spn+".pref");
            pref_out << oss.str() << endl;
            pref_out.close();
            config << "#preferred " << ProblemNumber << ".pref" << endl;
        }
        config << "#integrals " << ProblemNumber << ".intg" << endl;
        config << "#output " << ProblemNumber << ".tables" << endl;
        
        // *.intg
        ostringstream intg;
        intg << "{";
        for(int i=0; i<Integrals.nops(); i++) {
            intg << "{" << pn << "," << Integrals[i] << (i<Integrals.nops()-1 ? "}," : "}");
        }
        intg << "}" << endl;
        
        if(WorkingDir.length()<1) WorkingDir = to_string(getpid());
        if(!dir_exists(WorkingDir)) system(("mkdir -p "+WorkingDir).c_str());
        ofstream start_out(WorkingDir+"/"+spn+".start");
        start_out << sss << endl;
        start_out.close();
        
        ofstream config_out(WorkingDir+"/"+spn+".config");
        config_out << config.str() << endl;
        config_out.close();
        
        ofstream intg_out(WorkingDir+"/"+spn+".intg");
        intg_out << intg.str() << endl;
        intg_out.close();
    
    }
    
    /**
     * @brief Run FIRE reduction 
     */
    void FIRE::Run() {
        ostringstream cmd;
        cmd << "cd " << WorkingDir << " && $(which FIRE" << Version << ")";
        if(Version>5) cmd << " -silent";
        cmd << " -c " << ProblemNumber << " >/dev/null";
        system(cmd.str().c_str());
        system(("rm -rf "+WorkingDir+"/db"+to_string(ProblemNumber)).c_str());
        if(!file_exists(WorkingDir + "/" + to_string(ProblemNumber) + ".tables")) {
            system(cmd.str().c_str());
            system(("rm -rf "+WorkingDir+"/db"+to_string(ProblemNumber)).c_str());
        }
    }
    
    /**
     * @brief Import tables  
     */
    void FIRE::Import() {
        string spn = to_string(ProblemNumber);
        ifstream ifs(WorkingDir+"/"+spn+".tables");
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        
        string_replace_all(ostr, "\"", "");
        
        Parser tp;
        auto tp_lst = tp.Read(ostr);
        exmap id2F;
        for(auto item : tp_lst.op(1)) {
            if(!is_zero(item.op(1).op(0))) throw Error("FIRE::Reduce: pn is NOT 0.");
            id2F[item.op(0)] = F(ProblemNumber, item.op(1).op(1));
        }
        
        for(auto item : tp_lst.op(0)) {
            ex left = item.op(0).subs(id2F);
            ex right = 0;
            for(auto it : item.op(1)) {
                right += it.op(0).subs(id2F) * it.op(1);
            }
            right = right.subs(d==D); // replace d back to D
            if(is_zero(left-right)) MasterIntegrals.append(left);
            else Rules.append(left==right);
        }
        MasterIntegrals.sort();
        MasterIntegrals.unique();
     
        // handle Cuts not equal 1, using #preferred
        if(Cuts.nops()>0 && mi_pref.nops()<1) {
            lst cIntegrals;
            auto mis = MasterIntegrals;
            MasterIntegrals.remove_all();
            for(auto item : mis) {
                lst mi = ex_to<lst>(item.op(1));
                bool isOK = true;
                for(auto cx : Cuts) {
                    if(!is_zero(mi.op(ex_to<numeric>(cx).to_int()-1)-1)) {
                        isOK = false;
                        break;
                    }
                }
                if(isOK) MasterIntegrals.append(item);
                else {
                    cIntegrals.append(mi);
                    auto pi = mi;
                    vector<int> ipos, ineg;
                    for(int i=0; i<mi.nops(); i++) {
                        bool isCut = false;
                        for(auto cx : Cuts) {
                            if(is_zero(cx-1-i)) {
                                isCut = true;
                                break;
                            }
                        }
                        if(isCut) pi.let_op(i)=1;
                        else if(mi.op(i)<=0) ineg.push_back(i);
                        else ipos.push_back(i);
                    }
                    
                    lst mi_pref_tmp;
                    int max = 2;
                    ex total = pow(numeric(max), ineg.size());
                    for(numeric in=0; in<total; in++) {
                        auto cin = in;
                        auto pi2 = pi;
                        for(int i=0; i<ineg.size(); i++) {
                            int re = mod(cin,max).to_int();
                            pi2.let_op(ineg[i]) = ex(0)-re;
                            cin = (cin-re)/max;
                        }
                        mi_pref_tmp.append(pi2);
                    }
                    
                    int max2 = 2*max+1;
                    total = pow(numeric(max2), ipos.size());
                    for(auto pi : mi_pref_tmp) {
                        for(numeric in=0; in<total; in++) {
                            auto cin = in;
                            auto pi2 = pi;
                            for(int i=0; i<ipos.size(); i++) {
                                int re = mod(cin,max2).to_int();
                                pi2.let_op(ipos[i]) = re-max;
                                cin = (cin-re)/max2;
                            }
                            mi_pref.append(pi2);
                        }
                    }
                }
            }
            
            // Reduce again
            mi_pref.sort();
            mi_pref.unique();
            if(mi_pref.nops()>0) {
                FIRE fire;
                fire.Propagators = Propagators;
                fire.Internal = Internal;
                fire.External = External;
                fire.Replacements = Replacements;
                fire.Pairs = Pairs;
                fire.ProblemNumber = ProblemNumber;
                fire.Cuts = Cuts;
                fire.Integrals = cIntegrals;
                fire.mi_pref = mi_pref;
                mi_pref.remove_all();
                fire.WorkingDir = WorkingDir + "_C"+to_string(ProblemNumber);
                fire.Reduce();
                system(("rm -rf " + fire.WorkingDir).c_str());
                for(auto item : fire.MasterIntegrals) MasterIntegrals.append(item);
                auto rules = Rules;
                Rules.remove_all();
                for(auto r : rules) Rules.append(r.op(0)==subs_naive(r.op(1),fire.Rules));
                for(auto r : fire.Rules) Rules.append(r);
                MasterIntegrals.sort();
                MasterIntegrals.unique();
            }
        }
    }
    
}
