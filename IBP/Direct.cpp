/**
 * @file
 * @brief IBP with KIRA
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib {

    static matrix RowReduce(matrix mat) {
        static map<pid_t, Fermat> fermat_map;
        static int v_max = 0;

        auto pid = getpid();
        if((fermat_map.find(pid)==fermat_map.end())) { // init section
            fermat_map[pid].Init();
            v_max = 0;
        }
        Fermat &fermat = fermat_map[pid];
        
        lst rep_vs;
        ex tree = mat;
        for(const_preorder_iterator i = tree.preorder_begin(); i != tree.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e) || e.match(a(w))) rep_vs.append(e);
        }
        rep_vs.sort();
        rep_vs.unique();
        //sort_lst(rep_vs);

        exmap v2f;
        symtab st;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        int nrow = mat.rows();
        int ncol = mat.cols();
        
        ss << "Array m[" << nrow << "," << ncol+1 << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m]:=[(";
        for(int c=0; c<ncol; c++) {
            for(int r=0; r<nrow; r++) {
                ss << mat(r,c).subs(iEpsilon==0).subs(v2f) << ",";
            }
        }
        for(int r=0; r<nrow; r++) ss << "0,";
        ss << ")];" << endl;
        ss << "Redrowech([m]);" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");

        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![m" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
        
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("Direc::Export, last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        string_trim(ostr);
        
        ostr.erase(0, ostr.find(":=")+2);
        string_replace_all(ostr, "[", "{");
        string_replace_all(ostr, "]", "}");
        Parser fp(st);
        matrix mr(nrow, ncol);
        auto res = fp.Read(ostr);
        for(int r=0; r<nrow; r++) {
            auto cur = res.op(r);
            for(int c=0; c<ncol; c++) mr(r,c) = cur.op(c);
        }
        return mr;
    }

    /**
     * @brief Export input data for KIRA
     */
    void Direct::Export() {

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
            throw Error("Direct::Export: #(ISP) > #(Propagators).");
        }
        
        lst sp2s, s2sp, ss;
        int _pic=0;
        for(auto item : ISP) {
            _pic++;
            Symbol si("P"+to_string(_pic));
            ss.append(si);
            sp2s.append(w*item==w*si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst leqns;
        for(int i=0; i<ISP.nops(); i++) { // note NOT pdim
            auto eq = Propagators.op(i).expand().subs(iEpsilon==0); // drop iEpsilon
            eq = eq.subs(sp2s);
            eq = eq.subs(Replacements);
            if(eq.has(iWF(w))) throw Error("Direct::Export, iWF used in eq.");
            leqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(leqns, ss);
        if(s2p.nops() != ISP.nops()) throw Error("Direct::Export: lsolve failed.");
        
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
                ibp -= (a(i)+Shift[i]) * F(ns) * dp;
            }
            
            ibp = ibp * iep;
            ibp = ibp.expand();
            ibp = ibp.subs(sp2s);
            ibp = ibp.subs(Replacements);
            ibp = ibp.subs(s2p);
            
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
        
        // zero sector
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
                //Rlst.append(ns1);
            }
        }
        
        // Lee Zero Sector
        if(true) {
            //IsAlwaysZero = true;
            lst ns0;
            for(int i=0; i<pdim; i++) ns0.append(1);
            size_t tot = std::pow(2LL,pdim);
            for(size_t n=0; n<tot; n++) {
                int cn = n;
                lst ns1 = ns0;
                lst sector = ns0;
                for(int j=0; j<pdim; j++) {
                    if((cn%2)==1) { ns1.let_op(j) = -1; sector.let_op(j) = 0; }
                    cn /= 2;
                }
                //if(IsZero(sector)) Rlst.append(ns1);
                //else if(IsAlwaysZero) IsAlwaysZero = false;
            } 
            if(IsAlwaysZero) {
                lst ws;
                int pdim = Propagators.nops();
                for(int i=0; i<pdim; i++) ws.append(wild(i));
                Rules.append(F(ProblemNumber, ws)==0);
                return;
            }
        }
        
        // all sectors
        auto ibps_o = ibps;
        ibps.remove_all();
cout << ibps_o << endl;
exit(0);
        for(auto ibp : ibps_o) {
            for(int s=-3; s<=3; s++) {
                ibps.append(ibp.subs(a(w)==a(w)+s));
            }
        }
        
        size_t ltot = std::pow(3L,pdim); // see below
        auto ret = GiNaC_Parallel(ltot, [&](int idx)->ex {
            lst res;
            auto li = idx;
            auto cli = li;
            vector<int> csec;
            lst cibps = ibps;
            for(int i=0; i<pdim; i++) {
                auto im = (cli % 3)-1;
                if(im==0) cibps = ex_to<lst>(subs(cibps,a(i)==0));
                csec.push_back(im);
                cli /= 3;
            }
            cibps.sort();
            cibps.unique();
            map<int,exvector> sum_fs;
            
            exset fset;
            find(cibps,F(w),fset);

            for(auto fi : fset) {
                int sum = 0;
                auto ns = fi.op(0).subs(a(w)==0);
                for(int i=0; i<pdim; i++) {
                    auto ni = ns.op(i);
                    if(csec[i]==-1) sum -= ex2int(ni);
                    else if(csec[i]==1) sum += ex2int(ni);
                    else sum += abs(ex2int(ni));
                }
                sum_fs[sum].push_back(fi);
            }

            lst sum_tot;
            for(auto kv : sum_fs) sum_tot.append(lst{kv.first, kv.second.size()});
            sort_lst(sum_tot,false);

            matrix mat(cibps.nops(), fset.size());
            exvector fvec;
            int ccol = 0;
            for(auto st : sum_tot) {
                auto fs = sum_fs[ex2int(st.op(0))];
                sort_vec(fs);
                for(int cc=0; cc<fs.size(); cc++) {
                    fvec.push_back(fs[cc]);
                    int row = 0;
                    for(auto ibp : cibps) {
                        mat(row,ccol+cc) = ibp.coeff(fs[cc]);
                        row++;
                    }
                }
                ccol += fs.size();
            }

            auto mr = RowReduce(mat);
            
            for(int r=0; r<mat.rows(); r++) {
                int c1 = -1, c = -1;
                for(auto ti : sum_tot) {
                    int ct = ex2int(ti.op(1));
                    for(int k=0; k<ct; k++) {
                        c++;
                        if(mr(r,c).is_zero()) continue;
                        else if(c1==-1) c1=c;
                        else goto next_row;
                    }
                    if(c1!=-1) {
                        ex sol = 0;
                        for(int ci=c1+1; ci<mat.cols(); ci++) sol -= mr(r,ci)/mr(r,c1)*fvec[ci];
                        auto ns0 = fvec[c1].op(0).subs(a(w)==0);
                        exmap aSH;
                        for(int i=0; i<ns0.nops(); i++) {
                            if(csec[i]==0) continue;
                            if(!ns0.op(i).is_zero()) aSH[a(i)]=a(i)-ns0.op(i);
                        }
                        sol = sol.subs(aSH);

                        exset fs;
                        find(sol,F(w),fs);
                        vector<vector<int>> fns;
                        for(auto fi : fs) {
                            fi = fi.subs(a(w)==0);
                            vector<int> ns;
                            for(auto it : fi.op(0)) ns.push_back(0-ex2int(it));
                            fns.push_back(ns);
                        }
                        Condition cond;
                        for(int i=0; i<csec.size(); i++) {
                            if(csec[i]==-1) {
                                int min = 100000;
                                for(auto ns : fns) {
                                    if(ns[i]<min) min = ns[i];
                                }
                                cond.cs.push_back(make_pair(-1,min));
                            } else if(csec[i]==1) {
                                int max = -100000;
                                for(auto ns : fns) {
                                    if(ns[i]>max) max = ns[i];
                                }
                                cond.cs.push_back(make_pair(1,max));
                            } else {
                                cond.cs.push_back(make_pair(0,ex2int(ns0.op(i))));
                            }
                        }

                        res.append(lst{cond.cs2ex(),sol});
                        goto next_row;
                    }
                }
                next_row: ;
            }
            return res;
        }, "IBP");
        
        ConSolVec.clear();
        for(auto its : ret) {
        for(auto item : its) {
            Condition cond;
            cond.ex2cs(item.op(0));
            ConSolVec.push_back(make_pair(cond,item.op(1)));
        }}
        
    }
        
        
        
    
    
    /**
     * @brief invoke kira program for reduction
     */
    void Direct::Run() {
    
        // improve here, try by level
        exmap sols;
        auto intgs = Integrals;
        while(true) {
            exset fs;
            for(auto intg : intgs) {
                if(key_exists(sols, F(intg))) continue;
                for(auto cs : ConSolVec) {
                    if(cs.first.IsOK(intg)) {
                        exmap as;
                        for(int i=0; i<intg.nops(); i++) as[a(i)] = intg.op(i);
                        try { // improve here
                            //if(numer_denom(cs.second).op(1).subs(as).is_zero()) continue;
                            auto sol = cs.second.subs(as);
                            find(sol, F(w), fs);
                            sols[F(intg)] = sol;
                            break;
                        } catch(...) { }
                    }
                }
            }
            intgs.remove_all();
            for(auto fi : fs) intgs.append(fi.op(0));
            if(intgs.nops()<1) break;
        }
                
        // subs by level
        map<int,exmap> lvl_sol;
        lst lvls;
        for(auto sol : sols) {
            ex sum = 0;
            for(auto ni : sol.first.op(0)) sum += (ni<0 ? -ni : ni);
            lvl_sol[ex2int(sum)][sol.first] = sol.second;
            lvls.append(ex2int(sum));
        }
        lvls.sort();
        lvls.unique();
        //sort_lst(lvls);

        for(int i=0; i<lvls.nops(); i++) {
            int si = ex2int(lvls.op(i));
            for(auto kv : lvl_sol[si]) {
            for(int j=0; j<i; j++) {
                int sj = ex2int(lvls.op(j));
                lvl_sol[si][kv.first] = kv.second.subs(lvl_sol[sj]);
            }}
        }
                
        for(auto intg : Integrals) {
            ex fi = F(intg);
            ex sum = 0;
            for(auto ni : fi.op(0)) sum += (ni<0 ? -ni : ni);
            fi = fi.subs(lvl_sol[ex2int(sum)]);
            Rules.append(F(intg)==fi);
        }
        
    }

    /**
     * @brief import kira result
     */
    void Direct::Import() {
        
    }

}
