/**
 * @file
 * @brief IBP with Solver
 */

#include "BASIC.h"
#include "IBP.h"
#include <cmath>

namespace HepLib {

    // return MR & T, such that MR = T.M0
    static pair<matrix,matrix> RowReduce(matrix mat) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
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
        if(nrow>0) ss << "0";
        for(int r=1; r<nrow; r++) ss << ",0";
        ss << ")];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        bool sparse = false;
        if(nrow>1100) sparse = true;
        if(sparse) fermat.Execute("Sparse([m]);");
        fermat.Execute("Redrowech([m],[t]);");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        matrix mr(nrow, ncol);
        if(true) { // read matrix m
            if(sparse) ss << "![m]" << endl;
            else ss << "![m" << endl;
            auto ostr = fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            
            // make sure last char is 0
            if(ostr[ostr.length()-1]!='0') throw Error("RowReduce, last char is NOT 0.");
            ostr = ostr.substr(0, ostr.length()-1);
            string_trim(ostr);
            
            ostr.erase(0, ostr.find(":=")+2);
            size_t sn = ostr.length();
            char lc;
            for(size_t i=0; i<sn; i++) {
                char & c = ostr[i];
                if(c=='[') {
                    c = '{';
                    if(sparse && i>0 && lc=='}') ostr[i-1] = ',';
                } else if(c==']') c = '}';
                else if(sparse && (c==' '||c=='\t'||c=='\n'||c=='\r')) continue;
                lc = c;
            }
            
            Parser fp(st);
            auto res = fp.Read(ostr);
            if(sparse) {
                for(auto const & item : res) {
                    int r = -1;
                    for(auto const & it : item) {
                        if(r==-1) r = ex2int(it)-1;
                        else {
                            int c = ex2int(it.op(0))-1;
                            mr(r,c) = it.op(1);
                        }
                    }
                }
            } else {
                for(int r=0; r<nrow; r++) {
                    auto cur = res.op(r);
                    for(int c=0; c<ncol; c++) mr(r,c) = cur.op(c);
                }
            }
        }
        matrix mt(nrow, nrow);
        if(true) { // read matrix t
            if(sparse) ss << "![t]" << endl;
            else ss << "![t" << endl;
            auto ostr = fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            
            // make sure last char is 0
            if(ostr[ostr.length()-1]!='0') throw Error("Solver::Export, last char is NOT 0.");
            ostr = ostr.substr(0, ostr.length()-1);
            string_trim(ostr);
            
            ostr.erase(0, ostr.find(":=")+2);
            size_t sn = ostr.length();
            char lc;
            for(size_t i=0; i<sn; i++) {
                char & c = ostr[i];
                if(c=='[') {
                    c = '{';
                    if(sparse && i>0 && lc=='}') ostr[i-1] = ',';
                } else if(c==']') c = '}';
                else if(sparse && (c==' '||c=='\t'||c=='\n'||c=='\r')) continue;
                lc = c;
            }
            Parser fp(st);
            auto res = fp.Read(ostr);
            if(sparse) {
                for(auto const & item : res) {
                    int r = -1;
                    for(auto const & it : item) {
                        if(r==-1) r = ex2int(it)-1;
                        else {
                            int c = ex2int(it.op(0))-1;
                            mt(r,c) = it.op(1);
                        }
                    }
                }
            } else {
                for(int r=0; r<nrow; r++) {
                    auto cur = res.op(r);
                    for(int c=0; c<nrow; c++) mt(r,c) = cur.op(c);
                }
            }
        }
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m],[t]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        return make_pair(mr,mt);
    }
    
    // both row and column start from 0
    typedef map<int,map<int,ex>> SparseMatrix;
    static void RowReduce(SparseMatrix & smat, int pn) {
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        
        lst rep_vs;
        for(auto const & rv : smat) for(auto const & cv : rv.second) {
            ex tree = cv.second;
            for(const_preorder_iterator i = tree.preorder_begin(); i != tree.preorder_end(); ++i) {
                auto e = (*i);
                if(is_a<symbol>(e) || e.match(a(w))) rep_vs.append(e);
            }
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
        
        ostringstream oss;
        int nrow=-1, ncol=-1;
        int rcur = 0, rtot = smat.size();
        oss << "[";
        for(auto const & rv : smat) {
            rcur++;
            if(rv.second.size()<1) continue;
            int r = rv.first;
            if(r>nrow) nrow = r;
            oss << "[" << r+1 << ",";
            int ccur = 0, ctot = rv.second.size();
            for(auto const & cv : rv.second) {
                ccur++;
                int c = cv.first;
                if(c>ncol) ncol = c;
                oss << "[" << c+1 << "," << cv.second.subs(iEpsilon==0).subs(v2f) << "]";
                if(ccur!=ctot) oss << ",";
                else oss << "]";
            }
            if(rcur!=rtot) oss << "`" << endl;
            else oss << "];" << endl;
        }
        nrow = nrow+1;
        ncol = ncol+1;
        
        ss << "Array m[" << nrow << "," << ncol+1 << "] Sparse;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m]:=" << oss.str();
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        if(pn>0) fermat.Execute("Redrowech([m],,,"+to_string(pn)+");");
        else fermat.Execute("Redrowech([m]);");
        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        smat.clear();
        if(true) { // read matrix m
            ss << "![m]" << endl;
            auto ostr = fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            
            // make sure last char is 0
            if(ostr[ostr.length()-1]!='0') throw Error("RowReduce, last char is NOT 0.");
            ostr = ostr.substr(0, ostr.length()-1);
            string_trim(ostr);
            
            ostr.erase(0, ostr.find(":=")+2);
            size_t sn = ostr.length();
            char lc;
            for(size_t i=0; i<sn; i++) {
                char & c = ostr[i];
                if(c=='[') {
                    c = '{';
                    if(i>0 && lc=='}') ostr[i-1] = ',';
                } else if(c==']') c = '}';
                else if(c==' '||c=='\t'||c=='\n'||c=='\r') continue;
                lc = c;
            }
            
            Parser fp(st);
            auto res = fp.Read(ostr);
            for(auto const & item : res) {
                int r = -1;
                for(auto const & it : item) {
                    if(r==-1) r = ex2int(it)-1;
                    else {
                        int c = ex2int(it.op(0))-1;
                        smat[r][c] = it.op(1);
                    }
                }
            }
        }
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
    }
    
    void Solver::SolveSparse(const ex & sector, const map<int,int> & n2n) {
        int pdim = Propagators.nops();
        lst ibps = IBPs;
        
        // add more eqn to ibps
        if(true) {
            if(true) {
                int n = ibps.nops();
                for(int i=0; i<pdim; i++) {
                    for(int j1=0; j1<n; j1++) for(int j2=0; j2<n; j2++) {
                        ibps.append(ibps.op(i).subs(lst{a(j1)==a(j1)+1}).subs(lst{a(j2)==a(j2)-1}));
//                        ibps.append(ibps.op(i).subs(lst{a(j1)==a(j1)+1}).subs(lst{a(j2)==a(j2)+1}));
//                        ibps.append(ibps.op(i).subs(lst{a(j1)==a(j1)-1}).subs(lst{a(j2)==a(j2)-1}));
                    }
                }
            } else {
                int n = ibps.nops();
                for(int i=0; i<pdim; i++) {
                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+1));
//                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+2));
//                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+3));
//                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+4));
                }

                n = ibps.nops();
                for(int i=0; i<pdim; i++) {
                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-1));
//                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-2));
//                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-3));
//                    for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-4));
                }
            }
        } else {
            exset fs;
            find(ibps, F(w), fs);
            int n = ibps.nops();
            for(auto fi : fs) {
                exmap a2a;
                for(int i=0; i<fi.op(0).nops(); i++) a2a[a(i)] = fi.op(0).op(i);
                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a2a));
            }
            for(auto fi : fs) {
                exmap a2a;
                for(int i=0; i<fi.op(0).nops(); i++) a2a[a(i)] = fi.op(0).op(i)+1;
                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a2a));
            }
            for(auto fi : fs) {
                exmap a2a;
                for(int i=0; i<fi.op(0).nops(); i++) a2a[a(i)] = fi.op(0).op(i)-1;
                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a2a));
            }
        }
        
        // some a's replaced by integer
        vector<int> nfix(pdim);
        for(int i=0; i<pdim; i++) nfix[i] = 0;
        exmap a2n;
        for(auto const & kv : n2n) {
            nfix[kv.first] = 1;
            int k = kv.first;
            if(k<0) k = pdim + k;
            a2n[a(k)] = kv.second;
        }
        if(a2n.size()>0) {
            int nibps = ibps.nops();
            for(int i=0; i<nibps; i++) ibps.let_op(i) = ibps.op(i).subs(a2n);
        }
        
        lst cons;
        map<ex,lst,ex_is_less> cons_vec;
        while(ibps.nops()>0) {
            ibps.sort();
            ibps.unique();
            
            exset fset;
            find(ibps,F(w),fset);
            map<int,exvector> fmap_vec;
            for(auto fi : fset) { // sort fset
                int psum = 0, nsum = 0;
                auto ns = fi.op(0);
                for(int i=0; i<pdim; i++) {
                    auto ni = ns.op(i).subs(a(i)==0,nopat);
                    if(sector.op(i)==1) {
                        if(nfix[i]==1 && ni<=0) psum -= 100000;
                        else psum += ex2int(ni)-1; // -1 from FIRE
                    } else {
                        if(nfix[i]==1 && ni>0) nsum += ex2int(ni);
                        else nsum -= ex2int(ni);
                    }
                }
                int key = psum+nsum;
                fmap_vec[-key].push_back(fi); // - to reverse the order
            }
            
            map<ex,int,ex_is_less> f2n;
            int next_gi = 0;
            exvector fvec(fset.size());
            int nr = 1, nc = 0;
            for(auto const & ev : fmap_vec) {
                if(next_gi<1) next_gi = ev.second.size();
                for(auto const & v : ev.second) {
                    fvec[nc] = v;
                    f2n[v] = nc;
                    nc++;
                }
                nr++;
            }
            nr = ibps.nops();
            
            if(Verbose>0) {
                cout << "+----------------------------------------" << endl;
                cout << "| IBPs: " << nr << ", Fs: " << next_gi << "/" << nc << endl;
                cout << "+----------------------------------------" << endl;
            }
            
            SparseMatrix smat;
            for(int r=0; r<nr; r++) {
                auto ibp = ibps.op(r);
                auto cvs = collect_lst(ibp,F(w));
                for(auto & cv : cvs) {
                    int c = f2n[cv.op(1)];
                    ex cc = cv.op(0);
                    if(!cc.is_zero()) smat[r][c] = cc;
                }
            }
            RowReduce(smat,next_gi);
            
            ibps.remove_all();
            for(auto const & rv : smat) {
                int next_sol_ibp = -1; // -1:init, 0:next, 1:sol, 2:ibp
                ex sol = 0, ns, fc;
                exvector nv(pdim);
                exmap amap;
                for(auto const & cv : rv.second) {
                    if(next_sol_ibp==-1) {
                        if(cv.first<next_gi) {
                            next_sol_ibp = 0;
                            if(!is_zero(cv.second-1)) throw Error("Solver::SolveSparse, NOT 1.");
                            ns = fvec[cv.first].op(0).subs(a(w)==0);
                            for(int i=0; i<pdim; i++) {
                                if(nfix[i]==0) amap[a(i)] = a(i)-ns.op(i);
                            }
                            fc = fvec[cv.first].subs(amap,nopat);
                            for(int i=0; i<pdim; i++) {
                                if(nfix[i]==1) {
                                    nv[i] = ns.op(i);
                                    if(sector.op(i)==1 && nv[i]<1) goto next_row; // subsector
                                    else if(sector.op(i)!=1 && nv[i]>0) goto next_row; // out of this sector
                                } else nv[i] = (sector.op(i)==1 ? 1 : 0);
                            }
                            continue; // need continue
                        } else next_sol_ibp = 2; // ibp
                    }
                    
                    if(next_sol_ibp==0) {
                        if(cv.first<next_gi) goto next_row; // not a solution or ibp
                        next_sol_ibp = 1; // solution found
                    }
                    
                    if(true) {
                        ex f = fvec[cv.first];
                        ex nd = cv.second;
                        if(next_sol_ibp==1) {
                            f = f.subs(amap,nopat);
                            nd = nd.subs(amap,nopat);
                            ex num = factor_flint(nd.numer());
                            ex den = factor_flint(nd.denom());
                            ns = f.op(0).subs(a(w)==0); // reuse ns here
                            for(int i=0; i<pdim; i++) {
                                if(sector.op(i)!=1) {
                                    ex nsi = ns.op(i);
                                    if(nfix[i]==1) {
                                        if(nsi>0) goto next_row; // out of this sector
                                        else continue;
                                    }
                                    nsi = -nsi; // note '-' here, ni<=nsi
                                    // due to possible (ai+nsi-1)*F(...,ai-nsi,...), so nsi+1 may be OK
                                    while(true) {
                                        ex cci = num.subs(a(i)==1+nsi);
                                        if(!cci.is_zero() || nsi>=nv[i]) break;
                                        nsi++;
                                    }
                                    if(nsi<nv[i]) nv[i] = nsi;
                                }
                            }
                        
                            if(!is_a<mul>(den)) den = lst{ den };
                            for(const auto & item : den) {
                                ex it = item;
                                if(is_a<power>(it)) it = it.op(0);
                                if(!it.has(a(w)) || it.has(d)) continue;
                                exset as;
                                find(it,a(w),as);
                                if(as.size()!=1) {
                                    cout << "size is NOT 1: " << it << endl;
                                    throw Error("Solver::SolveSparse, error occured.");
                                }
                                ex aw = *(as.begin());
                                if(it.degree(aw)!=1) {
                                    cout << "degree is NOT 1: " << it << endl;
                                    throw Error("Solver::SolveSparse, error occured.");
                                }
                                ex a0 = -it.coeff(aw,0)/it.coeff(aw,1);
                                int an = ex2int(aw.op(0));
                                if(sector.op(an)!=1) {
                                    if(a0<=nv[an]) nv[an] = a0-1;
                                } else {
                                    if(a0>=nv[an]) nv[an] = a0+1;
                                }
                            }
                        }
                        sol -= nd * f;
                    }
                }
                
                if(next_sol_ibp==1) {
                    for(int i=0; i<pdim; i++) { // check
                        if(nfix[i]==1 && fc.op(0).op(i)!=nv[i]) {
                            cout << i << ": " << fc << endl << nv << endl;
                            throw Error("Solver::SolveSparse, error occured.");
                        } else if(nfix[i]==0 && fc.op(0).op(i)!=a(i)) {
                            cout << i << ": " << fc << endl << nv << endl;
                            throw Error("Solver::SolveSparse, error occured.");
                        }
                    }

                    if(true) {
                        int psum = 0, nsum = 0;
                        for(int i=0; i<pdim; i++) {
                            if(nfix[i]==1) continue;
                            if(sector.op(i)==1 && !nv[i].is_equal(1)) psum++;
                            if(sector.op(i)!=1 && !nv[i].is_zero()) nsum++;
                        }
                        if(psum<=1 && nsum<=1) {
                            ex con = vec2lst(nv);
                            cons.append(con);
                            cons_vec[con].append(sol);
                        }
                    }
                } else ibps.append(sol);
                
                next_row: ;
            }
cons.sort();
cons.unique();
sort_lst(cons,false);
cout << cons << endl << endl;
for(auto & item : cons) {
    if(item==lst{1,1,1,1,1,1,1,1,1,1,0,0}) {
        cout << endl << cons_vec[item] << endl << endl;
        exit(0);
    }
}
        }
        
        cons.sort();
        cons.unique();
        sort_lst(cons,false);

//        cout << endl;
//        for(auto & k : cons) {
//            cons_vec[k].sort();
//            cons_vec[k].unique();
//            cout << k << endl;
//            cout << cons_vec[k] << endl << endl;
//        }
//        cout << endl << cons << endl;
        

    }
    
    void Solver::IBP() {
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
            throw Error("Solver::IBP: #(ISP) > #(Propagators).");
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
            if(eq.has(iWF(w))) throw Error("Solver::IBP, iWF used in eq.");
            leqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(leqns, ss);
        if(s2p.nops() != ISP.nops()) throw Error("Solver::IBP, lsolve failed.");
        
        // 1st version
        if(true)
        if(DSP.nops()<1) {
            for(auto p1 : Internal)
            for(auto p2 : InExternal)
            DSP.append(lst{p1,p2});
        }
        
        // Lee version
        if(false)
        if(DSP.nops()<1) {
            for(int i=0; i<Internal.nops(); i++)
                DSP.append(lst{Internal.op(i), Internal.op(i==Internal.nops()-1 ? 0 : i+1)});
            for(auto p2 : InExternal) DSP.append(lst{Internal.op(0), p2});
            DSP.append(lst{Internal, Internal});
        }

        IBPs.remove_all();
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
            IBPs.append(res);
        }
    }
    
    void Solver::Solve(const ex & sector, const map<int,int> & n2n) {
        int pdim = Propagators.nops();
        lst ibps = IBPs;
        
        // add more eqn to ibps
        if(false) {
//            int n = ibps.nops();
//            for(int i=0; i<pdim; i++) {
//                for(int j1=0; j1<n; j1++) for(int j2=0; j2<n; j2++)
//                    ibps.append(ibps.op(i).subs(lst{a(j1)==a(j1)+1}).subs(lst{a(j2)==a(j2)-1}));
//            }
            
            int n = ibps.nops();
            for(int i=0; i<pdim; i++) {
                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+1));
//                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+2));
//                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+3));
//                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)+4));
            }

            //n = ibps.nops();
            for(int i=0; i<pdim; i++) {
                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-1));
//                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-2));
//                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-3));
//                for(int j=0; j<n; j++) ibps.append(ibps.op(j).subs(a(i)==a(i)-4));
            }
        } else {
            exset fs;
            find(ibps, F(w1,w2), fs);
cout << fs << endl;
exit(0);
        }
        
        // some a's replaced by integer
        vector<int> nfix(pdim);
        for(int i=0; i<pdim; i++) nfix[i] = 0;
        for(auto const & kv : n2n) nfix[kv.first] = 1;
        exmap a2n;
        for(auto const & kv : n2n) {
            int k = kv.first;
            if(k<0) k = pdim + k;
            a2n[a(k)] = kv.second;
        }
        int nibps = ibps.nops();
        for(int i=0; i<nibps; i++) ibps.let_op(i) = ibps.op(i).subs(a2n);
        
        lst cons;
        map<ex,lst,ex_is_less> cons_vec;
        while(ibps.nops()>0) {
            ibps.sort();
            ibps.unique();
            
            exvector fvec;
            exset fset;
            find(ibps,F(w),fset);
            int nmax = -10000000, pmax = -10000000;
            for(auto fi : fset) { // sort fset
                int psum = 0, nsum = 0;
                auto ns = fi.op(0);
                for(int i=0; i<pdim; i++) {
                    auto ni = ns.op(i).subs(a(i)==0,nopat);
                    if(sector.op(i)==1) {
                        if(nfix[i]==1 && ni<=0) psum -= 100000;
                        else psum += ex2int(ni)-1; // -1 from FIRE
                    } else {
                        if(nfix[i]==1 && ni>0) nsum += ex2int(ni);
                        else nsum -= ex2int(ni);
                    }
                }
                
                if(psum+nsum<pmax+nmax) continue;
                if(psum+nsum>pmax+nmax) {
                    pmax = psum;
                    nmax = nsum;
                    fvec.clear();
                }
                
                fvec.push_back(fi);
            }
            //sort_vec(fvec);

            int nr = ibps.nops();
            int nc = fvec.size();
            matrix mat(nr, nc);
            
            if(Verbose>0) {
                cout << "+----------------------------------------" << endl;
                cout << "| IBPs: " << nr << ", Fs: " << nc << endl;
                cout << "+----------------------------------------" << endl;
            }
            
            matrix bvec(nr,1);
            for(int r=0; r<nr; r++) {
                auto ibp = ibps.op(r);
                for(int c=0; c<nc; c++) {
                    mat(r,c) = ibp;
                    ibp = ibp.subs(fvec[c]==0,nopat);
                    mat(r,c) = (mat(r,c)-ibp).coeff(fvec[c]);
                }
                bvec(r,0) = ibp;
            }

            // note that: mr = mt.mul(mat)
            auto rt = RowReduce(mat);
            auto mr = rt.first;
            auto mt = rt.second;
            bvec = mt.mul(bvec);
            for(int r=0; r<nr; r++) bvec(r,0) = collect_ex(bvec(r,0),F(w),o_flintfD);
            
            //cout << "CHECK: 1=" << ex_to<matrix>(normal(mr.sub(mt.mul(mat)))).is_zero_matrix() << endl;
                
            auto ibp_con_sol_vec = GiNaC_Parallel(nr, [&](const int r)->ex {
                ex ibp=0, con=0;
                int cc = -1;
                for(int c=0; c<nc; c++) {
                    if(mr(r,c).is_zero()) continue;
                    if(mr(r,c)==1 && cc==-1) cc = c;
                    else goto next_row; // not a solution, skip
                }
                
                if(cc!=-1) { // a solution found N -> N-1, N-2, subsectors
                    auto ns = fvec[cc].op(0).subs(a(w)==0);
                    exmap amap;
                    for(int i=0; i<pdim; i++) {
                        if(nfix[i]==0) amap[a(i)] = a(i)-ns.op(i);
                    }
                    fvec[cc] = fvec[cc].subs(amap,nopat);
                    auto sol = bvec(r,0).subs(amap,nopat);
                    sol = normal_flint(sol,o_flintfD);

                    exvector nv(pdim);
                    for(int i=0; i<pdim; i++) {
                        if(nfix[i]==1) {
                            nv[i] = ns.op(i);
                            if(sector.op(i)==1 && nv[i]<1) goto next_row; // subsector
                            else if(sector.op(i)!=1 && nv[i]>0) goto next_row; // out of this sector
                        } else nv[i] = (sector.op(i)==1 ? 1 : 0);
                    }

                    ex num = sol.numer();
                    exset fset;
                    find(sol,F(w),fset);
                    for(auto const & f : fset) { // make sure not to increase the sector
                        ex cc = factor_flint(num.coeff(f));
                        ns = f.op(0).subs(a(w)==0); // reuse ns here
                        for(int i=0; i<pdim; i++) {
                            if(sector.op(i)!=1) {
                                ex nsi = ns.op(i);
                                if(nfix[i]==1) {
                                    if(nsi>0) goto next_row; // out of this sector
                                    else continue;
                                }
                                nsi = -nsi; // note '-' here, ni<=nsi
                                // due to possible (ai+nsi-1)*F(...,ai-nsi,...), so nsi+1 may be OK
                                while(true) {
                                    ex cci = cc.subs(a(i)==1+nsi);
                                    if(!cci.is_zero() || nsi>=nv[i]) break;
                                    nsi++;
                                }
                                if(nsi<nv[i]) nv[i] = nsi;
                            }
                        }
                    }
                    
                    ex den = sol.denom();
                    if(!is_a<mul>(den)) den = lst{den};
                    for(const auto & item : den) {
                        ex it = item;
                        if(is_a<power>(it)) it = it.op(0);
                        if(!it.has(a(w)) || it.has(d)) continue;
                        exset as;
                        find(it,a(w),as);
                        if(as.size()!=1) {
                            cout << "size is NOT 1: " << it << endl;
                            throw Error("Solver::Solve, error occured.");
                        }
                        ex aw = *(as.begin());
                        if(it.degree(aw)!=1) {
                            cout << "degree is NOT 1: " << it << endl;
                            throw Error("Solver::Solve, error occured.");
                        }
                        ex a0 = -it.coeff(aw,0)/it.coeff(aw,1);
                        int an = ex2int(aw.op(0));
                        if(sector.op(an)!=1) {
                            if(a0<=nv[an]) nv[an] = a0-1;
                        } else {
                            if(a0>=nv[an]) nv[an] = a0+1;
                        }
                    }
                    
                    for(int i=0; i<pdim; i++) { // check
                        if(nfix[i]==1 && fvec[cc].op(0).op(i)!=nv[i]) {
                            cout << fvec << endl << nv << endl;
                            throw Error("Solver::Solve, error occured.");
                        } else if(nfix[i]==0 && fvec[cc].op(0).op(i)!=a(i)) {
                            cout << fvec << endl << nv << endl;
                            throw Error("Solver::Solve, error occured.");
                        }
                    }
                    
                    //cout << "corner: " << nv << endl;
                    int psum = 0, nsum = 0;
                    for(int i=0; i<pdim; i++) {
                        if(nfix[i]==1) continue;
                        if(sector.op(i)==1 && !nv[i].is_equal(1)) psum++;
                        if(sector.op(i)!=1 && !nv[i].is_zero()) nsum++;
                    }
                    if(nsum<=1 && psum<=1) con = lst{ vec2lst(nv), sol };
                } else if(!bvec(r,0).is_zero()) ibp = bvec(r,0).numer();
                next_row: ;
                return lst{ibp, con};
            });
            
            ibps.remove_all();
            for(int r=0; r<nr; r++) {
                auto const & item = ibp_con_sol_vec[r];
                if(!item.op(0).is_zero()) ibps.append(item.op(0));
                if(is_a<lst>(item.op(1))) {
                    cons.append(item.op(1).op(0));
                    cons_vec[item.op(1).op(0)].append(item.op(1).op(1));
                }
            }
            
cons.sort();
cons.unique();
cout << cons << endl << endl;
        }
        
        cons.sort();
        cons.unique();
        sort_lst(cons,false);

        cout << endl;
        for(auto & k : cons) {
            cons_vec[k].sort();
            cons_vec[k].unique();
            cout << k << endl;
            cout << cons_vec[k] << endl << endl;
        }
        cout << endl << cons << endl;
    }

    /**
     * @brief Export input data for KIRA
     */
    void Solver::Export() {
        
    }

    
    /**
     * @brief invoke kira program for reduction
     */
    void Solver::Run() {
    
    }

    /**
     * @brief import kira result
     */
    void Solver::Import() {
        
    }

}
