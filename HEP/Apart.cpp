/**
 * @file
 * @brief Functions to perform partial fraction
 */

#include "HEP.h"

namespace HepLib {

    #ifndef DOXYGEN_SKIP
    
    unsigned ApartIR1_SERIAL::serial = GiNaC::function::register_new(function_options("ApartIR",1).do_not_evalf_params().overloaded(2));
    unsigned ApartIR2_SERIAL::serial = GiNaC::function::register_new(function_options("ApartIR",2).do_not_evalf_params().overloaded(2));
    
    #endif

    namespace {
        inline bool isOK(const ex &expr_in, const lst &vars) {
            ex item = expr_in;
            if(is_a<power>(item)) {
                if(!is_a<numeric>(item.op(1)) || !ex_to<numeric>(item.op(1)).is_integer()) return false;
                item = item.op(0);
            }
            exmap vs; symbol s;
            for(auto v : vars) vs[v] = s*v;
            return item.subs(vs,nopat).degree(s)<=1;
        }
    }
        
    /**
     * @brief convert ApartIR to ex
     * @param expr_in expression contains ApartIR
     * @return ApartIR converted into normal ex
     */
    ex ApartIR2ex(const ex & expr_in) {
        ex ret = expr_in;
        ret = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(e.match(ApartIR(1)) || e.match(ApartIR(1,w))) return 1;
            else if(!e.has(ApartIR(w1,w2))) return e;
            else if(e.match(ApartIR(w1, w2))) {
                ex vars = e.op(1);
                ex ret=1;
                matrix mat = ex_to<matrix>(e.op(0));
                for(int c=0; c<mat.cols(); c++) {
                    ex sum=0;
                    for(int r=0; r<mat.rows()-2; r++) sum += mat(r,c) * vars.op(r);
                    sum += mat(mat.rows()-2,c);
                    ret *= pow(sum, mat(mat.rows()-1,c));
                }
                return ret;
            } else return e.map(self);
        })(ret);
        return ret;
    }
    
    /**
     * @brief convert ApartIR to F(ps, ns), ns is like FIRE convention
     * @param expr_in expression contains ApartIR
     * @return ApartIR converted into F(ps, ns)
     */
    ex ApartIR2F(const ex & expr_in) {
        ex ret = expr_in;
        ret = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(e.match(ApartIR(1)) || e.match(ApartIR(1,w))) return 1;
            else if(!e.has(ApartIR(w1,w2))) return e;
            else if(e.match(ApartIR(w1, w2))) {
                ex vars = e.op(1);
                matrix mat = ex_to<matrix>(e.op(0));
                lst pns;
                for(int c=0; c<mat.cols(); c++) {
                    ex sum=0;
                    for(int r=0; r<mat.rows()-2; r++) sum += mat(r,c) * vars.op(r);
                    sum += mat(mat.rows()-2,c);
                    if(is_zero(mat(mat.rows()-1,c))) continue;
                    pns.append(lst{ sum, ex(0)-mat(mat.rows()-1,c) });
                }
                sort_lst(pns); // need sort_lst
                lst ps, ns;
                for(auto item : pns) {
                    ps.append(item.op(0));
                    ns.append(item.op(1));
                }
                return F(ps, ns);
            } else return e.map(self);
        })(ret);
        return ret;
    }
    
    inline ex air2pn(ex air) {
        matrix mat = ex_to<matrix>(air.op(0));
        auto vars = air.op(1);
        lst pns;
        for(int c=0; c<mat.cols(); c++) {
            ex sum=0;
            for(int r=0; r<mat.rows()-2; r++) sum += mat(r,c) * vars.op(r);
            sum += mat(mat.rows()-2,c);
            if(is_zero(mat(mat.rows()-1,c))) continue;
            pns.append(lst{ sum, mat(mat.rows()-1,c) });
        }
        sort_lst(pns); // need sort_lst
        lst ps, ns;
        for(auto item : pns) {
            ps.append(item.op(0));
            ns.append(item.op(1));
        }
        return lst{ps, ns};
    }
    
    inline ex pn2mat(ex ps, ex ns, ex vars) { // no need to use normal and removed
        lst pnlst;
        for(int i=0; i<ps.nops(); i++) pnlst.append(lst{ ps.op(i), ns.op(i) });
        sort_lst(pnlst); // need sort_lst
        int nrow=vars.nops(), ncol=pnlst.nops();
        exmap vars0;
        for(auto v : vars) vars0[v]=0;
        matrix mat(nrow+2, ncol);
        for(int c=0; c<ncol; c++) {
            ex pn = pnlst.op(c);
            mat(nrow+1,c) = normal(pn.op(1));
            auto tmp = pn.op(0);
            for(int r=0; r<nrow; r++) {
                mat(r,c) = (tmp.coeff(vars.op(r))); // remove exnormal
            }
            mat(nrow,c) = (tmp.subs(vars0,nopat)); // remove exnormal
        }
        return mat;
    }
    
    // e1, e2 are sorted list, check if e1 is a subset of e2 or not.
    inline bool is_subset(const ex & e1, const ex & e2) {
        int i1=0, i2=0;
        int t1 = e1.nops();
        int t2 = e2.nops();
        if(t1>t2) return false;
        while(i1<t1) {
            if(t1-i1>t2-i2) return false;
            if(e1.op(i1).is_equal(e2.op(i2))) {
                i1++; i2++; continue;
            } 
            i2++;
        }
        return true;
    }

    inline ex pn2p(ex pn) {
        lst l;
        auto p = pn.op(0);
        auto n = pn.op(1);
        for(int i=0; i<n.nops(); i++) {
            if(n.op(i).is_zero()) continue;
            l.append(p.op(i));
        }
        return l;
    }
    
    exmap ApartRules(const exvector &airs, bool irc) { // irc=true to include ApartIRC
        exmap rules;
        if(airs.size()<1) return rules;
        #define CHECK false
        int nlimit = 50;
        vector<exset> nps_set(airs[0].op(1).nops());
        for(auto air : airs) {
            auto pn = air2pn(air);
            auto l = pn2p(pn);
            nps_set[l.nops()-1].insert(l);
        }
        vector<exvector> nps(nps_set.size());
        for(int i=0; i<nps_set.size(); i++) nps[i] = exvector(nps_set[i].begin(), nps_set[i].end());
        nps_set.clear();
        
        exmap s2s;
        for(int i=0; i<nps.size()-1; i++) {
            auto psi = nps[i];
            if(psi.size()<1) continue;
            if(psi.size()<nlimit) {
                for(auto si : psi) {
                    for(int j=nps.size()-1; j>i; j--) {
                        auto psj = nps[j];
                        for(auto sj : psj) {
                            if(is_subset(si,sj)) { s2s[si] = sj; goto nextsi; }
                        }
                    }
                    nextsi: ;
                }
            } else {
                exvector psi_vec(psi.begin(), psi.end());
                auto ret = GiNaC_Parallel(psi_vec.size(), [&psi_vec,&nps,i](int idx)->ex{
                    auto si = psi_vec[idx];
                    for(int j=nps.size()-1; j>i; j--) {
                        auto psj = nps[j];
                        for(int jj=0; jj<psj.size(); jj++) {
                            //if(is_subset(si,sj)) return lst{si, sj};
                            auto sj = psj[jj];
                            if(is_subset(si,sj)) return lst{idx, j, jj};
                        }
                    }
                    return lst{ };
                }, "AR-"+to_string(nps.size()-1)+"-"+to_string(i+1));
                for(auto item : ret) if(item.nops()>1) {
                    int idx = ex_to<numeric>(item.op(0)).to_int();
                    int j = ex_to<numeric>(item.op(1)).to_int();
                    int jj = ex_to<numeric>(item.op(2)).to_int();
                    s2s[psi_vec[idx]] = nps[j][jj];
                }
            }
        }

        if(airs.size()<nlimit) {
            for(int k=0; k<airs.size(); k++) {
                auto pn = air2pn(airs[k]);
                auto p = pn.op(0);
                auto n = pn.op(1);
                exmap p2n;
                for(int i=0; i<p.nops(); i++) if(!is_zero(n.op(i))) p2n[p.op(i)] = n.op(i);
                auto pk = pn2p(pn);
                auto fi = s2s.find(pk);
                if(fi!=s2s.end()) {
                    auto pp = fi->second;
                    lst nn;
                    for(int i=0; i<pp.nops(); i++) nn.append(p2n[pp.op(i)]);
                    ex air = ApartIR(pn2mat(pp,nn,airs[k].op(1)),airs[k].op(1));
                    if(CHECK) { // check
                        ex eL=1, eR=1;
                        for(int i=0; i<p.nops(); i++) if(!n.op(i).is_zero()) eL *= pow(WF(p.op(i)), n.op(i));
                        for(int i=0; i<pp.nops(); i++) if(!nn.op(i).is_zero()) eR *= pow(WF(pp.op(i)), nn.op(i));
                        ex chk = normal(eL-eR);
                        if(!chk.is_zero()) {
                            cout << p << ", " << n << endl;
                            cout << pp << ", " << nn << endl;
                            throw Error("ApartRules: check failed!");
                        }
                    }
                    if(irc) air = ApartIRC(air);
                    rules[airs[k]] = air;
                } else if(irc) rules[airs[k]] = ApartIRC(airs[k]);
            }
        } else {
            auto ret = GiNaC_Parallel(airs.size(), [&airs,&s2s,irc](int k)->ex{
                auto pn = air2pn(airs[k]);
                auto p = pn.op(0);
                auto n = pn.op(1);
                exmap p2n;
                for(int i=0; i<p.nops(); i++) if(!is_zero(n.op(i))) p2n[p.op(i)] = n.op(i);
                auto pk = pn2p(pn);
                auto fi = s2s.find(pk);
                if(fi!=s2s.end()) {
                    auto pp = fi->second;
                    lst nn;
                    for(int i=0; i<pp.nops(); i++) nn.append(p2n[pp.op(i)]);
                    ex air = ApartIR(pn2mat(pp,nn,airs[k].op(1)),airs[k].op(1));
                    air = ApartIRC(air);
                    if(CHECK) { // check
                        ex eL=1, eR=1;
                        for(int i=0; i<p.nops(); i++) if(!n.op(i).is_zero()) eL *= pow(WF(p.op(i)), n.op(i));
                        for(int i=0; i<pp.nops(); i++) if(!nn.op(i).is_zero()) eR *= pow(WF(pp.op(i)), nn.op(i));
                        ex chk = normal(eL-eR);
                        if(!chk.is_zero()) {
                            cout << p << ", " << n << endl;
                            cout << pp << ", " << nn << endl;
                            throw Error("ApartRules: check failed!");
                        }
                    }
                    if(irc) air = ApartIRC(air);
                    return lst{ airs[k], air };
                } else if(irc) return lst{ airs[k], ApartIRC(airs[k]) };
                else return lst{ };
            }, "AR");
            //ReShare(ret,airs);
            for(auto lr : ret) if(lr.nops()>1) rules[lr.op(0)] = lr.op(1);
        }
        return rules;

    } 
    
    /**
     * @brief Apart on matrix
     * @param mat each column: [c1,...,cn,c0,n] -> (c1 x1+...+cn xn+c0)^n
     * @return sum of coefficient * ApartIR
     */
    ex Apart(const matrix & mat) {
        static exmap mat_cache;
        if(using_cache && cache_limit>0 && mat_cache.size() > cache_limit) mat_cache.clear();
        auto mat_itr = mat_cache.find(mat);
        if(mat_itr!=mat_cache.end()) return mat_itr->second;
        
        int nrow = mat.rows()-2;
        int ncol = mat.cols();
        lst null_vec;
        
        // null vector
        static exmap null_cache;
        if(cache_limit>0 && null_cache.size() > cache_limit) null_cache.clear();
        ex key = sub_matrix(mat,0,nrow,0,ncol);
        auto null_itr = null_cache.find(key);
        if(null_itr==null_cache.end()) {
            if(Apart_using_fermat) {
                Fermat &fermat = Fermat::get();
                int &v_max = fermat.vmax;
        
                lst rep_vs;
                ex tree = mat;
                for(const_preorder_iterator i = tree.preorder_begin(); i != tree.preorder_end(); ++i) {
                    auto e = (*i);
                    if(is_a<symbol>(e) || is_a<Pair>(e) || is_a<Eps>(e)) {
                        rep_vs.append(e);
                    }
                }
                rep_vs.sort();
                rep_vs.unique();
                
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
                
                ss << "Array m[" << nrow << "," << ncol+1 << "];" << endl;
                fermat.Execute(ss.str());
                ss.clear();
                ss.str("");
                
                ss << "[m]:=[(";
                for(int c=0; c<ncol; c++) {
                    for(int r=0; r<nrow; r++) {
                        ss << mat(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
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
                if(ostr[ostr.length()-1]!='0') throw Error("Apart: last char is NOT 0.");
                ostr = ostr.substr(0, ostr.length()-1);
                string_trim(ostr);
                
                ostr.erase(0, ostr.find(":=")+2);
                string_replace_all(ostr, "[", "{");
                string_replace_all(ostr, "]", "}");
                Parser fp(st);
                auto mat2 = fp.Read(ostr);

                exmap xs;
                for(int c=0; c<ncol; c++) xs[c] = iWF(c);
                for(int r=nrow-1; r>=0; r--) {
                    ex xadd = 0;
                    int pi=-1;
                    for(int c=0; c<ncol; c++) {
                        if(!is_zero(get_op(mat2,r,c)) && pi<0) pi = c;
                        else xadd -= get_op(mat2,r,c)*xs[c];
                    }
                    xs[pi] = xadd;
                }
                for(int c=0; c<ncol; c++) null_vec.append(xs[c]);
            } else {
                matrix v(ncol, 1);
                exmap sRepl;
                for (int c=0; c<ncol; c++) {
                    symbol t;
                    v(c,0) = t;
                    sRepl[t]=iWF(c);
                }
                // Solve M*V = 0
                matrix zero(nrow, 1);
                matrix s = ex_to<matrix>(key.subs(iEpsilon==0,nopat)).solve(v,zero);
                for(int r=0; r<ncol; r++) null_vec.append(s(r,0).subs(sRepl,nopat));
            }
            null_cache.insert({key,null_vec});
        } else {
            null_vec = ex_to<lst>(null_itr->second);
        }

        // check null & return ApartIR
        bool is_null = true;
        for(int c=0; c<ncol; c++) {
            if(!is_zero(null_vec.op(c))) {
                is_null = false;
                break;
            }
        }
        if(is_null) {
            ex res = ApartIR(mat);
            if(using_cache) mat_cache.insert({mat,res});
            return res;
        }
        
        // handle numerator
        int ni=-1;
        for(int c=0; c<ncol; c++) {
            if(mat(nrow+1,c)>0 && !is_zero(null_vec.op(c))) {
                ni = c;
                break;
            }
        }
        
        if(ni!=-1) {
            ex nvec;
            if(true) {
                exset wfs;
                find(null_vec.op(ni),iWF(w),wfs);
                if(wfs.size()<1) throw Error("Apart: something is wrong!");
                symbol s;
                
                int max = -1;
                for(auto wf : wfs) {
                    auto n1 = subs(null_vec, wf==s);
                    n1 = subs(n1, iWF(w)==0);
                    for(int i=0; i<3; i++) {
                        auto n2 = subs(n1, s==i);
                        if(!is_zero(n2.op(ni))) {
                            int nt = 0;
                            for(auto ii : n2) if(ii.is_zero()) nt++;
                            if(nt>max && nt!=ncol) { nvec = n2; max = nt; }
                            break;
                        }
                    }
                }
                if(nvec.has(iWF(w)) || is_zero(nvec.op(ni))) throw Error("Apart: iWF to int failed with "+ex2str(null_vec.op(ni)));
            }
            
            ex sol = 0;
            for(int c=0; c<ncol; c++) {
                if(c==ni) continue;
                sol -= nvec.op(c) * (iWF(c)-mat(nrow,c)); // iWF(c) refer to c-th column, and minus the const term
            }
            sol = sol/nvec.op(ni) + mat(nrow,ni); // last one: const term
            sol = collect_ex(pow(sol.subs(iEpsilon==0,nopat), mat(nrow+1,ni)), iWF(w)); // expand the numerator with power

            if(!is_a<add>(sol)) sol = lst{ sol };
            ex res = 0;
            for(auto item : sol) {
                int nzero = 0;
                matrix mat2(nrow+2, ncol-1);
                for(int c=0; c<ncol; c++) {
                    if(c==ni) continue;
                    int c2 = (c<ni ? c : c-1);
                    for(int r=0; r<nrow+1; r++) {
                        mat2(r,c2) = mat(r,c);
                    }
                    auto expn = mat(nrow+1,c) + item.degree(iWF(c));
                    if(is_zero(expn)) nzero++;
                    mat2(nrow+1,c2) = expn;
                }
                if(nzero>0) {
                    if((ncol-1-nzero)==0) {
                        res += ApartIR(1) * item.subs(iWF(w)==1);
                    } else {
                        matrix mat3(nrow+2, ncol-1-nzero);
                        int cc=0;
                        for(int c=0; c<ncol-1; c++) {
                            if(is_zero(mat2(nrow+1,c))) continue;
                            for(int r=0; r<nrow+2; r++) {
                                mat3(r,cc) = mat2(r,c);
                            }
                            cc++;
                        }
                        res += Apart(mat3) * item.subs(iWF(w)==1);
                    }
                } else {
                    res += Apart(mat2) * item.subs(iWF(w)==1);
                }
            }
            res = collect_ex(res,ApartIR(w));
            if(using_cache) mat_cache.insert({mat,res});
            return res;
        }
        
        // handle all denominators
        ex cres0 = 0;
        for(int c=0; c<ncol; c++) cres0 += mat(nrow,c)*null_vec.op(c);
        cres0 = cres0.subs(iEpsilon==0,nopat);
        
        ex nvec,cres;
        if(true) {
            exset wfs;
            find(lst{cres0, null_vec},iWF(w),wfs);
            if(wfs.size()<1) {
                cres = cres0;
                nvec = null_vec;
            } else {
                symbol s;
                int max = -1;
                bool c0 = true;
                for(auto wf : wfs) {
                    auto c1 = subs(cres0, wf==s);
                    auto n1 = subs(null_vec, wf==s);
                    c1 = subs(c1, iWF(w)==0);
                    n1 = subs(n1, iWF(w)==0);
                    for(int i=0; i<5; i++) {
                        auto c2 = subs(c1, s==i).normal();
                        auto n2 = subs(n1, s==i);
                        if(!c0 && c2.is_zero()) continue;
                        if(c0 && !c2.is_zero()) {
                            c0 = false;
                            nvec = n2; cres = c2;
                        } else {
                            int nt = 0;
                            for(auto ii : n2) if(ii.is_zero()) nt++;
                            if(nt>max && nt!=ncol) { nvec = n2; cres = c2; max = nt; }
                        }
                    }
                }
            }
        } // end if(true)

        if(nvec.has(iWF(w)) || cres.has(iWF(w))) {
            cout << cres0 << ", " << null_vec << endl;
            throw Error("Apart: iWF still left.");
        }
        
        // handle const is NOT zero
        if(!IsZero(cres)) {
            ex res=0;
            for(int c=0; c<ncol; c++) {
                if(is_zero(nvec.op(c))) continue;
                if(is_zero(mat(nrow+1,c)+1)) {
                    matrix mat2(nrow+2,ncol-1);
                    int ccc = 0;
                    for(int cc=0; cc<ncol; cc++) {
                        if(cc==c) continue;
                        for(int r=0; r<nrow+2; r++) mat2(r,ccc)=mat(r,cc);
                        ccc++;
                    }
                    res += Apart(mat2) * nvec.op(c);
                } else {
                    matrix mat2(mat);
                    mat2(nrow+1,c) = mat2(nrow+1,c)+1;
                    res += Apart(mat2) * nvec.op(c);
                }
            }
            res = res/cres;
            res = collect_ex(res,ApartIR(w));
            if(using_cache) mat_cache.insert({mat,res});
            return res;
        } else {
            int ni=-1;
            for(int c=0; c<ncol; c++) {
                if(!is_zero(nvec.op(c))) {
                    ni = c;
                    break;
                }
            }

            ex res=0;
            for(int c=0; c<ncol; c++) {
                if(is_zero(nvec.op(c)) || c==ni) continue;
                if(is_zero(mat(nrow+1,c)+1)) {
                    matrix mat2(nrow+2,ncol-1);
                    int ccc = 0;
                    for(int cc=0; cc<ncol; cc++) {
                        if(cc==c) continue;
                        for(int r=0; r<nrow+2; r++) mat2(r,ccc)=mat(r,cc);
                        ccc++;
                    }
                    int ni2 = ni>c ? ni-1 : ni;
                    mat2(nrow+1,ni2) = mat2(nrow+1,ni2)-1;
                    res -= Apart(mat2) * nvec.op(c);
                } else {
                    matrix mat2(mat);
                    mat2(nrow+1,c) = mat2(nrow+1,c)+1;
                    mat2(nrow+1,ni) = mat2(nrow+1,ni)-1;
                    res -= Apart(mat2) * nvec.op(c);
                }
            }
            res = res/nvec.op(ni);
            res = collect_ex(res,ApartIR(w));
            if(using_cache) mat_cache.insert({mat,res});
            return res;
        }
    }

    /**
     * @brief Apart on ex
     * @param expr_ino normal expresion, product of [ linear w.r.t. vars ]^n
     * @param vars_in independent variables
     * @param smap the sign map
     * @return sum of coefficient * ApartIR
     */
    ex Apart(const ex &expr_ino, const lst &vars_in, exmap smap) {
        // Apart on rational terms
        if(!is_a<lst>(expr_ino)) {
            auto cv_lst = collect_lst(expr_ino, vars_in);
            ex res = 0;
            for(auto item : cv_lst) res += item.op(0) * Apart(lst{item.op(1)}, vars_in, smap);
            res = collect_ex(res, ApartIR(w1,w2));

            // random check
            lst nlst;
            for(const_preorder_iterator i = res.preorder_begin(); i != res.preorder_end(); ++i) {
                auto e = (*i);
                if(is_a<symbol>(e) || is_a<Pair>(e)) nlst.append(e);
            }
            for(auto var : vars_in) nlst.append(var);
            nlst.sort();
            nlst.unique();
            exmap nrepl;
            auto pi = nextprime(3);
            for(auto ni : nlst) {
                pi = nextprime(pi+1);
                nrepl[ni] = ex(1)/pi;
            }
            nrepl[iEpsilon]=0;
            ex chk = ApartIR2ex(subs(res,nrepl))-subs(expr_ino,nrepl);
            chk = exnormal(chk);
            if(!is_zero(chk)) throw Error("Apart@1 random check Failed.");
            return res;
        }
        
        // Apart on monomial term
        if(expr_ino.nops()!=1) throw Error("Apart: wrong convention found!");
        ex expr_in = expr_ino.op(0);

        static exmap apart_cache;
        if(using_cache && cache_limit>0 && apart_cache.size() > cache_limit) apart_cache.clear();
        ex ckey = lst{expr_in, vars_in};
        auto itr = apart_cache.find(ckey);
        if(itr!=apart_cache.end()) return itr->second;
    
        exmap map1, map2;
        lst vars;
        for(int i=0; i<vars_in.nops(); i++) {
            auto v = vars_in.op(i);
            Symbol s("_apX"+to_string(i));
            map1[v]=s;
            map2[s]=v;
            vars.append(s);
        }
        exmap sgnmap;
        for(auto kv : smap) sgnmap[kv.first.subs(map1,nopat)] = kv.second.subs(map1,nopat);
        
        ex expr = expr_in.subs(map1,nopat);
        if(!is_a<mul>(expr)) expr = lst{expr};
        
        // check only, try normal_fermat if faild
        bool ok = true;
        for(auto item : expr) {
            bool has_var=false;
            for(auto v : vars) {
                if(item.has(v)) {
                    has_var=true;
                    break;
                }
            }
            if(has_var) {
                if(!isOK(item,vars)) {
                    ok = false;
                    break;
                }
            }
        }
        if(!ok) {
            expr = expr_in.subs(map1,nopat);
            expr = exnormal(expr,o_flintfD); // need option factor=true, factor denominator
            if(!is_a<mul>(expr)) expr = lst{expr};
        }
        
        lst pnlst;
        ex pref = 1;
        map<ex,int,ex_is_less> count_ip;
        exmap ie_map;
        for(auto item : expr) {
            bool has_var=false;
            for(auto v : vars) {
                if(item.has(v)) {
                    has_var=true;
                    break;
                }
            }
            if(has_var) {
                if(!isOK(item,vars)) {
                    cout << expr_in << endl;
                    cout << item << endl;
                    throw Error("Apart: item is not linear wrt vars.");
                }
                
                ex pc, nc;
                if(is_a<power>(item)) {
                    pc = item.op(0);
                    nc = item.op(1);
                } else {
                    pc = item;
                    nc = 1;
                }
                
                // consider sign
                bool has_sgn = false;
                // iEpsilon
                if(pc.has(iEpsilon)) { // iEpsilon first
                    ex si = 1; // default
                    auto itr = sgnmap.find(iEpsilon);
                    if(itr != sgnmap.end() && !is_zero(itr->second)) si = itr->second;
                    ex sign = si/pc.coeff(iEpsilon);
                    pref /= pow(sign, nc);
                    pc *= sign;
                    has_sgn = true;
                }
                // v from sgnmap
                if(!has_sgn) {
                    for(auto v : vars) {
                        auto cc = pc.coeff(v);
                        if(is_zero(cc) || !key_exists(sgnmap,v) || is_zero(sgnmap[v])) continue;
                        ex sign = sgnmap[v]/cc;
                        pref /= pow(sign, nc);
                        pc *= sign;
                        has_sgn = true;
                        break;
                    }
                }
                // v from vars
                if(!has_sgn) {
                    for(auto v : vars) {
                        if(key_exists(sgnmap,v)) continue; // already handled by sgnmap
                        auto cc = pc.coeff(v);
                        if(is_zero(cc) || !is_a<numeric>(cc)) continue;
                        ex sign = 1/cc;
                        pref /= pow(sign, nc);
                        pc *= sign;
                        break;
                    }
                }
                
                ex key = expand(pc.subs(iEpsilon==0,nopat));
                if(pc.has(iEpsilon)) {
                    if(ie_map.find(-key) != ie_map.end()) {
                        cout << expr_ino << endl;
                        cout << "Item 1: " << ie_map[-key].subs(map2,nopat) << endl;
                        cout << "Item 2: " << pc.subs(map2,nopat) << endl;
                        throw Error("iEpsilon Error: maybe pinch singularity?");
                    }
                    ie_map[key] = pc;
                }
                
                auto itr = count_ip.find(key);
                if(itr==count_ip.end()) itr = count_ip.find(-key);
                if(itr==count_ip.end()) count_ip[key] = 1;
                else itr->second++;
        
                pnlst.append(lst{ pc, nc });
            } else pref *= item;
        }
        // check whether needs to combine again
        bool needs_again = false;
        for(auto kv : count_ip) {
            if(kv.second>1) { needs_again = true; break; }
        }
        if(needs_again) {
            exmap imap;
            for(auto pn : pnlst) {
                auto key = pn.op(0);
                if(key.has(iEpsilon)) {
                    auto k = key.subs(iEpsilon==0,nopat);
                    if(imap.find(-k) != imap.end()) {
                        cout << "Item 1: " << imap[k] << endl;
                        cout << "Item 2: " << key << endl;
                        throw Error("iEpsilon Error: maybe pinch singularity?");
                    }
                    if(imap.find(k)==imap.end()) imap[k] = key;
                }
            }
            exmap p2n;
            for(auto pn : pnlst) {
                auto key = pn.op(0);
                auto itr = imap.find(key);
                if(itr!=imap.end()) key = itr->second;
                else {
                    itr = imap.find(-key);
                    if(itr!=imap.end()) {
                        key = itr->second;
                        pref *= pow(-1, pn.op(1));
                    }
                }
                p2n[key] = p2n[key] + pn.op(1);
            }
            pnlst.remove_all();
            for(auto kv : p2n) {
                if(is_zero(kv.second)) continue;
                auto k = kv.first;
                auto v = kv.second;
                if(is_a<numeric>(v) && ex_to<numeric>(v)>0) {
                    k = k.subs(iEpsilon==0,nopat);
                }
                pnlst.append(lst{k, v});
            }
        }
        sort_lst(pnlst); // need sort_lst
        
        if(pnlst.nops()==0) return apart_cache[ckey] = pref * ApartIR(1,vars_in);
        
        int nrow=vars.nops(), ncol=pnlst.nops();
        exmap vars0;
        for(auto v : vars) vars0[v]=0;
        matrix mat(nrow+2, ncol);
        for(int c=0; c<ncol; c++) {
            ex pn = pnlst.op(c);
            if(pn.op(0).is_zero() && pn.op(1).is_zero()) throw Error("Apart: 0^0 Found!");
            mat(nrow+1,c) = normal(pn.op(1));
            auto tmp = pn.op(0);
            for(int r=0; r<nrow; r++) {
                mat(r,c) = exnormal(tmp.coeff(vars.op(r)));
            }
            mat(nrow,c) = exnormal(tmp.subs(vars0,nopat));
        }

        ex ret = Apart(mat);

        auto cv_lst = collect_lst(ret,ApartIR(w));
        ret = 0;
        for(auto cv : cv_lst) {
            ret += pref * cv.op(0) * ApartIR(cv.op(1).op(0), vars).subs(map2,nopat);
        }

        return apart_cache[ckey] = ret;
    }
    
    /**
     * @brief Apart on ex
     * @param expr_ino input expression
     * @param loops list of loop Vector
     * @param extps list of external Vector
     * @param smap the sign map
     * @return sum of coefficient * ApartIR
     */
    ex Apart(const ex &expr_ino, const lst &loops, const lst & extps, exmap smap) {
        auto expr_in = expr_ino.subs(SP_map,nopat);
        auto expr = expr_in;
        
        lst sps;
        for(auto li : loops) {
            for(auto li2: loops) {
                auto item = SP(li, li2).subs(SP_map,nopat);
                if(is_a<Pair>(item)) sps.append(item);
            }
            for(auto ei: extps) {
                auto item = SP(li, ei).subs(SP_map,nopat);
                if(is_a<Pair>(item)) sps.append(item);
            }
        }
        sps.sort();
        sps.unique();
        
        auto cv_lst = collect_lst(expr, loops);
        ex res = 0;
        for(auto item : cv_lst) {
            res += item.op(0) * Apart(lst{item.op(1)}, sps, smap);
        }
        
        res = collect_ex(res, ApartIR(w1,w2));

        // random check
        lst nlst;
        for(const_preorder_iterator i = res.preorder_begin(); i != res.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e) || is_a<Pair>(e)) nlst.append(e);
        }
        nlst.sort();
        nlst.unique();
        exmap nrepl;
        auto pi = nextprime(3);
        for(auto ni : nlst) {
            pi = nextprime(pi+1);
            nrepl[ni] = 1/pi;
        }
        nrepl[iEpsilon]=0;
        ex chk = ApartIR2ex(subs(res,nrepl))-subs(expr_in,nrepl);
        chk = normal(chk);
        if(!is_zero(chk)) {
            throw Error("Apart@2 random check Failed.");
        }

        return res;
    }
    
    /**
     * @brief complete the ApartIR elements
     * @param expr_in input expression
     * @return ApartIR with complete matrix rank, ready for IBP reduction
     */
    ex ApartIRC(const ex & expr_in) {
        exmap cache;
        MapFunction map([&cache](const ex & e, MapFunction &self)->ex {
            if(!e.has(ApartIR(w1,w2))) return e;
            else if(e.match(ApartIR(w1,w2))) {
                auto i = cache.find(e);
                if(i!=cache.end()) return i->second;
                int n = e.op(1).nops();
                if((e.op(0)).is_equal(1)) {
                    matrix mat(n+2,n);
                    for(int r=0; r<n+2; r++) {
                        for(int c=0; c<n; c++) mat(r,c) = 0;
                    }
                    for(int i=0; i<n; i++) mat(i,i) = 1;
                    return cache[e] = ApartIR(mat, e.op(1));
                }
                if(!is_a<matrix>(e.op(0))) throw Error("ApartIRC: Not matrix : " + ex2str(e.op(0)));
                auto mat0 = ex_to<matrix>(e.op(0));
                matrix mat(n+2,n);
                int cc = mat0.cols();
                if(cc==n) mat=mat0;
                else {
                    // zero each element
                    for(int r=0; r<n+2; r++) {
                        for(int c=0; c<cc; c++) mat(r,c) = 0;
                    }
                    // n-row from mat0 to mat, note the last 2 rows still 0
                    for(int r=0; r<n; r++) {
                        for(int c=0; c<cc; c++) mat(r,c) = mat0(r,c);
                    }
                    for(int i=0; i<n; i++) {
                        mat(i,cc) = 1;
                        auto r = mat.rank();
                        if(r==n) break;
                        if(r==cc+1) cc++;
                        else mat(i,cc) = 0;
                    }
                    if(mat.rank()!=n) throw Error("ApartIRC failed, NOT full rank.");
                    // last 2 rows from mat0 to mat
                    for(int r=n; r<n+2; r++) {
                        for(int c=0; c<mat0.cols(); c++) mat(r,c) = mat0(r,c);
                    }
                }
                return cache[e] = ApartIR(mat, e.op(1));
            } else return e.map(self);
        });
        return map(expr_in);
    }
    

}


