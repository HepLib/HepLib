/**
 * @file
 * @brief Helpers to call FORM in FC
 */
 
#include "HEP.h"

namespace HepLib {

    //-----------------------------------------------------------
    // Extend Parser for form, copied from ginac/parser of GiNaC
    //-----------------------------------------------------------
    namespace {

        alignas(2) static ex SP_reader(const exvector& ev) {
            return SP(ev[0], ev[1]).subs(SP_map);
        }
        
        alignas(2) static ex LC_reader(const exvector& ev) {
            if(ev.size()==4) return (-I)*LC(ev[0], ev[1], ev[2], ev[3]);
            else return Eps(ev);
        }
        
        alignas(2) static ex SUNT_reader(const exvector& ev) {
            int n = ev.size();
            if(n==3) return SUNT(ev[0], ev[1], ev[2]);
            else if(n>3) {
                lst as;
                for(int i=0; i<n-2; i++) as.append(ev[i]);
                return SUNT(as,ev[n-2],ev[n-1]);
            } else throw Error("SUNT_reader: number of arguments less than 3.");
        }
        
        alignas(2) static ex TTR_reader(const exvector& ev) {
            if(ev.size()==1) return TTR(ev[0]);
            lst as;
            for(auto item : ev) as.append(item);
            return TTR(as);
        }
        
        alignas(2) static ex DGamma_reader(const exvector& ev) {
            int n = ev.size();
            if(n==1) return GAS(ev[0]);
            ex ret = 1;
            int rl = ex2int(ev[0]);
            for(int i=1; i<n; i++) ret = ret * GAS(ev[i],rl);
            return ret;
        }
        
        alignas(2) static ex gi_reader(const exvector& ev) {
            int n = ev.size();
            if(n==0) return GAS(1);
            else if(n==1) return GAS(1,ex2int(ev[0]));
            else throw Error("DGamma_reader: number of arguments more than 2.");
        }
        
        alignas(2) static ex Matrix_reader(const exvector& ev) {
            return Matrix(ev[0],ev[1],ev[2]);
        }
    }
    
    //-----------------------------------------------------------
    // Run Form Program
    //-----------------------------------------------------------
    namespace {
    
        struct mapGamma : public map_function {
        public:
            ex operator()(const ex &e) {
                if(is_a<DGamma>(e)) return DGamma(ex_to<DGamma>(e), gline);
                else if(!DGamma::has(e)) return e;
                else return e.map(*this);
            }
            mapGamma(int _gline) : gline(_gline) { };
        private:
            unsigned gline = 0;
        };
        
        struct mapTR : public map_function {
        public:
            ex operator()(const ex &e) {
                if (e.match(TR(w))) {
                    ex trs = e.op(0);
                    gline++;
                    trs = mapGamma(gline)(trs);
                    trs = DGamma(1, gline) * trs;
                    if(glmax<gline) {
                        glmax = gline;
                        if(glmax>128) throw Error("too large index with glmax>128.");
                    }
                    return TR(trs);
                } else if(is_a<add>(e)) {
                    ex res = 0;
                    unsigned gl = gline;
                    unsigned glmax = gline;
                    for(auto item : e) {
                        gline = gl;
                        res += (*this)(item);
                        if(glmax<gline) glmax = gline;
                    }
                    gline = glmax;
                    return res;
                } else return e.map(*this);
            }
            unsigned glmax = 0;
        private:
            unsigned gline = 0;
        };
        
        string init_script = R"EOF(
CFunction pow,sqrt,gamma,HF,Matrix,WF;
Tensor TTR(cyclic), f(antisymmetric), T, f4, colTp;
Symbols reX,I2R,NF,NA,d,I,Pi;
AutoDeclare Symbols gCF, trcN;
Dimension NA;
AutoDeclare Index colA;
Dimension NF;
AutoDeclare Index colF;

#procedure SUNTrace
Dimension NF;

repeat;
	id,once,TTR(?a) = T(?a,colF1,colF1);
	sum colF1;
	repeat;
		id,once,T(colA1?,colA2?,?a,colF1?,colF2?) = T(colA1,colF1,colF3)*T(colA2,?a,colF3,colF2);
		sum colF3;
        id,once,f4(colA1?,colA2?,colA3?,colA4?) = f(colA1,colA2,colA5) * f(colA5,colA3,colA4);
        sum colA5;
    endrepeat;
endrepeat;

#do colXi = 1,1
    if ( count(f,1) || match(T(colA1?,colF1?,colF2?)*T(colA1?,colF3?,colF4?)) ) redefine colXi "0";
    id,once,f(colA1?,colA2?,colA3?) = 1/I2R/i_*T(colA1,colF1,colF2)*T(colA2,colF2,colF3)*T(colA3,colF3,colF1)-1/I2R/i_*T(colA3,colF1,colF2)*T(colA2,colF2,colF3)*T(colA1,colF3,colF1);
    sum colF1,colF2,colF3;
    id T(colA1?,colF1?,colF2?)*T(colA1?,colF3?,colF4?) = colTp(colF1,colF2,colF3,colF4);
    #do colXj = 1,1
        if ( count(colTp,1) ) redefine colXj "0";
        .sort
        id,once,colTp(colF1?,colF2?,colF3?,colF4?) = I2R*(d_(colF1,colF4)*d_(colF2,colF3)-d_(colF1,colF2)*d_(colF3,colF4)/NF);
    #enddo
#enddo

repeat;
	id T(colA1?,?a,colF1?,colF2?)*T(colA2?,?b,colF2?,colF3?) = T(colA1,?a,colA2,?b,colF1,colF3);
endrepeat;
id	T(?a,colF1?,colF1?) = TTR(?a);
id	TTR(colA1?) = 0;
id	TTR(colA1?,colA2?) = I2R*d_(colA1,colA2);
.sort

#endprocedure
.global
        )EOF";
    }
    
    //-----------------------------------------------------------
    // form 
    //-----------------------------------------------------------
    namespace {
    ex runform(const ex &expr_in, int verb) {
        static map<pid_t, Form> form_map;
        auto pid = getpid();
        if((form_map.find(pid)==form_map.end())) { // init section
            ostringstream ss;
            ss << init_script << endl;
            form_map[pid].Init("form");
            form_map[pid].Execute(ss.str());
        }
        Form &fprc = form_map[pid];
        
        ex expr = expr_in.subs(SP_map);
        ex all_expr = expr;
        stringstream sss;
        FormFormat ids(sss);
        for(auto kv : SP_map) {
            ids << "id " << kv.first << "=" << kv.second << ";" << endl;
            all_expr += iWF(kv.first) * iWF(kv.second);
        }
        string idstr = sss.str();
        
        lst vec_lst, VD_lst, CF_lst, CA_lst, sym_lst;
        for(const_preorder_iterator i = all_expr.preorder_begin(); i != all_expr.preorder_end(); ++i) {
            if(is_a<Vector>(*i)) vec_lst.append(*i);
            else if(is_a<Index>(*i)) {
                if(ex_to<Index>(*i).type==Index::Type::VD) VD_lst.append(*i);
                else if(ex_to<Index>(*i).type==Index::Type::CF) CF_lst.append(*i);
                else if(ex_to<Index>(*i).type==Index::Type::CA) CA_lst.append(*i);
            } else if(is_a<symbol>(*i)) sym_lst.append(*i);
            else if(is_a<GiNaC::function>(*i)) {
                static vector<string> fun_vec = { 
                    "iWF", "WF", "TR", "sin", "cos", "HF", "TTR", "Matrix"
                };
                auto func = ex_to<GiNaC::function>(*i).get_name();
                bool ok = false;
                for(auto fi : fun_vec) {
                    if(fi==func) { ok = true; break; }
                }
                if(!ok) {
                    cout << (*i) << endl;
                    throw Error("runform: Functions not defined in FORM: "+func);
                }
            }
        }
        vec_lst.sort(); vec_lst.unique();
        VD_lst.sort(); VD_lst.unique();
        CF_lst.sort(); CF_lst.unique();
        CA_lst.sort(); CA_lst.unique();
        sym_lst.append(d);
        sym_lst.sort(); sym_lst.unique();
        
        stringstream ss;
        FormFormat ff(ss);
        ff << ".store" << endl;
        symtab st;
        if(sym_lst.nops()>0) {
            ff << "Symbols";
            for(auto sx : sym_lst) {
                auto s = ex_to<symbol>(sx);
                ff << " " << s.get_name();
                st[s.get_name()] = sx;
            }
            ff << ";" << endl;
        }
        if(vec_lst.nops()>0) {
            ff << "Vectors";
            for(auto vx : vec_lst) {
                auto v = ex_to<Vector>(vx);
                ff << " " << v;
                st[v.name.get_name()] = v;
                for(auto vvx : vec_lst) {
                    auto vv = ex_to<Vector>(vvx);
                    st[v.name.get_name()+"__"+vv.name.get_name()] = SP(v,vv).subs(SP_map);
                }
            }
            ff << ";" << endl;
        }
        if(VD_lst.nops()>0) {
            if(form_using_dim4) ff << "Dimension 4;" << endl;
            else ff << "Dimension d;" << endl;
            ff << "Indices";
            for(auto ix : VD_lst) {
                auto i = ex_to<Index>(ix);
                ff << " " << i;
                st[i.name.get_name()] = i;
                auto itr = Index::Dimension.find(ix);
                if(itr!=Index::Dimension.end()) ff << "=" << itr->second;
            }
            ff << ";" << endl;
        }
        if(CF_lst.nops()>0) {
            ff << "Dimension NF;" << endl;
            ff << "Indices";
            for(auto ix : CF_lst) {
                auto i = ex_to<Index>(ix);
                ff << " " << i;
                st[i.name.get_name()] = i;
                auto itr = Index::Dimension.find(ix);
                if(itr!=Index::Dimension.end()) ff << "=" << itr->second;
            }
            ff << ";" << endl;
        }
        if(CA_lst.nops()>0) {
            ff << "Dimension NA;" << endl;
            ff << "Indices";
            for(auto ix : CA_lst) {
                auto i = ex_to<Index>(ix);
                ff << " " << i;
                st[i.name.get_name()] = i;
                auto itr = Index::Dimension.find(ix);
                if(itr!=Index::Dimension.end()) ff << "=" << itr->second;
            }
            ff << ";" << endl;
        }
        ff << ".global" << endl;
        ff << endl;
        
        // trace and contract
        bool islst = is_a<lst>(expr);
        lst expr_lst;
        if(islst) expr_lst = ex_to<lst>(expr);
        else expr_lst.append(expr);
        auto total = expr_lst.nops();
        
        string ostr;
        int gid = 1;
        ostr = "{";
        int c = 1;
        int kid = 0;
        map<ex,int,ex_is_less> e2i_map;
        for(auto it : expr_lst) {
            ff << ".store" << endl;
            ex item = it;
            item = MapFunction([](const ex & e, MapFunction &self)->ex{
                if(!e.has(TR(w))) return e;
                else if(e.match(TR(w))) {
                    ex nd = e.op(0).normal().numer_denom();
                    return TR(nd.op(0))/nd.op(1);
                } else return e.map(self);
            })(item);
            item = normal(item); // TODO: normal
            // pull out global common factor
            item = collect_common_factors(item);
            item = CoPat(item,[](const ex &e)->bool{return Index::has(e) || DGamma::has(e) || Eps::has(e);});
            auto ckey = item.op(0);
            if(!is_a<numeric>(ckey)) {
                int ckid;
                if(e2i_map.find(ckey)==e2i_map.end()) {
                    kid++;
                    e2i_map[ckey] = kid;
                    ckid = kid;
                } else ckid = e2i_map[ckey];
                string gCF = "gCF" + to_string(ckid);
                st[gCF] = item.op(0);
                item = Symbol(gCF) * item.op(1);
            } else {
                item = ckey * item.op(1);
            }
            
            // pull out color factor
            auto cv_lst = collect_lst(item, [](const ex &e)->bool{return Index::hasc(e);});
            item=0;
            exvector color_vec;
            for(int i=0; i<cv_lst.nops(); i++) {
                auto it = cv_lst.op(i);
                auto cc = it.op(0);
                auto vv = it.op(1);
                color_vec.push_back(vv);
                item += cc * Symbol("[cl"+to_string(i)+"]");
            }
            
            ostringstream coss;
            for(int i=0; i<color_vec.size(); i++) {
                coss << " [cl" << i << "]";
                ff << "L [cl" << i << "]=" << color_vec[i] << ";" << endl;
                ff << ".sort" << endl;
                ff << "#call SUNTrace" << endl;
                ff << ".sort" << endl;
                if(form_using_su3) {
                    ff << "id NF^reX?=3^reX;" << endl;
                    ff << "id I2R=1/2;" << endl;
                    ff << "id NA^reX?=8^reX;" << endl;
                    ff << ".sort" << endl;
                }
            }
            ff << endl;
            
            // methods to handle TR objects
            auto tr_mode = form_trace_mode;
            exset trs;
            find(item,TR(w),trs);
            if(form_trace_mode==form_trace_auto) {
                if(trs.size()<2) tr_mode = form_trace_all;
                else tr_mode = form_trace_each_each;
            }
            if(tr_mode==form_trace_all) {
                mapTR tr;
                item = tr(item);
                ff << "L [o]=" << item << ";" << endl;
                ff << ".sort" << endl;
                ff << "#do i = 1,1" << endl;
                ff << "id once HF(reX?)=reX;" << endl;
                ff << idstr;
                ff << "if(count(HF,1)>0) redefine i \"0\";" << endl;
                ff << ".sort" << endl;
                ff << "#enddo" << endl;
                ff << "contract 0;" << endl;
                ff << ".sort" << endl;
                ff << idstr << ".sort" << endl;
                ff << endl;
                for(int gl=1; gl<=tr.glmax; gl++) {
                    if(form_using_gamma5) {
                        ff << ".sort" << endl;
                        if(form_using_dim4) ff << "Dimension 4;" << endl;
                        else ff << "Dimension d;" << endl;
                        ff << "Indices [g5_i1], [g5_i2], [g5_i3], [g5_i4];" << endl;
                        ff << "id g_(" << gl << ",5_) = e_([g5_i1],[g5_i2],[g5_i3],[g5_i4])*g_(" << gl << ",[g5_i1],[g5_i2],[g5_i3],[g5_i4])/24;" << endl;
                        ff << "sum [g5_i1],[g5_i2],[g5_i3],[g5_i4];" << endl;
                        ff << ".sort" << endl;
                    }
                    if(form_using_dim4) ff << "trace4 " << gl << ";" << endl;
                    else ff << "tracen " << gl << ";" << endl;
                    ff << ".sort" << endl;
                    ff << idstr << ".sort" << endl;
                }
            } else if(tr_mode==form_trace_each_all) {
                exmap tr2v;
                int trn=0;
                exvector trvec;
                for(auto tr : trs) {
                    tr2v[tr] = Symbol("[tr"+to_string(trn)+"]");
                    auto tri = mapGamma(gid)(tr.op(0));
                    trvec.push_back(tri);
                    trn++;
                }
                item = item.subs(tr2v);
                for(int i=0; i<trvec.size(); i++) {
                    ff << "L [tr" << i << "]=" << trvec[i] << ";" << endl;
                    if(form_using_gamma5) {
                        ff << ".sort" << endl;
                        if(form_using_dim4) ff << "Dimension 4;" << endl;
                        else ff << "Dimension d;" << endl;
                        ff << "Indices [g5_i1], [g5_i2], [g5_i3], [g5_i4];" << endl;
                        ff << "id g_(" << gid << ",5_) = e_([g5_i1],[g5_i2],[g5_i3],[g5_i4])*g_(" << gid << ",[g5_i1],[g5_i2],[g5_i3],[g5_i4])/24;" << endl;
                        ff << "sum [g5_i1],[g5_i2],[g5_i3],[g5_i4];" << endl;
                    }
                    if(form_using_dim4) ff << "trace4 " << gid << ";" << endl;
                    else ff << "tracen " << gid << ";" << endl;
                    ff << ".sort" << endl;
                    ff << idstr << ".sort" << endl;
                }
                ff << "L [o]=" << item << ";" << endl;
                ff << ".sort" << endl;
                ff << "#do i = 1,1" << endl;
                ff << "id once HF(reX?)=reX;" << endl;
                ff << idstr;
                ff << "if(count(HF,1)>0) redefine i \"0\";" << endl;
                ff << ".sort" << endl;
                ff << "#enddo" << endl;
                ff << idstr << ".sort" << endl;
                ff << endl;
            } else if(tr_mode==form_trace_each_each) {
                exmap tr2v;
                int trn=0;
                exvector trvec;
                for(auto tr : trs) {
                    tr2v[tr] = Symbol("trcN"+to_string(trn));
                    auto tri = mapGamma(gid)(tr.op(0));
                    trvec.push_back(tri);
                    trn++;
                }
                item = item.subs(tr2v);

                if(color_vec.size()>0) ff << "drop" << coss.str() << ";" << endl;
                ff << "L [o]=" << item << ";" << endl;
                ff << ".sort" << endl;
                ff << "#do i = 1,1" << endl;
                ff << "id once HF(reX?)=reX;" << endl;
                ff << idstr;
                ff << "if(count(HF,1)>0) redefine i \"0\";" << endl;
                ff << ".sort" << endl;
                ff << "#enddo" << endl;
                ff << endl;

                for(int i=0; i<trvec.size(); i++) {
                    ff << "L [tr" << i << "]=" << trvec[i] << ";" << endl;
                    if(form_using_gamma5) {
                        ff << ".sort" << endl;
                        if(form_using_dim4) ff << "Dimension 4;" << endl;
                        else ff << "Dimension d;" << endl;
                        ff << "Indices [g5_i1], [g5_i2], [g5_i3], [g5_i4];" << endl;
                        ff << "id g_(" << gid << ",5_) = e_([g5_i1],[g5_i2],[g5_i3],[g5_i4])*g_(" << gid << ",[g5_i1],[g5_i2],[g5_i3],[g5_i4])/24;" << endl;
                        ff << "sum [g5_i1],[g5_i2],[g5_i3],[g5_i4];" << endl;
                    }
                    if(form_using_dim4) ff << "trace4 " << gid << ";" << endl;
                    else ff << "tracen " << gid << ";" << endl;
                    ff << ".sort" << endl;
                    ff << idstr << ".sort" << endl;
                    ff << "drop [tr" << i << "];" << endl;
                    ff << "id trcN" << i << "=[tr" << i << "];" << endl;
                    ff << ".sort" << endl;
                    ff << idstr << ".sort" << endl;
                    ff << endl;
                }
            } else {
                throw Error("runform: unsupported form_trace_mode = " + to_string(form_trace_mode));
            }
            
            ff << "contract 0;" << endl;
            ff << ".sort" << endl;
            ff << idstr << ".sort" << endl;
            
            if(verb==1) {
                cout << "\r                                     \r";
                cout << PRE << "\\--Form Script @ " << c << " / " << total << flush;
            } else if(verb>1) {
                cout << "--------------------------------------" << endl;
                cout << "Form Script @ " << c << " / " << total << endl;
                cout << "--------------------------------------" << endl;
                cout << ss.str() << endl;
            }
            
            auto script = ss.str();
            string_replace_all(script, "sin(", "sin_(");
            string_replace_all(script, "cos(", "cos_(");
            try {
                auto otmp = fprc.Execute(script);
                if(verb>2) {
                    cout << "--------------------------------------" << endl;
                    cout << "Form Output @ " << c << " / " << total << endl;
                    cout << "--------------------------------------" << endl;
                    cout << otmp << endl;
                }
            
                ostr += otmp;
                ss.clear();
                ss.str("");
                
                if(c<total) ostr += ",";
                c++;
            } catch(Error& err) {
                form_map.erase(pid);
                throw;
            }
        }
        if(verb==1) cout << endl;
        ostr += "}";
                
        string_replace_all(ostr, "[", "(");
        string_replace_all(ostr, "]", ")");
        string_replace_all(ostr, "\\\n", "");
        string_replace_all(ostr, " ","");
        for(auto v : vec_lst) {
            string pat(ex_to<Vector>(v).name.get_name());
            string from = pat+"(";
            string to = "d_("+pat+",";
            string_replace_all(ostr, from, to);
            from = pat+".";
            to = pat+"__";
            string_replace_all(ostr, from, to);
        }
        
        string_replace_all(ostr, "d_(", "SP(");
        string_replace_all(ostr, "e_(", "LC(");
        string_replace_all(ostr, "sin_(", "sin(");
        string_replace_all(ostr, "cos_(", "cos(");
        string_replace_all(ostr, "g_(", "DG(");
        string_replace_all(ostr, ",5_", ",5");
        string_replace_all(ostr, ",6_", ",6");
        string_replace_all(ostr, ",7_", ",7");
        string_replace_all(ostr, "gi_(", "GI(");
        string_replace_all(ostr, "TTR(", "TTRX(");
        
        st["I2R"] = TF;
        st["NA"] = NA;
        st["NF"] = NF;
        st["I"] = I;
        st["i_"] = I;

        Parser fp(st);
        fp.FTable[make_pair("SP", 2)] = SP_reader;
        fp.FTable[make_pair("LC", 4)] = LC_reader;
        fp.FTable[make_pair("Matrix", 3)] = Matrix_reader;
        for(int i=1; i<30; i++) fp.FTable[make_pair("T", i)] = SUNT_reader;
        for(int i=1; i<30; i++) fp.FTable[make_pair("TTRX", i)] = TTR_reader;
        for(int i=1; i<30; i++) fp.FTable[make_pair("DG", i)] = DGamma_reader;
        fp.FTable[make_pair("GI", 0)] = gi_reader;
        fp.FTable[make_pair("GI", 1)] = gi_reader;
        ex ret = fp.Read(ostr);
        if(!islst) ret = ret.op(0);
        return ret;
    }}
    
    /**
     * @brief evalulate expr in form program, see also the form_trace_mode and form_expand_mode
     * @param expr the input expression
     * @param verb for verb output
     * @return result with index contract, trace performed, etc.
     */
    ex form(const ex & iexpr, int verb) {
        ex expr = iexpr;
        if(!is_a<lst>(expr)) {
            static MapFunction mf([](const ex & e, MapFunction & self)->ex{
                if(!e.has(TR(w))) return e;
                else if(e.match(TR(w))) {
                    if(is_a<mul>(e.op(0))) {
                        ex c = 1, v = 1;
                        for(auto const & ei : e.op(0)) {
                            if(DGamma::has(ei)) v *= ei;
                            else c *= ei;
                        }
                        return c*TR(v);
                    } else return e;
                } else return e.map(self);
            });
            expr = mf(expr);
        }
        if(form_expand_mode==form_expand_none || is_a<lst>(expr)) return runform(expr, verb);
        else if(form_expand_mode==form_expand_tr) {
            auto cv_lst = collect_lst(expr.subs(SP_map), TR(w));
            lst to_lst;
            for(auto cv : cv_lst) to_lst.append(cv.op(0)*cv.op(1));
            lst out_lst = ex_to<lst>(runform(to_lst, verb));
            
            ex ret = 0;
            for(int i=0; i<cv_lst.nops(); i++) ret += out_lst.op(i);
            
            return ret.subs(SP_map);
        } else if(form_expand_mode==form_expand_ci) {
            auto cv_lst = collect_lst(expr.subs(SP_map), [](const ex & e)->bool {
                return e.has(TR(w)) || Index::hasc(e);
            });
            lst to_lst;
            for(auto cv : cv_lst) to_lst.append(cv.op(0)*cv.op(1));
            lst out_lst = ex_to<lst>(runform(to_lst, verb));
            
            ex ret = 0;
            for(int i=0; i<cv_lst.nops(); i++) ret += out_lst.op(i);
            
            return ret.subs(SP_map);
        } else if(form_expand_mode==form_expand_li) {
            auto cv_lst = collect_lst(expr.subs(SP_map), [](const ex & e)->bool {
                return e.has(TR(w)) || Index::hasv(e);
            });
            lst to_lst;
            for(auto cv : cv_lst) to_lst.append(cv.op(0)*cv.op(1));
            lst out_lst = ex_to<lst>(runform(to_lst, verb));
            
            ex ret = 0;
            for(int i=0; i<cv_lst.nops(); i++) ret += out_lst.op(i);
            
            return ret.subs(SP_map);
        } else if(form_expand_mode==form_expand_all) {
            auto cv_lst = collect_lst(expr.subs(SP_map), [](const ex & e)->bool {
                return e.has(TR(w)) || SUNT::has(e) || SUNF::has(e) || SUNF4::has(e) || Index::has(e) || DGamma::has(e);
            });
            lst to_lst;
            for(auto cv : cv_lst) to_lst.append(cv.op(1));
            lst out_lst = ex_to<lst>(runform(to_lst, verb));
            
            ex ret = 0;
            for(int i=0; i<cv_lst.nops(); i++) ret += cv_lst.op(i).op(0) * out_lst.op(i);
            
            return ret.subs(SP_map);
        } else throw Error("form: unsupported form_expand_mode = " + to_string(form_expand_mode));
        return 0;
    }

    /**
     * @brief make the charge conjugate operaton, M -> C^{-1} . M^T . C w.r.t. a Matrix object
     * @param expr the input expression
     * @return returned charge conjugated expression
     */
    ex charge_conjugate(const ex & expr) {
        if(expr.has(Matrix(w1,w2,w3))) throw Error("charge_conjugate: Matrix found.");
        if(!DGamma::has(expr)) return expr;
        if(is_a<DGamma>(expr)) {
            DGamma g = ex_to<DGamma>(expr);
            if(is_a<Vector>(g.pi) || is_a<Index>(g.pi)) return -expr;
            else if(is_zero(g.pi-1) || is_zero(g.pi-5)) return expr;
        }
        
        if(is_a<add>(expr)) {
            ex ret = 0;
            for(auto item : expr) ret += charge_conjugate(item);
            return ret;
        } else if(is_a<mul>(expr)) {
            ex ret = 1;
            for(auto item : expr) ret *= charge_conjugate(item);
            return ret;
        } else if(is_a<ncmul>(expr)) {
            ex ret = 1;
            int n = expr.nops();
            for(int i=n-1; i>-1; i--) ret *= charge_conjugate(expr.op(i));
            return ret;
        }
        cout << DGamma::has(expr) << " : " << expr << endl;
        throw Error("charge_conjugate: unexpected region.");
        return 0;
    }

}

