/**
 * @file
 * @brief Helpers to call FORM in FC
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "FC.h"

namespace HepLib::FC {

    //-----------------------------------------------------------
    // Extend Parser for form, copied from ginac/parser of GiNaC
    //-----------------------------------------------------------
    namespace {

        alignas(2) static ex SP_reader(const exvector& ev) {
            return SP(ev[0], ev[1]).subs(SP_map);
        }
        
        alignas(2) static ex LC_reader(const exvector& ev) {
            return LC(ev[0], ev[1], ev[2], ev[3]);
        }
        
        alignas(2) static ex SUNT_reader(const exvector& ev) {
            return SUNT(ex_to<Index>(ev[0]), ex_to<Index>(ev[1]), ex_to<Index>(ev[2]));
        }
    }
    
    //-----------------------------------------------------------
    // Run Form Program
    //-----------------------------------------------------------
    namespace {
    
        struct mapGamma : public map_function {
        public:
            ex operator()(const ex &e) {
                if(is_a<DiracGamma>(e)) return DiracGamma(ex_to<DiracGamma>(e), gline);
                else if(!DiracGamma::has(e)) return e;
                else return e.map(*this);
            }
            mapGamma(int _gline) : gline(_gline) { };
        private:
            unsigned gline = 0;
        };
        
        struct mapTR : public map_function {
        public:
            ex operator()(const ex &e) {
                if (!e.has(TR(w))) return e;
                else if (e.match(TR(w))) {
                    ex gs = e.op(0);
                    gline++;
                    gs = mapGamma(gline)(gs);
                    gs = DiracGamma(1, gline) * gs;
                    if(glmax<gline) {
                        glmax = gline;
                        if(glmax>128) throw Error("too large index with glmax>128.");
                    }
                    return TR(gs);
                } else if(is_a<add>(e)) {
                    ex res = 0;
                    unsigned gl = gline;
                    for(auto item : e) {
                        gline = gl;
                        res += (*this)(item);
                    }
                    return res;
                } else return e.map(*this);
            }
            unsigned glmax = 0;
        private:
            unsigned gline = 0;
        };
        
        // SU(N) : Phys.Rev., D14, 1536 (1976)
        string init_script = R"EOF(
CFunction pow,sqrt,gamma;
Tensor T,f(antisymmetric);
Tensor colTp;
Symbols gCF,I2R,NF,NA,D,I,Pi;
Dimension NA;
AutoDeclare Index colAj;
Dimension NF;
AutoDeclare Index colFi;

#procedure SUNTrace
Dimension NF;
#do colXi = 1,1
    if ( count(f,1) || match(T(colFi1?,colFi2?,colAj1?)*T(colFi3?,colFi4?,colAj1?)) ) redefine colXi "0";
    id,once,f(colAj1?,colAj2?,colAj3?) = 1/I2R/i_*T(colFi1,colFi2,colAj1)*T(colFi2,colFi3,colAj2)*T(colFi3,colFi1,colAj3)-1/I2R/i_*T(colFi1,colFi2,colAj3)*T(colFi2,colFi3,colAj2)*T(colFi3,colFi1,colAj1);
    sum colFi1,colFi2,colFi3;
    id T(colFi1?,colFi2?,colAj1?)*T(colFi3?,colFi4?,colAj1?) = colTp(colFi1,colFi2,colFi3,colFi4);
    #do colXj = 1,1
        if ( count(colTp,1) ) redefine colXj "0";
        .sort
        id,once,colTp(colFi1?,colFi2?,colFi3?,colFi4?) = I2R*(d_(colFi1,colFi4)*d_(colFi2,colFi3)-d_(colFi1,colFi2)*d_(colFi3,colFi4)/NF);
    #enddo
    #do colXk = 1,1
        if ( match(T(colFi1?,colFi1?,colAj1?)) ) redefine colXk "0";
        .sort
        id,once,T(colFi1?,colFi1?,colAj1?) = 0;
    #enddo
#enddo

#endprocedure
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
        
        ex expr = expr_in.subs(SP_map).subs(HF(w)==w);
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
                    "iWF", "TR", "sin", "cos"
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
        sym_lst.append(D);
        sym_lst.sort(); sym_lst.unique();
        
        stringstream ss;
        FormFormat ff(ss);
        symtab st;
        if(vec_lst.nops()>0) {
            ff << "Vectors";
            for(auto vx : vec_lst) {
                auto v = ex_to<Vector>(vx);
                ff << " " << v;
                st[v.name.get_name()] = v;
                for(auto vvx : vec_lst) {
                    auto vv = ex_to<Vector>(vvx);
                    st[v.name.get_name()+"_"+vv.name.get_name()] = SP(v,vv).subs(SP_map);
                }
            }
            ff << ";" << endl;
        }
        if(VD_lst.nops()>0) {
            ff << "Dimension D;" << endl;
            ff << "Indices";
            for(auto ix : VD_lst) {
                auto i = ex_to<Index>(ix);
                ff << " " << i;
                st[i.name.get_name()] = i;
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
            }
            ff << ";" << endl;
        }
        if(sym_lst.nops()>0) {
            ff << "Symbols";
            for(auto sx : sym_lst) {
                auto s = ex_to<symbol>(sx);
                ff << " " << s.get_name();
                st[s.get_name()] = sx;
            }
            ff << ";" << endl;
        }
        
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
        for(auto it : expr_lst) {
            ex item = it;
            // pull out global common factor
            item = collect_common_factors(item);
            item = CoPat(item,[](const ex &e)->bool{return Index::has(e) || DiracGamma::has(e);});
            st["gCF"] = item.op(0);
            item = Symbol("gCF") * item.op(1);
            // pull out color factor
            auto cv_lst = mma_collect_lst(item, [](const ex &e)->bool{return Index::hasc(e);});
            item=0;
            exvector color_vec;
            for(int i=0; i<cv_lst.nops(); i++) {
                auto it = cv_lst.op(i);
                auto cc = it.op(0);
                auto vv = it.op(1);
                color_vec.push_back(vv);
                item += cc * Symbol("[cl"+to_string(i)+"]");
            }
            for(int i=0; i<color_vec.size(); i++) {
                ff << "L [cl" << i << "]=" << color_vec[i] << ";" << endl;
                ff << ".sort" << endl;
                ff << "#call SUNTrace" << endl;
                ff << ".sort" << endl;
            }
            
            // two method to handle TR objects
            if(trace_method==1) {
                mapTR tr;
                item = tr(item);
                ff << "L [o]=" << item << ";" << endl;
                ff << ".sort" << endl;
                ff << idstr << ".sort" << endl;
                for(int gl=1; gl<=tr.glmax; gl++) {
                    ff << "tracen " << gl << ";" << endl;
                    ff << ".sort" << endl;
                    ff << idstr << ".sort" << endl;
                }
            } else if(trace_method==2) {
                exset trs;
                find(item,TR(w),trs);
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
                    ff << "tracen " << gid << ";" << endl;
                    ff << ".sort" << endl;
                    ff << idstr << ".sort" << endl;
                }
                ff << "L [o]=" << item << ";" << endl;
                ff << ".sort" << endl;
            }
            ff << idstr << ".sort" << endl;
            ff << "contract 0;" << endl;
            ff << ".sort" << endl;
            ff << idstr << ".sort" << endl;
            
            if(verb==1) {
                cout << "\r                                     \r";
                cout << "  \\--Form Script @ " << c << " / " << total << flush;
            } else if(verb>1) {
                cout << "--------------------------------------" << endl;
                cout << "Form Script @ " << c << " / " << total << endl;
                cout << "--------------------------------------" << endl;
                cout << ss.str() << endl;
            }
            
            auto script = ss.str();
            string_replace_all(script, "sin(", "sin_(");
            string_replace_all(script, "cos(", "cos_(");
            auto otmp = fprc.Execute(script);
        
            if(verb>2) {
                cout << "--------------------------------------" << endl;
                cout << "Form Output @" << c << " / " << total << endl;
                cout << "--------------------------------------" << endl;
                cout << otmp << endl;
            }
        
            ostr += otmp;
            ss.clear();
            ss.str("");
            
            if(c<total) ostr += ",";
            c++;
        }
        if(verb==1) cout << endl;
        ostr += "}";
        
        string_replace_all(ostr, "[", "(");
        string_replace_all(ostr, "]", ")");
        for(auto v : vec_lst) {
            string pat(ex_to<Vector>(v).name.get_name());
            string from = pat+"(";
            string to = "SP("+pat+",";
            string_replace_all(ostr, from, to);
            from = pat+".";
            to = pat+"_";
            string_replace_all(ostr, from, to);
        }
        
        string_replace_all(ostr, "d_(", "SP(");
        string_replace_all(ostr, "e_(", "LC(");
        string_replace_all(ostr, "sin_", "sin");
        string_replace_all(ostr, "cos_", "cos");
        
        st["I2R"] = ex(1)/2;
        st["NA"] = NA;
        st["NF"] = NF;
        st["I"] = I;
        st["i_"] = I;

        Parser fp(st);
        fp.FTable[make_pair("SP", 2)] = SP_reader;
        fp.FTable[make_pair("LC", 4)] = LC_reader;
        fp.FTable[make_pair("T", 3)] = SUNT_reader;
        ex ret = fp.Read(ostr);
        if(!islst) ret = ret.op(0);
        return ret;
    }}
    
    /**
     * @brief evalulate expr in form program
     * @param expr the input expression
     * @param all true for sending all into form, otherwise will use mma_collect w.r.t. Index/DiracGamma
     * @return result with index contract, trace performed, etc.
     */
    ex form(const ex &expr, bool all, int verb) {
        if(all || is_a<lst>(expr)) return runform(expr, verb);
        auto cv_lst = mma_collect_lst(expr.subs(SP_map), [](const ex & e)->bool {
            return e.has(TR(w)) || SUNT::has(e) || SUNF::has(e) || Index::has(e) || DiracGamma::has(e);
        });
        lst to_lst;
        for(auto cv : cv_lst) to_lst.append(cv.op(1));
        lst out_lst = ex_to<lst>(runform(to_lst, verb));
        
        ex ret = 0;
        for(int i=0; i<cv_lst.nops(); i++) ret += cv_lst.op(i).op(0) * out_lst.op(i);
        
        return ret.subs(SP_map);
    }


}

