#include "FC.h"

namespace HepLib::FC {
    
    //-----------------------------------------------------------
    // Extend Parser for form, copied from ginac/parser of GiNaC
    //-----------------------------------------------------------
    namespace {

        alignas(2) static ex SP_reader(const exvector& ev) {
            return SP(ev[0], ev[1]).subs(sp_map);
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
                if (e.match(TR(w))) {
                    ex gs = e.op(0);
                    gline++;
                    ss << "tracen " << gline << ";" << endl;
                    gs = mapGamma(gline)(gs);
                    gs = DiracGamma(1, gline) * gs;
                    return TR(gs);
                } else if (!e.has(TR(w))) return e;
                else return e.map(*this);
            }
            stringstream ss;
        private:
            unsigned gline = 0;
        };
        
        void string_replace_all(string &str, const string &from, const string &to) {
            size_t start_pos = 0;
            while((start_pos = str.find(from, start_pos)) != string::npos) {
                str.replace(start_pos, from.length(), to);
                start_pos += to.length();
            }
        }
        
    }
    
    //-----------------------------------------------------------
    // form 
    //-----------------------------------------------------------
    static HepLib::Form form_proc;
    static map<pid_t, bool> init_map;
    
    namespace {
    ex runform(const ex &expr_in, bool verb) {
        if((init_map.find(PID)==init_map.end()) || !init_map[PID]) { // init section
            ostringstream ss;
            ss << "CFunction pow,sqrt,Pi,gamma;" << endl;
            ss << "Symbols D, NF, NA, I;" << endl;
            ss << "#include- " << InstallPrefix <<  "/include/SUN.h" << endl;
            form_proc.Init("form");
            form_proc.Execute(ss.str());
            init_map[PID] = true;
        }
        
        ex expr = expr_in.subs(sp_map);
        lst vec_lst, VD_lst, CF_lst, CA_lst, sym_lst;
        for(const_preorder_iterator i = expr.preorder_begin(); i != expr.preorder_end(); ++i) {
            if(is_a<Vector>(*i)) vec_lst.append(*i);
            else if(is_a<Index>(*i)) {
                if(ex_to<Index>(*i).type==Index::Type::VD) VD_lst.append(*i);
                else if(ex_to<Index>(*i).type==Index::Type::CF) CF_lst.append(*i);
                else if(ex_to<Index>(*i).type==Index::Type::CA) CA_lst.append(*i);
            } else if(is_a<symbol>(*i)) sym_lst.append(*i);
            else if(is_a<GiNaC::function>(*i)) {
                static vector<string> fun_vec = { 
                    "TR", "sin", "cos"
                };
                auto func = ex_to<GiNaC::function>(*i).get_name();
                bool ok = false;
                for(auto fi : fun_vec) {
                    if(fi==func) { ok = true; break; }
                }
                if(!ok) throw Error("runform: Functions not defined in FORM: "+func);
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
                    st[v.name.get_name()+"_"+vv.name.get_name()] = SP(v,vv).subs(sp_map);
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
        string ostr;
        if(is_a<lst>(expr)) {
            auto total = expr.nops();
            ostr = "{";
            int c = 1;
            for(auto item : expr) {
                mapTR tr;
                item = tr(item);
                ff << "L [o]=" << item << ";" << endl;
                ff << "contract 0;" << endl;
                ff << "#call SUNTrace" << endl;
                ff << tr.ss.str();
                ff << "contract 0;" << endl;
                
                if(verb) {
                    cout << "--------------------------------------" << endl;
                    cout << "Form Script @" << c-1 << " / " << total-1 << endl;
                    cout << "--------------------------------------" << endl;
                    cout << ss.str() << endl;
                }
                
                auto script = ss.str();
                string_replace_all(script, "sin(", "sin_(");
                string_replace_all(script, "cos(", "cos_(");
                auto otmp = form_proc.Execute(script);
            
                if(verb) {
                    cout << "--------------------------------------" << endl;
                    cout << "Form Output @" << c-1 << " / " << total-1 << endl;
                    cout << "--------------------------------------" << endl;
                    cout << otmp << endl;
                }
            
                ostr += otmp;
                ss.clear();
                ss.str("");
                
                if(c<expr.nops()) ostr += ",";
                c++;
            }
            ostr += "}";
        } else {
            mapTR tr;
            expr = tr(expr);
            ff << "L [o]=" << expr << ";" << endl;
            ff << "contract 0;" << endl;
            ff << "#call SUNTrace" << endl;
            ff << tr.ss.str();
            ff << "contract 0;" << endl;
            
            if(verb) {
                cout << "--------------------------------------" << endl;
                cout << "Form Script:" << endl;
                cout << "--------------------------------------" << endl;
                cout << ss.str() << endl;
            }
            auto script = ss.str();
            string_replace_all(script, "sin(", "sin_(");
            string_replace_all(script, "cos(", "cos_(");
            ostr = form_proc.Execute(script);
            ss.clear();
            ss.str("");
            
            if(verb) {
                cout << "--------------------------------------" << endl;
                cout << "Form Output:" << endl;
                cout << "--------------------------------------" << endl;
                cout << ostr << endl;
            }
        }
        
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
        return ret;
    }}
    
    /**
     * @brief evalulate expr in form program
     * @param expr the input expression
     * @param all true for sending all into form, otherwise will use mma_collect w.r.t. Index/DiracGamma
     * @return 返回说明
     * @retval result with index contract, trace performed, etc.
     */
    ex form(const ex &expr, bool all, bool verb) {
        if(all || is_a<lst>(expr)) return runform(expr, verb);
        
        if(expr.has(coVF(w))) throw Error("form error: expr has coVF already.");
        auto ret = mma_collect(expr.subs(sp_map), [](const ex & e)->bool {
            return Index::has(e) || DiracGamma::has(e);
        },false,true);
        
        lst to_lst;
        int current = 0;
        ret = MapFunction([&](const ex & e, MapFunction &self)->ex{
            if(e.match(coVF(w))) {
                to_lst.append(e.op(0));
                return coVF(current++);
            } else if (!e.has(coVF(w))) return e;
            else return e.map(self);
        })(ret);
        
        lst out_lst = ex_to<lst>(runform(to_lst, verb));
        ret = MapFunction([&](const ex & e, MapFunction &self)->ex{
            if(e.match(coVF(w))) {
                return out_lst.op(ex_to<numeric>(e.op(0)).to_int());
            } else if (!e.has(coVF(w))) return e;
            else return e.map(self);
        })(ret);
        return ret.subs(sp_map);
    }


}

