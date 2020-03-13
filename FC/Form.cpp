#include "FC.h"

namespace HepLib::FC {
    
    //-----------------------------------------------------------
    // Extend Parser for form, copied from ginac/parser of GiNaC
    //-----------------------------------------------------------
    namespace {
    
        auto init_first = get_default_reader(); // IMPORTANT: initialize first
    
        ex SP_reader(const exvector& ev) {
            return SP(ev[0], ev[1]);
        }
        
        ex LC_reader(const exvector& ev) {
            return LC(ev[0], ev[1], ev[2], ev[3]);
        }
        
        const prototype_table& form_func_reader() {
            using std::make_pair;
            static bool initialized = false;
            static prototype_table dr;
            if(!initialized) {
                dr = get_default_reader();
                dr[make_pair("SP", 2)] = SP_reader;
                dr[make_pair("LC", 4)] = LC_reader;
                initialized = true;
            }
            return dr;
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
                } else {
                    return e.map(*this);
                }
            }
            stringstream ss;
        private:
            unsigned gline = 0;
        };
        
        void replace_all(string &str, const string &from, const string &to) {
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
    ex form(ex expr, bool verb) {
        lst vec_lst, VD_lst, CF_lst, CA_lst, sym_lst;
        for(const_preorder_iterator i = expr.preorder_begin(); i != expr.preorder_end(); ++i) {
            if(is_a<Vector>(*i)) vec_lst.append(*i);
            else if(is_a<Index>(*i)) {
                if(ex_to<Index>(*i).IndexType==Index::Type::VD) VD_lst.append(*i);
                else if(ex_to<Index>(*i).IndexType==Index::Type::CF) CF_lst.append(*i);
                else if(ex_to<Index>(*i).IndexType==Index::Type::CA) CA_lst.append(*i);
            } else if(is_a<symbol>(*i)) sym_lst.append(*i);
        }
        vec_lst.sort(); vec_lst.unique();
        VD_lst.sort(); VD_lst.unique();
        CF_lst.sort(); CF_lst.unique();
        CA_lst.sort(); CA_lst.unique();
        sym_lst.append(D);
        sym_lst.sort(); sym_lst.unique();
        
        stringstream ss;
        FormFormat ff(ss);
        ff << "#include- " << InstallPrefix <<  "/include/color.h" << endl;
        ff << "Symbols D, NR, NA;" << endl;
        symtab st;
        if(vec_lst.nops()>0) {
            ff << "Vectors";
            for(auto vx : vec_lst) {
                auto v = ex_to<Vector>(vx);
                ff << " " << v;
                st[v.name.get_name()] = v;
                for(auto vvx : vec_lst) {
                    auto vv = ex_to<Vector>(vvx);
                    st[v.name.get_name()+"_"+vv.name.get_name()] = v(vv);
                }
                for(auto ix : VD_lst) {
                    auto i = ex_to<Index>(ix);
                    st[v.name.get_name()+"_"+i.name.get_name()] = v(i);
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
                for(auto jx : VD_lst) {
                    auto j = ex_to<Index>(jx);
                    st[i.name.get_name()+"_"+j.name.get_name()] = i(j);
                }
            }
            ff << ";" << endl;
        }
        if(CF_lst.nops()>0) {
            ff << "Dimension NR;" << endl;
            ff << "Indices";
            for(auto ix : CF_lst) {
                auto i = ex_to<Index>(ix);
                ff << " " << i;
                st[i.name.get_name()] = i;
                for(auto jx : CF_lst) {
                    auto j = ex_to<Index>(jx);
                    st[i.name.get_name()+"_"+j.name.get_name()] = i(j);
                }
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
                for(auto jx : CA_lst) {
                    auto j = ex_to<Index>(jx);
                    st[i.name.get_name()+"_"+j.name.get_name()] = i(j);
                }
            }
            ff << ";" << endl;
        }
        if(sym_lst.nops()>0) {
            ff << "Symbols";
            for(auto sx : sym_lst) {
                auto s = ex_to<symbol>(sx);
                ff << " " << s;
                st[s.get_name()] = s;
            }
            ff << ";" << endl;
        }
        
        // trace and contract
        mapTR tr;
        expr = tr(expr);
        ff << "L [o]=" << expr << ";" << endl;
        ff << "contract 0;" << endl;
        ff << "#call docolor" << endl;
        ff << tr.ss.str();
        ff << "contract 0;" << endl;
        
        if(verb) {
            cout << "--------------------------------------" << endl;
            cout << "Form Script:" << endl;
            cout << "--------------------------------------" << endl;
            cout << ss.str() << endl;
        }

        HepLib::Form form;
        form.Init(string(InstallPrefix)+string("/bin/form"));
        auto ostr = form.Execute(ss.str());
        form.Exit();
        
        if(verb) {
            cout << "--------------------------------------" << endl;
            cout << "Form Output:" << endl;
            cout << "--------------------------------------" << endl;
            cout << ostr << endl;
        }

        for(auto v : vec_lst) {
            string pat(ex_to<Vector>(v).name.get_name());
            string from = pat+"(";
            string to = "("+pat+"_";
            replace_all(ostr, from, to);
            from = pat+".";
            to = pat+"_";
            replace_all(ostr, from, to);
        }
        
        replace_all(ostr, "d_(", "SP(");
        replace_all(ostr, "e_(", "LC(");
        
        st["I2R"] = ex(1)/2;
        st["cA"] = CA;
        st["cR"] = CF;
        st["NA"] = NA;
        st["NR"] = NF;
        
        parser form_parser(st, false, form_func_reader());
        ex ret = form_parser(ostr);
        return ret;
    }


}

