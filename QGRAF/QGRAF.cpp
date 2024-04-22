/**
 * @file
 * @brief QGRAF related
 */

#include "QGRAF.h"

namespace HepLib::QGRAF {

    namespace {
        string di_("di");
        string li_("li");
        string ti_("ti");
        string fi_("fi");
        string ci_("ci");
        string ai_("ai");
        string rdi_("rdi");
        string rli_("rli");
        string rti_("rti");
        string rfi_("rfi");
        string rci_("rci");
        string rai_("rai");
        inline string n2s(ex fn) {
            int n = ex_to<numeric>(fn).to_int();
            return (n<0 ? "m" : "") + to_string(abs(n));
        }
    }

    REGISTER_FUNCTION(Propagator, do_not_evalf_params())
    REGISTER_FUNCTION(InField, do_not_evalf_params())
    REGISTER_FUNCTION(OutField, do_not_evalf_params())
    
    #ifndef DOXYGEN_SKIP
    
    unsigned Field2_SERIAL::serial = GiNaC::function::register_new(function_options("Field",2).do_not_evalf_params().overloaded(2));
    unsigned Field3_SERIAL::serial = GiNaC::function::register_new(function_options("Field",3).do_not_evalf_params().overloaded(2));
    unsigned Vertex2_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",2).do_not_evalf_params().overloaded(5));
    unsigned Vertex3_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",3).do_not_evalf_params().overloaded(5));
    unsigned Vertex4_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",4).do_not_evalf_params().overloaded(5));
    unsigned Vertex5_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",5).do_not_evalf_params().overloaded(5));
    unsigned Vertex6_SERIAL::serial = GiNaC::function::register_new(function_options("Vertex",6).do_not_evalf_params().overloaded(5));
    
    #endif

    Index DI(ex fn) { return Index(di_+n2s(fn),Index::Type::VD); }
    Index LI(ex fn) { return Index(li_+n2s(fn),Index::Type::VD); }
    Index TI(ex fn) { return Index(ti_+n2s(fn),Index::Type::CF); }
    Index FI(ex fn) { return Index(fi_+n2s(fn),Index::Type::CF); }
    Index CI(ex fn) { return Index(ci_+n2s(fn),Index::Type::CA); }
    Index AI(ex fn) { return Index(ai_+n2s(fn),Index::Type::CA); }
    Index RDI(ex fn) { return Index(rdi_+n2s(fn),Index::Type::VD); }
    Index RLI(ex fn) { return Index(rli_+n2s(fn),Index::Type::VD); }
    Index RTI(ex fn) { return Index(rti_+n2s(fn),Index::Type::CF); }
    Index RFI(ex fn) { return Index(rfi_+n2s(fn),Index::Type::CF); }
    Index RCI(ex fn) { return Index(rci_+n2s(fn),Index::Type::CA); }
    Index RAI(ex fn) { return Index(rai_+n2s(fn),Index::Type::CA); }
    
    /**
     * @brief generte the Amplitudes
     * @param st symtab for the parser, no need for Symbol, usually used for momentum Vector
     * @return Amplitudes for current QGRAF object
     */
    lst Process::Amplitudes(symtab st) {
        auto rc = system("rm -f qgraf.dat qgraf.out qgraf.sty qgraf.mod");
        std::ofstream style;
        style.open("qgraf.sty", ios::out);
        style << Style << endl;
        style.close();
        
        std::ofstream model;
        style.open("qgraf.mod", ios::out);
        style << Model << endl;
        style.close();
        
        std::ofstream ofs;
        ofs.open("qgraf.dat", ios::out);
        ofs << "output='qgraf.out';" << endl;
        ofs << "style='qgraf.sty';" << endl;
        ofs << "model='qgraf.mod';" << endl;
        ofs << "in=" << In << ";" << endl;
        ofs << "out=" << Out << ";" << endl;
        ofs << "loops=" << Loops << ";" << endl;
        ofs << "loop_momentum=" << LoopPrefix << ";" << endl;
        ofs << "options=" << Options << ";" << endl;
        for(auto vs : Others) ofs << vs << ";" << endl;
        ofs.close();
        
        if(Debug) rc = system((InstallPrefix+"/bin/qgraf").c_str());
        else rc = system((InstallPrefix+"/bin/qgraf > /dev/null").c_str());
        
        ifstream ifs("qgraf.out");
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        
        if(!Debug) {
            if(access("qgraf.dat",F_OK)!=-1) remove("qgraf.dat");
            if(access("qgraf.mod",F_OK)!=-1) remove("qgraf.mod");
            if(access("qgraf.sty",F_OK)!=-1) remove("qgraf.sty");
            if(access("qgraf.out",F_OK)!=-1) remove("qgraf.out");
        }
        
        const char* rm_chars = " \t\v\r\n,";
        if(!ostr.empty()) {
            ostr.erase(0, ostr.find_first_not_of(rm_chars));
            ostr.erase(ostr.find_last_not_of(rm_chars)+1);
        }
        ostr = "{" + ostr + "}";

        Parser amp(st);
        ex ret = amp.Read(ostr);
        
        return ex_to<lst>(ret);
    }
    
    /**
     * @brief generate the topological lines for the amplitude
     * @param amp the element from Amplitudes
     * @return topoligical lines, {{vid1,fs1}, {vid2,fs2}, lst{f1,f2}/lst{exf},mom}
     */
    lst TopoLines(const ex & amp) {
        map<ex,int,ex_is_less> v2id, fid2vid;
        map<int,ex> vid2fs; // fileds in the vertex
        int cid = 0;
        lst lines;
        for(const_preorder_iterator i = amp.preorder_begin(); i != amp.preorder_end(); ++i) {
            ex e = (*i);
            if(isFunction(e,"OutField")) {
                lines.append(lst{ e.op(1), iWF(e.op(1)), lst{e.op(0)}, e.op(2) });
            } else if(isFunction(e,"InField")) {
                lines.append(lst{ iWF(e.op(1)), e.op(1), lst{e.op(0)}, e.op(2) });
            } else if(isFunction(e, "Propagator")) {
                lines.append(lst{ iWF(e.op(0).op(1)), iWF(e.op(1).op(1)), lst{e.op(0).op(0), e.op(1).op(0)}, e.op(2) });
            } else if(isFunction(e, "Vertex")) {
                auto itr = v2id.find(e);
                if(itr==v2id.end()) {
                    cid++;
                    v2id.emplace(make_pair(e,cid));
                    lst fs;
                    for(auto f : e) {
                        ex fid = f.op(1);
                        fid2vid[fid] = cid;
                        fs.append(f.op(0));
                    }
                    vid2fs[cid] = fs;
                }
            }
        }
        lines = ex_to<lst>(MapFunction([&vid2fs,&fid2vid](const ex &e, MapFunction &self)->ex{
            if(!e.has(iWF(w))) return e;
            else if(e.match(iWF(w))) {
                auto vid = fid2vid[e.op(0)];
                return lst{vid, vid2fs[vid]};
            } else return e.map(self);
        })(lines));
        
        return lines.sort();
    }
    
    /**
     * @brief generate Feynman diagrams for the amplitudes, in PDF format
     * @param amps refers to Amplitudes
     * @param fn the filename of the PDF
     * @return nonthing, check pdf file
     */
    void DrawPDF(const lst & amps, string fn) {
        int id=0, rc;
        exvector amp_vec;
        for(auto item : amps) amp_vec.push_back(item);
        string tex_path = to_string(getpid()) + "_TeX/";
        if(!dir_exists(tex_path)) rc = system(("mkdir -p "+tex_path).c_str());
        int limit = 300;
        
        GiNaC_Parallel(amp_vec.size(), [&amp_vec,tex_path](int idx)->ex {
            auto amp = amp_vec[idx];
            ofstream out(tex_path+to_string(idx)+".tex");
            out << "\\documentclass[tikz]{standalone}" << endl;
            out << "\\usepackage{tikz-feynman}" << endl;
            out << "\\tikzfeynmanset{compat=1.1.0}" << endl;
            out << "\\begin{document}" << endl;
            out << "\\feynmandiagram{" << endl;
            auto lines = TopoLines(amp);

            exmap bend_map;
            std::map<ex,int,ex_is_less> vtex_map; // vertex option, only once
            for(auto l : lines) {
                lst ll = lst{l.op(0), l.op(1)};
                bool isExt = (is_a<numeric>(ll.op(0)) && ll.op(0)<0) || (is_a<numeric>(ll.op(1)) && ll.op(1)<0);
                ll.sort();
                bend_map[ll] = bend_map[ll] + 1;
                
                auto fidL = (is_a<numeric>(l.op(1)) ? l.op(1) : l.op(1).op(0));
                out << "\"" << fidL << "\"";
                if(fidL<0) {
                    out << "[particle=";
                    if(InOutTeX[fidL].length()>0) out << InOutTeX[fidL];
                    else out << fidL;
                    out << "]";
                } else if(vtex_map[l.op(1).op(1)]==0) {
                    out << VerTeX[l.op(1).op(1)];
                    vtex_map[l.op(1).op(1)]=1;
                }

                out << "  --[";
                auto f = l.op(2).op(0);
                if(LineTeX[f].length()>0) {
                    if(!isExt) out << LineTeX[f];
                    else {
                        auto cpos = LineTeX[f].find(", edge");
                        if(cpos>0) out << LineTeX[f].substr(0,cpos);
                        else out << LineTeX[f];
                    }
                }
                if(bend_map[ll]>2) out << ",half right";
                else if(bend_map[ll]>1) out << ",half left";
                if(is_zero(l.op(0)-l.op(1))) out << ",loop,distance=2cm";
                out << "]";
                
                auto fidR = (is_a<numeric>(l.op(0)) ? l.op(0) : l.op(0).op(0));
                out << "  \"" << fidR << "\"";
                if(fidR<0) {
                    out << "[particle=";
                    if(InOutTeX[fidR].length()>0) out << InOutTeX[fidR];
                    else out<< fidR;
                    out << "]";
                } else if(vtex_map[l.op(0).op(1)]==0) {
                    out << VerTeX[l.op(0).op(1)];
                    vtex_map[l.op(0).op(1)]=1;
                }
                out << ";" << endl;
            }
            out << "};" << endl;
            out << "\\end{document}" << endl;
            out.close();
            auto rc = system(("cd "+tex_path+" && echo X | lualatex " + to_string(idx) + " 1>/dev/null").c_str());
            return 0;
        }, "TeX");
        
        ofstream out(tex_path+"diagram.tex");
        out << "\\let\\mypdfximage\\pdfximage" << endl;
        out << "\\def\\pdfximage{\\immediate\\mypdfximage}" << endl;
        out << "\\documentclass{standalone}" << endl;
        out << "\\usepackage{graphicx}" << endl;
        out << "\\usepackage{adjustbox}" << endl;
        out << "\\begin{document}" << endl;
        out << "\\begin{adjustbox}{valign=T,width=\\textwidth}" << endl;
        out << "\\begin{tabular}{|cc|cc|cc|cc|}" << endl;
        out << "\\hline" << endl;
        int total = amps.nops();
        int namps = total;
        if((total%4)!=0) total = (total/4+1)*4;
        for(int i=0 ; i<total; i++) {
            
            if((i!=0) && (i+1!=total) && (i%limit)==0) {
                out << "\\end{tabular}" << endl << endl;
                out << "\\begin{tabular}{|cc|cc|cc|cc|}" << endl;
                out << "\\hline" << endl;
            }
            
            out << "{\\tiny " << i+1 << "}&" << endl;
            if(i<namps) {
                out << "\\includegraphics[keepaspectratio,";
                out << "height=0.22\\textwidth,";
                out << "width=0.22\\textwidth]";
                out << "{"<<i<<".pdf}" << endl;
            }
            if((i+1)%4==0) out << "\\\\ \\hline";
            else out << "&";
        }
        out << "\\end{tabular}" << endl;
        out << "\\end{adjustbox}" << endl;
        out << "\\end{document}" << endl;
        out.close();
        if(Debug) rc = system(("cd "+tex_path+" && pdflatex diagram && mv diagram.pdf ../"+fn).c_str());
        else rc = system(("cd "+tex_path+" && echo X | pdflatex diagram 1>/dev/null && mv diagram.pdf ../"+fn).c_str());
        if(!Debug) rc = system(("rm -r "+tex_path).c_str());
    }
    
    /**
     * @brief cut the amplitude, and return the connected parts, may have many different cuts
     * @param amp one of the Amplitudes
     * @param prop the line type to be cutted, can be lst{g,g} or lst{ lst{g,g}, lst{A,A} }
     * @param n the number of lines to be cutted
     * @return vector of lst, each element in the vector, is actually a lst of lst, different connectted parts
     */
    vector<lst> ShrinkCut(ex amp, lst prop, int n) {
        vector<lst> ret;
        if(prop.nops()<1) throw Error("ShrinkCut: no cut provided!");
        auto tls = TopoLines(amp);
        vector<int> cls_vec;
        for(int i=0; i<tls.nops(); i++) {
            auto pi = tls.op(i).op(2);
            if(pi.nops()<2) continue;
            if(!is_a<lst>(prop.op(0))) { // prop : lst{ g, g }
                if(is_zero(pi.op(0)-prop.op(0)) && is_zero(pi.op(1)-prop.op(1))) cls_vec.push_back(i);
            } else {
                for(auto iprop : prop) {
                    if(is_zero(pi.op(0)-iprop.op(0)) && is_zero(pi.op(1)-iprop.op(1))) cls_vec.push_back(i);
                }
            }
        }
        if(cls_vec.size()<n) return ret;
        
        Combinations(cls_vec.size(), n, [n,&ret,cls_vec,tls](const int * is)->void {
            int cls[n];
            for(int i=0; i<n; i++) cls[i] = cls_vec[is[i]];
            
            // cut each line into 2 half-lines, labeled with 0
            auto tls2 = tls;
            for(auto ci : cls) {
                auto ol = tls2.op(ci);
                tls2.let_op(ci) = lst{ol.op(0), 0, lst{ol.op(2).op(0)}, ol.op(3)};
                tls2.append(lst{0, ol.op(1), lst{ol.op(2).op(1)}, ol.op(3)} );
            }
            
            // shrink internal lines
            int last = 0;
            int ntls2 = tls2.nops();
            while(true) {
                ex lp = 0;
                for(int i=last; i<ntls2; i++) {
                    auto li = tls2.op(i);
                    if(is_zero(li) || li.op(2).nops()<2) continue;
                    last = i;
                    lp = li;
                    tls2.let_op(last) = 0;
                    break;
                }
                if(is_zero(lp)) break;
                
                for(int i=0; i<ntls2; i++) {
                    if(is_zero(tls2.op(i))) continue;
                    if(is_zero(tls2.op(i).op(0)-lp.op(0))) tls2.let_op(i).let_op(0) = lp.op(1);
                    if(is_zero(tls2.op(i).op(1)-lp.op(0))) tls2.let_op(i).let_op(1) = lp.op(1);
                }
            }

            // final connected parts
            map<int, lst> con_map;
            for(auto li : tls2) {
                if(is_zero(li)) continue;
                ex key, val;
                ex fiL = li.op(0);
                if(is_a<lst>(fiL)) fiL = fiL.op(0);
                ex fiR = li.op(1);
                if(is_a<lst>(fiR)) fiR = fiR.op(0);
                if(fiL>0 && fiR<0) {
                    con_map[ex_to<numeric>(fiL).to_int()].append(fiR);
                } else if(fiR>0 && fiL<0) {
                    con_map[ex_to<numeric>(fiR).to_int()].append(fiL);
                } else {
                    if(fiL>0 && is_zero(fiR)) key = fiL;
                    else if(fiR>0 && is_zero(fiL)) key = fiR;
                    else throw Error("ShrinkCut: unexpcected point reached.");
                    val = li.op(2).op(0);
                    con_map[ex_to<numeric>(key).to_int()].append(val);
                }
            }
            
            lst item;
            for(auto kv : con_map) item.append(kv.second.sort());
            ret.push_back(item);
        });
        
        return ret;
    }
    
    /**
     * @brief check a amplitude has a loop w.r.t. propapagtor
     * @param amp one of the Amplitudes
     * @param prop the line type in the loop
     * @return true if the corresponding loop exsits
     */
    bool HasLoop(ex amp, lst prop) {
        auto tls = TopoLines(amp);
        int ntls = tls.nops();
        // shrink internal lines of prop
        int last = 0;
        while(true) {
            ex lp = 0;
            for(int i=last; i<ntls; i++) {
                auto li = tls.op(i);
                if(is_zero(li) || li.op(2).nops()<2) continue; // external line
                auto cpi = li.op(2);
                bool m1 = !(is_zero(cpi.op(0)-prop.op(0)) && is_zero(cpi.op(1)-prop.op(1)));
                bool m2 = !(is_zero(cpi.op(0)-prop.op(1)) && is_zero(cpi.op(1)-prop.op(0)));
                if(m1 && m2) continue; // other propagators
                if(is_zero(li.op(0)-li.op(1))) return true; // a loop found
                last = i;
                lp = li;
                tls.let_op(last) = 0;
                break;
            }
            if(is_zero(lp)) break;
            
            for(int i=0; i<ntls; i++) {
                if(is_zero(tls.op(i))) continue;
                if(is_zero(tls.op(i).op(0)-lp.op(0))) tls.let_op(i).let_op(0) = lp.op(1);
                if(is_zero(tls.op(i).op(1)-lp.op(0))) tls.let_op(i).let_op(1) = lp.op(1);
            }
        }
        return false;
    }
    
    /**
     * @brief Propagator for quark
     * @param e expression with head of Propagator
     * @param m the quark mass
     * @return the quark propagator
     */
    ex QuarkPropagator(const ex & e, const ex & m) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto mom = e.op(2);
        return I * SP(TI(fi1),TI(fi2)) * Matrix(GAS(mom)+GAS(1)*m, DI(fi1),DI(fi2)) / (SP(mom)-m*m);
    }
    
    /**
     * @brief Propagator for lepton
     * @param e expression with head of Propagator
     * @param m the quark mass
     * @return the lepton propagator
     */
    ex LeptonPropagator(const ex & e, const ex & m) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto mom = e.op(2);
        return I * Matrix(GAS(mom)+GAS(1)*m, DI(fi1),DI(fi2)) / (SP(mom)-m*m);
    }
    
    /**
     * @brief Propagator for gluon
     * @param e expression with head of Propagator
     * @param xi R-xi gauge parameter
     * @return the gloun propagator for R-xi gauge
     */
    ex GluonPropagator(const ex & e, const ex & xi) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto mom = e.op(2);
        return (-I) * SP(CI(fi1),CI(fi2)) * ( SP(LI(fi1),LI(fi2)) / SP(mom) - (1-xi)*SP(mom,LI(fi1))*SP(mom,LI(fi2))/pow(SP(mom),2) );
    }
        
    /**
     * @brief Propagator for gluon ghost
     * @param e expression with head of Propagator
     * @param xi the result wil not depend on xi, since null gluon mass
     * @return the gluon ghost propagator
     */
    ex GluonGhostPropagator(const ex & e, const ex & xi, const ex & eta_G) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto mom = e.op(2);
        return I * eta_G * SP(CI(fi1),CI(fi2)) / SP(mom);
    }
    
    /**
     * @brief Propagator for A/Z/W boson
     * @param e expression with head of Propagator
     * @param m the boson mass, 0/MZ/MW
     * @param xi R-xi gauge parameter
     * @return the A/Z/W boson propagator
     */
    ex AZWPropagator(const ex & e, const ex & m, const ex & xi) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto mom = e.op(2);
        return -I/(SP(mom)-m*m) * (SP(LI(fi1),LI(fi2))-(1-xi)*SP(LI(fi1))*SP(LI(fi2))/(SP(mom)-xi*m*m));
    }
    
    /**
     * @brief Propagator for A/Z/W ghost
     * @param e expression with head of Propagator
     * @param m mass for A/Z/W -> 0/MZ/MW
     * @param xi R-xi gauge parameter
     * @return the A/Z/W ghost propagator @ R-xi gauge
     */
    ex AZWGhostPropagator(const ex & e, const ex & m, const ex & xi, const ex & eta_G) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto mom = e.op(2);
        return eta_G * I / (SP(mom)-xi*m*m);
    }
    
    /**
     * @brief Propagator for scalar
     * @param e expression with head of Propagator
     * @param m mass for the particle
     * @param xi R-xi gauge parameter, xi=1 for higgs
     * @return the scalar propagator
     */
    ex ScalarPropagator(const ex & e, const ex & m, const ex & xi) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto mom = e.op(2);
        return I / (SP(mom)-xi*m*m);
    }
    
    /**
     * @brief q-qbar-g vertex
     * @param e expression with head of Vertex
     * @param eta_s +1/-1 the sign convention
     * @return the q-qbar-g vertex
     */
    ex QuarkGluonVertex(const ex & e, const ex & eta_s) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto fi3 = e.op(2).op(1);
        return -I*eta_s*gs*Matrix(GAS(LI(fi3)),DI(fi1),DI(fi2))*SUNT(CI(fi3),TI(fi1),TI(fi2));
    }

    /**
     * @brief g-g-g vertex
     * @param e expression with head of Vertex
     * @param eta_s +1/-1 the sign convention
     * @return the g-g-g vertex
     */
    ex Gluon3Vertex(const ex & e, const ex & eta_s) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto fi3 = e.op(2).op(1);
        auto mom1 = e.op(0).op(2);
        auto mom2 = e.op(1).op(2);
        auto mom3 = e.op(2).op(2);
        return -eta_s*gs*SUNF(CI(fi1),CI(fi2),CI(fi3))*(
            SP(mom1-mom2,LI(fi3))*SP(LI(fi1),LI(fi2)) +
            SP(mom2-mom3,LI(fi1))*SP(LI(fi2),LI(fi3)) +
            SP(mom3-mom1,LI(fi2))*SP(LI(fi3),LI(fi1))
        );
    }
    
    /**
     * @brief g-g-g-g vertex
     * @param e expression with head of Vertex
     * @param eta_s +1/-1 the sign convention
     * @return the g-g-g-g vertex
     */
    ex Gluon4Vertex(const ex & e, const ex & eta_s) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto fi3 = e.op(2).op(1);
        auto fi4 = e.op(3).op(1);
        return -I*gs*gs*(
            SUNF4(CI(fi1),CI(fi2),CI(fi3),CI(fi4))*(SP(LI(fi1),LI(fi3))*SP(LI(fi2),LI(fi4))-SP(LI(fi1),LI(fi4))*SP(LI(fi2),LI(fi3))) +
            SUNF4(CI(fi1),CI(fi3),CI(fi2),CI(fi4))*(SP(LI(fi1),LI(fi2))*SP(LI(fi3),LI(fi4))-SP(LI(fi1),LI(fi4))*SP(LI(fi2),LI(fi3))) +
            SUNF4(CI(fi1),CI(fi4),CI(fi2),CI(fi3))*(SP(LI(fi1),LI(fi2))*SP(LI(fi4),LI(fi3))-SP(LI(fi1),LI(fi3))*SP(LI(fi4),LI(fi2)))
        );
    }
    
    /**
     * @brief ghbar-gh-g vertex
     * @param e expression with head of Vertex
     * @param eta_s +1/-1 the sign convention
     * @param eta_G +1/-1 the sign convention
     * @return the ghbar-gh-g vertex
     */
    ex GhostGluonVertex(const ex & e, const ex & eta_s, const ex & eta_G) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto fi3 = e.op(2).op(1);
        auto mom1 = e.op(0).op(2);
        return eta_s*eta_G*gs*SUNF(CI(fi1),CI(fi2),CI(fi3))*SP(mom1,LI(fi3));
    }

    /**
     * @brief qbar-q-A vertex
     * @param e expression with head of Vertex
     * @param eq the charge
     * @param eta_e +1/-1 the sign convention
     * @return the q-qbar-A vertex
     */
    ex QuarkAVertex(const ex & e, const ex & eq, const ex & eta_e) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto fi3 = e.op(2).op(1);
        return -I*eta_e*eq*Matrix(GAS(LI(fi3)),DI(fi1),DI(fi2))*SP(TI(fi1),TI(fi2));
    }
    
    /**
     * @brief l-lbar-A vertex
     * @param e expression with head of Vertex
     * @param eq the charge
     * @param eta_e +1/-1 the sign convention
     * @return the q-qbar-A vertex
     */
    ex LeptonAVertex(const ex & e, const ex & eq, const ex & eta_e) {
        auto fi1 = e.op(0).op(1);
        auto fi2 = e.op(1).op(1);
        auto fi3 = e.op(2).op(1);
        return -I*eta_e*eq*Matrix(GAS(LI(fi3)),DI(fi1),DI(fi2));
    }

    /**
     * @brief Change Index from left to right, only affect li/di/ci/ti, external index start with dim/lim/cim/tim will also be changed if all=true
     * @param e input expression 
     * @param all if true external index will also be changed
     * @return Index changed
     */
     ex IndexL2R(ex e, bool all) {
        static MapFunction map([all](const ex &e, MapFunction &self)->ex {
            if(!Index::has(e)) return e;
            else if(is_a<Index>(e)) {
                auto idx = ex_to<Index>(e);
                auto nstr = idx.name.get_name();
                if(!all && nstr.rfind("lim",0)==0) return e;
                else if(!all && nstr.rfind("dim",0)==0) return e;
                else if(!all && nstr.rfind("cim",0)==0) return e;
                else if(!all && nstr.rfind("tim",0)==0) return e;
                else if(nstr.rfind("li",0)==0 || nstr.rfind("di",0)==0 || nstr.rfind("ci",0)==0 || nstr.rfind("ti",0)==0) return Index("r"+nstr, idx.type);
                else return e;
            }
            else return e.map(self);
        });
        return map(e);
     }
     
     ex IndexCC(ex e, bool all) {
        return IndexL2R(e,all);
     }
     
    /**
     * @brief polarization sum for gluon 
     * @param qi gluon qgraf index
     * @param color true for QCD, false for QED
     * @return g^{i ir} delta^{i ir}
     */
    ex GluonSumL(int qi, bool color) {
        if(color) return -SP(LI(qi), RLI(qi)) * SP(CI(qi), RCI(qi));
        else return -SP(LI(qi), RLI(qi));
    }
    
    /**
     * @brief polarization sum for gluon
     * @param qi gluon qgraf index
     * @param color true for QCD, false for QED
     * @return g^{i ir} delta^{i ir}
     */
    ex GluonSumR(int qi, bool color) {
        if(color) return -SP(RLI(qi),LI(qi)) * SP(RCI(qi),CI(qi));
        else return -SP(RLI(qi),LI(qi));
    }

    /**
     * @brief polarization sum for quark 
     * @param qi quark qgraf index
     * @param p anti-quark momentum vector
     * @param m anti-quark mass
     * @param color true for QCD, false for QED
     * @return Quark summation
     */
    ex QuarkSumL(int qi, ex p, ex m, bool color) {
        if(color) return Matrix(GAS(p)+m*GAS(1), RDI(qi), DI(qi)) * SP(RTI(qi), TI(qi));
        else return Matrix(GAS(p)+m*GAS(1), RDI(qi), DI(qi));
    }
    
    /**
     * @brief polarization sum for quark
     * @param qi quark qgraf index
     * @param p anti-quark momentum vector
     * @param m anti-quark mass
     * @param color true for QCD, false for QED
     * @return Quark summation
     */
    ex QuarkSumR(int qi, ex p, ex m, bool color) {
        if(color) return Matrix(GAS(p)+m*GAS(1), DI(qi), RDI(qi)) * SP(TI(qi),RTI(qi));
        else return Matrix(GAS(p)+m*GAS(1), DI(qi), RDI(qi));
    }

    /**
     * @brief polarization sum for anti-quark 
     * @param qi anti-quark qgraf index
     * @param p anti-quark momentum vector
     * @param m anti-quark mass
     * @param color true for QCD, false for QED
     * @return anti-Quark summation
     */
    ex AntiQuarkSumL(int qi, ex p, ex m, bool color) {
        if(color) return Matrix(GAS(p)-m*GAS(1), DI(qi), RDI(qi)) * SP(TI(qi), RTI(qi));
        else return Matrix(GAS(p)-m*GAS(1), DI(qi), RDI(qi));
    }
    
    /**
     * @brief polarization sum for anti-quark
     * @param qi anti-quark qgraf index
     * @param p anti-quark momentum vector
     * @param m anti-quark mass
     * @param color true for QCD, false for QED
     * @return anti-Quark summation
     */
    ex AntiQuarkSumR(int qi, ex p, ex m, bool color) {
        if(color) return Matrix(GAS(p)-m*GAS(1), RDI(qi), DI(qi)) * SP(RTI(qi), TI(qi));
        else return Matrix(GAS(p)-m*GAS(1), RDI(qi), DI(qi));
    }
    
    /**
     * @brief polarization sum for ghost 
     * @param qi ghost qgraf index
     * @return g^{i ir} delta^{i ir}
     */
    ex GhostSumL(int qi) {
        return SP(CI(qi), RCI(qi));
    }
    
    /**
     * @brief polarization sum for ghost
     * @param qi ghost qgraf index
     * @return g^{i ir} delta^{i ir}
     */
    ex GhostSumR(int qi) {
        return SP(RCI(qi),CI(qi));
    }
    
    /**
     * @brief polarization sum for anti-ghost, note that we will add extra - here
     * @param qi anti-ghost qgraf index
     * @return g^{i ir} delta^{i ir}
     */
    ex AntiGhostSumL(int qi) {
        return -SP(CI(qi), RCI(qi));
    }
    
    /**
     * @brief polarization sum for anti-ghost, note that we will add extra - here
     * @param qi anti-ghost qgraf index
     * @return g^{i ir} delta^{i ir}
     */
    ex AntiGhostSumR(int qi) {
        return -SP(RCI(qi),CI(qi));
    }
    
    /**
     * @brief polarization sum for total angular momentum
     * @param qi qgraf index
     * @param p the total momentum
     * @return -g^{qi, rqi} + p^qi p^rqi/p.p
     */
    ex J1SumL(int qi, ex p) {
        if(is_zero(p)) return -SP(LI(qi),RLI(qi));
        else return -SP(LI(qi),RLI(qi)) + SP(p,LI(qi)) * SP(p,RLI(qi)) / SP(p);
    }
    
    /**
     * @brief polarization sum for total angular momentum
     * @param qi qgraf index
     * @param p the momentum
     * @return -g^{qi, rqi} + p^qi p^rqi/p.p
     */
    ex J1SumR(int qi, ex p) {
        if(is_zero(p)) return -SP(RLI(qi),LI(qi));
        else return -SP(RLI(qi),LI(qi)) + SP(p,RLI(qi)) * SP(p,LI(qi)) / SP(p);
    }
    
    // https://arxiv.org/abs/1209.6213v2
    lst Models::FeynRulesSM(const lst & amps, const ex & xi) {
        lst ret;
        for(auto item : amps) ret.append(FeynRulesSM(item,xi));
        return ret;
    }
    ex Models::FeynRulesSM(const ex & amp, const ex & xi) {
        if(is_a<lst>(amp)) return FeynRulesSM(ex_to<lst>(amp), xi);
        
        static Symbol CW("CW"); // cos(theta)
        static Symbol SW("SW"); // sin(theta)
        //static Symbol C2W("C2W"); // cos(2theta)
        static ex CW2 = CW*CW;
        static ex SW2 = SW*SW;
        static ex C2W = CW2-SW2;
        
        static Symbol U("U"); // U-quark
        static Symbol Ubar("Ubar"); // anti U-quark
        static Symbol D("D"); // D-quark
        static Symbol Dbar("Dbar"); // anti D-quark
        static Symbol C("C"); // C-quark
        static Symbol Cbar("Cbar"); // anti C-quark
        static Symbol S("S"); // S-quark
        static Symbol Sbar("Sbar"); // anti S-quark
        static Symbol T("T"); // T-quark
        static Symbol Tbar("Tbar"); // anti T-quark
        static Symbol B("B"); // B-quark
        static Symbol Bbar("Bbar"); // anti B-quark
        
        static Symbol g("g"); // gluon
        static Symbol gh("gh"); // gluon ghost
        static Symbol ghbar("ghbar"); // anti gluon ghost
        
        static Symbol A("A"); // photon
        static Symbol Wm("Wm"); // W-
        static Symbol Wp("Wp"); // W+
        static Symbol Z("Z"); // Z
        
        static Symbol ghA("ghA"); // photon ghost
        static Symbol ghAbar("ghAbar"); // anti photon ghost
        static Symbol ghWm("ghWm"); // W- ghost
        static Symbol ghWmbar("ghWmbar"); // anti W- ghost
        static Symbol ghWp("ghWp"); // W+ ghost
        static Symbol ghWpbar("ghWpbar"); // anti W+ ghost
        static Symbol ghZ("ghZ"); // Z ghost
        static Symbol ghZbar("ghZbar"); // anti Z ghost

        static Symbol em("em"); // e-
        static Symbol ep("ep"); // e+
        static Symbol ne("ne"); // e-neutrino
        static Symbol nebar("nebar"); // anti e-neutrino
        static Symbol mum("mum"); // mu-
        static Symbol mup("mup"); // mu+
        static Symbol nmu("nmu"); // mu-neutrino
        static Symbol nmubar("nmubar"); // anti mu-neutrino
        static Symbol taum("taum"); // tau-
        static Symbol taup("taup"); // tau+
        static Symbol ntau("ntau"); // tau-neutrino
        static Symbol ntaubar("ntaubar"); // anti tau-neutrino
        
        static Symbol chi("chi"); // Z goldstone
        static Symbol phim("phim"); // W- goldstone
        static Symbol phip("phip"); // W+ goldstone
        static Symbol H("H"); // higgs
        
        static auto M = [](const string & si)->ex {
            return Symbol("M"+si);
        };
        static ex MW = M("W");
        static ex MZ = M("Z");
        static ex MH = M("H");
        
        static Symbol EL("e"); // charge e
        
        static ex eta_s = -1;
        // Peskin's book convention
        static ex eta = -1;
        static ex eta_prime = -1;
        static ex eta_Z = 1;
        static ex eta_theta = 1;
        static ex eta_Y = 1;
        static ex eta_e = -1;
        static ex eta_G = 1;
        
        static ex GEW = eta_e*eta*eta_theta * EL/SW;
        static ex EL2 = EL*EL;
        static ex MW2 = MW*MW;
        static ex MZ2 = MZ*MZ;
        static ex MH2 = MH*MH;
        static ex GEW2 = GEW*GEW;
        static ex sqrt2 = sqrt(ex(2));
        
        map<string,ex> T3;
        T3["U"] = T3["C"] = T3["T"] = ex(1)/2;
        T3["D"] = T3["S"] = T3["B"] = -ex(1)/2;
        T3["ne"] = T3["nmu"] = T3["ntau"] = ex(1)/2;
        T3["em"] = T3["mum"] = T3["taum"] = -ex(1)/2;
        T3["ep"] = T3["mup"] = T3["taup"] = -ex(1)/2; // for sure
        
        map<string,ex> Q;
        Q["U"] = Q["C"] = Q["T"] = ex(2)/3;
        Q["D"] = Q["S"] = Q["B"] = -ex(1)/3;
        Q["ne"] = Q["nmu"] = Q["ntau"] = 0;
        Q["em"] = Q["mum"] = Q["taum"] = -1;
        Q["ep"] = Q["mup"] = Q["taup"] = -1; // for sure
        
        auto is_qbar_q = [&](const ex pi1, const ex pi2)->bool {
             return (pi1==Ubar && pi2==U) || (pi1==Dbar && pi2==D) || (pi1==Cbar && pi2==C) || (pi1==Sbar && pi2==S) || (pi1==Tbar && pi2==T) || (pi1==Bbar && pi2==B);
        };
        auto is_lbar_l = [&](const ex pi1, const ex pi2)->bool {
             return (pi1==ep && pi2==em) || (pi1==mup && pi2==mum) || (pi1==taup && pi2==taum);
        };
        auto is_lbar_n = [&](const ex pi1, const ex pi2)->bool {
             return (pi1==ep && pi2==ne) || (pi1==mup && pi2==nmu) || (pi1==taup && pi2==ntau);
        };
        auto is_nbar_n = [&](const ex pi1, const ex pi2)->bool {
            return (pi1==nebar && pi2==ne) || (pi1==nmubar && pi2==nmu) || (pi1==ntaubar && pi2==ntau);
        };
        auto is_nbar_l = [&](const ex pi1, const ex pi2)->bool {
            return (pi1==nebar && pi2==em) || (pi1==nmubar && pi2==mum) || (pi1==ntaubar && pi2==taum);
        };
        auto is_uqbar = [&](const ex pi)->bool {
            return (pi==Ubar) || (pi==Cbar) || (pi==Tbar);
        };
        auto is_uq = [&](const ex pi)->bool {
            return (pi==U) || (pi==C) || (pi==T);
        };
        auto is_dqbar = [&](const ex pi)->bool {
            return (pi==Dbar) || (pi==Sbar) || (pi==Bbar);
        };
        auto is_dq = [&](const ex pi)->bool {
            return (pi==D) || (pi==S) || (pi==B);
        };
        auto is_uqbar_dq = [&](const ex pi1, const ex pi2)->bool {
            return is_uqbar(pi1) && is_dq(pi2);
        };
        auto is_dqbar_uq = [&](const ex pi1, const ex pi2)->bool {
            return is_dqbar(pi1) && is_uq(pi2);
        };
        auto CKM = [&](const string & si1, const string & si2)->ex {
            return Symbol("V"+si1+si2);
        };
        
        auto gfV = [&](const string & si)->ex {
            return T3[si]/2-Q[si]*SW2;
        };
        auto gfA = [&](const string & si)->ex {
            return T3[si]/2;
        };
        
        auto fr = MapFunction([&](const ex &e, MapFunction &self)->ex {
            if(isFunction(e,"OutField") || isFunction(e,"InField")) return 1;
            else if(isFunction(e, "Propagator")) {
                auto pi = e.op(0).op(0);
                auto si = ex2str(pi);
                if(pi==U || pi==D || pi==C || pi==S || pi==T || pi==B) {
                    return QuarkPropagator(e, M(si));
                } else if(pi==g) {
                    return GluonPropagator(e, xi);
                } else if(pi==gh) {
                    return GluonGhostPropagator(e, xi, eta_G);
                } else if(pi==A) {
                    return AZWPropagator(e, 0, xi);
                } else if(pi==Z) {
                    return AZWPropagator(e, MZ, xi);
                } else if(pi==Wm) {
                    return AZWPropagator(e, MW, xi);
                } else if(pi==ghA) {
                    return AZWGhostPropagator(e, 0, xi, eta_G);
                } else if(pi==ghZ) {
                    return AZWGhostPropagator(e, MZ, xi, eta_G);
                } else if(pi==ghWm || pi==ghWp) {
                    return AZWGhostPropagator(e, MW, xi, eta_G);
                } else if(pi==em || pi==mum || pi==taum) {
                    return LeptonPropagator(e, M(si));
                } else if(pi==ne || pi==nmu || pi==ntau) {
                    return LeptonPropagator(e, 0);
                } else if(pi==H) {
                    return ScalarPropagator(e, MH, 1);
                } else if(pi==chi) {
                    return ScalarPropagator(e, MZ, xi);
                } else if(pi==phim) {
                    return ScalarPropagator(e, MW, xi);
                } else {
                    cout << endl << e << endl;
                    throw Error("Propagator Not defined!");
                }
            } else if(isFunction(e, "Vertex")) {
                auto pi1 = e.op(0).op(0);
                auto pi2 = e.op(1).op(0);
                auto pi3 = e.op(2).op(0);
                auto fi1 = e.op(0).op(1);
                auto fi2 = e.op(1).op(1);
                auto fi3 = e.op(2).op(1);
                auto mom1 = e.op(0).op(2);
                auto mom2 = e.op(1).op(2);
                auto mom3 = e.op(2).op(2);
                int nn = e.nops();
                
                auto si1 = ex2str(pi1);
                auto si2 = ex2str(pi2);
                string_replace_all(si1, "bar", "");
                string_replace_all(si2, "bar", "");
                auto si = si2;

                if(nn==3) {
                    if(pi1==ghbar && pi2==gh && pi3==g) {
                        // ghbar-gh-g
                        return GhostGluonVertex(e, eta_s, eta_G);
                    } else if(pi1==g && pi2==g && pi3==g) {
                        // g^3
                        return Gluon3Vertex(e, eta_s);
                    } else if((pi3==g) && is_qbar_q(pi1, pi2)) {
                        // qbar-q-g
                        return QuarkGluonVertex(e, eta_s);
                    } else if((pi3==A) && is_qbar_q(pi1, pi2)) {
                        // qbar-q-A
                        return QuarkAVertex(e, EL*Q[si], eta_e);
                    } else if((pi3==A) && is_lbar_l(pi1, pi2)) {
                        // lbar-l-A
                        return LeptonAVertex(e, EL*Q[si], eta_e);
                    } else if((pi3==Z) && is_qbar_q(pi1, pi2)) {
                        // qbar-q-Z
                        return -I*eta*eta_Z*GEW/CW*(gfV(si)*Matrix(GAS(LI(fi3)),DI(fi1),DI(fi2))-gfA(si)*Matrix(GAS(LI(fi3))*GAS(5),DI(fi1),DI(fi2)))*SP(TI(fi1),TI(fi2));
                    } else if((pi3==Z) && (is_lbar_l(pi1, pi2) || is_nbar_n(pi1, pi2))) {
                        // lbar-l-Z & nbar-n-Z
                        return -I*eta*eta_Z*GEW/CW*(gfV(si)*Matrix(GAS(LI(fi3)),DI(fi1),DI(fi2))-gfA(si)*Matrix(GAS(LI(fi3))*GAS(5),DI(fi1),DI(fi2)));
                    } else if(pi1==Wp && pi2==Wm && pi3==A) {
                        // Wp-Wm-A
                        return -I*eta_e*EL*(SP(LI(fi1),LI(fi2))*SP(mom2-mom1,LI(fi3))+SP(LI(fi2),LI(fi3))*SP(mom3-mom2,LI(fi1))+SP(LI(fi3),LI(fi1))*SP(mom1-mom3,LI(fi2)));
                    } else if(pi1==Wp && pi2==Wm && pi3==Z) {
                        // Wp-Wm-Z
                        return -I*eta*eta_Z*GEW*CW*(SP(LI(fi1),LI(fi2))*SP(mom2-mom1,LI(fi3))+SP(LI(fi2),LI(fi3))*SP(mom3-mom2,LI(fi1))+SP(LI(fi3),LI(fi1))*SP(mom1-mom3,LI(fi2)));
                    } else if(is_uqbar(pi1) && is_dq(pi2) && pi3==Wp) {
                        // uqbar-dq-Wp
                        auto ckm = CKM(si1, si2);
                        return -I*eta*GEW/sqrt2*Matrix(GAS(LI(fi3))-GAS(LI(fi3))*GAS(5),DI(fi1),DI(fi2))/2*ckm*SP(TI(fi1),TI(fi2));
                    } else if(is_dqbar(pi1) && is_uq(pi2) && pi3==Wm) {
                        // dqbar-uq-Wm
                        auto ckm = CKM(si1, si2);
                        return -I*eta*GEW/sqrt2*Matrix(GAS(LI(fi3))-GAS(LI(fi3))*GAS(5),DI(fi1),DI(fi2))/2*ckm*SP(TI(fi1),TI(fi2));
                    } else if( (is_nbar_l(pi1, pi2) && pi3==Wp) || (is_lbar_n(pi1, pi2) && pi3==Wm) ) {
                        // nbar-l-Wp & lbar-n-Wm
                        return -I*eta*GEW/sqrt2*Matrix(GAS(LI(fi3))-GAS(LI(fi3))*GAS(5),DI(fi1),DI(fi2))/2;
                    } else if(is_qbar_q(pi1, pi2) && (pi3==H)) {
                        // qbar-q-H
                        return -I*GEW/2*M(si)/MW*Matrix(GAS(1),DI(fi1),DI(fi2))*SP(TI(fi1),TI(fi2));
                    } else if( (is_lbar_l(pi1, pi2) || is_nbar_n(pi1, pi2)) && (pi3==H)) {
                        // lbar-l-H & nbar-n-H
                        return -I*GEW/2*M(si)/MW*Matrix(GAS(1),DI(fi1),DI(fi2));
                    } else if(is_qbar_q(pi1, pi2) && (pi3==chi)) {
                        // qbar-q-chi
                        return -GEW*T3[si]*M(si)/MW*Matrix(GAS(5),DI(fi1),DI(fi2))*SP(TI(fi1),TI(fi2));
                    } else if( (is_lbar_l(pi1, pi2) || is_nbar_n(pi1, pi2)) && (pi3==chi)) {
                        // lbar-l-chi & nbar-n-chi
                        return -GEW*T3[si]*M(si)/MW*Matrix(GAS(5),DI(fi1),DI(fi2));
                    } else if(is_uqbar_dq(pi1, pi2) && pi3==phip) {
                        // uqbar-dq-phip
                        auto ckm = CKM(si1,si2);
                        return I*GEW/sqrt2*(M(si1)/MW*Matrix(GAS(1)-GAS(5),DI(fi1),DI(fi2))/2-M(si2)/MW*Matrix(GAS(1)+GAS(5),DI(fi1),DI(fi2))/2)*ckm*SP(TI(fi1),TI(fi2));
                    } else if(is_dqbar_uq(pi1, pi2) && pi3==phim) {
                        // dqbar-uq-phim
                        auto ckm = CKM(si1,si2);
                        return I*GEW/sqrt2*(M(si2)/MW*Matrix(GAS(1)+GAS(5),DI(fi1),DI(fi2))/2-M(si1)/MW*Matrix(GAS(1)-GAS(5),DI(fi1),DI(fi2))/2)*ckm*SP(TI(fi1),TI(fi2));
                    } else if(is_nbar_l(pi1, pi2) && pi3==phip) {
                        // nbar-l-phip
                        return -I*GEW/sqrt2*M(si2)/MW*Matrix(GAS(1)+GAS(5),DI(fi1),DI(fi2))/2; // Mnl = 0
                    } else if(is_lbar_n(pi1, pi2) && pi3==phim) {
                        // lbar-n-phim
                        return -I*GEW/sqrt2*M(si1)/MW*Matrix(GAS(1)-GAS(5),DI(fi1),DI(fi2))/2; // Mnl = 0
                    } else if(pi1==A && pi2==phip && pi3==phim) {
                        // A-phip-phim
                        return -I*eta_e*EL*SP(mom2-mom3,LI(fi1));
                    } else if(pi1==Z && pi2==phip && pi3==phim) {
                        // Z-phip-phim
                        return -I*eta*eta_Z*GEW*C2W/(2*CW)*SP(mom2-mom3,LI(fi1));
                    } else if(pi1==Wp && pi2==phim && pi3==H) {
                        // Wp-phim-H
                        return I/2*eta*GEW*SP(mom2-mom3,LI(fi1));
                    } else if(pi1==Wm && pi2==phip && pi3==H) {
                        // Wm-phip-H
                        return -I/2*eta*GEW*SP(mom2-mom3,LI(fi1));
                    } else if( ((pi1==Wp && pi2==phim) || (pi1==Wm && pi2==phip)) && pi3==chi) {
                        // Wp-phim-chi & Wm-phip-chi
                        return -eta*GEW/2*SP(mom2-mom3,LI(fi1));
                    } else if(pi1==Z && pi2==chi && pi3==H) {
                        // Z-chi-H
                        return -eta*eta_Z*GEW/(2*CW)*SP(mom2-mom3,LI(fi1));
                    } else if( ((pi1==phim && pi2==Wp) || (pi1==phip && pi2==Wm)) && pi3==A) {
                        // phim-Wp-A & phip-Wm-A
                        return I*eta_e*eta*EL*MW*SP(LI(fi2),LI(fi3));
                    } else if( ((pi1==phim && pi2==Wp) || (pi1==phip && pi2==Wm)) && pi3==Z) {
                        // phim-Wp-Z & phip-Wm-Z
                        return -I*eta_Z*GEW*MZ*SW2*SP(LI(fi2),LI(fi3));
                    } else if(pi1==H && pi2==Wp && pi3==Wm) {
                        // H-Wp-Wm
                        return I*GEW*MW*SP(LI(fi2),LI(fi3));
                    } else if(pi1==H && pi2==Z && pi3==Z) {
                        // H-Z-Z
                        return I*GEW/CW*MZ*SP(LI(fi2),LI(fi3));
                    } else if(pi1==H && pi2==phim && pi3==phip) {
                        // H-phim-phip
                        return -I/2*GEW*MH2/MW;
                    } else if(pi1==H && pi2==H && pi3==H) {
                        // H-H-H
                        return -3/ex(2)*I*GEW*MH2/MW;
                    } else if(pi1==H && pi2==chi && pi3==chi) {
                        // H-chi-chi
                        return -I/2*GEW*MH2/MW;
                    } else if(pi1==ghWpbar && pi2==ghWp && pi3==A) {
                        // ghWpbar-ghWp-A
                        return I*eta_G*eta_e*EL*SP(mom1,LI(fi3));
                    } else if(pi1==ghWmbar && pi2==ghWm && pi3==A) {
                        // ghWmbar-ghWm-A
                        return -I*eta_G*eta_e*EL*SP(mom1,LI(fi3));
                    } else if(pi1==ghWpbar && pi2==ghWp && pi3==Z) {
                        // ghWpbar-ghWp-Z
                        return I*eta_G*eta*eta_Z*GEW*CW*SP(mom1,LI(fi3));
                    } else if(pi1==ghWmbar && pi2==ghWm && pi3==Z) {
                        // ghWmbar-ghWm-Z
                        return -I*eta_G*eta*eta_Z*GEW*CW*SP(mom1,LI(fi3));
                    } else if(pi1==ghWpbar && pi2==ghZ && pi3==Wp) {
                        // ghWpbar-ghZ-Wp
                        return -I*eta_G*eta*eta_Z*GEW*CW*SP(mom1,LI(fi3));
                    } else if(pi1==ghWmbar && pi2==ghZ && pi3==Wm) {
                        // ghWmbar-ghZ-Wm
                        return I*eta_G*eta*eta_Z*GEW*CW*SP(mom1,LI(fi3));
                    } else if(pi1==ghWpbar && pi2==ghA && pi3==Wp) {
                        // ghWpbar-ghA-Wp
                        return -I*eta_G*eta_e*EL*SP(mom1,LI(fi3));
                    } else if(pi1==ghWmbar && pi2==ghA && pi3==Wm) {
                        // ghWmbar-ghA-Wm
                        return I*eta_G*eta_e*EL*SP(mom1,LI(fi3));
                    } else if(pi1==ghZbar && pi2==ghWp && pi3==Wm) {
                        // ghZbar-ghWp-Wm
                        return -I*eta_G*eta*GEW*CW*SP(mom1,LI(fi3));
                    } else if(pi1==ghZbar && pi2==ghWm && pi3==Wp) {
                        // ghZbar-ghWm-Wp
                        return I*eta_G*eta*GEW*CW*SP(mom1,LI(fi3));
                    } else if(pi1==ghAbar && pi2==ghWp && pi3==Wm) {
                        // ghAbar-ghWp-Wm
                        return -I*eta_G*eta_e*EL*SP(mom1,LI(fi3));
                    } else if(pi1==ghAbar && pi2==ghWm && pi3==Wp) {
                        // ghAbar-ghWm-Wp
                        return I*eta_G*eta_e*EL*SP(mom1,LI(fi3));
                    } else if(pi1==ghWpbar && pi2==ghWp && pi3==chi) {
                        // ghWpbar-ghWp-chi
                        return eta_G*GEW/2*xi*MW;
                    } else if(pi1==ghWmbar && pi2==ghWm && pi3==chi) {
                        // ghWmbar-ghWm-chi
                        return -eta_G*GEW/2*xi*MW;
                    } else if( ((pi1==ghWpbar && pi2==ghWp) ||(pi1==ghWmbar && pi2==ghWm)) && pi3==H) {
                        // ghWpbar-ghWp-H & ghWmbar-ghWm-H
                        return -I/2*eta_G*GEW*xi*MW;
                    } else if(pi1==ghZbar && pi2==ghZ && pi3==H) {
                        // ghZbar-ghZ-H
                        return -eta_G*I*GEW/(2*CW)*xi*MZ;
                    } else if( pi1==ghZbar && ((pi2==ghWp && pi3==phim) || (pi2==ghWm && pi3==phip)) ) {
                        // ghZbar-ghWp-phim & ghZbar-ghWm-phip
                        return I/2*eta_G*eta_Z*GEW*xi*MZ;
                    } else if( ((pi1==ghWpbar && pi3==phip) || (pi1==ghWmbar && pi3==phim)) && pi2==ghZ) {
                        // ghWpbar-ghZ-phip & ghWmbar-ghZ-phim
                        return -I*eta_G*eta_Z*GEW*C2W/(2*CW)*xi*MW;
                    } else if( ((pi1==ghWpbar && pi3==phip) || (pi1==ghWmbar && pi3==phim)) && pi2==ghA) {
                        // ghWpbar-ghA-phip & ghWmbar-ghA-phim
                        return -I*eta_G*eta_e*eta*EL*xi*MW;
                    } else {
                        cout << endl << e << endl;
                        throw Error("Vertex Not defined!");
                    }
                } else if(nn==4) {
                    auto pi4 = e.op(3).op(0);
                    auto fi4 = e.op(3).op(1);
                    auto mom4 = e.op(3).op(2);
                    
                    if(pi1==g && pi2==g && pi3==g && pi4==g) {
                        // g^4
                        return Gluon4Vertex(e, eta_s);
                    } else if(pi1==Wp && pi2==Wm && pi3==A && pi4==A) {
                        // Wp-Wm-A-A
                        return -I*EL2*(2*SP(LI(fi1),LI(fi2))*SP(LI(fi3),LI(fi4))-SP(LI(fi1),LI(fi3))*SP(LI(fi2),LI(fi4))-SP(LI(fi1),LI(fi4))*SP(LI(fi2),LI(fi3)));
                    } else if(pi1==Wp && pi2==Wm && pi3==Z && pi4==Z) {
                        // Wp-Wm-Z-Z
                        return -I*GEW2*CW2*(2*SP(LI(fi1),LI(fi2))*SP(LI(fi3),LI(fi4))-SP(LI(fi1),LI(fi3))*SP(LI(fi2),LI(fi4))-SP(LI(fi1),LI(fi4))*SP(LI(fi2),LI(fi3)));
                    } else if(pi1==Wp && pi2==Wm && pi3==A && pi4==Z) {
                        // Wp-Wm-A-Z
                        return -I*eta_e*eta*eta_Z*EL*GEW*CW*(2*SP(LI(fi1),LI(fi2))*SP(LI(fi3),LI(fi4))-SP(LI(fi1),LI(fi3))*SP(LI(fi2),LI(fi4))-SP(LI(fi1),LI(fi4))*SP(LI(fi2),LI(fi3)));
                    } else if(pi1==Wp && pi2==Wp && pi3==Wm && pi4==Wm) {
                        // Wp-Wp-Wm-Wm
                        return I*GEW2*(2*SP(LI(fi1),LI(fi2))*SP(LI(fi3),LI(fi4))-SP(LI(fi1),LI(fi3))*SP(LI(fi2),LI(fi4))-SP(LI(fi1),LI(fi4))*SP(LI(fi2),LI(fi3)));
                    } else if(pi1==Wp && pi2==Wm && pi3==H && pi4==H) {
                        // Wp-Wm-H-H
                        return I/2*GEW2*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Wp && pi2==Wm && pi3==chi && pi4==chi) {
                        // Wp-Wm-chi-chi
                        return I/2*GEW2*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Z && pi2==Z && pi3==H && pi4==H) {
                        // Z-Z-H-H
                        return I/2*GEW2/CW2*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Z && pi2==Z && pi3==chi && pi4==chi) {
                        // Z-Z-chi-chi
                        return I/2*GEW2/CW2*SP(LI(fi1),LI(fi2));
                    } else if(pi1==A && pi2==A && pi3==phip && pi4==phim) {
                        // A-A-phip-phim
                        return 2*I*EL2*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Z && pi2==Z && pi3==phip && pi4==phim) {
                        // Z-Z-phip-phim
                        return I/2*pow(GEW*C2W/CW,2)*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Wp && pi2==Wm && pi3==phip && pi4==phim) {
                        // Wp-Wm-phip-phim
                        return I/2*GEW2*SP(LI(fi1),LI(fi2));
                    } else if( ((pi1==Wp && pi3==phim) || (pi1==Wm && pi3==phip)) && pi2==Z && pi4==H) {
                        // Wp-Z-phim-H & Wm-Z-phip-H
                        return -I*eta_Z*GEW2*SW2/(2*CW)*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Wm && pi2==Z && pi3==phip && pi4==chi) {
                        // Wm-Z-phip-chi
                        return -eta_Z*GEW2*SW2/(2*CW)*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Wp && pi2==Z && pi3==phim && pi4==chi) {
                        // Wp-Z-phim-chi
                        return eta_Z*GEW2*SW2/(2*CW)*SP(LI(fi1),LI(fi2));
                    } else if( ((pi1==Wm && pi3==phip) || (pi1==Wp && pi3==phim)) && pi2==A && pi4==H) {
                        // Wm-A-phip-H & Wp-A-phim-H
                        return I/2*eta_e*eta*EL*GEW*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Wp && pi2==A && pi3==phim && pi4==chi) {
                        // Wp-A-phim-chi
                        return -1/ex(2)*eta_e*eta*EL*GEW*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Wm && pi2==A && pi3==phip && pi4==chi) {
                        // Wm-A-phip-chi
                        return 1/ex(2)*eta_e*eta*EL*GEW*SP(LI(fi1),LI(fi2));
                    } else if(pi1==Z && pi2==A && pi3==phip && pi4==phim) {
                        // Z-A-phip-phim
                        return I*eta_e*eta*eta_Z*EL*GEW*C2W/CW*SP(LI(fi1),LI(fi2));
                    } else if(pi1==phim && pi2==phip && pi3==phim && pi4==phip) {
                        // phim-phip-phim-phip
                        return -I/2*GEW2*MH2/MW2;
                    } else if(pi1==H && pi2==H && pi3==phim && pi4==phip) {
                        // H-H-phim-phip
                        return -I/4*GEW2*MH2/MW2;
                    } else if(pi1==chi && pi2==chi && pi3==phim && pi4==phip) {
                        // chi-chi-phim-phip
                        return -I/4*GEW2*MH2/MW2;
                    } else if(pi1==H && pi2==H && pi3==H && pi4==H) {
                        // H-H-H-H
                        return -3/ex(4)*I*GEW2*MH2/MW2;
                    } else if(pi1==H && pi2==H && pi3==chi && pi4==chi) {
                        // H-H-chi-chi
                        return -I/4*GEW2*MH2/MW2;
                    } else if(pi1==chi && pi2==chi && pi3==chi && pi4==chi) {
                        // chi-chi-chi-chi
                        return -3/ex(4)*I*GEW2*MH2/MW2;
                    } else {
                        cout << endl << e << endl;
                        throw Error("Vertex Not defined!");
                    }
                } else {
                    cout << endl << e << endl;
                    throw Error("Vertex Not defined!");
                }
            } else return e.map(self);
            return e;
        });
        return fr(amp);
    }
        
        
    lst Models::FeynRulesQCD(const lst & amps, const ex & xi) {
        lst ret;
        for(auto item : amps) ret.append(FeynRulesQCD(item,xi));
        return ret;
    }
    ex Models::FeynRulesQCD(const ex & amp, const ex & xi) {
        if(is_a<lst>(amp)) return FeynRulesQCD(ex_to<lst>(amp), xi);
        
        static Symbol A("A"), Q("Q"), Qbar("Qbar"), q("q"), qbar("qbar"), g("g"), gh("gh"), ghbar("ghbar");
        static Symbol m("m"), eq("eq"), eQ("eQ");
        
        auto fr = MapFunction([&](const ex &e, MapFunction &self)->ex {
            if(isFunction(e,"OutField") || isFunction(e,"InField")) return 1;
            else if(isFunction(e, "Propagator")) {
                if(e.op(0).op(0)==q) {
                    return QuarkPropagator(e, 0);
                } else if(e.op(0).op(0)==Q) {
                    return QuarkPropagator(e, m);
                } else if(e.op(0).op(0)==g) {
                    return GluonPropagator(e, xi);
                } else if(e.op(0).op(0)==gh) {
                    return GluonGhostPropagator(e, xi);
                } else {
                    cout << "expr: " << e << endl;
                    throw Error("Propagator Not Found!");
                }
            } else if(isFunction(e, "Vertex")) {
                if(e.nops()==3 && ((e.op(0).op(0)==qbar && e.op(1).op(0)==q) || (e.op(0).op(0)==Qbar && e.op(1).op(0)==Q)) && (e.op(2).op(0)==g || e.op(2).op(0)==A) ) {
                    // qbar-q-g or Qbar-Q-g or g -> A
                    if(e.op(2).op(0)==g) return QuarkGluonVertex(e);
                    else {
                        auto fi1 = e.op(0).op(1);
                        auto fi2 = e.op(1).op(1);
                        auto fi3 = e.op(2).op(1);
                        if(e.op(1).op(0)==q) QuarkAVertex(e, eq);
                        else if(e.op(1).op(0)==Q) QuarkAVertex(e, eQ);
                        else throw Error("Vertex Error.");
                    }
                } else if(e.nops()==3 && e.op(0).op(0)==ghbar && e.op(1).op(0)==gh) {
                    // ghbar-gh-g
                    return GhostGluonVertex(e);
                } else if(e.nops()==3 && e.op(0).op(0)==g && e.op(1).op(0)==g && e.op(2).op(0)==g) {
                    // g^3
                    return Gluon3Vertex(e);
                } else if(e.nops()==4 && e.op(0).op(0)==g && e.op(1).op(0)==g && e.op(2).op(0)==g && e.op(3).op(0)==g) {
                    // g^4
                    return Gluon4Vertex(e);
                } else {
                    cout << "expr: " << e << endl;
                    throw Error("Vertex Not Found!");
                }
            } else return e.map(self);
            return e;
        });
        return fr(amp);
    }
    
}

