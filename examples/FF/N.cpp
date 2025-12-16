#include "HepLib.h"
using namespace HepLib;
using namespace SD;

Symbol kp("kp"), pp("pp"), m("m");
Symbol p("p"), n("n"), k1("k1"), k2("k2"), k3("k3"), q1("q1"), q2("q2");
Symbol K1("K1"), K2("K2"), K3("K3");

string cm_path = "cm"; // input directory
string SD_path = "SD"; // save directory
int loops = 1;
int tloops = 1;
ex SF = 1;

ex zz = 45/ex(100);
bool use_eps = true;
bool isCut = true;

void extra_info() {
    if(Verbose>0) cout << "  z = " << zz << endl;
}

void ExportNull(const string & prefix);
void Prepare(int idx) {

    extra_info();
    
    Symbol::set(NA,8);
    Symbol::set(NF,3);
    Symbol::set(kp,1);
    Symbol::set(m,1);
    Symbol::set(d,4-2*ep);
    Symbol::set(gs,sqrt(4*Pi));
    Symbol::set(Symbol("z"),zz);
    
    ex m2 = m*m;
    ex k1p = z(1)*(1-zz)*kp;
    ex k1m = pow(K1,2)/(2*k1p);
    ex k2p = z(2)*(1-zz)*kp;
    ex k2m = pow(K2,2)/(2*k2p);
    
    ex pp = zz/2*kp;
    ex pm = m2/(2*pp);
    ex ppv = pp;
    
    Symbol::set("pp",pp);
    
    ex cf = file2ex(cm_path + "/" + to_string(idx) + ".txt");
    
    auto cv_lst = collect_lst(cf,F(w1,w2));
    if(cv_lst.nops()!=1) throw Error("something is wrong, not the format: c * F(w1,w2).");

    auto pns = cv_lst.op(0).op(1);
    auto pref = cv_lst.op(0).op(0);
    
    FeynmanParameter fp;
    fp.LoopMomenta = lst{ };
    fp.tLoopMomenta = lst{ };
    if(loops>0) fp.LoopMomenta.append(q1);
    if(loops>1) fp.LoopMomenta.append(q2);
    if(tloops>0) fp.tLoopMomenta.append(K1);
    if(tloops>1) fp.tLoopMomenta.append(K2);
    
    fp.lReplacement[p*p] = m2;
    fp.lReplacement[n*n] = 0;
    fp.lReplacement[p*n] = pp;
    fp.lReplacement[k1*k1] = 0;
    fp.lReplacement[k2*k2] = 0;
    fp.lReplacement[n*k1] = k1p;
    fp.lReplacement[n*k2] = k2p;
    fp.lReplacement[p*k1] = k1p*pm+k1m*pp;
    fp.lReplacement[p*k2] = k2p*pm+k2m*pp;
    fp.lReplacement[k1*k2] = k1p*k2m+k2p*k1m-K1*K2;
            
    fp.tReplacement[p*p] = m2;
    fp.tReplacement[n*n] = 0;
    fp.tReplacement[p*n] = pp;
    
    if(tloops>0) fp.nReplacement[z(w)] = ex(1)/tloops;
    
    fp.Prefactor = pow(2*Pi, (2*ep-4)*loops) * pref;
    for(int i=0; i<pns.op(0).nops(); i++) {
        if(isCut && i<tloops+1) {
            if(pns.op(1).op(i)!=1) throw Error("Cut error");
            continue;
        }
        if(is_zero(pns.op(1).op(i))) continue;
        fp.Propagator.append(pns.op(0).op(i).expand());
        fp.Exponent.append(pns.op(1).op(i));
    }
    
    ex M = 2*m;
    ex zPrefactor = QCD::FF::zIntFactor(tloops, SF, kp, zz, M);
    fp.Prefactor = fp.Prefactor * zPrefactor;
    
    SecDec work;
    work.CheckEnd = true;
    work.eps_lst = str2lst("{{eps,0},{ep,0}}");
    work.Initialize(fp);
    
    if(tloops>0)
    for(auto &fe : work.FunExp) {
        int xn = pns.op(0).nops();
        exmap z2x;
        lst zs;
        ex zFactor = 1;
        for(int i=1; i<=tloops; i++) {
            auto xi = x(xn+i-1);
            if(fe.has(xi)) {
                cout << "fe = " << fe << endl;
                throw Error("Check failed: fe has xi.");
            }
            z2x[z(i)] = xi;
            zs.append(xi);
            zFactor /= xi;
        }
        let_op_append(fe, 2, zs);
        
        fe.let_op(0).let_op(0) = fe.op(0).op(0) * zFactor;
        fe.let_op(0) = subs(fe.op(0), z2x);
        
        auto tmp = Factor(fe.op(0).op(0)).subs(Symbol::vmap);

        if(tmp.has(x(w)) && is_a<mul>(tmp)) {
            ex rem = 1;
            for(auto item : tmp) {
                if(!item.has(x(w))) {
                    rem *= item;
                } else if(item.match(pow(w0, w1))) {
                    let_op_append(fe, 0, item.op(0));
                    let_op_append(fe, 1, item.op(1));
                } else {
                    if(!item.is_polynomial(zs)) {
                        cout << item << endl;
                        throw Error("Factor failed, item is not a polynormial w.r.t zs");
                    }
                    let_op_append(fe, 0, item);
                    let_op_append(fe, 1, 1);
                }
            }
            fe.let_op(0).let_op(0) = rem;
        } else if(tmp.has(x(w)) && tmp.match(pow(w0, w1))) {
            fe.let_op(0).let_op(0) = 1;
            let_op_append(fe, 0, tmp.op(0));
            let_op_append(fe, 1, tmp.op(1));
        } 
        
        if(use_eps) {
            for(auto zi : zs) {
                let_op_append(fe, 0, zi);
                let_op_append(fe, 1, eps);
            }
        }

    }
    
    work.Normalizes();
    work.XReOrders();
    work.Normalizes();
    work.Scalelesses();
    work.ChengWu();
    work.RemoveDeltas();
    work.SDPrepares();
    work.EpsExpands();
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    work.CIPrepares(ikey.str());

    if(work.IsZero) {
        ostringstream ifn;
        ifn << SD_path << "/" << idx;
        ExportNull(ifn.str());
        return;
    }
}

void Contour(int idx) {
    extra_info();
    
    SecDec work;
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    work.Contours(ikey.str());
}

ex Integrate(int idx, int ii = -1) {
    extra_info();
    
    SecDec work;
    work.eps_lst = str2lst("{{eps,0},{ep,0}}");
                
    work.EpsAbs = 1E-3;
    work.LambdaSplit = 5;
    work.CTTryPTS = 100000;
    work.CTryL = 1;
    work.CTryR = 1;
    work.CTryM = 1;
    work.CTLaMax = 1;
    work.CTryRRatio = 1;
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
        
    auto intor = new HCubature();
    intor->MPXLimit = 1E-3;
    intor->QXLimit = 5E-1;
    
    //auto intor = new HCubatureMP();
    //intor->MPDigits = 100;
    
    intor->RunPTS = 100000;
    intor->RunMAX = 10;
    //intor->use_last = true;
    work.Integrator = intor;
    
    work.Integrates(ikey.str().c_str(), "", ii);
    
    return work.ResultError;
}

ex ReIntegrate(int idx, qREAL err) {
    extra_info();
    
    SecDec work;
    
    work.eps_lst = str2lst("{{eps,0},{ep,0}}");
                
    work.EpsAbs = 1E-4;
    work.LambdaSplit = 5;
    work.CTTryPTS = 100000;
    work.CTryL = 1;
    work.CTryR = 1;
    work.CTryM = 1;
    work.CTLaMax = 1;
    work.CTryRRatio = 1;
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
        
    auto intor = new HCubatureMP();
    intor->MPDigits = 100;
    
    //auto intor = new HCubature();
    //intor->MPXLimit = 1E-2;
    //intor->QXLimit = 5E-1;
    //intor->MPXDim = 1;
    //intor->QXDim = 2;
    
    intor->RunPTS = 100000;
    intor->RunMAX = 30;
    //intor->use_last = true;
    work.Integrator = intor;
    
    work.ReIntegrates(ikey.str().c_str(), "", err);
    
    return work.ResultError;
}

int main(int argc, char** argv) {

    //Debug = true;
    
    const char* arg_a = NULL;
    const char* arg_n = NULL;
    const char* arg_i = NULL;
    const char* arg_v = "0";
    const char* arg_e = "0";
    bool is_o = false;
    bool to_skip = false;
    const char* arg_o = "";
    const char* arg_l = "0";
    
    string exe_name(argv[0]);
    // handle options
    for (int opt; (opt = getopt(argc, argv, "a:n:i:o:v:e:z:se:v:")) != -1;) {
        switch (opt) {
            case 'a': arg_a = optarg; break;
            case 'n': arg_n = optarg; break;
            case 'i': arg_i = optarg; break;
            case 'o': arg_o = optarg; is_o = true; break;
            case 'v': arg_v = optarg; break;
            case 'e': arg_e = optarg; break;
            case 's': to_skip = true; break;
            case 'z': zz = stoi(optarg)/ex(100); SD_path = SD_path + "/" + optarg;  break;
            default:
                cout << "supported options: -a A -n N -i I -o O -v V -e E -z Z -s -v V" << endl;
                cout << "A: can be p(prepare), c(contour), i(integrate), r(result), a(analyse)" << endl;
                cout << "N: only apply A on prefix N." << endl;
                cout << "I: only apply A(=i) on part I in prefix N." << endl;
                cout << "O: output filename." << endl;
                cout << "V: integer value, Verbose = V" << endl;
                cout << "E: print when error > E." << endl;
                cout << "Z: integer value, the actual z = Z/100." << endl;
                cout << "-s: if <n>.res.gar exists then skip." << endl;
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
            
    ostringstream cmd;
    cmd << "mkdir -p " << SD_path;
    system(cmd.str().c_str());
    cmd.clear();
    cmd.str("");
    cmd << "ls " << cm_path << "|wc -l";
    auto owc = RunOS(cmd.str().c_str());
    int nmi;
    stringstream ss;
    ss << owc;
    ss >> nmi;
    
    int in = 0;
    if(arg_n != NULL) in = stoi(arg_n);
    int ii = -1;
    if(arg_n != NULL && arg_i != NULL) ii = stoi(arg_i);
    Verbose = stoi(arg_v);
    numeric ee(arg_e);
    
    auto is_to_skip = [&](int i) -> bool {
        return file_exists(SD_path+"/"+to_string(i)+".res.gar") || file_exists(SD_path+"/"+to_string(i)+".null");
    };
    
    // call Prepare(index)
    if(arg_a == NULL || !strcmp(arg_a, "p")) {
        vector<int> fid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            if(to_skip && is_to_skip(i)) continue;
            try {
                cout << "Preparing " << i << "/" << nmi << endl;
                Prepare(i);
            } catch(exception& e) {
                cout << e.what() << endl;
                fid.push_back(i);
            } catch(...) {
                cout << "other uncatch error" << endl;
                fid.push_back(i);
            }
            cout << endl;
        }
        
        if(fid.size()>0) {
            cout << "Failed IDs: {";
            for(auto item : fid) cout << item << ", ";
            cout << "}" << endl;
            return 0;
        }
    }

    // call Contour(index)
    if(arg_a == NULL || !strcmp(arg_a, "c") || !strcmp(arg_a, "ci")) {
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            if(to_skip && is_to_skip(i)) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << "Contouring - " << i << "/" << nmi << endl;
            Contour(i);
            cout << endl;
        }
    }
    
    // call Integrate(index, sub_index)
    if(arg_a == NULL || !strcmp(arg_a, "i") || !strcmp(arg_a, "ci")) {
        ex fRes = 0;
        vector<int> nid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            if(to_skip && is_to_skip(i)) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << "Integrating - " << i << "/" << nmi << endl;
            auto res = Integrate(i, ii);
            fRes += res;
            cout << endl;
            if(res.has(SD::NaN)) nid.push_back(i);
            else if(ee>0 && VEMaxErr(res)>ee) {
                cout << "n = " << i << endl;
                cout << res << endl << endl;
            }
        }
        
        if(nid.size()>0) {
            cout << "NaN IDs: {";
            for(auto item : nid) cout << item << ", ";
            cout << "}" << endl;
            return 0;
        }
    }
    
    // call ReIntegrate(index, err)
    if(arg_a != NULL && !strcmp(arg_a, "er") && ee>0) {
        ex fRes = 0;
        vector<int> nid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            if(to_skip && is_to_skip(i)) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << "ReIntegrating - " << i << "/" << nmi << endl;
            auto res = ReIntegrate(i, ex2q(ee));
            fRes += res;
            cout << endl;
            if(res.has(SD::NaN)) nid.push_back(i);
            else if(ee>0 && VEMaxErr(res)>ee) {
                cout << "n = " << i << endl;
                cout << res << endl << endl;
            }
        }
        
        if(nid.size()>0) {
            cout << "NaN IDs: {";
            for(auto item : nid) cout << item << ", ";
            cout << "}" << endl;
            return 0;
        }
    }
    
    // show the result
    if(arg_a == NULL || !strcmp(arg_a, "r")) {
        ex fRes = 0;
        vector<int> nid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) {
                ostringstream ifn;
                ifn << SD_path << "/" << i << ".null";
                if(!file_exists(ifn.str().c_str())) {
                    cout << RED << "File NOT Found: " << ifn.str() << RESET << endl;
                    exit(1);
                }
                continue;
            }
            
            ss.clear();
            ss.str("");
            ss << SD_path << "/" << i << ".res.gar";
            if(!file_exists(ss.str().c_str())) {
                cout << RED << "File NOT Found: " << ss.str() << RESET << endl;
                exit(1);
            }
            auto res = garRead(ss.str());
            fRes += res;
            if(res.has(SD::NaN)) nid.push_back(i);            
            else if(ee>0 && VEMaxErr(res)>ee) {
                cout << "n = " << i << endl;
                cout << res << endl << endl;
            }
        }
        
        if(nid.size()>0) {
            cout << "NaN IDs: {";
            for(auto item : nid) cout << item << ", ";
            cout << "}" << endl;
            return 0;
        }
        
        fRes = VESimplify(fRes);
        extra_info();
        cout << endl << "Final Result";
        if(in>0) cout << "[ n=" << in << " ]";
        cout << ":" << endl;
        cout << "------" << endl;
        cout << fRes << endl;
        cout << "------" << endl;
        cout << VEResult(fRes) << endl << endl << endl;
        
        if(is_o) {
            ofstream ofs;
            system(("mkdir -p $(dirname '"+ string(arg_o) +"')").c_str());
            ofs.open(arg_o, ios::out);
            if (!ofs) throw runtime_error("failed to open final out file!");
            ofs << fRes << endl;
            ofs.close();
        }
    }
    
    // analyse result
    if(arg_a != NULL && !strcmp(arg_a, "a")) {
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".res.gar";
            map<string, ex> gar;
            garRead(ss.str(), gar);
            auto relst = ex_to<lst>(gar["relst"]);
            int max_index;
            ex max_re = -1;
            for(int ii=0; ii<relst.nops(); ii++) {
                auto tmp = VEMaxErr(relst.op(ii));
                if(abs(tmp)>max_re) {
                    max_re = abs(tmp);
                    max_index = ii;
                }
                if(ee>0 && tmp>ee) {
                    cout << exe_name << " -n " << i << " -a i -i " << (ii+1) << " -v " << Verbose << endl;
                    cout << "# " << tmp << endl;
                    cout << endl;
                }
            }
            cout << "# [" << i << "] Max Index: " << max_index+1 << " / " << relst.nops() << endl;
            cout << "# " << max_re << endl;
            cout << endl << endl;
        }
    }

    return 0;
}

void ExportNull(const string & prefix) {
    ostringstream cmd;
    cmd << "rm -f " << prefix << "[-.]*";
    system(cmd.str().c_str());
    
    ofstream ofs;
    ofs.open(prefix+".null", ios::out);
    if (!ofs) throw runtime_error("failed to open final null file!");
    ofs << "Result is Null." << endl;
    ofs.close();
}
