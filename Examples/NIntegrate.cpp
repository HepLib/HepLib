#include "FC.h"
#include "SD.h"

using namespace HepLib;
using namespace FC;
using namespace Quarkonium;
using namespace SD;

string cm_path = "cm"; // input directory
string SD_path = "SD"; // save directory

Symbol z1("z1"), z2("z2");
Symbol p("p"), q1("q1"), k1s("k1s"), k2s("k2s");
int epN = 0;

void ExportNull(const string & fn) {
    ofstream ofs;
    ofs.open(fn, ios::out);
    if (!ofs) throw runtime_error("failed to open final null file!");
    ofs << "Result is Null." << endl;
    ofs.close();
}

void Prepare(int idx) {
    
    ex cf = file2ex(cm_path + "/" + to_string(idx) + ".txt");
    cf = mma_collect(cf,F(w),true,true);
    if(!is_a<mul>(cf)) throw Error("something is wrong.");
    auto ps = cf.subs(lst{coVF(w)==w,coCF(w)==1}).op(0);
    auto pref = cf.subs(lst{coVF(w)==1,coCF(w)==w});
    exmap ps_map;
    for(auto item : ps) ps_map[item] += 1;
    
    pref *= pow(1-z1,-ep) * pow(z1,1-2*ep) * pow(1-z2,-ep) * pow(z2,-ep) / (3456*pow(Pi,3));
    pref *= (243*pow(2,2+2*ep)*pow(Pi,1+ep)*pow(16*(-9+pow(Pi,2))+ep*(84-87*pow(Pi,2)+672*zeta(3)),-1)*tgamma(1-ep))/5;
    
    FeynmanParameter fp;
    for(auto kv : ps_map) {
        fp.Propagators.append(kv.first);
        fp.Exponents.append(kv.second);
    }
    fp.LoopMomenta = lst{ q1 };
    fp.lReplacements[p*p] = 1;
    fp.lReplacements[k1s*k1s] = 0;
    fp.lReplacements[k2s*k2s] = 0;
    fp.lReplacements[k1s*p] = z1;
    fp.lReplacements[k2s*p] = 1-z1*z2;
    fp.lReplacements[k1s*k2s] = 2*(z1-z1*z2);
    
    SecDec work;       
    work.epN = epN;
    work.CheckEnd = true;
    
    work.Initialize(fp);
    if(work.IsZero) {
        ostringstream ifn;
        ifn << SD_path << "/" << idx << ".null";
        ExportNull(ifn.str().c_str());
        ostringstream cmd;
        cmd << "rm " << SD_path << "/" << idx << "[-.]*";
        system(cmd.str().c_str());
        return;
    }
    
    for(auto &fe : work.FunExp) {
        int xn = fe.op(2).op(0).nops();
        exmap z2x;
        lst zs = {x(xn+1),x(xn+2)};
        z2x[z1]=x(xn+1);
        z2x[z2]=x(xn+2);
        
        fe.let_op(0).let_op(0) = fe.op(0).op(0) * pref;
        fe.let_op(0) = subs(fe.op(0), z2x);
        
        auto tmp = Factor(fe.op(0).op(0));

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
        
    }
    
    work.Normalizes();
    work.XReOrders();
    work.Normalizes();
    work.Scalelesses();
    work.ChengWu();
    work.RemoveDeltas();
    work.SDPrepares();
    work.EpsEpExpands();
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    work.CIPrepares(ikey.str());

    if(work.IsZero) {
        ostringstream ifn;
        ifn << SD_path << "/" << idx << ".null";
        ExportNull(ifn.str());
    }
}

void Contour(int idx) {
    SecDec work;    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    work.Contours(ikey.str());
}

ex Integrate(int idx, int ii = -1) {
    SecDec work;
    work.epN = epN;
    
    work.use_ErrMin = false;
    ErrMin::err_min = 1E-3;
    ErrMin::MaxRND = 20;
    ErrMin::hjRHO = 0.5;
            
    work.EpsAbs = 1E-5;
    work.RunPTS = 1000000;
    work.RunMAX = 20;
    work.LambdaSplit = 5;
    work.TryPTS = 1000000;
    work.CTryLeft = 1;
    work.CTryRight = 1;
    work.CTry = 1;
    work.ReIm = 1;
    work.CTryRightRatio = 2;
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
        
    auto intor = new HCubature();
    intor->MPXLimit = 5E-3;
    intor->QXLimit = 5E-1;
    //intor->use_last = true;
    work.Integrator = intor;
    
    work.Integrates(ikey.str().c_str(), "", ii);
    
    return work.ResultError;
}

int main(int argc, char** argv) {

    //SD::debug = true;
    char* dc = getenv("SD_DLCLOSE");
    if(dc!=NULL && !strcmp(dc, "0")) SecDec::use_dlclose = false;
    
    const char* arg_a = NULL;
    const char* arg_n = NULL;
    const char* arg_i = NULL;
    const char* arg_v = "1";
    const char* arg_e = "1E-2";
    bool is_o = false;
    const char* arg_o = "";
    
    // handle options
    for (int opt; (opt = getopt(argc, argv, "a:n:i:o:v:e:")) != -1;) {
        switch (opt) {
            case 'a': arg_a = optarg; break;
            case 'n': arg_n = optarg; break;
            case 'i': arg_i = optarg; break;
            case 'o': arg_o = optarg; is_o = true; break;
            case 'v': arg_v = optarg; break;
            case 'e': arg_e = optarg; break;
            default:
                cout << "supported options: -a A -n N -i I -o O -v V -e E" << endl;
                cout << "A: can be p(prepare), c(contour), i(integrate), r(result), ar(analyse result)." << endl;
                cout << "N: only apply A on prefix N." << endl;
                cout << "I: only apply A(=i) on part I in prefix N." << endl;
                cout << "O: output filename." << endl;
                cout << "V: verbose level." << endl;
                cout << "E: print when error > E for specific N." << endl;
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
    
    ostringstream cmd;
    cmd << "mkdir -p " << SD_path;
    if(!dir_exists(SD_path)) system(cmd.str().c_str());
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
    
    // call Prepare(index)
    if(arg_a == NULL || !strcmp(arg_a, "p")) {
        vector<int> fid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            try {
                cout << WHITE << "-" << i << "/" << nmi << " :" << RESET << endl;
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
    
    // call Prepare(index)
    if(arg_a == NULL || !strcmp(arg_a, "c") || !strcmp(arg_a, "ci")) {
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << WHITE << "-" << i << "/" << nmi << " :" << RESET << endl;
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
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << WHITE << "-" << i << "/" << nmi << " :" << RESET << endl;
            auto res = Integrate(i, ii);
            fRes += res;
            cout << endl;
            if(res.has(SD::NaN)) nid.push_back(i);
            
            auto err = res.subs(lst{VE(w1,w2)==w2, CV(w1,w2)==w2});
            if(abs(err.coeff(ep,epN).coeff(eps,0))>ee) {
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
                    cout << Color_Error << "File NOT Found: " << ifn.str() << RESET << endl;
                    exit(1);
                }
                continue;
            }
            
            ss.clear();
            ss.str("");
            ss << SD_path << "/" << i << ".res.gar";
            if(!file_exists(ss.str().c_str())) {
                cout << Color_Error << "File NOT Found: " << ss.str() << RESET << endl;
                exit(1);
            }
            auto res = garRead(ss.str());
            fRes += res;
            if(res.has(SD::NaN)) nid.push_back(i);
            
            auto err = res.subs(lst{VE(w1,w2)==w2,CV(w1,w2)==w2});
            if(abs(err.coeff(ep,epN).coeff(eps,0))>ee) {
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
        
        fRes = VESimplify(fRes, epN);
        cout << endl << "Final Result";
        if(in>0) cout << " [ n=" << in << " ]";
        cout << ":" << endl;
        cout << "------" << endl;
        cout << fRes << endl;
        cout << "------" << endl;
        cout << VEResult(fRes) << endl << endl << endl;
        
        if(is_o) {
            ofstream ofs;
            ofs.open(arg_o, ios::out);
            if (!ofs) throw runtime_error("failed to open final out file!");
            ofs << fRes << endl;
            ofs.close();
        }
    }
    
    // analyse result
    if(arg_a != NULL && (!strcmp(arg_a, "ar") || !strcmp(arg_a, "ae"))) {
        if(in<1) {
            cout << "n should be setted with -a ar option!" << endl;
            return 1;
        }
        
        bool ae = !strcmp(arg_a, "ae");
        stringstream ss;
        ss << SD_path << "/" << in << ".res.gar";
        map<string, ex> gar;
        garRead(ss.str(), gar);
        auto relst = ex_to<lst>(gar["relst"]);
        int max_index;
        ex max_re = -1;
        for(int i=0; i<relst.nops(); i++) {
            auto tmp = relst.op(i).subs(lst{VE(w1,w2)==w2,CV(w1,w2)==w2});
            tmp = tmp.subs(VE(w1, w2)==w2);
            tmp = tmp.expand().coeff(ep, epN).coeff(eps,0).evalf();
            if(abs(tmp)>max_re) {
                max_re = abs(tmp);
                max_index = i;
            }
            if(ae && abs(tmp)>ee) {
                cout << "# " << VESimplify(relst.op(i).expand().coeff(ep, epN)) << endl;
                cout << "./NIntegrate -n " << arg_n << " -a i -i " << (i+1) << endl;
            }
        }
        cout << "# Max Index: " << max_index+1 << " / " << relst.nops() << endl;
        cout << "# " << VESimplify(relst.op(max_index).expand().coeff(ep, epN)) << endl;
    }

    return 0;
}
