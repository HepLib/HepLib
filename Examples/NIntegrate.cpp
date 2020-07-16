#include "FC.h"
#include "SD.h"

using namespace HepLib;
using namespace FC;
using namespace Quarkonium;
using namespace SD;

string cm_path = "cm"; // input directory
string SD_path = "SD"; // save directory

numeric sN("111936/1000"), mN("168/100");
Symbol s("s"), m("m"), p("p"), q1("q1"), q2("q2"), k("k");
Symbol Chi0("Chi0"), Chi1("Chi1"), Chi2("Chi2");
Symbol zm("zm");
int loops;
int epN = 0;

void ExportNull(const string & prefix);
void Prepare(int idx) {
    
    Symbol::Assign(NA,8);
    Symbol::Assign(NF,3);
    Symbol::Assign(s,sN);
    Symbol::Assign(m,mN);
    Symbol::Assign(D,4-2*ep);
    Symbol::Assign(gs,sqrt(4*Pi)); // as to ZZ with asBare
    
    if(loops==0) Symbol::Assign(zm,RC::Zm(m,2)-1);
    else if(loops==1) Symbol::Assign(zm,RC::Zm(m,1)-1);
    else if(loops==2) Symbol::Assign(zm,RC::Zm(m,0)-1);
    
    ex ZZ = 1;
    if(loops==0) ZZ *= pow(RC::Z2("Q",m,2),2) * pow(RC::asBare(2),0);
    else if(loops==1) ZZ *= pow(RC::Z2("Q",m,1),2) * pow(RC::asBare(1),1);
    else if(loops==2) ZZ *= pow(RC::Z2("Q",m,0),2) * pow(RC::asBare(0),2);
    
    ex cf = file2ex(cm_path + "/" + to_string(idx) + ".txt");
    
    auto cv_lst = mma_collect_lst(cf,F(w1,w2));
    if(cv_lst.nops()!=1) throw Error("something is wrong.");

    auto pns = cv_lst.op(0).op(1);
    auto pref = cv_lst.op(0).op(0);
    pref *= ZZ; // Renormalization Constant
    pref = pref.subs(Symbol::AssignMap).subs(Symbol::AssignMap);
    pref = mma_series(pref,as,2);
    
    if(is_zero(pns-1)) {
        auto res = mma_series(pref,ep,epN);
        cout << res << endl;
        ostringstream ifn;
        ifn << SD_path << "/" << idx;
        ExportNull(ifn.str());
        return;
    }
    
    FeynmanParameter fp;
    if(loops==1) fp.LoopMomenta = lst{q1};
    else if(loops==2) fp.LoopMomenta = lst{q1, q2};
    
    for(int i=0; i<pns.op(0).nops(); i++) {
        fp.Propagators.append(pns.op(0).op(i).expand());
        fp.Exponents.append(pns.op(1).op(i));
    }
    fp.Prefactor = pref * pow(2*Pi, (2*ep-4)*loops);
    fp.lReplacements[p*p] = m*m;
    fp.lReplacements[k*k] = 0;
    fp.lReplacements[p*k] = s/4-m*m;
    
    SecDec work;       
    work.epN = epN;
    work.CheckEnd = true;

    work.Initialize(fp);
    if(work.IsZero) {
        ostringstream ifn;
        ifn << SD_path << "/" << idx;
        ExportNull(ifn.str());
        return;
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
            
    work.EpsAbs = 1E-3;
    work.RunPTS = 300000;
    work.RunMAX = 10;
    work.LambdaSplit = 5;
    work.TryPTS = 300000;
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
    const char* arg_v = "0";
    const char* arg_e = "1";
    bool is_o = false;
    const char* arg_o = "";
    const char* arg_l = "0";
    
    // handle options
    for (int opt; (opt = getopt(argc, argv, "a:n:i:o:v:e:l:")) != -1;) {
        switch (opt) {
            case 'a': arg_a = optarg; break;
            case 'n': arg_n = optarg; break;
            case 'i': arg_i = optarg; break;
            case 'o': arg_o = optarg; is_o = true; break;
            case 'v': arg_v = optarg; break;
            case 'e': arg_e = optarg; break;
            case 'l': arg_l = optarg; break;
            default:
                cout << "supported options: -a A -n N -i I -o O -v V -e E -l L" << endl;
                cout << "A: can be p(prepare), c(contour), i(integrate), r(result), a(analyse)" << endl;
                cout << "N: only apply A on prefix N." << endl;
                cout << "I: only apply A(=i) on part I in prefix N." << endl;
                cout << "O: output filename." << endl;
                cout << "V: verbose level." << endl;
                cout << "E: print when error > E." << endl;
                cout << "L: the loops: 0, 1 or 2." << endl;
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
    
    loops = stoi(arg_l);
    cm_path = cm_path + to_string(loops);
    SD_path = SD_path + to_string(loops);
    
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
    
    // call Prepare(index)
    if(arg_a == NULL || !strcmp(arg_a, "p")) {
        vector<int> fid;
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
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

    // call Prepare(index)
    if(arg_a == NULL || !strcmp(arg_a, "c") || !strcmp(arg_a, "ci")) {
        for(int i=1; i<=nmi; i++) {
            if(in>0 && i!=in) continue;
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
            stringstream ss;
            ss << SD_path << "/" << i << ".ci.gar";
            if(!file_exists(ss.str().c_str())) continue;
            cout << "Integrating - " << i << "/" << nmi << endl;
            auto res = Integrate(i, ii);
            fRes += res;
            cout << endl;
            if(res.has(SD::NaN)) nid.push_back(i);
            else if(VEMaxErr(res)>ee) {
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
                    cout << ErrColor << "File NOT Found: " << ifn.str() << RESET << endl;
                    exit(1);
                }
                continue;
            }
            
            ss.clear();
            ss.str("");
            ss << SD_path << "/" << i << ".res.gar";
            if(!file_exists(ss.str().c_str())) {
                cout << ErrColor << "File NOT Found: " << ss.str() << RESET << endl;
                exit(1);
            }
            auto res = garRead(ss.str());
            fRes += res;
            if(res.has(SD::NaN)) nid.push_back(i);            
            else if(VEMaxErr(res)>ee) {
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
        if(in>0) cout << "[ n=" << in << " ]";
        cout << ":" << endl;
        cout << "------" << endl;
        cout << fRes << endl;
        cout << "------" << endl;
        cout << VEResult(fRes) << endl << endl << endl;
        
        if(is_o) {
            ofstream ofs;
            ofs.open(arg_o, ios::out);
            if (!ofs) throw Error("failed to open final out file!");
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
                    cout << "./NIntegrate -n " << i << " -a i -i " << (ii+1) << endl;
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
    if (!ofs) throw Error("failed to open final null file!");
    ofs << "Result is Null." << endl;
    ofs.close();
}
