// A Template for Sector Decomposition

#include "HepLib.h"

using namespace HepLib;
using namespace SD;

string cm_path = "cm"; // input directory
string SD_path = "SD"; // save directory
int epN = 0;
int epsN = 0;

void ExportNull(const string & prefix);
int CheckNull(const string & prefix);
void Prepare(int idx) {

    XIntegrand xint;
    xint.Functions = lst{ 1, x(0)+x(1) };
    xint.Exponents = lst{ 1, -1+ep };
    
    SecDec work;
    work.Initialize(xint);

    work.Normalizes();
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
        ifn << SD_path << "/" << idx;
        ExportNull(ifn.str());
        return;
    }
    
}

void Contour(int idx) {
    static SecDec work;
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
    work.Contours(ikey.str());
}

ex Integrate(int idx, int ni = -1) {
    SecDec work;
    work.epN = epN;
    work.epsN = epsN;
    
    work.use_ErrMin = false;
    ErrMin::err_min = 1E-3;
    ErrMin::MaxRND = 20;
    ErrMin::hjRHO = 0.5;
            
    work.EpsAbs = 1E-3;
    work.RunPTS = 1000000;
    work.RunMAX = 5;
    work.LambdaSplit = 5;
    work.TryPTS = 1000000;
    work.CTryLeft = 1;
    work.CTryRight = 1;
    work.CTry = 1;
    //work.ReIm = 1;
    work.CTryRightRatio = 2;
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
        
    auto intor = new HCubature();
    intor->MPXLimit = 5E-12;
    intor->QXLimit = 5E-8;
    intor->use_last = false;
    intor->MPXDim = 0;
    work.Integrator = intor;
    
    work.Integrates(ikey.str(), "", ni);
    
    return work.ResultError;
}

ex ReIntegrate(int idx, qREAL err) {
    SecDec work;
    work.epN = epN;
    work.epsN = epsN;
    
    work.use_ErrMin = false;
    ErrMin::err_min = 1E-3;
    ErrMin::MaxRND = 20;
    ErrMin::hjRHO = 0.5;
            
    work.EpsAbs = 1E-4;
    work.RunPTS = 1000000;
    work.RunMAX = 10;
    work.LambdaSplit = 5;
    work.TryPTS = 1000000;
    work.CTryLeft = 1;
    work.CTryRight = 1;
    work.CTry = 1;
    //work.ReIm = 1;
    work.CTryRightRatio = 2;
    
    ostringstream ikey;
    ikey << SD_path << "/" << idx;
        
    auto intor = new HCubature();
    intor->MPXLimit = 5E-5;
    intor->QXLimit = 5E-3;
    intor->use_last = true;
    work.Integrator = intor;
    
    work.ReIntegrates(ikey.str(), "", err);
    
    return work.ResultError;
}

int main(int argc, char** argv) {

    string exe = argv[0];
    
    string arg_a = "";
    string arg_n = "";
    string arg_i = "";
    string arg_v = "0";
    string arg_e = "-1";
    bool is_o = false;
    string arg_o = "";
    
    // handle options
    for (int opt; (opt = getopt(argc, argv, "a:n:i:o:v:e:dz:")) != -1;) {
        switch (opt) {
            case 'a': arg_a = optarg; break;
            case 'n': arg_n = optarg; break;
            case 'i': arg_i = optarg; break;
            case 'o': arg_o = optarg; is_o = true; break;
            case 'v': arg_v = optarg; break;
            case 'e': arg_e = optarg; break;
            case 'd': Debug = true; break;
            default:
                cout << "supported options: -a A -n N -i I -o O -v V -e E" << endl;
                cout << "A: can be p/c/i/r OR e/em/er" << endl;
                cout << "    p: only run Prepare precedure." << endl;
                cout << "    c: only run Contour precedure." << endl;
                cout << "    i: only run Integrate precedure." << endl;
                cout << "    r: only show the result." << endl;
                cout << "    e: to show the index with (err > E)." << endl;
                cout << "    em: only show the index with MAX error." << endl;
                cout << "    er: ReIntegrate for (error > E)." << endl;
                cout << "N: only apply A on prefix N." << endl;
                cout << "I: only apply A(=i) on part I in prefix N." << endl;
                cout << "O: output filename." << endl;
                cout << "V: verbose level." << endl;
                cout << "E: the Error limit." << endl;
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
    int nmi = 1;
    if(dir_exists(cm_path)) {
        cmd << "ls " << cm_path << "|wc -l";
        auto owc = RunOS(cmd.str().c_str());
        stringstream ss;
        ss << owc;
        ss >> nmi;
    }
    
    int in = 0;
    if(arg_n!="") in = stoi(arg_n);
    int ni = -1;
    if(arg_n!="" && arg_i!="") ni = stoi(arg_i);
    Verbose = stoi(arg_v);
    numeric ee(arg_e.c_str());
    
    //-------------------------------------------------------------
    vector<int> fid;
    vector<int> nid;
    ex fRes = 0;
    for(int n=1; n<=nmi; n++) {
        if(in>0 && n!=in) continue;
        try {
        
            // call Prepare(index)
            if(arg_a=="" || arg_a=="p") {
                cout << "Preparing " << n << "/" << nmi << endl;
                Prepare(n);
                cout << endl;
            }
            
            // call Contour(index)
            if(arg_a=="" || (arg_a=="c")) {
                int chk_null = CheckNull(SD_path+"/"+to_string(n));
                if(0==chk_null) throw Error("CheckNull failed.");
                if(2==chk_null) continue;
                cout << "Contouring - " << n << "/" << nmi << endl;
                Contour(n);
                cout << endl;
            }
            
            // call Integrate(index, sub_index)
            if(arg_a=="" || (arg_a=="i")) {
                int chk_null = CheckNull(SD_path+"/"+to_string(n));
                if(0==chk_null) throw Error("CheckNull failed.");
                if(2==chk_null) continue;
                cout << "Integrating - " << n << "/" << nmi << endl;
                auto res = Integrate(n, ni);
                cout << endl;
            }
            
            // call ReIntegrate(index, err)
            if(arg_a=="er" && ee>0) {
                int chk_null = CheckNull(SD_path+"/"+to_string(n));
                if(0==chk_null) throw Error("CheckNull failed.");
                if(2==chk_null) continue;
                cout << "ReIntegrating - " << n << "/" << nmi << endl;
                auto res = ReIntegrate(n, ex2q(ee));
                cout << endl;
            }
            
            // show the result
            if(arg_a=="" || (arg_a=="r")) {
                int chk_null = CheckNull(SD_path+"/"+to_string(n));
                if(0==chk_null) throw Error("CheckNull failed.");
                if(2==chk_null) continue;
                
                stringstream ss;
                ss << SD_path << "/" << n << ".res.gar";
                if(!file_exists(ss.str())) {
                    cout << RED << "File NOT Found: " << ss.str() << RESET << endl;
                    exit(1);
                }
                auto res = garRead(ss.str());
                fRes += res;
                if(res.has(SD::NaN)) nid.push_back(n);
                else if(ee>0 && VEMaxErr(res)>ee) {
                    cout << "n=" << n << ", err_max=" << VEMaxErr(res) << endl << endl;
                }
            }
            
            // analyse result
            if((arg_a=="e") || (arg_a=="em")) {
                int chk_null = CheckNull(SD_path+"/"+to_string(n));
                if(0==chk_null) throw Error("CheckNull failed.");
                if(2==chk_null) continue;
                
                stringstream ss;
                ss << SD_path << "/" << n << ".res.gar";
                map<string, ex> gar;
                garRead(ss.str(), gar);
                auto relst = ex_to<lst>(gar["relst"]);
                int max_index;
                ex max_re = -1;
                for(int i=0; i<relst.nops(); i++) {
                    auto tmp = VEMaxErr(relst.op(i));
                    if(abs(tmp)>max_re) {
                        max_re = abs(tmp);
                        max_index = i;
                    }
                    if(ee>0 && abs(tmp)>ee) {
                        if(arg_a=="e") {
                            cout << exe << " -n " << n << " -a i -i " << (i+1) << endl;
                            cout << "# " << VEMaxErr(relst.op(i)) << endl;
                            cout << endl;
                        }
                    }
                }
                if(arg_a=="em") {
                    cout << "# [" << n << "] Max Index: " << max_index+1 << " / " << relst.nops() << endl;
                    cout << "# " << max_re << endl;
                    cout << endl << endl;
                }
            }
            
        } catch(exception& e) {
            cout << e.what() << endl;
            fid.push_back(n);
        } catch(...) {
            cout << "other uncatch error" << endl;
            fid.push_back(n);
        }
    }
    
    if(fid.size()>0) {
        cout << "Failed IDs: {";
        for(auto item : fid) cout << item << ", ";
        cout << "}" << endl;
        return 0;
    }
    if(nid.size()>0) {
        cout << "NaN IDs: {";
        for(auto item : fid) cout << item << ", ";
        cout << "}" << endl;
        return 0;
    }
    
    if(arg_a=="" || arg_a=="r") {
        fRes = VESimplify(fRes, epN, epsN);
        cout << endl << "Final Result:";
        if(in>0) cout << "[ n=" << in << " ]";
        cout << ":" << endl;
        cout << "------" << endl;
        cout << fRes << endl;
        cout << "------" << endl;
        cout << VEResult(fRes) << endl << endl << endl;
    }
    
    if(is_o) {
        ofstream ofs;
        ofs.open(arg_o, ios::out);
        if (!ofs) throw runtime_error("failed to open final out file!");
        ofs << fRes << endl;
        ofs.close();
    }
    //-------------------------------------------------------------
    
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

int CheckNull(const string & prefix) {
    stringstream ss;
    ss << prefix << ".ci.gar";
    if(!file_exists(ss.str())) {
        ostringstream ifn;
        ifn << prefix << ".null";
        if(!file_exists(ifn.str())) return 0;
        return 2;
    } else return 1;
}
