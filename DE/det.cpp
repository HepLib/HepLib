#include "functions.cpp"

int main(int argc, char *argv[]) {
    const char *x_name = "x";
    const char *ep_name = "ep";
    const char *m_path = "mat.m";
    const char *f_path = "fmat.m";
    const char *n_path = "nmat.m";
    const char *s_path = "smat.m";
    const char *t_path = "";
    
    numeric ep_ex;
    bool is_ep = false;
    
    for (int opt; (opt = getopt(argc, argv, "x:e:m:f:n:s:t:N:")) != -1;) {
        switch (opt) {
            case 'x': x_name = optarg; break;
            case 'e': ep_name = optarg; break;
            case 'm': m_path = optarg; break;
            case 'f': f_path = optarg; break;
            case 'n': n_path = optarg; break;
            case 's': s_path = optarg; break;
            case 't': t_path = optarg; break;
            case 'N': ep_ex = numeric(optarg); is_ep = true; break;
            default: return 1;
        }
    }
    argc -= optind;
    argv += optind;
        
    if((argc == 1) && !strcmp(argv[0], "fuchsify")) {
        symbol x(x_name), ep(ep_name);
        parser reader(symtab({{ep_name, ep}, {x_name, x}}));
        matrix mat = load_matrix(m_path, reader);
        if(is_ep) matrix_map_inplace(mat, [&](auto &&e) { return e.subs(ep == ep_ex); });

        cout << "---------------------------------------" << endl;
        cout << mat.rows() << " x " << mat.cols() << " matrix loaded from " << m_path << endl;
        
        auto mt = fuchsify(mat, x);
        save_matrix(f_path, mt.first);
        if(strlen(t_path)>0) save_matrix(t_path, mt.second);
        cout << "---------------------------------------" << endl;
        cout << "f-matrix saved to " << f_path << endl;
        if(strlen(t_path)>0) cout << "t-matrix saved to " << t_path << endl;
    } else if((argc == 1) && !strcmp(argv[0], "normalize")) {
        symbol x(x_name), ep(ep_name);
        parser reader(symtab({{ep_name, ep}, {x_name, x}}));
        matrix mat = load_matrix(f_path, reader);
        if(is_ep) matrix_map_inplace(mat, [&](auto &&e) { return e.subs(ep == ep_ex); });
        
        cout << "---------------------------------------" << endl;
        cout << mat.rows() << " x " << mat.cols() << " matrix loaded from " << f_path << endl;
        
        auto mt = normalize(mat, x);
        save_matrix(n_path, mt.first);
        if(strlen(t_path)>0) save_matrix(t_path, mt.second);
        cout << "---------------------------------------" << endl;
        cout << "n-matrix saved to " << n_path << endl;
        if(strlen(t_path)>0) cout << "t-matrix saved to " << t_path << endl;
    } else if ((argc == 1 || argc == 2) && !strcmp(argv[0], "shearing")) {
        symbol x(x_name), ep(ep_name);
        parser reader(symtab({{ep_name, ep}, {x_name, x}}));
        matrix mat = load_matrix(f_path, reader);
        if(is_ep) matrix_map_inplace(mat, [&](auto &&e) { return e.subs(ep == ep_ex); });
        int epN = -19790923;
        if(argc == 2) epN = stoi(argv[1]);

        cout << "---------------------------------------" << endl;
        cout << mat.rows() << " x " << mat.cols() << " matrix loaded from " << f_path << endl;

        auto mt = shearing(mat, x, ep, epN);
        save_matrix(n_path, mt.first);
        if(strlen(t_path)>0) save_matrix(t_path, mt.second);
        cout << "---------------------------------------" << endl;
        cout << "n-matrix saved to " << n_path << endl;
        if(strlen(t_path)>0) cout << "t-matrix saved to " << t_path << endl;
    } else if ((argc > 1) && !strcmp(argv[0], "dess")) {
        Digits = 100;
        symbol x(x_name), ep(ep_name);
        parser reader(symtab({{ep_name, ep}, {x_name, x}}));
        
        int xN = 10, epN = -19790923;
        ex x0 = x;
        if(argc>1) xN = stoi(argv[1]);
        if(argc>2) epN = stoi(argv[2]);
        if(argc>3) x0 = numeric(argv[3]);
            
        matrix mat = load_matrix(n_path, reader);
        
        cout << "---------------------------------------" << endl;
        cout << mat.rows() << " x " << mat.cols() << " matrix loaded from " << f_path << endl;
        
        auto smat = dess(mat, x, ep, x0, xN, epN);
        cout << "---------------------------------------" << endl;
        save_matrix(s_path, smat);
        cout << "s-matrix saved to " << s_path << endl;
        
    } else if ((argc == 2) && !strcmp(argv[0], "show")) {
        symbol x(x_name), ep(ep_name);
        parser reader(symtab({{ep_name, ep}, {x_name, x}}));
        m_path = argv[1];
        matrix mat = load_matrix(m_path, reader);
        if(is_ep) matrix_map_inplace(mat, [&](auto &&e) { return e.subs(ep == ep_ex); });

        cout << "---------------------------------------" << endl;
        cout << mat.rows() << " x " << mat.cols() << " matrix loaded from " << m_path << endl;
        cout << "---------------------------------------" << endl;
        
        int pr = prank(mat, x);
        matrix a0 = a0_matrix(mat, x, pr);
        a0 = normal(a0);
        
        cout << "matrix information:" << endl;
        cout << "  Poincare Rank: " << pr << endl;
        
        auto ev_map = eigenvalues(a0);
        vector<exvector> ev_groups;
        for(auto &kv : ev_map) {
            auto ev = kv.first;
            bool in_g = false;
            for(auto &gi : ev_groups) {
                if(in_g) break;
                for(auto &ev_i : gi) {
                    auto diff_ev = normal(ev_i-ev);
                    if(is_a<numeric>(diff_ev) && ex_to<numeric>(diff_ev).is_integer()) {
                        in_g = true;
                        gi.push_back(ev);
                        break;
                    }
                }
            }
            if(!in_g) {
                exvector gi_new;
                gi_new.push_back(ev);
                ev_groups.push_back(gi_new);
            }
        }
        
        cout << "  eigen values of A0 summary:" << endl;
        for(auto &gi : ev_groups) {
            sort(gi.begin(),gi.end(),[&](const auto &a, const auto &b){
                return normal(a-b).info(info_flags::positive);
            });
            ostringstream ostr;
            cout << "    ";
            for(auto &ev_i : gi) {
                ostr << ev_i << "[" << ev_map[ev_i] << "],  ";
            }
            string str = ostr.str();
            cout << str.substr(0,str.size()-3) << endl;
        }
        cout << "---------------------------------------" << endl;
        
        
        cout << "  jordan block of A0:" << endl;
        ostringstream ostr;
        int pn = 0;
        for(auto kv : jordan(a0).second) {
            pn++;
            ostr << kv.first << "[" << kv.second << "],  ";
            if(pn>=4) {
                string str = ostr.str();
                ostr.str("");
                ostr.clear();
                cout << "    " << str.substr(0,str.size()-3) << endl;
                pn = 0;
            }
        }
        
    } else if(false) {
        cout << "---------------------------------------" << endl;
        cout << "missing or invalid arguments!" << endl;
    } else {
    
        cout << "det fuchsify / det shearing / det dess xN epN x0" << endl;
        
    }
    
    
    cout << "\r" << flush;
    cout << "---------------------------------------" << endl;

    return 0;
}






