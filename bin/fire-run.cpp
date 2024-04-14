#include "HepLib.h"
#include <getopt.h>

using namespace HepLib;

int main(int argc, char **argv) {
    
    string arg_host;
    string arg_port;
    string arg_exe;
    string arg_options;
    string arg_suffix;
    bool arg_nolog = false;
    
    { // handle long options
        int opt;
        int option_index = 0;
        const char *optstring = "";
        static struct option long_options[] = {
            { "h", required_argument, 0, 0 },
            { "host", required_argument, 0, 0 },
            { "p", required_argument, 0, 0 },
            { "port", required_argument, 0, 0 },
            { "o", required_argument, 0, 0 },
            { "options", required_argument, 0, 0 },
            { "e", required_argument, 0, 0 },
            { "exe", required_argument, 0, 0 },
            { "s", required_argument, 0, 0},
            { "suffix", required_argument, 0, 0},
            { "L", no_argument, 0, 0},
            { "nolog", no_argument, 0, 0},
            { "help", no_argument, 0, 0},
            {0, 0, 0, 0}
        };
        while((opt =getopt_long_only(argc, argv, optstring, long_options, &option_index))!= -1) {
            switch (option_index) {
                case 0:
                case 1:
                    arg_host = optarg;
                    break;
                case 2:
                case 3:
                    arg_port = optarg;
                    break;
                case 4:
                case 5:
                    arg_options = optarg;
                    break;
                case 6:
                case 7:
                    arg_exe = optarg;
                    break;
                case 8:
                case 9:
                    arg_suffix = optarg;
                    break;
                case 10:
                case 11:
                    arg_nolog = true;
                    break;
                default:
                    cout << "fire-run, a wrapper to run FIRE" << endl;
                    cout << "Supported Options:" << endl;
                    cout << "  --host/-h: server host/ip." << endl;
                    cout << "  --port/-p: server port." << endl;
                    cout << "  --exe/-e: path to FIRE execute." << endl;
                    cout << "  --options/-o: options passed to FIRE command." << endl;
                    cout << "  --suffix/-s: suffix to log/db/tables filename, also passed to FIRE command." << endl;
                    cout << "  --nolog/-L: disable generating log files." << endl;
                    cout << "  --help: print this information." << endl;
                    exit(1);
            }
        }
        argc -= optind;
        argv += optind;
    }
    
    if(arg_host=="" || arg_port=="") {
        cout << "host and port must be set!" << endl;
        return 0;
    }
    
    string sip = arg_host;
    string sport = arg_port;
    string exe = arg_exe;
    string options = "";
    if(arg_options!="") options = " " + arg_options;
    string suffix = "";
    if(arg_suffix!="") {
        options += " -suffix " + arg_suffix;
        suffix = "-" + arg_suffix;
    }
    bool nolog = arg_nolog;
        
    while(true) {
          
        std::string si = Server::Next(sip, sport);
        if(!file_exists(si+".start")) continue;
        if(file_exists(si+suffix+".log")) continue;
        if(file_exists(si+suffix+".tables")) continue;
        system(("touch "+si+suffix+".log").c_str());
        system(("echo $HOSTNAME > "+si+suffix+".log").c_str());
        
        if(!file_exists(exe)) {
            cout << "exe: " << exe << endl;
            cout << "FIRE NOT FOUND on HOST: " << RunOS("hostname") << endl;
            break;
        }
        
        if(nolog) system((exe+" -c "+si+options).c_str());
        else system((exe+" -c "+si+options+" >> "+si+suffix+".log").c_str());
        system(("dbdir=$(cat "+si+".config | grep database | sed s/'#database '//);test -d $dbdir"+suffix+" && rm -rf $dbdir"+suffix).c_str());
        
        if(file_exists(si+suffix+".log")) remove((si+suffix+".log").c_str());
    }

    return 0;
}
