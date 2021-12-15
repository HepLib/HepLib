#include "BASIC.h"
#include <netdb.h>
#include <getopt.h>

using namespace std;
using namespace HepLib;
#define MAXSIZE 4096

int main(int argc, char** argv) {
    
    string arg_p = "8899";
    string arg_s = "localhost";
    
    // handle long options
    int opt;
    int digit_optind = 0;
    int option_index = 0;
    const char *string = "";
    static struct option long_options[] = {
        { "port", required_argument, NULL, 'p' },
        { "server", required_argument, NULL, 's' },
        { "help", no_argument, NULL, 'h' },
        {NULL, 0, NULL, 0}
    };
    while((opt =getopt_long_only(argc,argv,string,long_options,&option_index))!= -1) {
        switch (opt) {
            case 'p': arg_p = optarg; break;
            case 's': arg_s = optarg; break;
            default:
                cout << "A ParRun Client" << endl;
                cout << "Supported Options:" << endl;
                cout << "  --port: server port @server/@client." << endl;
                cout << "  --server: server ip or hostname @client." << endl;
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
    int port = stoi(arg_p);
    
    
    std::string sip = arg_s;
    if(sip.length()<1) {
        cout << "server ip or hostname is required for client mode." << endl;
        exit(1);
    }
    
    while(true) {
        int sockfd, n,rec_len;  
        char recvline[MAXSIZE], sendline[MAXSIZE];  
        char buf[MAXSIZE];  
        struct sockaddr_in servaddr;  
            
        if( (sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {  
            cout << "create socket error(" << errno << "): " << strerror(errno) << endl;  
            exit(1);  
        }  
      
        memset(&servaddr, 0, sizeof(servaddr));  
        servaddr.sin_family = AF_INET;  
        servaddr.sin_port = htons(port);  
        struct hostent *hext;
        if ( (hext = gethostbyname(sip.c_str())) == NULL ) {
            cout << "gethostbyname error for " << sip << endl;  
            exit(1);
        }
        memcpy(&servaddr.sin_addr, hext->h_addr_list[0], hext->h_length);
       
        if( connect(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr)) < 0) {  
            cout << "connect error(" << errno << "): " << strerror(errno) << endl;  
            exit(1);  
        }  
         
        if((rec_len = recv(sockfd, buf, MAXSIZE,0)) == -1) {  
           perror("recv error: ");  
           exit(1);  
        }  
        
        buf[rec_len]  = '\0';  
        std::string data = buf;
        close(sockfd);  
        
        if(file_exists(data+".log")) continue;
        // check .exit
        auto pid = getpid();
        if(file_exists(to_string(pid)+".exit")) exit(0);
        
//        std::string cmd = arg_c;
//        string_replace_all(cmd, "[i]", data);
//        system(cmd.c_str());
        
//        void* main_module = dlopen(sofn.str().c_str(), RTLD_NOW);
//        if(main_module == nullptr) {
//            main_module = dlopen(("./"+sofn.str()).c_str(), RTLD_NOW);
//            if(main_module == nullptr) {
//                cout << "dlerror(): " << dlerror() << endl;
//                throw Error("Integrates: could not open main module!");
//            }
//        }
        
//        fp = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
//        if(fp==NULL) {
//            cout << "dlerror(): " << dlerror() << endl;
//            throw Error("Integrates: fp==NULL");
//        }
            
    }    
    
    return 0;
}
