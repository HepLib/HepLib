#include "BASIC.h"
#include <netdb.h>
#include <getopt.h>
#include <dlfcn.h>

using namespace std;
using namespace HepLib;
#define MAXSIZE 65536

int main(int argc, char** argv) {

    if(!file_exists("PRC.gar")) {
        cout << "No PRC.gar Found!" << endl;
        return 0;
    }
    
    string sip = "localhost";
    int port;
    string sofn, func;
    
    if(true) { // scope for ar/dict
        archive ar;
        ifstream in("PRC.gar");
        in >> ar;
        in.close();
        
        map<string, ex> dict;
        for(int i=0; i<ar.num_expressions(); i++) {
            string name;
            ex res = ar.unarchive_ex(GiNaC_archive_Symbols, name, i);
            dict[name] = res;
        }
        
        port = stoi(ex_to<Symbol>(dict["port"]).get_name());
        sofn = ex_to<Symbol>(dict["so"]).get_name();
        func = ex_to<Symbol>(dict["fun"]).get_name();
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
            close(sockfd); 
            cout << "gethostbyname error for " << sip << endl;  
            exit(1);
        }
        memcpy(&servaddr.sin_addr, hext->h_addr_list[0], hext->h_length);
       
        if( connect(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr)) < 0) { 
            close(sockfd); 
            cout << "connect error(" << errno << "): " << strerror(errno) << endl;  
            if(errno != 54) sleep(3);
            continue;
        }
        
        if((rec_len = recv(sockfd, buf, MAXSIZE,0)) == -1) {
            close(sockfd); 
            cout << "recv error(" << errno << "): " << strerror(errno) << endl;  
            sleep(3);
            continue;
        }
        
        buf[rec_len]  = '\0';  
        std::string id = buf;
        close(sockfd);  
        
        void* module = dlopen(sofn.c_str(), RTLD_NOW);
        if(module == nullptr) {
            module = dlopen(sofn.c_str(), RTLD_NOW);
            if(module == nullptr) {
                cout << "dlerror(): " << dlerror() << endl;
                throw Error("Integrates: could not open main module!");
            }
        }
        
        auto run = (RUN)dlsym(module, func.c_str());
        if(run==NULL) {
            cout << "dlerror(): " << dlerror() << endl;
            throw Error("PRC: dlsym NOT Found!");
        }
        run(id);
        
    } 
    
    return 0;
}
