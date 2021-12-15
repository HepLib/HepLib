#include "BASIC.h"
#include <netdb.h>
#include <getopt.h>

#define MAXSIZE 4096

namespace HepLib {

    string cpp = "c++";
    
    exvector ParRun(lst elst, exvector input, string code_block) {
        
        auto pid = getpid();
        ostringstream cmd;
        cmd << "mkdir -p " << pid;
        system(cmd.str().c_str());
        cmd.clear();
        cmd.str("");
        
        cmd << pid << "/ParRun.cpp";
        ofstream out(cmd.str());
        cmd.clear();
        cmd.str("");
        
        out << "#include \"HepLib.h\"" << endl;
        out << "ex ParRunFunc(lst input, int idx) {" << endl;
        out << code_block << endl;
        out << "}" << endl;
        out.close();
        
        #ifdef _USE_FLOAT128
        cmd << cpp << " " << LIB_FLAGS <<  " -Wl,-rpath,. -rdynamic -fPIC -D_USE_FLOAT128 -shared -lHepLib -lquadmath -lmpfr -lgmp " << " -o " << pid << "/ParRun.so " << pid << "/ParRun.cpp";
        cmd << " -lHepLib -lquadmath -lmpfr -lgmp";
        #else
        cmd << cpp << " " << LIB_FLAGS <<  " -Wl,-rpath,. -rdynamic -fPIC -shared -lHepLib -lmpfr -lgmp " << " -o " << pid << "/ParRun.so " << pid << "/ParRun.cpp";
        cmd << " -lHepLib -lmpfr -lgmp";
        #endif
        system(cmd.str().c_str());
        cmd.clear();
        cmd.str("");
    
        int socket_fd;
        struct sockaddr_in servaddr;  
        char buff[MAXSIZE];  
        
        int port = 8890;
        int total = input.size();
        int round = 3;
        
        if( (socket_fd = socket(AF_INET, SOCK_STREAM, 0)) == -1 ) {  
            cout << "create socket error(" << errno << "): " << strerror(errno) << endl;  
            exit(1);  
        }  
        
        memset(&servaddr, 0, sizeof(servaddr));  
        servaddr.sin_family = AF_INET;  
        servaddr.sin_addr.s_addr = htonl(INADDR_ANY);  
        servaddr.sin_port = htons(port); 
        
        if( bind(socket_fd, (struct sockaddr*)&servaddr, sizeof(servaddr)) == -1) {  
            cout << "bind socket error(" << errno << "): " << strerror(errno) << endl;  
            exit(1);  
        }  
        
        if( listen(socket_fd, 10) == -1) {  
            cout << "listen socket error(" << errno << "): " << strerror(errno) << endl;  
            exit(1);  
        }
        
        cout << endl << "Started @ " << now() << endl;
        cout << "  Server Port: " << port << endl;
        
        for(int r=0; r<round; r++) {
            for(int c=0; c<total; c++) {
                auto current = c;
                cout << "\r                                     \r";
                cout << "  Server: " << current << " / " << (total-1) << " @ " << now(false) << flush;
                
                if(!file_exists(to_string(current)+".log")) {
                    int connect_fd;
                    if( (connect_fd = accept(socket_fd, (struct sockaddr*)NULL, NULL)) == -1) {
                        cout << "accept socket error(" << errno << "): " << strerror(errno) << endl;
                        continue;
                    }
                    
                    cmd.clear();
                    cmd.str("");
                    cmd << pid << " " << current;
                    string data = cmd.str();
                    if(send(connect_fd, data.c_str(),data.length(),0) == -1) perror("send error");
                    close(connect_fd);
                }
                
            }
            cout << endl;
        }
        
        cout << "Finished @ " << now() << endl << endl;
        close(socket_fd);
        exit(0);
        
    }
    
}
