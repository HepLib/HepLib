#include "BASIC.h"
#include <netdb.h>
#include <getopt.h>

#define MAXSIZE 4096

namespace HepLib {

    string cpp = "c++";
    
    exvector ParRun(int total, string dir) {
        
        int socket_fd;
        struct sockaddr_in servaddr;  
        char buff[MAXSIZE];  
        
        int port = 8890;
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
                    
                    string data = dir+"/"+to_string(current);
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
