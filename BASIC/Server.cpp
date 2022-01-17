#include "BASIC.h"
#include <netdb.h>
#include <sstream>

namespace HepLib {

    string Server::Next(string sip, string sport) {
        static const int MAXSIZE = 4096;
        auto port = stoi(sport);
        
        int sockfd, n,rec_len;  
        char recvline[MAXSIZE], sendline[MAXSIZE];  
        char buf[MAXSIZE];  
        struct sockaddr_in servaddr;  
            
        if( (sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) { 
            ostringstream oss;
            oss << "create socket error(" << errno << "): " << strerror(errno);  
            throw Error(oss.str());
        }
      
        memset(&servaddr, 0, sizeof(servaddr));  
        servaddr.sin_family = AF_INET;  
        servaddr.sin_port = htons(port);  
        struct hostent *hext;
        if ( (hext = gethostbyname(sip.c_str())) == NULL ) {
            ostringstream oss;
            oss << "gethostbyname error for " << sip << endl;  
            throw Error(oss.str());
        }
        memcpy(&servaddr.sin_addr, hext->h_addr_list[0], hext->h_length);
       
        if( connect(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr)) < 0) {  
            ostringstream oss;
            oss << "connect error(" << errno << "): " << strerror(errno) << endl;  
            throw Error(oss.str());
        }  
         
        if((rec_len = recv(sockfd, buf, MAXSIZE,0)) == -1) {  
            ostringstream oss;
            oss << "recv error(" << errno << "): " << strerror(errno) << endl;  
            throw Error(oss.str());
        }  
        
        buf[rec_len]  = '\0';  
        string ret = buf;
        close(sockfd);  
        return ret;
    }
    
    void Server::Start() {
        
        int port = Port;
        auto total = Total;
        auto round = Round;
        
        int socket_fd;
        struct sockaddr_in servaddr;  
                
        if( (socket_fd = socket(AF_INET, SOCK_STREAM, 0)) == -1 ) {  
            throw Error(string("Server: ") + strerror(errno));  
        } 
        
        memset(&servaddr, 0, sizeof(servaddr));  
        servaddr.sin_family = AF_INET;  
        servaddr.sin_addr.s_addr = htonl(INADDR_ANY);  
        servaddr.sin_port = (port>0 ? htons(port) : 0); 
        
        if( bind(socket_fd, (struct sockaddr*)&servaddr, sizeof(servaddr)) == -1) {  
            throw Error(string("Server: ") + strerror(errno));  
        } 
        if(port<=0) {
            struct sockaddr_in connAddr;
            socklen_t len = sizeof(connAddr);
            if(getsockname(socket_fd, (sockaddr*)&connAddr, &len) !=0) throw Error("Server: retrived Port failed!");
            Port = port = ntohs(connAddr.sin_port);
        }
        
        if( listen(socket_fd, 10) == -1) throw Error(string("Server: ") + strerror(errno));  
        map<string,ex> dict;
        char hostname[1024];
        gethostname(hostname, 1024);
        dict["port"] = Port;
        dict["DL"] = Symbol(DL);
        dict["FUNC"] = Symbol(FUNC);
        garWrite("PRC.gar", dict);
        
        if(Verbose>0) cout << endl << "Started @ " << now() << endl;        
        for(int r=0; r<round; r++) {
            for(int c=0; c<total; c++) {
                auto current = c;
                if(Verbose>1) {
                    cout << "\r                                     \r";
                    cout << "  Server[" << port << "]: " << current << " / " << (total-1) << " @ " << now(false) << flush;
                }
                
                string skip = Skip;
                string_replace_all(skip, "[ID]", to_string(current));
                if(file_exists(skip)) continue;
                
                int connect_fd;
                if( (connect_fd = accept(socket_fd, (struct sockaddr*)NULL, NULL)) == -1) {
                    cout << "Server: " << strerror(errno) << endl;
                    continue;
                }
                struct linger so_linger;
                so_linger.l_onoff = 1; 
                so_linger.l_linger = 0; 
                setsockopt(connect_fd, SOL_SOCKET, SO_LINGER, &so_linger, sizeof so_linger); 
                
                string data = to_string(current); // data = ID
                if(send(connect_fd, data.c_str(), data.length(),0) == -1) {
                    if(Verbose>1) cout << "Server: " << strerror(errno) << endl;
                }
                close(connect_fd);
            }
            if(Verbose>1) cout << endl;
        }
        if(Verbose>0) cout << "Finished @ " << now() << endl << endl;
        close(socket_fd);
    }
    
}
