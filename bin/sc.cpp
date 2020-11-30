#include "BASIC.h"
#include <netdb.h>
#include <getopt.h>

using namespace std;
using namespace HepLib;
#define MAXSIZE 4096

int main(int argc, char** argv) {
    
    string arg_t = "-1";
    string arg_p = "8899";
    string arg_s = "localhost";
    string arg_c = "echo [i]";
    string arg_r = "3";
    
    // handle long options
    int opt;
    int digit_optind = 0;
    int option_index = 0;
    const char *string = "";
    static struct option long_options[] = {
        { "total", required_argument, NULL, 't' },
        { "port", required_argument, NULL, 'p' },
        { "server", required_argument, NULL, 's' },
        { "command", required_argument, NULL, 'c' },
        { "round", required_argument, NULL, 'r'},
        { "help", no_argument, NULL, 'h' },
        {NULL, 0, NULL, 0}
    };
    while((opt =getopt_long_only(argc,argv,string,long_options,&option_index))!= -1) {
        switch (opt) {
            case 't': arg_t = optarg; break;
            case 'p': arg_p = optarg; break;
            case 's': arg_s = optarg; break;
            case 'c': arg_c = optarg; break;
            case 'r': arg_r = optarg; break;
            default:
                cout << "A simple Server/Client, bypass with [i].log" << endl;
                cout << "Supported Options:" << endl;
                cout << "  --total: total elements @server." << endl;
                cout << "  --port: server port @server/@client." << endl;
                cout << "  --server: server ip or hostname @client." << endl;
                cout << "  --command: command with [i] replaced @client." << endl;
                cout << "  --round: round to be cycled @server." << endl;
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
    int port = stoi(arg_p);
    int total = stoi(arg_t);
    int round = stoi(arg_r);
    
    if(total > 0) {
        int socket_fd;
        struct sockaddr_in servaddr;  
        char buff[MAXSIZE];  
        int n;  
        
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
                    
                    std::string data = to_string(current);
                    if(send(connect_fd, data.c_str(),data.length(),0) == -1) perror("send error");
                    close(connect_fd);
                }
                
                // check .exit
                auto pid = getpid();
                if(file_exists(to_string(pid)+".exit")) {
                    cout << "Exit @ " << now() << endl << endl;
                    close(socket_fd);
                    exit(0);
                }
            }
            cout << endl;
        }
        
        cout << "Finished @ " << now() << endl << endl;
        close(socket_fd);
        exit(0);
    } else {
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
            
            std::string cmd = arg_c;
            string_replace_all(cmd, "[i]", data);
            system(cmd.c_str());
        }
        exit(0); 
    }
    
    return 0;
}
