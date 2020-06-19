#include "Basic.h"
#include <netdb.h>

using namespace std;
using namespace HepLib;
#define MAXSIZE 4096

int main(int argc, char** argv) {
    
    string arg_m = "s";
    string arg_p = "8899";
    string arg_s = "localhost";
    string arg_c = "echo [i]";
    string cm_path = "cm";
    string sd_path = "SD";
    bool ignore = false;
    
    // handle options
    for (int opt; (opt = getopt(argc, argv, "m:p:s:c:i")) != -1;) {
        switch (opt) {
            case 'm': arg_m = optarg; break;
            case 'p': arg_p = optarg; break;
            case 's': arg_s = optarg; break;
            case 'c': arg_c = optarg; break;
            case 'i': ignore=true; break;
            default:
                printf("supported options: -m A -p P -s S -c C.\n");
                printf("M: s(server), c(client).\n");
                printf("P: server port.\n");
                printf("S: server ip or hostname.\n");
                printf("C: command on client, [i] will be replaced.\n");
                printf("-i: to ignore file_exist check for res/null.\n");
                exit(1);
        }
    }
    argc -= optind;
    argv += optind;
    int port = stoi(arg_p);
    
    if(arg_m == "s") {
        int socket_fd, connect_fd;  
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
        
        ostringstream cmd;
        cmd << "ls " << cm_path << "|wc -l";
        auto owc = RunOS(cmd.str().c_str());
        int total = stoi(owc);
        int current = 1;
         
        while(true) {  
            if(current>total) current=1;
            cout << "\r                               \r";
            cout << "Server: " << current << " / " << total << " @ " << now(false) << flush;
            if( (connect_fd = accept(socket_fd, (struct sockaddr*)NULL, NULL)) == -1) {  
                cout << "accept socket error(" << errno << "): " << strerror(errno) << endl;  
                continue;  
            }  
            
            bool second_run = false;
            while(file_exists(to_string(current)+".log") || (!ignore && file_exists(sd_path + "/" + to_string(current) + ".res.gar") || file_exists(sd_path + "/" + to_string(current) + ".null"))) {
                if(current>=total) {
                    current=0;
                    if(second_run) break;
                    second_run = true;
                }
                current++;
            }
            
            #pragma omp parallel 
            #pragma omp for firstprivate(connect_fd) nowait
            for(int i=0; i<2; i++) {
                if(omp_get_thread_num()!=0) {
                    string data = to_string(current);
                    if(send(connect_fd, data.c_str(),data.length(),0) == -1) perror("send error");  
                    close(connect_fd);  
                }
            } 
            
            wait(NULL); 
            close(connect_fd);  
            if(current==0) break;
            current++;
        }
        
        cout << " @ " << now(false) << endl;
        close(socket_fd);
        exit(0);
    } else {
        string sip = arg_s;
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
            string data = buf;  
            close(sockfd);  
            
            if(data=="0") break;
            
            string cmd = arg_c;
            cmd += " > [i].log;rm [i].log";
            string_replace_all(cmd, "[i]", data);
            system(cmd.c_str());
        }
        exit(0); 
    }
    

    return 0;
}
