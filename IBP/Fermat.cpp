#include "IBP.h"

namespace HepLib {

void RunFermat(char const * ifn) {

    int P2C[2], C2P[2];
    pipe(P2C);
    pipe(C2P);
    
    auto pid = fork();
    if (pid == 0) { // child process
        close(P2C[1]);
        close(C2P[0]);
        dup2(C2P[1], 1);
        close(C2P[1]);
        dup2(P2C[0], 0);
        close(P2C[0]);
        system("./ferm64/fer64");
        exit(0);
    } else { // parent process
        close(P2C[0]);
        close(C2P[1]);
        auto to=fdopen(P2C[1],"w");
        auto from=fdopen(C2P[0],"r");
        fprintf(to, "&(R=%s);\n", ifn);
        fflush(to);
        
        char lc = ' ';
        bool first = true;
        while(true) {
            auto c = fgetc(from);
            //putchar(c);
            if(c==EOF) break;
            if(lc=='\n' && c=='>') {
                if(!first) break;
                first = false;
            }
            lc = c;
        }
        fclose(from);
        fclose(to);
        kill(pid,SIGKILL);
        
        close(P2C[1]);
        close(C2P[0]);
    }
}


}
