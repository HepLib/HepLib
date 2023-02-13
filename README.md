### Install `HepLib` and its related requirements

1. **Download** `install.sh`
```
wget --no-check-certificate -qO install.sh 'https://heplib.github.io/install.sh'
chmod a+x install.sh
```

2. **Run** `install.sh` to install `HepLib`
```
INSTALL_PATH=/path/heplib_to_be_installed jn=8 ./install.sh
```

### **Try** a simple example
1. Prepare a `C++` file named `trace.cpp` with the following content
```cpp
#include "HepLib.h"
using namespace HepLib;
int main(int argc, char** argv) {
    Index mu("mu"), nu("nu");
    Vector p1("p1"), p2("p2");
    Symbol m("m");
    //note GAS(1) in gline, corresponds to the identity matrix
    ex gline = GAS(p1)*GAS(mu)*(GAS(p2)+m*GAS(1))*GAS(mu);
    ex trace = form(TR(gline));
    cout << trace << endl;
    return 0;
}
```

2. Compile `trace.cpp` using `heplib++` from `HepLib`
```bash
INSTALL_PATH/bin/heplib++ -o trace trace.cpp
```

3. Run `trace` to get the trace for ${\rm Tr}[\gamma\cdot p_1 \ \gamma^\mu\ (\gamma\cdot p_2+m)\ \gamma_\mu]$
```bash
./trace
# -4*d*p2.p1+8*p2.p1
``` 
