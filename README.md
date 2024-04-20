### **How** to Install `HepLib`
- Download and extract the `all-in-one` archive
```bash
wget https://heplib.github.io/download/Install.tar.gz
tar xfv Install.tar.gz
```

- Install using `install.sh`
```bash
cd Install
INSTALL_PATH=<Path to Install> jn=8 ./install.sh
# INSTALL_PATH=<Path to Install> jn=8 ./install-M1.sh # Apple Silicon Chip
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



### Copyright 


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/.
