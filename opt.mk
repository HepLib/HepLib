osName := $(shell uname -s)


ifeq ($(osName),Darwin)
opt = -w -std=gnu++14 -fPIC -I ../include -L .. -L /usr/local/feng/lib -I /usr/local/feng/include
else
opt = -w -std=gnu++14 -fPIC -I ../include -L .. -I. -L $$HOME/usr/local/lib -L $$HOME/usr/local/lib64 -I $$HOME/usr/local/include -L $$HOME/glibc/lib -L /WORK/app/glibc/2.14/lib -L $$HOME/usr/lib -I $$HOME/usr/include
endif



ifeq ($(osName),Darwin)
prefix = /usr/local/feng
else
prefix = $$HOME/usr/local
endif



ifeq ($(osName),Darwin)
WSTP_PATH = /Applications/Mathematica.app/Contents/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions
opt := $(opt) -L $(WSTP_PATH) -lc++ -lWSTPi4 -framework Foundation
else
WSTP_PATH = /usr/local/Wolfram/Mathematica/12.0/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions
opt := $(opt) -L $(WSTP_PATH) -lm -lpthread -lrt -lstdc++ -ldl -luuid -lWSTP64i4
endif
#disable WSTP
#WSTP_PATH = 



