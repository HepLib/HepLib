osName := $(shell uname -s)

ifeq ($(osName),Darwin)
opt = -w -std=gnu++14 -I ../include -fPIC -L /usr/local/feng/lib -I /usr/local/feng/include
else
opt = -w -std=gnu++14 -I ../include -fPIC -I. -L $$HOME/usr/local/lib -I $$HOME/usr/local/include
endif

ifeq ($(osName),Darwin)
prefix = /usr/local/feng
else
prefix = $$HOME/usr/local
endif





