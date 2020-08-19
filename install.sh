#!/bin/bash

if [ -z $prefix ]; then
    export prefix=/usr/local
fi

export jn=8
export install_gmp='yes'
export install_mpfr='yes'
export install_cln='yes'
export install_ginac='yes'
export install_cuba='yes'
export install_minuit2='yes'
export install_qhull='yes'
export install_heplib='yes'

export CWD=$PWD

# install GMP
if [ $install_gmp == 'yes' ]; then
    export pkg="gmp-6.2.0"
    if [ ! -f $pkg.tar.gz ]; then
        wget --no-check-certificate https://gmplib.org/download/gmp/$pkg.tar.gz
    fi
    tar zxf $pkg.tar.gz
    cd $pkg
    ./configure --prefix=$prefix
    make -j $jn
    make install
    cd $CWD
    rm -rf $pkg
fi

# install MPFR
if [ $install_mpfr == 'yes' ]; then
    export pkg="mpfr-4.1.0"
    if [ ! -f $pkg.tar.gz ]; then
        wget --no-check-certificate https://www.mpfr.org/mpfr-current/$pkg.tar.gz
    fi
    tar zxf $pkg.tar.gz
    cd $pkg
    if [ $install_gmp == 'yes' ]; then
        ./configure --prefix=$prefix --with-gmp=$prefix --enable-float128 --enable-thread-safe
    else
        ./configure --prefix=$prefix --enable-float128 --enable-thread-safe
    fi
    make -j $jn
    make install
    cd $CWD
    rm -rf $pkg
fi

# install CLN
if [ $install_cln == 'yes' ]; then
    export pkg="cln-1.3.6"
    if [ ! -f $pkg.tar.bz2 ]; then
        wget --no-check-certificate https://www.ginac.de/CLN/$pkg.tar.bz2
    fi
    tar jxf $pkg.tar.bz2
    cd $pkg
    if [ $install_gmp == 'yes' ]; then
        ./configure --prefix=$prefix --with-gmp=$prefix
    else
        ./configure --prefix=$prefix
    fi
    make -j $jn
    make install
    cd $CWD
    rm -rf $pkg
fi

# install GiNaC
if [ $install_ginac == 'yes' ]; then
    export pkg="ginac-1.7.11"
    if [ ! -f $pkg.tar.bz2 ]; then
        wget --no-check-certificate https://www.ginac.de/$pkg.tar.bz2
    fi
    tar jxf $pkg.tar.bz2
    cd $pkg
    ./configure --prefix=$prefix PKG_CONFIG_PATH=$prefix/lib/pkgconfig
    make -j $jn
    make install
    cd $CWD
    rm -rf $pkg
fi

# install CUBA
if [ $install_cuba == 'yes' ]; then
    export pkg="Cuba-4.2"
    if [ ! -f $pkg.tar.gz ]; then
        wget --no-check-certificate http://www.feynarts.de/cuba/$pkg.tar.gz
    fi
    tar zxf $pkg.tar.gz
    cd $pkg
    ./configure --prefix=$prefix CFLAGS="-fPIC" CXXFLAGS="-fPIC"
    make
    make install
    cd $CWD
    rm -rf $pkg
fi

# install ROOT::Minuit2
if [ $install_minuit2 == 'yes' ]; then
    export pkg="Minuit2-5.34.14"
    if [ ! -f $pkg.tar.gz ]; then
        wget --no-check-certificate http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/$pkg.tar.gz
    fi
    tar zxf $pkg.tar.gz
    cd $pkg
    ./configure --prefix=$prefix
    make -j $jn
    make install
    cd $CWD
    rm -rf $pkg
fi

# install QHull
if [ $install_qhull == 'yes' ]; then
    export pkg="qhull-2020.2"
    if [ ! -f $pkg.zip ]; then
        wget --no-check-certificate http://www.qhull.org/download/$pkg.zip
    fi
    unzip $pkg.zip
    cd $pkg
    cp Makefile Makefile.bak
    cat Makefile.bak | sed "s/\/usr\/local/\$\$prefix/g" > Makefile
    make
    make install
    cd $CWD
    rm -rf $pkg
fi

# install HepLib
if [ $install_heplib == 'yes' ]; then
    echo "install HepLib"
fi
