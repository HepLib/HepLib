CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)
SET(CMAKE_Fortran_COMPILER gfortran)

#-----------------------------------------------
# Configuration Here
#-----------------------------------------------
# set the install prefix
IF(APPLE)
    SET(PREFIX_DIR "/usr/local/feng")
ELSE()
    SET(PREFIX_DIR "$ENV{HOME}/usr/local")
ENDIF()

# set the "-I" information for gcc/g++, ignored while not exists
SET(INC_DIRS
    "/usr/local/feng/include"
    "/dev/shm/usr/local/include"
    "/dev/shm/usr/include"
    "$ENV{HOME}/usr/local/include"
    "$ENV{HOME}/usr/include"
)

# set the "-L" information for gcc/g++, ignored while not exists
SET(LIB_DIRS
    "/usr/local/feng/lib"
    "/dev/shm/usr/local/lib"
    "/dev/shm/usr/lib"
    "/dev/shm/glibc/lib"
    "$ENV{HOME}/usr/local/lib"
    "$ENV{HOME}/usr/lib"
    "$ENV{HOME}/glibc/lib"
)

# Other things:
# -MPFR: need to be compiled with option: --enable-float128
# -Cuba: need to be compile with CFLAGS="-fPIC" CXXFLAGS="-fPIC"
#
#-----------------------------------------------
# Configure End
#-----------------------------------------------

INCLUDE_DIRECTORIES(BEFORE include)
FOREACH(dir ${INC_DIRS})
IF(EXISTS ${dir})
    INCLUDE_DIRECTORIES(${dir})
ENDIF()
ENDFOREACH()

FOREACH(dir ${LIB_DIRS})
IF(EXISTS ${dir})
    LINK_DIRECTORIES(${dir})
ENDIF()
ENDFOREACH()

LINK_LIBRARIES(
    ginac cln Minuit2 mpfr qhullstatic quadmath gomp cubaq HepLib
)

# OpenMP
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu++14")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")

PROJECT(HepLib_Test VERSION 1.0 LANGUAGES C CXX Fortran)

AUX_SOURCE_DIRECTORY(. SRCS)
FOREACH(cpp ${SRCS})
    GET_FILENAME_COMPONENT(exe ${cpp} NAME_WE)
    ADD_EXECUTABLE(${exe} ${cpp})
ENDFOREACH()

