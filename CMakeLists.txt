CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)
SET(CMAKE_Fortran_COMPILER gfortran)

#-----------------------------------------------
# Configure
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
    "$ENV{HOME}/usr/local/include"
    "$ENV{HOME}/usr/include"
)

# set the "-L" information for gcc/g++, ignored while not exists
SET(LIB_DIRS
    "/usr/local/feng/lib"
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

SET(CMAKE_INSTALL_PREFIX "${PREFIX_DIR}" CACHE PATH "install prefix" FORCE)

CONFIGURE_FILE(include/Init.cpp "${CMAKE_CURRENT_BINARY_DIR}/OTHER/Init.cpp" @ONLY)
CONFIGURE_FILE(include/HepLib.pc "${CMAKE_CURRENT_BINARY_DIR}/OTHER/HepLib.pc" @ONLY)

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
    ginac cln Minuit2 mpfr qhullstatic quadmath gomp cubaq
)

#-----------------------------------------------
PROJECT(HepLib VERSION 1.0 LANGUAGES C CXX Fortran)
AUX_SOURCE_DIRECTORY(Basic BASIC_SRCS)
AUX_SOURCE_DIRECTORY(SD SD_SRCS)
AUX_SOURCE_DIRECTORY(FC FC_SRCS)
AUX_SOURCE_DIRECTORY(IBP IBP_SRCS)
AUX_SOURCE_DIRECTORY("${CMAKE_CURRENT_BINARY_DIR}/OTHER" OTHER)
#-----------------------------------------------
# check WSTP
#-----------------------------------------------
IF(APPLE)
    FIND_LIBRARY(HAS_WSTP WSTPi4 ${LINK_DIRECTORIES})
    IF(HAS_WSTP)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lc++ -lWSTPi4 -framework Foundation")
    ELSE()
        STRING(REPLACE "Basic/WSKernel.cpp" "" BASIC_SRCS "${BASIC_SRCS}")
    ENDIF(HAS_WSTP)
ELSE()
    FIND_LIBRARY(HAS_WSTP WSTP64i4 ${LINK_DIRECTORIES})
    IF(HAS_WSTP)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -lpthread -lrt -lstdc++ -ldl -luuid -lWSTP64i4")
    ELSE()
        STRING(REPLACE "Basic/WSKernel.cpp" "" BASIC_SRCS "${BASIC_SRCS}")
    ENDIF(HAS_WSTP)
ENDIF(APPLE)

#-----------------------------------------------
# install library HepLib
#-----------------------------------------------
ADD_LIBRARY(HepLib SHARED ${BASIC_SRCS} ${SD_SRCS} ${FC_SRCS} ${IBP_SRCS} ${OTHER})
INSTALL(TARGETS HepLib DESTINATION lib)

#-----------------------------------------------
# install include/*.h
#-----------------------------------------------
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.h")

#-----------------------------------------------
# install bin/*
#-----------------------------------------------
FILE(GLOB cpp_files "bin/*.cpp")
FOREACH(exe_cpp ${cpp_files})
    GET_FILENAME_COMPONENT(exe ${exe_cpp} NAME_WE)
    ADD_EXECUTABLE(${exe} bin/${exe}.cpp)
    TARGET_LINK_LIBRARIES(${exe} HepLib)
    INSTALL(TARGETS ${exe} DESTINATION bin)
ENDFOREACH()

ADD_EXECUTABLE(qgraf "bin/qgraf-3.1.4.f")
INSTALL(TARGETS qgraf DESTINATION bin)
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/OTHER/HepLib.pc" DESTINATION "lib/pkgconfig")

# OpenMP
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

MESSAGE(STATUS "Install Prefix: ${CMAKE_INSTALL_PREFIX}")

