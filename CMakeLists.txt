# usage:
# cmake -DCMAKE_INSTALL_PREFIX=p0 -DINC_PATH="p1;p2" -DLIB_PATH="p3;p4" path
#  p0 is the directory where HepLib will be installed
#  if the required libraries[ CLN/GiNaC/Flint etc] are not installed to p0,
#  one need to set p1/p3 the corresponding include and lib path
# e.g. on MacOS-M1
# cmake -DCMAKE_INSTALL_PREFIX=$HOME/usr -DINC_PATH=/opt/homebrew/include -DLIB_PATH=/opt/homebrew/lib ..

CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)
SET(CMAKE_FORTRAN_COMPILER gfortran)

#-----------------------------------------------
# qgraf
#-----------------------------------------------
ADD_EXECUTABLE(qgraf "QGRAF/qgraf-3.1.4.f")

#-----------------------------------------------
# INC_PATH & LIB_PATH
#-----------------------------------------------
SET(INC_PATH "${CMAKE_INSTALL_PREFIX}/include;${INC_PATH}")
LIST(REMOVE_DUPLICATES INC_PATH)
SET(INC_FLAGS "")
SET(INC_FLAGS_2 "")
FOREACH(dir ${INC_PATH})
IF(EXISTS ${dir})
    INCLUDE_DIRECTORIES(${dir})
    SET(INC_FLAGS "${INC_FLAGS} -I'${dir}'")
    SET(INC_FLAGS_2 "${INC_FLAGS_2};-I'${dir}'")
ENDIF()
ENDFOREACH()
STRING(FIND "${INC_FLAGS}" "${CMAKE_INSTALL_PREFIX}/include" INC_POS)
IF( "${INC_POS}" STREQUAL "-1" )
    SET(INC_FLAGS "${INC_FLAGS} -I'${CMAKE_INSTALL_PREFIX}/include'")
ENDIF()
STRING(FIND "${INC_FLAGS_2}" "${CMAKE_INSTALL_PREFIX}/include" INC_POS_2)
IF( "${INC_POS_2}" STREQUAL "-1" )
    SET(INC_FLAGS_2 "${INC_FLAGS_2};-I'${CMAKE_INSTALL_PREFIX}/include'")
ENDIF()

SET(LIB_PATH "${CMAKE_INSTALL_PREFIX}/lib;${LIB_PATH}")
LIST(REMOVE_DUPLICATES LIB_PATH)
SET(LIB_FLAGS "-Wl,-rpath,.")
FOREACH(dir ${LIB_PATH})
IF(EXISTS ${dir})
    LINK_DIRECTORIES(${dir})
    SET(LIB_FLAGS "${LIB_FLAGS} -Wl,-rpath,'${dir}' -L'${dir}'")
ENDIF()
ENDFOREACH()
STRING(FIND "${LIB_FLAGS}" "${CMAKE_INSTALL_PREFIX}/lib" LIB_POS)
IF( "${LIB_POS}" STREQUAL "-1" )
    SET(LIB_FLAGS "${LIB_FLAGS} -Wl,-rpath,'${CMAKE_INSTALL_PREFIX}/lib' -L'${CMAKE_INSTALL_PREFIX}/lib'")
ENDIF()

#-----------------------------------------------
# External Libraries
#-----------------------------------------------
LINK_LIBRARIES(ginac cln mpfr qhullstatic_r quadmath gomp dl flint gmp jemalloc)
#LINK_LIBRARIES(ginac cln mpfr qhullstatic_r quadmath gomp dl flint gmp)

#-----------------------------------------------
# HepLib - Include & Configure
#-----------------------------------------------
INCLUDE_DIRECTORIES(BEFORE MB)
INCLUDE_DIRECTORIES(BEFORE IBP)
INCLUDE_DIRECTORIES(BEFORE SD)
INCLUDE_DIRECTORIES(BEFORE DE)
INCLUDE_DIRECTORIES(BEFORE AMF)
INCLUDE_DIRECTORIES(BEFORE HEP)
INCLUDE_DIRECTORIES(BEFORE QGRAF)
INCLUDE_DIRECTORIES(BEFORE QCD)
#INCLUDE_DIRECTORIES(BEFORE EX)
INCLUDE_DIRECTORIES(BEFORE BASIC)
INCLUDE_DIRECTORIES(BEFORE INC)

CONFIGURE_FILE(INC/Init.cpp "${CMAKE_CURRENT_BINARY_DIR}/OINC/Init.cpp" @ONLY)
CONFIGURE_FILE(INC/heplib++ "${CMAKE_CURRENT_BINARY_DIR}/OINC/heplib++" @ONLY)
 
#-----------------------------------------------
PROJECT(HepLib VERSION 1.5 LANGUAGES C CXX Fortran)
AUX_SOURCE_DIRECTORY(BASIC BASIC_SRCS)
AUX_SOURCE_DIRECTORY(SD SD_SRCS)
AUX_SOURCE_DIRECTORY(DE DE_SRCS)
AUX_SOURCE_DIRECTORY(AMF AMF_SRCS)
AUX_SOURCE_DIRECTORY(HEP HEP_SRCS)
AUX_SOURCE_DIRECTORY(QGRAF QGRAF_SRCS)
AUX_SOURCE_DIRECTORY(QCD QCD_SRCS)
AUX_SOURCE_DIRECTORY(IBP IBP_SRCS)
AUX_SOURCE_DIRECTORY(EX EX_SRCS)
AUX_SOURCE_DIRECTORY("${CMAKE_CURRENT_BINARY_DIR}/OINC" OINC_SRCS)

#-----------------------------------------------
# install library HepLib
#-----------------------------------------------
ADD_LIBRARY(HepLib SHARED ${BASIC_SRCS} ${HEP_SRCS} ${QGRAF_SRCS} ${QCD_SRCS} ${IBP_SRCS} ${SD_SRCS} ${DE_SRCS} ${AMF_SRCS} ${EX_SRCS} ${OINC_SRCS} ${MPQ_SRCS})
SET_TARGET_PROPERTIES(HepLib PROPERTIES INSTALL_RPATH "$ORIGIN;${LIB_PATH}")
INSTALL(TARGETS HepLib DESTINATION lib)

#-----------------------------------------------
# install header
#-----------------------------------------------
INSTALL(FILES "${CMAKE_SOURCE_DIR}/INC/HepLib.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/BASIC/BASIC.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/BASIC/Symbols.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/BASIC/PoolPipe.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/MB/MB.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/IBP/IBP.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/SD/SD.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/DE/DE.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/DE/exFlint.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/AMF/AMF.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/SD/NFunctions.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/SD/mpreal.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/HEP/HEP.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/QGRAF/QGRAF.h" DESTINATION include)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/QCD/QCD.h" DESTINATION include)
#INSTALL(FILES "${CMAKE_SOURCE_DIR}/EX/HepLibW.h" DESTINATION include)
#INSTALL(FILES "${CMAKE_SOURCE_DIR}/EX/HepLib.i" DESTINATION include)
# backup the source
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/BASIC" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/bin" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/DE" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/AMF" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/doc" DESTINATION src)
#INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/EX" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/HEP" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/IBP" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/INC" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/MB" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/QCD" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/QGRAF" DESTINATION src)
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/SD" DESTINATION src)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/CMakeLists.txt" DESTINATION src)

#-----------------------------------------------
# install bin/*
#-----------------------------------------------
FILE(GLOB cpp_files "bin/*.cpp")
FOREACH(exe_cpp ${cpp_files})
    GET_FILENAME_COMPONENT(exe ${exe_cpp} NAME_WE)
    ADD_EXECUTABLE(${exe} bin/${exe}.cpp)
    SET_TARGET_PROPERTIES(${exe} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib;${LIB_PATH}")
    TARGET_LINK_LIBRARIES(${exe} HepLib)
    INSTALL(TARGETS ${exe} DESTINATION bin)
ENDFOREACH()

INSTALL(TARGETS qgraf DESTINATION bin)
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/OINC/heplib++" DESTINATION "bin"
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)

#-----------------------------------------------
# OpenMP support
#-----------------------------------------------
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -fopenmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -fopenmp")

#-----------------------------------------------
# Command Line Tools 15.0 - bug work around
#-----------------------------------------------
# if you get an error message like
# 'libunwind: _Unwind_GetTextRelBase - _Unwind_GetTextRelBase() not implemented'
# Please uncomment the following line
#if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ld_classic")
#endif()


#-----------------------------------------------
# PreCompiled Headers
#-----------------------------------------------
ADD_CUSTOM_TARGET(
    pch
    COMMAND ${CMAKE_CXX_COMPILER} -c -w ${INC_FLAGS_2} -x c++-header ${CMAKE_INSTALL_PREFIX}/include/HepLib.h
    COMMAND ${CMAKE_CXX_COMPILER} -c -w ${INC_FLAGS_2} -x c++-header ${CMAKE_INSTALL_PREFIX}/include/NFunctions.h
    COMMENT "Compiling headers ..."
)
