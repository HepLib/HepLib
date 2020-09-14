# usage 
# cmake -DCMAKE_INSTALL_PREFIX=p0 -DINC_PATH="p1;p2" -DLIB_PATH="p3;p4" path

CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)
SET(CMAKE_FORTRAN_COMPILER gfortran)

#-----------------------------------------------
# INC_PATH & LIB_PATH
#-----------------------------------------------
SET(INC_PATH "${INC_PATH};${CMAKE_INSTALL_PREFIX}/include")
LIST(REMOVE_DUPLICATES INC_PATH)
SET(INC_FLAGS "")
FOREACH(dir ${INC_PATH})
IF(EXISTS ${dir})
    INCLUDE_DIRECTORIES(${dir})
    SET(INC_FLAGS "${INC_FLAGS} -I'${dir}'")
ENDIF()
ENDFOREACH()
STRING(FIND "${INC_FLAGS}" "${CMAKE_INSTALL_PREFIX}/include" INC_POS)
IF( "${INC_POS}" STREQUAL "-1" )
    SET(INC_FLAGS "${INC_FLAGS} -I'${CMAKE_INSTALL_PREFIX}/include'")
ENDIF()

SET(LIB_PATH "${LIB_PATH};${CMAKE_INSTALL_PREFIX}/lib")
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

INCLUDE_DIRECTORIES(BEFORE include)
CONFIGURE_FILE(include/Init.cpp "${CMAKE_CURRENT_BINARY_DIR}/OTHER/Init.cpp" @ONLY)
CONFIGURE_FILE(include/HepLib.pc "${CMAKE_CURRENT_BINARY_DIR}/OTHER/HepLib.pc" @ONLY)
CONFIGURE_FILE(include/heplib++ "${CMAKE_CURRENT_BINARY_DIR}/OTHER/heplib++" @ONLY)

#-----------------------------------------------
# External Libraries
#-----------------------------------------------
LINK_LIBRARIES(ginac cln Minuit2 mpfr qhullstatic quadmath gomp cubaq)

#-----------------------------------------------
PROJECT(HepLib VERSION 1.0 LANGUAGES C CXX Fortran)
AUX_SOURCE_DIRECTORY(Basic BASIC_SRCS)
AUX_SOURCE_DIRECTORY(SD SD_SRCS)
AUX_SOURCE_DIRECTORY(FC FC_SRCS)
AUX_SOURCE_DIRECTORY(IBP IBP_SRCS)
AUX_SOURCE_DIRECTORY(GiNaC GiNaC_SRCS)
AUX_SOURCE_DIRECTORY("${CMAKE_CURRENT_BINARY_DIR}/OTHER" OTHER)

#-----------------------------------------------
# WSTP - Wolfram Symbolic Transfer Protocol
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
ADD_LIBRARY(HepLib SHARED ${BASIC_SRCS} ${SD_SRCS} ${FC_SRCS} ${IBP_SRCS} ${GiNaC_SRCS} ${OTHER})
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
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/OTHER/heplib++" DESTINATION "bin" 
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)

# OpenMP
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=gnu++14")
