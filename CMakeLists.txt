CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

#-----------------------------------------------
# Configuration Here
#-----------------------------------------------
IF(APPLE)
    SET(CMAKE_INSTALL_PREFIX "/usr/local/feng")
ELSE()
    SET(CMAKE_INSTALL_PREFIX "$ENV{HOME}/usr/local")
ENDIF()

SET(INC_DIRS
    "/usr/local/feng/include"
    $ENV{HOME}/usr/local/include
    $ENV{HOME}/usr/include
)

SET(LIB_DIRS
    "/usr/local/feng/lib"
    $ENV{HOME}/usr/local/lib
    $ENV{HOME}/usr/lib
)

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
    ginac cln Minuit2 mpfr qhullstatic quadmath gomp cubaq
)

PROJECT(HepLib VERSION 1.0 LANGUAGES C CXX)
AUX_SOURCE_DIRECTORY(ExGiNaC ExGiNaC_SRCS)
AUX_SOURCE_DIRECTORY(SD SD_SRCS)

#-----------------------------------------------
# check WSTP
#-----------------------------------------------
IF(APPLE)
    FIND_LIBRARY(HAS_WSTP WSTPi4 ${LINK_DIRECTORIES})
    IF(HAS_WSTP)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lc++ -lWSTPi4 -framework Foundation")
    ELSE()
        STRING(REPLACE "ExGiNaC/WSKernel.cpp" "" ExGiNaC_SRCS "${ExGiNaC_SRCS}")
    ENDIF(HAS_WSTP)
ELSE()
    FIND_LIBRARY(HAS_WSTP WSTP64i4 ${LINK_DIRECTORIES})
    IF(HAS_WSTP)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -lpthread -lrt -lstdc++ -ldl -luuid -lWSTP64i4")
    ELSE()
        STRING(REPLACE "ExGiNaC/WSKernel.cpp" "" ExGiNaC_SRCS "${ExGiNaC_SRCS}")
    ENDIF(HAS_WSTP)
ENDIF(APPLE)

#-----------------------------------------------
# install *.h #
#-----------------------------------------------
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include
DESTINATION include
FILES_MATCHING PATTERN "*.h")

#-----------------------------------------------
# install library HepLib #
#-----------------------------------------------
ADD_LIBRARY(HepLib SHARED ${ExGiNaC_SRCS} ${SD_SRCS})
INSTALL(TARGETS HepLib DESTINATION lib)

#-----------------------------------------------
# install bin/*
#-----------------------------------------------
FILE(GLOB_RECURSE cpp_files "bin/*.cpp")
FOREACH(exe_cpp ${cpp_files})
    GET_FILENAME_COMPONENT(exe ${exe_cpp} NAME_WE)
    ADD_EXECUTABLE(${exe} bin/${exe}.cpp)
    TARGET_LINK_LIBRARIES(${exe} HepLib)
    INSTALL(TARGETS ${exe} DESTINATION bin)
ENDFOREACH()

MESSAGE(STATUS "Install Prefix: ${CMAKE_INSTALL_PREFIX}")

