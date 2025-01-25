# Install script for directory: /home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/ggrbb/code/pbrt-v3/tools/src/ext/openexr/IlmBase/Iex/libIex.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexBaseExc.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexMathExc.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexThrowErrnoExc.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexErrnoExc.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexMacros.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/Iex.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexNamespace.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexExport.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Iex/IexForward.h"
    )
endif()

