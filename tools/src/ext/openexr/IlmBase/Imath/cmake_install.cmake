# Install script for directory: /home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/ggrbb/code/pbrt-v3/tools/src/ext/openexr/IlmBase/Imath/libImath.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathBoxAlgo.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathBox.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathColorAlgo.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathColor.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathEuler.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathExc.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathExport.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathForward.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathFrame.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathFrustum.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathFrustumTest.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathFun.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathGL.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathGLU.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathHalfLimits.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathInt64.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathInterval.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathLimits.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathLineAlgo.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathLine.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathMath.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathMatrixAlgo.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathMatrix.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathNamespace.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathPlane.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathPlatform.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathQuat.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathRandom.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathRoots.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathShear.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathSphere.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathVecAlgo.h"
    "/home/ggrbb/code/pbrt-v3/src/ext/openexr/IlmBase/Imath/ImathVec.h"
    )
endif()

