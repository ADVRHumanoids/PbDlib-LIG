CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(pbdlib)

INCLUDE(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
IF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")
ELSEIF(COMPILER_SUPPORTS_CXX0X)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ELSE()
    MESSAGE(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
ENDIF()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(MacroInstallLib)

FIND_PACKAGE(Armadillo 5.4 REQUIRED)
SET(CMAKE_MAJOR_VERSION 1)
SET(CMAKE_MINOR_VERSION 0)
SET(CMAKE_PATCH_VERSION 0)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

ADD_LIBRARY(pbd SHARED
		src/datapoints.cpp
		src/demonstration.cpp
		src/tpdemonstration.cpp
		src/taskparameters.cpp
		src/mvn.cpp
		src/gmm.cpp
		src/gmr.cpp
		src/tpgmm.cpp
		src/lqr.cpp
		src/hmm.cpp
		src/hsmm.cpp
		src/adhsmm.cpp
		src/quaternion.cpp
		src/trajgmm.cpp
		src/mpc.cpp
		src/trajMPC.cpp
)

#Uncomment for debugging
#add_definitions(-g)

################################################################################
### samples
################################################################################
OPTION(PBDLIB_BUILD_TEST "Build PbDLib examples" 1)
IF (PBDLIB_BUILD_TEST)
  ADD_SUBDIRECTORY(examples)
ENDIF ()

TARGET_LINK_LIBRARIES(pbd ${ARMADILLO_LIBRARIES})

# Copy data/ to build so that the examples can load from that location
FILE(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/examples/)

################################################################################
### install
################################################################################
SET(ROOT_INSTALL_LIBDIR lib)
SET(ROOT_INSTALL_INCLUDEDIR include)
FILE(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/pbdlib/*.h)

INSTALL(TARGETS pbd DESTINATION ${ROOT_INSTALL_LIBDIR})	
INSTALL(FILES ${headers} DESTINATION ${ROOT_INSTALL_INCLUDEDIR}/pbdlib/)

library_install(pbd 1 0 0)


