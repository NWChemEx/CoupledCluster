cmake_minimum_required(VERSION ${CMAKE_VERSION})
project(CoupledCluster VERSION ${PROJECT_VERSION} LANGUAGES CXX)
include(TargetMacros)

#Add the current directory's sources to the list
set(CoupledCluster_SRCS CoupledCluster.cpp)

#Add the current directory's header files to the list
set(CoupledCluster_INCLUDES )

set(CoupledCluster_CFLAGS )
set(CoupledCluster_LFLAGS )

cmsb_add_library(coupledcluster CoupledCluster_SRCS CoupledCluster_INCLUDES CoupledCluster_CFLAGS CoupledCluster_LFLAGS)
target_include_directories(coupledcluster PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)


