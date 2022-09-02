cmake_minimum_required(VERSION ${CMAKE_VERSION})
project(CoupledCluster-Methods VERSION 0.0.0 LANGUAGES CXX)
find_package(CoupledCluster REQUIRED)
include(TargetMacros)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/util/external)
include(${CMAKE_CURRENT_LIST_DIR}/scf/scf.cmake)
 include(${CMAKE_CURRENT_LIST_DIR}/cc/cc.cmake)
