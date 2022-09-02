cmake_minimum_required(VERSION 3.22.0)

include(FetchContent)

project(CoupledCluster VERSION 0.0.0 LANGUAGES C CXX Fortran)
include(${PROJECT_SOURCE_DIR}/cmake/check_user_env.cmake)

find_package( CMakeBuild QUIET )

if( NOT ${CMakeBuild_FOUND} )
  if (NOT CMSB_TAG)
    set(CMSB_TAG 858b1b20df53b779bdbcd1b4cf7f674d9b7077a8)
    if (ENABLE_DEV_MODE)
      set(CMSB_TAG main)
    endif()
  endif()

  FetchContent_Declare(
    CMakeBuild
    GIT_REPOSITORY https://github.com/NWChemEx-Project/CMakeBuild.git
    GIT_TAG ${CMSB_TAG}
  )

  FetchContent_MakeAvailable( CMakeBuild )
  set(CMSB_MACROS ${CMakeBuild_SOURCE_DIR}/cmake/macros)

  FetchContent_MakeAvailable( linalg-cmake-modules )
  set(CMSB_LAM_PATH ${linalg-cmake-modules_SOURCE_DIR})
endif()

set(CMSB_PROJECTS CoupledCluster)

set(CoupledCluster_DEPENDENCIES TAMM)

set(BUILD_TESTS OFF)
set(BUILD_LIBINT ON)

set(CoupledCluster_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CoupledCluster)
set(CoupledCluster_METHODS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/methods/)

if(NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

set(CMSB_BASISSET_DIR ${CMAKE_CURRENT_SOURCE_DIR}/methods/basis)
build_cmsb_module(${CMAKE_CURRENT_LIST_DIR})

