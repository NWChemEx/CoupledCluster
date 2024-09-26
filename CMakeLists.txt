# Copyright 2024 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.14)
set(VERSION 1.0.0) #TODO: get from git
project(coupledcluster VERSION "${VERSION}" LANGUAGES CXX)

include(FetchContent)
enable_language(C)
FetchContent_Declare(
    nwx_cmake
    GIT_REPOSITORY https://github.com/NWChemEx/NWXCMake
)
FetchContent_MakeAvailable(nwx_cmake)
list(APPEND CMAKE_MODULE_PATH "${nwx_cmake_SOURCE_DIR}/cmake")

set(
    CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
    CACHE STRING "" FORCE
)

include(nwx_versions)
include(get_cmaize)
include(nwx_cxx_api_docs)

set(CoupledCluster_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(CoupledCluster_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
nwx_cxx_api_docs("${CoupledCluster_SOURCE_DIR}" "${CoupledCluster_INCLUDE_DIR}")

set(DEPENDENCIES "")

cmaize_option_list(
    BUILD_TESTING OFF "Should we build the tests?"
    BUILD_PYBIND11_PYBINDINGS ON "Build pybind11 python3 bindings?"
    ENABLE_EXPERIMENTAL_FEATURES OFF "Build features which are not 1.0-ed yet?"
)

cmaize_find_or_build_dependency(
    scf
    URL github.com/NWChemEx/SCF
    BUILD_TARGET scf
    FIND_TARGET nwx::scf
    CMAKE_ARGS BUILD_TESTING=OFF
)
list(APPEND DEPENDENCIES scf)

cmaize_find_or_build_dependency(
    tamm
    URL github.com/NWChemEx-Project/TAMM
    VERSION main
    BUILD_TARGET tamm
    FIND_TARGET tamm::tamm
    CMAKE_ARGS MODULES="CC;DFT"
)
list(APPEND DEPENDENCIES tamm)

cmaize_find_or_build_dependency(
    exachem
    URL github.com/ExaChem/exachem
    VERSION main
    BUILD_TARGET exachem
    FIND_TARGET exachem::exachem
    CMAKE_ARGS MODULES="CC;DFT"
)
list(APPEND DEPENDENCIES exachem)

cmaize_find_or_build_dependency(
    Libint2
    URL github.com/evaleev/libint
    VERSION 2.7.2
    BUILD_TARGET int2
    FIND_TARGET Libint2::int2
)
list(APPEND DEPENDENCIES Libint2)

cmaize_add_library(
    coupledcluster
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/cc"
    INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include/cc"
    DEPENDS "${DEPENDENCIES}"
)

include(nwx_pybind11)
nwx_add_pybind11_module(
    ${PROJECT_NAME}
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/python"
    DEPENDS "${PROJECT_NAME}"
)

cmaize_find_or_build_dependency(
    chemcache 
    URL github.com/NWChemEx/ChemCache
    BUILD_TARGET chemcache
    FIND_TARGET nwx::chemcache
    CMAKE_ARGS BUILD_TESTING=OFF
               BUILD_PYBIND11_PYBINDINGS=${BUILD_PYBIND11_PYBINDINGS}
)

cmaize_add_executable(
    cc_driver
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/examples/driver"
    DEPENDS coupledcluster chemcache
)

if("${BUILD_TESTING}") 
    include(CTest)
    set(PYTHON_TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/tests/python")
    set(CXX_TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/tests/cxx")

    # Need to call with MPI, so we have to deconstruct the cmaize_add_test and
    # nwx_pybind11_tests methods
    cmaize_find_or_build_dependency(
        Catch2
        URL github.com/catchorg/Catch2
        BUILD_TARGET Catch2
        FIND_TARGET Catch2::Catch2
        VERSION v3.6.0
    )
    cmaize_add_executable(
        test_coupledcluster
        SOURCE_DIR "${CXX_TEST_DIR}"
        INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/src/cc"
        DEPENDS Catch2 coupledcluster chemcache
    )
    add_test(
        NAME "test_coupledcluster"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
                "${CMAKE_BINARY_DIR}/test_coupledcluster"
    )

    add_test(
        NAME py_test_unit_${PROJECT_NAME}
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
                "${Python_EXECUTABLE}" 
                "${PYTHON_TEST_DIR}/unit_tests/test_coupledcluster.py"
    )
    nwx_python_path(
        TEST_PYTHONPATH
        SUBMODULES chemcache simde chemist pluginplay parallelzone
    )
    set_tests_properties(
        py_test_unit_${PROJECT_NAME}
        PROPERTIES ENVIRONMENT "${TEST_PYTHONPATH}"
    )
endif()