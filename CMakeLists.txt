cmake_minimum_required(VERSION 3.14)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" VERSION)
project(ccsd VERSION "${VERSION}" LANGUAGES CXX)

set(
        CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
        CACHE STRING "" FORCE
)

include(get_cpp)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("${CMAKE_CURRENT_SOURCE_DIR}/ccsd")

### Options ###
option(BUILD_TESTING "Should we build the tests?" OFF)
option(BUILD_INTEGRATION_TESTING "Should we build the integration tests?" OFF)
option(BUILD_PYBINDINGS "Should we build Python3 bindings?" OFF)

cpp_find_or_build_dependency(
        simde
        URL github.com/NWChemEx-Project/SimDE
        PRIVATE TRUE
        BUILD_TARGET simde
        FIND_TARGET nwx::simde
        CMAKE_ARGS BUILD_TESTING=OFF
)

cpp_find_or_build_dependency(
         tamm
         URL github.com/NWChemEx-Project/TAMM
         VERSION main
         BUILD_TARGET tamm
         FIND_TARGET tamm::tamm
         CMAKE_ARGS BUILD_METHODS=OFF
)

cpp_add_library(
        coupledcluster
        SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/cc"
        INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include/cc"
        DEPENDS nwx::simde tamm::tamm
)

if("${BUILD_TESTING}")
    cpp_find_or_build_dependency(
            Catch2
            VERSION v2.x
            URL github.com/catchorg/Catch2
            BUILD_TARGET Catch2
            VERSION v2.13.3
            FIND_TARGET Catch2::Catch2
    )
    cpp_find_or_build_dependency(
            mokup
            URL github.com/NWChemEX-Project/Mokup
            PRIVATE
            BUILD_TARGET mokup
            FIND_TARGET nwx::mokup
    )
    cpp_add_tests(
            test_coupledcluster
            SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/tests"
            INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include/cc"
            DEPENDS Catch2::Catch2 coupledcluster nwx::mokup
    )
endif()

if("${BUILD_INTEGRATION_TESTING}")
    cpp_find_or_build_dependency(
            Catch2
            VERSION v2.x
            URL github.com/catchorg/Catch2
            BUILD_TARGET Catch2
            VERSION v2.13.3
            FIND_TARGET Catch2::Catch2
    )
    cpp_find_or_build_dependency(
            integrals
            URL github.com/NWChemEx-Project/Integrals
            PRIVATE TRUE
            BUILD_TARGET integrals
            FIND_TARGET nwx::integrals
            CMAKE_ARGS BUILD_TESTING=OFF
    )
    cpp_find_or_build_dependency(
            scf
            URL github.com/NWChemEx-Project/SCF
            PRIVATE TRUE
            BUILD_TARGET scf
            FIND_TARGET nwx::scf
            CMAKE_ARGS BUILD_TESTING=OFF
    )
    cpp_add_tests(
            integration_coupledcluster
            SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/integration_tests"
            INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/ccsd"
            DEPENDS Catch2::Catch2 coupledcluster nwx::integrals nwx::scf
    )
endif()
