# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.13)

project(SEALExamples VERSION 4.1.2 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find BLAS library
find_package(BLAS REQUIRED)

# If not called from root CMakeLists.txt
if(NOT DEFINED SEAL_BUILD_EXAMPLES)
    set(SEAL_BUILD_EXAMPLES ON)

    set(SEAL_DIR "../seallibs/lib/cmake/SEAL-4.1")
    
    # Import Microsoft SEAL
    find_package(SEAL 4.1.2 EXACT REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()

if(SEAL_BUILD_EXAMPLES)
    add_executable(test_example)

    target_sources(test_example
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/test.cpp
            ${CMAKE_CURRENT_LIST_DIR}/utils.cpp
    )

    if(TARGET SEAL::seal)
        target_link_libraries(test_example PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(test_example PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

    # Add new test executable for CPMM
    add_executable(test_cpmm)

    target_sources(test_cpmm
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/cpmm.cpp
            ${CMAKE_CURRENT_LIST_DIR}/utils.cpp
    )

    if(TARGET SEAL::seal)
        target_link_libraries(test_cpmm PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(test_cpmm PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

    if(BLAS_FOUND)
        target_link_libraries(test_cpmm PRIVATE ${BLAS_LIBRARIES})
        if (BLAS_INCLUDE_DIRS)
            target_include_directories(test_cpmm PRIVATE ${BLAS_INCLUDE_DIRS})
        endif()
    endif()
endif()
