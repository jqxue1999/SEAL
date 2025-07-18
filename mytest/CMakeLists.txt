# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.13)

project(SEALExamples VERSION 4.1.2 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find BLAS library
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)

# Enable AVX2 optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma -march=native -mtune=native")
message(STATUS "AVX2 optimization enabled: -mavx2 -mfma -march=native -mtune=native")

# Find nlohmann/json library
find_package(nlohmann_json 3.2.0 REQUIRED)
message(STATUS "Found nlohmann/json")

# If not called from root CMakeLists.txt
if(NOT DEFINED SEAL_BUILD_EXAMPLES)
    set(SEAL_BUILD_EXAMPLES ON)

    set(SEAL_DIR "../seallibs/lib/cmake/SEAL-4.1")
    
    # Import Microsoft SEAL
    find_package(SEAL 4.1.2 EXACT REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()

if(SEAL_BUILD_EXAMPLES)
    list(APPEND CMAKE_PREFIX_PATH "/home/jiaq/Research/SEAL/benchmarklibs")
    find_package(benchmark REQUIRED)

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

    # # 新的统一测试程序
    # add_executable(example)

    # target_sources(example
    #     PRIVATE
    #         ${CMAKE_CURRENT_LIST_DIR}/example.cpp
    #         ${CMAKE_CURRENT_LIST_DIR}/cpmm.cpp
    #         ${CMAKE_CURRENT_LIST_DIR}/ccmm.cpp
    #         ${CMAKE_CURRENT_LIST_DIR}/digits.cpp
    #         ${CMAKE_CURRENT_LIST_DIR}/cpscale.cpp
    #         ${CMAKE_CURRENT_LIST_DIR}/utils.cpp
    # )

    # if(TARGET SEAL::seal)
    #     target_link_libraries(example PRIVATE SEAL::seal)
    # elseif(TARGET SEAL::seal_shared)
    #     target_link_libraries(example PRIVATE SEAL::seal_shared)
    # else()
    #     message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    # endif()

    # if(BLAS_FOUND)
    #     target_link_libraries(example PRIVATE ${BLAS_LIBRARIES})
    #     if (BLAS_INCLUDE_DIRS)
    #         target_include_directories(example PRIVATE ${BLAS_INCLUDE_DIRS})
    #     endif()
    # endif()
    
    # # Link nlohmann/json library
    # target_link_libraries(example PRIVATE nlohmann_json::nlohmann_json)
    # target_link_libraries(example PRIVATE OpenMP::OpenMP_CXX)

    # Benchmark for plain_blas.cpp
    add_executable(plain_blas_bench)
    target_sources(plain_blas_bench
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/plain_blas_bench.cpp
            ${CMAKE_CURRENT_LIST_DIR}/plain_blas.cpp
    )
    target_link_libraries(plain_blas_bench PRIVATE benchmark::benchmark pthread OpenMP::OpenMP_CXX BLAS::BLAS)

    # Benchmark for plain_omp.cpp
    add_executable(plain_omp_bench)
    target_sources(plain_omp_bench
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/plain_omp_bench.cpp
            ${CMAKE_CURRENT_LIST_DIR}/plain_omp.cpp
    )
    target_link_libraries(plain_omp_bench PRIVATE benchmark::benchmark pthread OpenMP::OpenMP_CXX)

    # Benchmark for coefficients_seal.cpp
    add_executable(coefficients_seal_bench)

    target_sources(coefficients_seal_bench
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/coefficients_seal_bench.cpp
            ${CMAKE_CURRENT_LIST_DIR}/coefficients_seal.cpp
            ${CMAKE_CURRENT_LIST_DIR}/cpmm.cpp
            ${CMAKE_CURRENT_LIST_DIR}/coeff_modulus_config.cpp
    )

    if(TARGET SEAL::seal)
        target_link_libraries(coefficients_seal_bench PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(coefficients_seal_bench PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

    if(BLAS_FOUND)
        target_link_libraries(coefficients_seal_bench PRIVATE ${BLAS_LIBRARIES})
        if (BLAS_INCLUDE_DIRS)
            target_include_directories(coefficients_seal_bench PRIVATE ${BLAS_INCLUDE_DIRS})
        endif()
    endif()
    
    # Link nlohmann/json library
    target_link_libraries(coefficients_seal_bench PRIVATE nlohmann_json::nlohmann_json)
    target_link_libraries(coefficients_seal_bench PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(coefficients_seal_bench PRIVATE benchmark::benchmark pthread)

    # Benchmark for coefficients_blas.cpp
    add_executable(coefficients_blas_bench)
    target_sources(coefficients_blas_bench
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/coefficients_blas_bench.cpp
            ${CMAKE_CURRENT_LIST_DIR}/coefficients_blas.cpp
            ${CMAKE_CURRENT_LIST_DIR}/cpmm.cpp
            ${CMAKE_CURRENT_LIST_DIR}/cpscale.cpp
            ${CMAKE_CURRENT_LIST_DIR}/coeff_modulus_config.cpp
    )

    if(TARGET SEAL::seal)
        target_link_libraries(coefficients_blas_bench PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(coefficients_blas_bench PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

    if(BLAS_FOUND)
        target_link_libraries(coefficients_blas_bench PRIVATE ${BLAS_LIBRARIES})
        if (BLAS_INCLUDE_DIRS)
            target_include_directories(coefficients_blas_bench PRIVATE ${BLAS_INCLUDE_DIRS})
        endif()
    endif()
    
    # Link nlohmann/json library
    target_link_libraries(coefficients_blas_bench PRIVATE nlohmann_json::nlohmann_json)
    target_link_libraries(coefficients_blas_bench PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(coefficients_blas_bench PRIVATE benchmark::benchmark pthread)

    # Benchmark for digits_bench.cpp
    add_executable(digits_bench)
    target_sources(digits_bench
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/digits_bench.cpp
            ${CMAKE_CURRENT_LIST_DIR}/digits.cpp
            ${CMAKE_CURRENT_LIST_DIR}/cpmm.cpp
            ${CMAKE_CURRENT_LIST_DIR}/coeff_modulus_config.cpp
    )

    if(TARGET SEAL::seal)
        target_link_libraries(digits_bench PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(digits_bench PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

    if(BLAS_FOUND)
        target_link_libraries(digits_bench PRIVATE ${BLAS_LIBRARIES})
        if (BLAS_INCLUDE_DIRS)
            target_include_directories(digits_bench PRIVATE ${BLAS_INCLUDE_DIRS})
        endif()
    endif()
    
    # Link nlohmann/json library
    target_link_libraries(digits_bench PRIVATE nlohmann_json::nlohmann_json)
    target_link_libraries(digits_bench PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(digits_bench PRIVATE benchmark::benchmark pthread)

    # Benchmark for baseline_bench.cpp
    add_executable(baseline_bench)
    target_sources(baseline_bench
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/baseline_bench.cpp
            ${CMAKE_CURRENT_LIST_DIR}/baseline.cpp
            ${CMAKE_CURRENT_LIST_DIR}/coeff_modulus_config.cpp
    )

    if(TARGET SEAL::seal)
        target_link_libraries(baseline_bench PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(baseline_bench PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

    if(BLAS_FOUND)
        target_link_libraries(baseline_bench PRIVATE ${BLAS_LIBRARIES})
        if (BLAS_INCLUDE_DIRS)
            target_include_directories(baseline_bench PRIVATE ${BLAS_INCLUDE_DIRS})
        endif()
    endif()
    
    # Link nlohmann/json library
    target_link_libraries(baseline_bench PRIVATE nlohmann_json::nlohmann_json)
    target_link_libraries(baseline_bench PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(baseline_bench PRIVATE benchmark::benchmark pthread)



endif()
