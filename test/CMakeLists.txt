cmake_minimum_required(VERSION 3.2)
include(ExternalProject)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

project(iitii)
INCLUDE(ExternalProject)
include_directories(${PROJECT_SOURCE_DIR}/..)

# After build, double-check external/src/htslib-stamp/htslib-build-out.log
# to make sure CFLAGS were set as expected
ExternalProject_Add(htslib
    URL https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND bash -c "autoreconf && CFLAGS='-gdwarf -O3 -DNDEBUG -fvisibility=hidden' ./configure --enable-libcurl=no"
    BUILD_COMMAND bash -c "make -n && make -j$(nproc)"
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_BUILD ON
  )
ExternalProject_Get_Property(htslib source_dir)
set(HTSLIB_INCLUDE_DIR ${source_dir}/htslib)
include_directories(${HTSLIB_INCLUDE_DIR})
ExternalProject_Get_Property(htslib binary_dir)
add_library(libhts STATIC IMPORTED)
set(LIBHTS_A ${binary_dir}/libhts.a)
set_property(TARGET libhts PROPERTY IMPORTED_LOCATION ${LIBHTS_A})

ExternalProject_Add(CTPL
    URL https://github.com/vit-vit/CTPL/archive/v.0.0.2.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(CTPL source_dir)
set(CTPL_INCLUDE_DIR ${source_dir})
include_directories(${CTPL_INCLUDE_DIR})

ExternalProject_Add(catch
  URL https://github.com/philsquared/Catch/archive/v2.9.1.zip
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property(catch source_dir)
include_directories(${source_dir}/single_include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17 -Wall -Werror=return-type -Wno-reorder -fdiagnostics-color=auto -ffast-math")

add_executable(test_iitii util.h test_iitii.cc)
add_dependencies(test_iitii catch htslib)
target_link_libraries(test_iitii libhts libz.a libbz2.a liblzma.a libdeflate.a)

add_executable(ideal_benchmark util.h ideal_benchmark.cc)
add_dependencies(ideal_benchmark htslib)
target_link_libraries(ideal_benchmark libhts libz.a libbz2.a liblzma.a libdeflate.a)

add_executable(gnomad_benchmark util.h gnomad_benchmark.cc)
add_dependencies(gnomad_benchmark htslib)
target_link_libraries(gnomad_benchmark libhts libz.a libbz2.a liblzma.a libdeflate.a)

include(CTest)
add_test(NAME unit_tests COMMAND bash -c "LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so ./test_iitii -d yes")
