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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17 -Wall -Werror=return-type -Wno-reorder -fdiagnostics-color=auto")

add_executable(test_iitii test_iitii.cc)
add_dependencies(test_iitii catch)

include(CTest)
add_test(NAME unit_tests COMMAND ./test_iitii -d yes)
