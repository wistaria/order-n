cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(order_n NONE)

set(ALPS_HOME_ENV $ENV{ALPS_HOME})
if(ALPS_ROOT_DIR OR ALPS_HOME_ENV)
  find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME}
    NO_SYSTEM_ENVIRONMENT_PATH)
  message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
  include(${ALPS_USE_FILE})
else(ALPS_ROOT_DIR OR ALPS_HOME_ENV)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
  find_package(Boost REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
  include(add_alps_test)
endif(ALPS_ROOT_DIR OR ALPS_HOME_ENV)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# enable C and C++ compilers
enable_language(C CXX)

enable_testing()
include(CTest)

include_directories(include)
add_subdirectory(infinite-pmj)
