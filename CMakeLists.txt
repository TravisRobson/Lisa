
cmake_minimum_required( VERSION 3.14 )

project( Lisa VERSION 0.1.0 LANGUAGES C CXX )


option(ENABLE_LISA_OPENMP  "Use OpenMP"             OFF)
option(ENABLE_LISA_PROFILE "Enable LISA profiling"  ON)
option(ENABLE_LISA_TESTS   "Enable LISA tests"      OFF)
option(BUILD_SHARED_LIBS   "Build shared libraries" OFF) # Default, for Python ctypes interface.
option(LISA_USE_GSL        "Build with GSL"         ON)
option(BUILD_LEGACY_LISA   "Build legacy code"      ON) # \todo I will want to turn this off 
                                                        # when I put on GitHub for copyright reasons.

if (ENABLE_LISA_OPENMP)
 find_package(OpenMP REQUIRED)
endif()

if (LISA_USE_GSL)
  find_package(GSL REQUIRED)
endif()


if (BUILD_SHARED_LIBS)
  set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()


if (ENABLE_LISA_TESTS)
  include(CTest)
  enable_testing()
endif()


# specify C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# \todo create better way to manage compile flags
# \todo set C flags?
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow")
# -Werror -Wunused-parameter

add_library(lisa "")

add_subdirectory(src)


