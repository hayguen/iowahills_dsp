cmake_minimum_required(VERSION 2.8)
project(IOWA_HILLS_FILTERS)

# installs Iowa Hills filters to CMAKE_INSTALL_PREFIX

option(USE_DEBUG_ASAN      "use GCC's address sanitizer?" OFF)
option(DISABLE_LINK_WITH_M "Disables linking with m library to build with clangCL from MSVC" OFF)

# C90 requires the gcc extensions for function attributes like always_inline
# C99 provides the function attributes: no gcc extensions required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# populate what to install
set(INSTALL_TARGETS "")
set(INSTALL_HEADERS "")


if (USE_DEBUG_ASAN)
  set(ASANLIB "asan")
else()
  set(ASANLIB "")
endif()


########################################################################
# select the release build type by default to get optimization flags
########################################################################
if ( NOT CMAKE_BUILD_TYPE )
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)

if ( CMAKE_C_COMPILER_ID MATCHES "MSVC" )
  # using Visual Studio C++
  message(STATUS "INFO: detected MSVC: will not link math lib m")
  set(MATHLIB "")
else()
  if(DISABLE_LINK_WITH_M)
  else()
    message(STATUS "INFO: detected NO MSVC: ${CMAKE_C_COMPILER_ID}: will link math lib m")
    set(MATHLIB "m")
  endif()
endif()

set(STDCXXLIB "")
if (MINGW)
  set(STDCXXLIB "stdc++")
endif()


set( HDRS
  include/iowahills/CplxDMath.hpp
  include/iowahills/FFTCode.h
  include/iowahills/FIRFilterCode.h
  include/iowahills/FreqSamplingCode.h
  include/iowahills/IIRFilterCode.h
  include/iowahills/LowPassPrototypes.h
  include/iowahills/LowPassRoots.h
  include/iowahills/NewParksMcClellan.h
  include/iowahills/PFiftyOneRevE.h
  include/iowahills/QuadRootsRevH.h
)

set( SRCS
  src/FFTCode.cpp
  src/FIRFilterCode.cpp
  src/FreqSamplingCode.cpp
  src/IIRFilterCode.cpp
  src/LowPassPrototypes.cpp
  src/LowPassRoots.cpp
  src/NewParksMcClellan.cpp
  src/PFiftyOneRevE.cpp
  src/QuadRootsRevH.cpp
  ${HDRS}
)

set(INSTALL_HEADERS ${INSTALL_HEADERS} ${HDRS})

######################################################

add_library(iowa_hills_dsp STATIC ${SRCS})
target_compile_definitions(iowa_hills_dsp PRIVATE _USE_MATH_DEFINES)
if (USE_DEBUG_ASAN)
  target_compile_options(iowa_hills_dsp PRIVATE "-fsanitize=address")
endif()

target_compile_options(iowa_hills_dsp PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -Wextra -pedantic>)
target_link_libraries( iowa_hills_dsp ${MATHLIB} )

target_include_directories(iowa_hills_dsp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(INSTALL_TARGETS ${INSTALL_TARGETS} iowa_hills_dsp)

install( TARGETS ${INSTALL_TARGETS}  DESTINATION lib)
install( FILES   ${INSTALL_HEADERS}  DESTINATION include/iowahills)

add_custom_target(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_SOURCE_DIR}/uninstall.cmake"
)


#######################################################
# add_subdirectory(examples) ? -> No!
# Have subdirectory as standalone test - without having this CMakeLists.txt
######################################################

add_executable(example examples/FilterKitMain.cpp examples/FilterKitMain.h)
target_compile_definitions(example PRIVATE _USE_MATH_DEFINES)
if (USE_DEBUG_ASAN)
  target_compile_options(example PRIVATE "-fsanitize=address")
endif()
target_compile_options(example PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -Wextra -pedantic>)
target_link_libraries( example iowa_hills_dsp ${MATHLIB} ${ASANLIB} )
