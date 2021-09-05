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
  include/iowahills/limits.h
  include/iowahills/CplxDMath.hpp
  include/iowahills/windowing.h
  include/iowahills/fft.h
  include/iowahills/goertzel.h
  include/iowahills/fir.h
  include/iowahills/freq_sampling.h
  include/iowahills/iir.h
  include/iowahills/lowpass_prototypes.h
  include/iowahills/lowpass_roots.h
  include/iowahills/parks_mcclellan.h
  include/iowahills/PFiftyOneRevE.h
  include/iowahills/QuadRootsRevH.h
)

set( SRCS
  src/non_throwing_vector.h
  src/fft.cpp
  src/goertzel.cpp
  src/windowing.cpp
  src/fir.cpp
  src/freq_sampling.cpp
  src/iir.cpp
  src/lowpass_prototypes.cpp
  src/lowpass_roots.cpp
  src/parks_mcclellan.cpp
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
target_compile_options(iowa_hills_dsp PRIVATE $<$<C_COMPILER_ID:Clang>:-Wall -Wextra -pedantic>)
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
target_compile_options(example PRIVATE $<$<C_COMPILER_ID:Clang>:-Wall -Wextra -pedantic>)

target_link_libraries( example iowa_hills_dsp ${MATHLIB} ${ASANLIB} )


#######################################################
# tools
######################################################

add_executable(pmord tools/pmord_tool.cpp)
target_compile_definitions(pmord PRIVATE _USE_MATH_DEFINES)
if (USE_DEBUG_ASAN)
  target_compile_options(pmord PRIVATE "-fsanitize=address")
endif()
target_compile_options(pmord PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -Wextra -pedantic>)
target_compile_options(pmord PRIVATE $<$<C_COMPILER_ID:Clang>:-Wall -Wextra -pedantic>)
target_link_libraries( pmord iowa_hills_dsp ${MATHLIB} ${ASANLIB} )


#######################################################
# tests
######################################################

add_executable(test_fft tests/test_fft.cpp)
target_compile_definitions(test_fft PRIVATE _USE_MATH_DEFINES)
if (USE_DEBUG_ASAN)
  target_compile_options(test_fft PRIVATE "-fsanitize=address")
endif()
target_compile_options(test_fft PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -Wextra -pedantic>)
target_compile_options(test_fft PRIVATE $<$<C_COMPILER_ID:Clang>:-Wall -Wextra -pedantic>)
target_link_libraries( test_fft iowa_hills_dsp ${MATHLIB} ${ASANLIB} )


enable_testing()

add_test(NAME test_fft
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_fft"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)


