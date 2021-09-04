cmake_minimum_required(VERSION 2.8)
project(IOWA_HILLS_FILTERKIT)

option(USE_DEBUG_ASAN      "use GCC's address sanitizer?" OFF)
option(DISABLE_LINK_WITH_M "Disables linking with m library to build with clangCL from MSVC" OFF)

# C90 requires the gcc extensions for function attributes like always_inline
# C99 provides the function attributes: no gcc extensions required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (USE_DEBUG_ASAN)
  set(ASANLIB "asan")
else()
  set(ASANLIB "")
endif()


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


######################################################

add_executable(example FilterKitMain.cpp FilterKitMain.h)
target_compile_definitions(example PRIVATE _USE_MATH_DEFINES)
if (USE_DEBUG_ASAN)
  target_compile_options(example PRIVATE "-fsanitize=address")
endif()

target_compile_options(example PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -Wextra -pedantic>)
target_link_libraries( example iowa_hills_dsp ${MATHLIB} ${ASANLIB} )

