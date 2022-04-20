cmake_minimum_required(VERSION 3.20)

project(blom-benchmarks VERSION 1.0.0 LANGUAGES C Fortran)

#option(WITH_ACC "compile with OpenACC" OFF)
#option(WITH_ACC_CPU "compile with OpenACC multicore" OFF)
# https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_%3CLANG%3E_COMPILER_ID
add_library(options_flags INTERFACE)
target_compile_options(options_flags INTERFACE
  # or use $<COMPILE_LANG_AND_ID:language,compiler_ids>
  $<$<Fortran_COMPILER_ID:Intel>:
    -O3
   >

   $<$<Fortran_COMPILER_ID:GNU>:
    -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none -ffree-line-length-none
    -ftree-vectorize
    -ftree-vectorizer-verbose=5
    $<$<VERSION_GREATER_EQUAL:${CMAKE_Fortran_COMPILER_VERSION},10>:
      -fallow-argument-mismatch # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
     >
   >

  $<$<Fortran_COMPILER_ID:PGI>:
  -fast -fastsse -O3 -Minline -Munroll
    $<$<BOOL:${WITH_ACC}>:
      -acc -ta=tesla:cc80
      -Minfo=accel
     >
    $<$<BOOL:${WITH_ACC_CPU}>:
      -acc -ta=multicore
      -Minfo=accel
     >
   >
  )

target_link_options(options_flags INTERFACE
  $<$<Fortran_COMPILER_ID:PGI>:
    $<$<BOOL:${WITH_ACC}>:
      -acc -ta=tesla:cc80
    >
    $<$<BOOL:${WITH_ACC_CPU}>:
      -acc -ta=multicore
    >
  >
)

add_library(loop_utils OBJECT wallclock.c wallclock.F90)

add_executable(advect_pbmin            $<TARGET_OBJECTS:loop_utils> advect_pbmin.F90)
add_executable(advect_eitvel           $<TARGET_OBJECTS:loop_utils> advect_eitvel.F90)

target_link_libraries(advect_pbmin     options_flags)
target_link_libraries(advect_eitvel    options_flags)

