cmake_minimum_required(VERSION 3.20)

project(blom-benchmarks VERSION 1.0.0 LANGUAGES C Fortran)

option(WITH_FLAGS_AVX2 "flag" ON)
option(WITH_FLAGS_NOPREC_DIV "flag" ON)
option(WITH_FLAGS_UNROLL "flag" ON)
option(WITH_FLAGS_PREFETCH "flag" ON)

option(ADVECT_EITVEL_OPT1 "enable opt1 within advect_eitvel" ON)
# removed option(ADVECT_EITVEL_OPT2 "enable opt2 within advect_eitvel" ON)
option(ADVECT_EITVEL_OPT3 "enable opt3 within advect_eitvel" ON)
option(WITH_ACC "compile with OpenACC" OFF)
#option(WITH_ACC_CPU "compile with OpenACC multicore" OFF)
# https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_%3CLANG%3E_COMPILER_ID
add_library(options_flags INTERFACE)
target_compile_definitions(options_flags INTERFACE
    $<$<BOOL:${ADVECT_EITVEL_OPT1}>:ADVECT_EITVEL_OPT1>
    $<$<BOOL:${ADVECT_EITVEL_OPT3}>:ADVECT_EITVEL_OPT3>
  )
target_compile_options(options_flags INTERFACE
  # or use $<COMPILE_LANG_AND_ID:language,compiler_ids>
  -g
  $<$<Fortran_COMPILER_ID:Intel>:
    -O3
    # https://developer.amd.com/wp-content/resources/Compiler%20Options%20Quick%20Ref%20Guide%20for%20AMD%20EPYC%207xx2%20Series%20Processors.pdf
    $<$<BOOL:${WITH_FLAGS_AVX2}>:     -march=core-avx2 >
    $<$<BOOL:${WITH_FLAGS_NOPREC_DIV}>:-no-prec-div >
    $<$<BOOL:${WITH_FLAGS_UNROLL}>:   -unroll-aggressive >
    $<$<BOOL:${WITH_FLAGS_PREFETCH}>: -qopt-prefetch >

    -qopenmp
    -qopt-report-phase=vec,loop -qopt-report=5
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
  $<$<Fortran_COMPILER_ID:Intel>:
    -qopenmp
    >
  $<$<Fortran_COMPILER_ID:PGI>:
    $<$<BOOL:${WITH_ACC}>:
      -acc -ta=tesla:cc80
    >
    $<$<BOOL:${WITH_ACC_CPU}>:
      -acc -ta=multicore
    >
  >
)

add_library(loop_utils OBJECT wallclock.c wallclock.F90 random.F90)

add_executable(advect_pbmin            $<TARGET_OBJECTS:loop_utils> advect_pbmin.F90)
add_executable(advect_eitvel           $<TARGET_OBJECTS:loop_utils> advect_eitvel.F90)
add_executable(advect_eitvel_opt1      $<TARGET_OBJECTS:loop_utils> advect_eitvel_opt1.F90)
add_executable(remap_zero              $<TARGET_OBJECTS:loop_utils> remap_zero.F90)
add_executable(remap_limited_grad      $<TARGET_OBJECTS:loop_utils> remap_limited_grad.F90)
add_executable(remap_vel_u             $<TARGET_OBJECTS:loop_utils> remap_vel_u.F90)

target_link_libraries(advect_pbmin       options_flags)
target_link_libraries(advect_eitvel      options_flags)
target_link_libraries(advect_eitvel_opt1 options_flags)
target_link_libraries(remap_zero         options_flags)
target_link_libraries(remap_limited_grad options_flags)
target_link_libraries(remap_vel_u        options_flags)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/Fypp.cmake)
set(FYPP_FLAGS )
if(WITH_ACC)
	set(FYPP_FLAGS ${FYPP_FLAGS} --define=LOOP_ACC)
endif()
PreProcessFyppFiles(
  advect_pbmin_opt1.fypp
)

add_executable(advect_pbmin_opt1            $<TARGET_OBJECTS:loop_utils> advect_pbmin_opt1.F90)
target_link_libraries(advect_pbmin_opt1     options_flags)
