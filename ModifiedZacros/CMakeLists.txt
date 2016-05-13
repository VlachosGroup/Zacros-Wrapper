CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
project ("Graph Theoretical Kinetic Monte-Carlo" Fortran)
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
option(doopenmp    "Will compile with openmp, if found." on)


set(F90_SOURCES compat.f90 constants_module.f90 energetics_handle_module.f90
                energetics_setup_module.f90 error_module.f90
                gauss_quadrature_module.f90 graph_functions_module.f90
                heap_functions_module.f90 
                kmc_simulation_handle_module.f90 lattice_handle_module.f90
                lattice_setup_module.f90 mechanism_setup_module.f90
                mt19937_module.f90 output_handle_module.f90
                parser_module.f90 random_deviates_module.f90
                rates_handle_module.f90 restarts_module.f90
                sampling_handle_module.f90 simulation_setup_module.f90
                state_setup_module.f90)

   
set(FREELINE_SOURCES error_module.f90 simulation_setup_module.f90
                     energetics_setup_module.f90 mechanism_setup_module.f90
                     lattice_setup_module.f90 state_setup_module.f90
                     output_handle_module.f90 restarts_module.f90)

enable_language(Fortran)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  SET_SOURCE_FILES_PROPERTIES( ${FREELINE_SOURCES} PROPERTIES COMPILE_FLAGS
    -ffree-line-length-0 )
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  set (CMAKE_Fortran_FLAGS "-kind=byte -ieee=full")
endif()


if(doopenmp)
  find_package(OpenMP_Fortran)
  if(OpenMP_Fortran_FLAGS)
    set(OPENMP_SOURCES rates_handle_module.f90 zacros_main.f90
                       lattice_handle_module.f90 kmc_simulation_handle_module.f90)
    set_source_files_properties( ${OPENMP_SOURCES} 
                                 PROPERTIES 
                                 COMPILE_FLAGS ${OpenMP_Fortran_FLAGS} )     
  else(OpenMP_Fortran_FLAGS)
    message("No OpenMP flag found.")
    set(doopenmp FALSE)
  endif(OpenMP_Fortran_FLAGS)
endif(doopenmp)


add_library(zacros STATIC ${F90_SOURCES})

add_executable(zacros.x zacros_main.f90)
target_link_libraries(zacros.x zacros)
if(doopenmp)
  set_target_properties(zacros.x PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS})
endif()

find_package(PythonInterp)

if(${PYTHONINTERP_FOUND} STREQUAL "FALSE")
  message("No python intepreter found. Cannot run tests.")
  set(dotests FALSE)
else(${PYTHONINTERP_FOUND} STREQUAL "FALSE")
  set(dotests TRUE)
  ENABLE_TESTING()
  add_subdirectory(tests)
endif(${PYTHONINTERP_FOUND} STREQUAL "FALSE")

