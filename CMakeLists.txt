cmake_minimum_required(VERSION 3.17.2...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

find_package(LAPACK)
find_package(BLAS)
if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
  OUTPUT ac0_libmodule.c ac0_lib-f2pywrappers2.f90
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/pyscf/cas_ac0/accas_lib.f90"
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          "${CMAKE_CURRENT_SOURCE_DIR}/src/pyscf/cas_ac0/accas_lib.f90" -m ac0_lib --lower)

python_add_library(ac0_lib MODULE "${CMAKE_CURRENT_BINARY_DIR}/ac0_libmodule.c"
                    "${CMAKE_CURRENT_BINARY_DIR}/ac0_lib-f2pywrappers2.f90"
                    "${CMAKE_CURRENT_SOURCE_DIR}/src/pyscf/cas_ac0/accas_lib.f90" WITH_SOABI)
target_link_libraries(ac0_lib PRIVATE fortranobject ${lapackblas_libraries})

install(TARGETS ac0_lib DESTINATION pyscf/cas_ac0)
