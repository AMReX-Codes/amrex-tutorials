amrex_tutorial_add(ExternalMPI SOURCES main.cpp INPUTS inputs REQUIRES MPI
                   TEST_ARGS inputs_ctest)

# The program only reads the inputs file given as the first argument, so the
# test arguments (floating-point traps) are appended to a copy of it.
if(TARGET Base_ExternalMPI)
    file(READ inputs _inputs)
    list(JOIN AMReXTutorials_TEST_ARGS "\n" _test_args)
    file(CONFIGURE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/inputs_ctest
         CONTENT "${_inputs}${_test_args}\n" @ONLY)
endif()
