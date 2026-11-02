# Two programs that exchange MultiFab data through amrex::MPMD. They are
# tested together below, not separately.
amrex_tutorial_add(MPMD_Main  SOURCES main.cpp  INPUTS inputs REQUIRES MPI NO_TEST)
amrex_tutorial_add(MPMD_Agent SOURCES agent.cpp INPUTS inputs REQUIRES MPI NO_TEST)

if(NOT TARGET Coupling_MPMD_Main)
    return()
endif()

file(COPY agent.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_package(MPI REQUIRED COMPONENTS CXX)
set(_main ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:Coupling_MPMD_Main> inputs)
list(APPEND _main ${AMReXTutorials_TEST_ARGS})
set(_agent_args inputs ${AMReXTutorials_TEST_ARGS})

# Some launchers crash on a multi-program command line (Open MPI 5 with
# Ubuntu's PMIx 5.0.9 segfaults in PMIx_Setenv). Probe with a no-op program
# and, if disabling the PMIx "ompi" plugin avoids the crash, run the tests
# with it disabled.
set(_probe ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
    ${CMAKE_COMMAND} -E true : ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_COMMAND} -E true)
set(_test_env)
execute_process(COMMAND ${_probe} RESULT_VARIABLE _probe_rc
                OUTPUT_QUIET ERROR_QUIET TIMEOUT 30)
if(NOT _probe_rc EQUAL 0)
    execute_process(COMMAND ${CMAKE_COMMAND} -E env PMIX_MCA_pmdl=^ompi ${_probe}
                    RESULT_VARIABLE _probe_rc OUTPUT_QUIET ERROR_QUIET TIMEOUT 30)
    if(_probe_rc EQUAL 0)
        message(STATUS "Coupling_MPMD: launcher needs PMIX_MCA_pmdl=^ompi for MPMD runs")
        set(_test_env PMIX_MCA_pmdl=^ompi)
    endif()
endif()

# mpiexec -n 2 MPMD_Main inputs : -n 2 MPMD_Agent inputs
add_test(NAME Coupling_MPMD
    COMMAND ${_main} : ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:Coupling_MPMD_Agent> ${_agent_args}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(_test_env)
    set_tests_properties(Coupling_MPMD PROPERTIES ENVIRONMENT "${_test_env}")
endif()

# mpiexec -n 2 MPMD_Main inputs : -n 2 python3 agent.py inputs
# The agent exchanges raw data, so pyAMReX and AMReX must have the same precision.
if(AMReXTutorials_PYTHON)
    if(AMReX_PRECISION STREQUAL "DOUBLE")
        find_package(Python3 REQUIRED COMPONENTS Interpreter)
        add_test(NAME Coupling_MPMD_Python
            COMMAND ${_main} : ${MPIEXEC_NUMPROC_FLAG} 2 ${Python3_EXECUTABLE} agent.py inputs
                    dim=${AMReXTutorials_DIM}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        if(_test_env)
            set_tests_properties(Coupling_MPMD_Python PROPERTIES ENVIRONMENT "${_test_env}")
        endif()
    else()
        message(STATUS "Skipping Coupling_MPMD_Python: pyAMReX uses double precision")
    endif()
endif()
