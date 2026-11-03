# Two programs that exchange MultiFab data through amrex::MPMD. They are
# tested together below, not separately.
amrex_tutorial_add(MPMD_Main  SOURCES main.cpp  INPUTS inputs REQUIRES MPI NO_TEST)
amrex_tutorial_add(MPMD_Agent SOURCES agent.cpp INPUTS inputs REQUIRES MPI NO_TEST)

if(NOT TARGET Coupling_MPMD_Main)
    return()
endif()

configure_file(agent.py ${CMAKE_CURRENT_BINARY_DIR}/agent.py COPYONLY)

find_package(MPI REQUIRED COMPONENTS CXX)
set(_main ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:Coupling_MPMD_Main> ${MPIEXEC_POSTFLAGS} inputs)
list(APPEND _main ${AMReXTutorials_TEST_ARGS})
set(_agent_args inputs ${AMReXTutorials_TEST_ARGS})

# If a multi-program launch fails, run the tests with the PMIx "ompi" plugin
# disabled when that makes the launch work. The result is cached, so the
# launcher is only probed on the first configure.
if(NOT DEFINED CACHE{AMReXTutorials_MPMD_TEST_ENV})
    set(_probe ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
        ${CMAKE_COMMAND} -E true : ${MPIEXEC_NUMPROC_FLAG} 1 ${CMAKE_COMMAND} -E true)
    set(_env)
    execute_process(COMMAND ${_probe} RESULT_VARIABLE _probe_rc
                    OUTPUT_QUIET ERROR_QUIET TIMEOUT 30)
    if(NOT _probe_rc EQUAL 0)
        execute_process(COMMAND ${CMAKE_COMMAND} -E env PMIX_MCA_pmdl=^ompi ${_probe}
                        RESULT_VARIABLE _probe_rc OUTPUT_QUIET ERROR_QUIET TIMEOUT 30)
        if(_probe_rc EQUAL 0)
            message(STATUS "Coupling_MPMD: launcher needs PMIX_MCA_pmdl=^ompi for MPMD runs")
            set(_env PMIX_MCA_pmdl=^ompi)
        endif()
    endif()
    set(AMReXTutorials_MPMD_TEST_ENV "${_env}" CACHE INTERNAL
        "Environment for the MPMD tests")
endif()
set(_test_env ${AMReXTutorials_MPMD_TEST_ENV})

# mpiexec -n 2 MPMD_Main inputs : -n 2 MPMD_Agent inputs
add_test(NAME Coupling_MPMD
    COMMAND ${_main} : ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:Coupling_MPMD_Agent>
            ${MPIEXEC_POSTFLAGS} ${_agent_args}
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
            COMMAND ${_main} : ${MPIEXEC_NUMPROC_FLAG} 2 ${Python3_EXECUTABLE}
                    ${MPIEXEC_POSTFLAGS} agent.py inputs
                    dim=${AMReXTutorials_DIM} ${AMReXTutorials_GPU_TEST_ARGS}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        if(_test_env)
            set_tests_properties(Coupling_MPMD_Python PROPERTIES ENVIRONMENT "${_test_env}")
        endif()
    else()
        message(STATUS "Skipping Coupling_MPMD_Python: pyAMReX uses double precision")
    endif()
endif()
