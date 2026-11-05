# amrex_tutorial_add(<name>
#     SOURCES <files>...
#     [INPUTS <files>...]          # copied to the build directory
#     [REQUIRES <features>...]     # EB FFT PARTICLES LSOLVERS LSOLVERS_EM LSOLVERS_INCFLO
#                                  # AMRLEVEL MPI GPU NOGPU SUNDIALS
#     [DIMS <dims>...]             # supported AMReX_SPACEDIM values, default 1 2 3
#     [TEST_ARGS <args>...]        # default: the first INPUTS file
#     [NPROCS <n>]                 # MPI ranks for the test, default 2
#     [NO_TEST])
#
# Creates the executable target <Topic>_<name>, where <Topic> is the name of
# the parent directory, and registers a CTest test that runs it.  The tutorial
# is skipped if a required feature or the current dimension is unavailable.

# The tests do not trap floating-point exceptions (amrex.fpe_trap_*):
# vectorized code raises them in inactive SIMD lanes, see Base/SIMD.

# The MPI launcher for the tests (MPIEXEC_EXECUTABLE and friends)
if(AMReX_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

function(_amrex_tutorial_has_feature _feature _out)
    set(_has OFF)
    if(_feature STREQUAL "EB")
        set(_has ${AMReX_EB})
    elseif(_feature STREQUAL "FFT")
        set(_has ${AMReX_FFT})
    elseif(_feature STREQUAL "PARTICLES")
        set(_has ${AMReX_PARTICLES})
    elseif(_feature STREQUAL "LSOLVERS")
        set(_has ${AMReX_LINEAR_SOLVERS})
    elseif(_feature STREQUAL "LSOLVERS_EM")
        set(_has ${AMReX_LINEAR_SOLVERS_EM})
    elseif(_feature STREQUAL "LSOLVERS_INCFLO")
        set(_has ${AMReX_LINEAR_SOLVERS_INCFLO})
    elseif(_feature STREQUAL "AMRLEVEL")
        set(_has ${AMReX_AMRLEVEL})
    elseif(_feature STREQUAL "MPI")
        set(_has ${AMReX_MPI})
    elseif(_feature STREQUAL "SUNDIALS")
        set(_has ${AMReX_SUNDIALS})
    elseif(_feature STREQUAL "GPU")
        if(NOT AMReX_GPU_BACKEND STREQUAL "NONE")
            set(_has ON)
        endif()
    elseif(_feature STREQUAL "NOGPU")
        if(AMReX_GPU_BACKEND STREQUAL "NONE")
            set(_has ON)
        endif()
    else()
        message(FATAL_ERROR "amrex_tutorial_add: unknown feature '${_feature}'")
    endif()
    if(_has)
        set(${_out} ON PARENT_SCOPE)
    else()
        set(${_out} OFF PARENT_SCOPE)
    endif()
endfunction()

function(amrex_tutorial_add _name)
    cmake_parse_arguments(PARSE_ARGV 1 _arg
        "NO_TEST" "NPROCS" "SOURCES;INPUTS;REQUIRES;DIMS;TEST_ARGS")

    cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH _parent)
    cmake_path(GET _parent FILENAME _topic)
    set(_target ${_topic}_${_name})

    if(NOT _arg_DIMS)
        set(_arg_DIMS 1 2 3)
    endif()
    if(NOT AMReXTutorials_DIM IN_LIST _arg_DIMS)
        message(STATUS "Skipping ${_target}: not available in ${AMReXTutorials_DIM}D")
        return()
    endif()
    foreach(_feature IN LISTS _arg_REQUIRES)
        _amrex_tutorial_has_feature(${_feature} _has)
        if(NOT _has)
            message(STATUS "Skipping ${_target}: requires ${_feature}")
            return()
        endif()
    endforeach()

    add_executable(${_target} ${_arg_SOURCES})
    target_include_directories(${_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(${_target} PRIVATE AMReX::amrex_${AMReXTutorials_DIM}d)
    set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${_name})
    if(AMReX_GPU_BACKEND STREQUAL "CUDA")
        setup_target_for_cuda_compilation(${_target})
    endif()

    # configure_file re-copies an input file when it changes.
    foreach(_input IN LISTS _arg_INPUTS)
        cmake_path(GET _input FILENAME _input_name)
        configure_file(${_input} ${CMAKE_CURRENT_BINARY_DIR}/${_input_name} COPYONLY)
    endforeach()

    if(_arg_NO_TEST)
        return()
    endif()

    if(DEFINED _arg_TEST_ARGS)
        set(_test_args ${_arg_TEST_ARGS})
    elseif(_arg_INPUTS)
        list(GET _arg_INPUTS 0 _first_input)
        cmake_path(GET _first_input FILENAME _test_args)
    endif()

    set(_command $<TARGET_FILE:${_target}> ${_test_args})
    if(AMReX_MPI)
        if(NOT _arg_NPROCS)
            set(_arg_NPROCS 2)
        endif()
        set(_command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_arg_NPROCS}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${_target}> ${MPIEXEC_POSTFLAGS} ${_test_args})
    endif()
    add_test(NAME ${_target} COMMAND ${_command}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()
