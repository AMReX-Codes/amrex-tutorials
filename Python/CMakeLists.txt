# pyAMReX tutorials are Python scripts. With -DAMReXTutorials_PYTHON=ON, they
# are run as CTest tests with the Python3 interpreter, which must be able to
# import amrex (pyAMReX), e.g., -DPython3_EXECUTABLE=/path/to/venv/bin/python.
if(NOT AMReXTutorials_PYTHON)
    return()
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# The tests get the GPU arguments but no floating-point traps, which NumPy
# does not expect.
foreach(_name IN ITEMS MultiFab HeatEquation)
    set(_dir ${CMAKE_CURRENT_BINARY_DIR}/${_name})
    # configure_file copies the files again when they change.
    foreach(_file IN ITEMS main.py inputs)
        configure_file(${_name}/${_file} ${_dir}/${_file} COPYONLY)
    endforeach()
    set(_command ${Python3_EXECUTABLE} main.py inputs ${AMReXTutorials_GPU_TEST_ARGS})
    if(AMReX_MPI)
        find_package(MPI REQUIRED COMPONENTS CXX)
        set(_command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} ${_command} ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(NAME Python_${_name} COMMAND ${_command} WORKING_DIRECTORY ${_dir})
endforeach()
