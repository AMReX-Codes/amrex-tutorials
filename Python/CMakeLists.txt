# pyAMReX tutorials are Python scripts. With -DAMReXTutorials_PYTHON=ON, they
# are run as CTest tests with the Python3 interpreter, which must be able to
# import amrex (pyAMReX), e.g., -DPython3_EXECUTABLE=/path/to/venv/bin/python.
if(NOT AMReXTutorials_PYTHON)
    return()
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter)

foreach(_name IN ITEMS MultiFab HeatEquation)
    set(_dir ${CMAKE_CURRENT_BINARY_DIR}/${_name})
    file(COPY ${_name}/main.py ${_name}/inputs DESTINATION ${_dir})
    set(_command ${Python3_EXECUTABLE} main.py inputs)
    if(AMReX_MPI)
        find_package(MPI REQUIRED COMPONENTS CXX)
        set(_command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} ${_command} ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(NAME Python_${_name} COMMAND ${_command} WORKING_DIRECTORY ${_dir})
endforeach()
