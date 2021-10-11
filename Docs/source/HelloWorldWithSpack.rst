.. role:: cpp(code)
   :language: c++

.. _hello_world_with_spack:


Build HelloWorld with SPACK
===========================

.. admonition:: **Time to Complete**: 2 mins
   :class: warning

   **GOALS:**
     - Install AMReX with SPACK
     - Build with CMake
     - Run HelloWorld

This tutorial will walk through the steps involved for install AMReX with
SPACK. The source code of this example can be found  at ``amrex-tutorials/GuidedTutorials/HelloWorld/``
and is shown below.

SPACK
-----

SPACK is a package manager for HPC. To install SPACK on your system follow the
directions here..


::

   spack install amrex %gcc@9.3.0 +shared dimensions="2"

This command will tell SPACE to install AMReX with the shared library option enabled.


.. table:: Available build options

   +----------------+------------------------------------------+------------------+
   | Option         | Description                              | Values           |
   +================+==========================================+==================+
   | dimensions     | Dimensionality                           | 3\*, 2           |
   +----------------+------------------------------------------+------------------+
   | shared         | Build Shared Library                     | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | mpi            | Build with MPI support                   | True\*, False    |
   +----------------+------------------------------------------+------------------+
   | openmp         | Build with OpenMP support                | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | precision      | Real Precision                           | single, double\* |
   +----------------+------------------------------------------+------------------+
   | eb             | Build Embedded Boundary Classes          | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | fortran        | Build Fortran API                        | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | linear_solvers | Build Linear Solvers                     | True\*, False    |
   +----------------+------------------------------------------+------------------+
   | amrdata        | Build Data Services                      | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | particles      | Build Particle Classes                   | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | plotfile_tools | Build Plotfile Tools (i.e. fcompare)     | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | hdf5           | Enable HDF5-basedd I/O                   | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | hypre          | Enable Hypre Interfacces                 | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | petsc          | Enable PETSc Interfacces                 | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | sundials       | Enable SUNDIALS Interfacces              | True, False\*    |
   +----------------+------------------------------------------+------------------+
   | pic            | Enable PIC                               | True, False\*    |
   +----------------+------------------------------------------+------------------+



Build HelloWorld
----------------

::

   spack install amrex %gcc@9.3.0

This command will install AMReX compiled with gcc version 9.3.0.
To configure our build of ``HelloWorld`` we will
need to know where the files were installed. To do this, type:

::

  spack find -px amrex

Note the path of the install directory. At the end, add ``/lib/cmake/AMReX``.

Navigate to ``GuidedTutorials/HelloWorld/``. Then create and enter a build directory.

::

  mkdir build
  cd build


Next call CMake within the build directory. It's necessary to specific the
location of the configuration file to ensure the version of AMReX installed by
SPACK is used. If CMake cannot find the installed package, it will error or attempt
to download its own copy of AMReX from Git. To build with the install instance of
AMReX type the command:

::

  cmake .. -DAMReX_ROOT=<Spack Parent Dir>/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/amrex-20.11-6mdqeinhhi5ynvq6r6ywt5c77qlc6lfx/lib/cmake/AMReX

This will configure the files for building. Please note that the exact location
of your installation will differ. The example above is only illustrative.
To build type:

::

  cmake --build . -j2

This tells CMake to build the executable using the configuration
files in the current directory. If the build is successful
you will find a ``HelloWorld`` executable in the current directory.

To run ``HelloWorld`` with 4 MPI processes type:

::

  mpiexec -n 4 ./HelloWorld





To build the ExampleCodes:




