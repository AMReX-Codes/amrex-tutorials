.. role:: cpp(code)
   :language: c++

.. _hello_world_with_spack:


Build HelloWorld with Spack
===========================

.. admonition:: **Time to Complete**: 2 mins
   :class: warning

   **GOALS:**
     - Install AMReX with Spack
     - Build with CMake
     - Run HelloWorld

This tutorial will walk through the steps involved for install AMReX with
Spack. The source code of this example can be found  at ``amrex-tutorials/GuidedTutorials/HelloWorld/``
and is shown below.

Spack
-----

Spack is a package manager for HPC. To install Spack on your system follow the
directions:   

Once Spack is installed,  the Spack environment

Install Spack: 
Integrate Spack into your shell:

After these steps are completed, you can type spack commands at the shell prompt. 
For example, the command,

::

   spack install amrex %gcc@9.3.0 +shared dimensions="2"

will tell Spack to install AMReX with the shared library option enabled, 
and AMReX_SPACEDIM=2.


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

Step 1: Install AMReX via Spack
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

   spack install amrex %gcc@9.3.0

This command will install AMReX compiled with gcc version 9.3.0.
To configure our build of ``HelloWorld`` we will
need to know where the files were installed. To do this, type:

::

  spack location -i amrex

Note the path of the install directory. At the end, add ``/lib/cmake/AMReX``.


.. note::
   Spack's install paths can be long. A nifty way to deal with this ist to store
   the install directory as a variable. To do this, use the command: 
   ::

      AMREX_INSTALLDIR=$(spack location -i amrex)

   The path stored in this variable can then be referrenced later in the install
   process.


Step 2
~~~~~~

Navigate to ``GuidedTutorials/HelloWorld/``. Then create and enter a build directory.

::

  mkdir build
  cd build


Next call CMake within the build directory. It's necessary to specific the
location of the configuration file to ensure the version of AMReX installed by
Spack is used. If CMake cannot find the installed package, it will error or attempt
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




