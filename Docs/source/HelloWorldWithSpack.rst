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




.. highlight:: c++

::

     #include <AMReX.H>
     #include <AMReX_Print.H>

     int main(int argc, char* argv[])
     {
         amrex::Initialize(argc,argv);
         {
             amrex::Print() << "Hello world from AMReX version "
                            << amrex::Version() << "\n";
         }
         amrex::Finalize();
     }

The main body of this short example contains three statements.  Usually the
first and last statements for the :cpp:`int main(...)` function of every
program should be calling :cpp:`amrex::Initialize` and :cpp:`amrex::Finalize`,
respectively. The second statement calls :cpp:`amrex::Print` to print out a
string that includes the AMReX version returned by the :cpp:`amrex::Version`
function. Finally, the third statement calls :cpp:`amrex::Finalize` to clean up
data structures that are necessary for proper AMReX operation.

Notice the braces placed between :cpp:`amrex::Initialize` and
:cpp:`amrex::Finalize`. It is considered a good programming practice to insert
these braces such that it is guaranteed that anything executed in the code is
done after AMReX has been initialized, and before AMReX is finalized.

The example code includes two AMReX header files. Note that the name
of all AMReX header files starts with ``AMReX_`` (or just AMReX in the case of
AMReX.H). All AMReX C++ functions are in the :cpp:`amrex` namespace.

Building the Code with GNU Make
-------------------------------

You build the code in the ``amrex-tutorials/GuidedTutorials/HelloWorld/`` directory.
Typing ``make`` will start the compilation process and result in an executable
named ``main3d.gnu.DEBUG.ex``. The name shows that the GNU compiler with debug
options set by AMReX is used.  It also shows that the executable is built for
3D. Although this simple example code is dimension independent, dimensionality
does matter for all non-trivial examples. The build process can be adjusted by
modifying the ``amrex-tutorials/GuidedTutorials/HelloWorld/GNUmakefile`` file.  More
details on how to build AMReX can be found in :ref:`Chap:BuildingAMReX`.

Running the Code
----------------

The example code can be run as follows,

.. highlight:: console

::

      ./main3d.gnu.DEBUG.ex

The result may look like,

.. highlight:: console

::

      AMReX (17.05-30-g5775aed933c4-dirty) initialized
      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty
      AMReX (17.05-30-g5775aed933c4-dirty) finalized

The version string means the current commit 5775aed933c4 (note that the first
letter g in g577.. is not part of the hash) is based on 17.05 with 30
additional commits and the AMReX work tree is dirty (i.e. there are uncommitted
changes).

In the GNU make file, ``GNUmakefile``,  there are compilation options for DEBUG mode (less optimized
code with more error checking), dimensionality, compiler type, and flags to
enable MPI and/or OpenMP parallelism.  If there are multiple instances of a
parameter, the last instance takes precedence.

Parallelization
---------------

Now let's build with MPI by typing ``make USE_MPI=TRUE`` (alternatively you can
set ``USE_MPI=TRUE`` in the GNUmakefile). This should make an executable named
``main3d.gnu.DEBUG.MPI.ex``. Note MPI in the file name. You can then run,

.. highlight:: console

::

      mpiexec -n 4 ./main3d.gnu.DEBUG.MPI.ex amrex.v=1

The result may look like,

.. highlight:: console

::

      MPI initialized with 4 MPI processes
      AMReX (17.05-30-g5775aed933c4-dirty) initialized
      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty
      AMReX (17.05-30-g5775aed933c4-dirty) finalized

If the compilation fails, you are referred to :ref:`Chap:BuildingAMReX` for
more details on how to configure the build system.  The *optional* command line
argument ``amrex.v=1`` sets the AMReX verbosity level
to 1 to print the number of MPI processes used.  The default verbosity
level is 1, and you can pass ``amrex.v=0`` to turn it off.
More details on how runtime parameters are handled can be found in
section :ref:`sec:basics:parmparse`.

If you want to build with OpenMP, type make ``USE_OMP=TRUE``.  This should make
an executable named ``main3d.gnu.DEBUG.OMP.ex``. Note OMP in the file name.
Make sure the ``OMP_NUM_THREADS`` environment variable is set on your system.
You can then run,

.. highlight:: console

::

      OMP_NUM_THREADS=4 ./main3d.gnu.DEBUG.OMP.ex

The result may look like,

.. highlight:: console

::

      OMP initialized with 4 OMP threads
      AMReX (17.05-30-g5775aed933c4-dirty) initialized
      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty
      AMReX (17.05-30-g5775aed933c4-dirty) finalized

Note that you can build with both ``USE_MPI=TRUE`` and ``USE_OMP=TRUE``.  You
can then run,

.. highlight:: console

::

      OMP_NUM_THREADS=4 mpiexec -n 2 ./main3d.gnu.DEBUG.MPI.OMP.ex

The result may look like,

.. highlight:: console

::

      MPI initialized with 2 MPI processes
      OMP initialized with 4 OMP threads
      AMReX (17.05-30-g5775aed933c4-dirty) initialized
      Hello world from AMReX version 17.05-30-g5775aed933c4-dirty
      AMReX (17.05-30-g5775aed933c4-dirty) finalized
