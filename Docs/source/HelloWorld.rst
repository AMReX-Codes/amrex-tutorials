.. role:: cpp(code)
   :language: c++

.. _guided_hello_world:


Hello World
====================

The source code of this example is at ``amrex-tutorials/GuidedTutorials/Basic/HelloWorld_C/``
and is also shown below.

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

Building the Code
-----------------

You build the code in the ``amrex-tutorials/GuidedTutorials/Basic/HelloWorld_C/`` directory.
Typing ``make`` will start the compilation process and result in an executable
named ``main3d.gnu.DEBUG.ex``. The name shows that the GNU compiler with debug
options set by AMReX is used.  It also shows that the executable is built for
3D. Although this simple example code is dimension independent, dimensionality
does matter for all non-trivial examples. The build process can be adjusted by
modifying the ``amrex-tutorials/GuidedTutorials/Basic/HelloWorld_C/GNUmakefile`` file.  More
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

In the GNUmakefile there are compilation options for DEBUG mode (less optimized
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
