.. _HYPRE:

HYPRE for ABecLaplacian_C Tutorial
==================================

The following directions explain how to configure, build and run
the ``ABecLaplacian_C`` example using the external solver library HYPRE.
HYPRE is a library of scalable linear solvers and multigrid methods. For
more information about HYPRE see their website_. More information on
setting HYPRE options within AMReX is available in the User's Guide in
the `External Solvers`_ section.

.. _website: https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods

.. _`External Solvers`: https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#external-solvers

Building with HYPRE via GNUMake
-------------------------------

#. Clone, configure, make install and set the ``HYPRE_DIR`` environment variable:

   .. code-block:: bash

      git clone https://github.com/hypre-space/hypre.git
      cd hypre/src
      ./configure
      make install
      export HYPRE_DIR=/path_to_hypre_dir/hypre/src/hypre

   .. note::

      If HYPRE fails the configure step, it may be necessary to manually specify
      several options. This can be done during the configure step. For example,
      one might replace :code:`./configure` with,

      .. code-block:: bash

         ./configure CXX=CC CC=cc FC=ftn --with-MPI


#. Navigate to the location of the ABecLaplacian_C example,
   ``amrex-tutorials/ExampleCodes/LinearSolvers/ABecLaplacian_C``.

#. Make the executable by typing,

   .. code-block:: bash

      make -j4 USE_HYPRE=TRUE

   This will compile the code with HYPRE enabled using 4 processes.

#. To run the HYPRE enabled example, use the ``inputs.hypre`` file
   as input.

   .. code-block:: bash

      ./ABecLaplacian_C inputs.hypre

#. To verify the HYPRE solver has been used, look for the line,

   .. code-block:: literal

      HYPRE BoomerAMG: Num. iterations = 4; Relative residual = 3.92236062e-05

   in the output. For this example, the code should quickly reach the end with
   only a handful of iterations at each step. The code successfully completes
   with the AMReX finalized output.


Building with HYPRE via CMake
-----------------------------

#. Follow the directions above for cloning, configuring, make installing
   and setting the environment variable.

#. Next navigate to the ExampleCodes parent directory,
   ``/amrex-tutorials/ExampleCodes/``. CMake is currently
   configured to compile multiple examples at once using the file ``CMakeLists.txt``.
   Create and enter a build directory to hold the configuration and compiled files.

#. From the build directory, call CMake with the following configuration:

   .. code-block:: bash

      cmake .. -DAMReX_HYPRE=ON \
      -DHYPRE_LIBRARIES=${HYPRE_DIR}\lib\libHYPRE.a \
      -DHYPRE_INCLUDE_DIRS=${HYPRE_DIR}\include \
      -DAMReX_LINEAR_SOLVERS=ON \
      -DAMReX_FORTRAN=ON

   The first option tells CMake to enable the HYPRE interface in AMReX. The
   following two options specify the location of the HYPRE library we
   compiled in the first step, and the directory of HYPRE's required files
   that we also installed in the first step. In this example, the
   environment variable ``HYPRE_DIR`` is used to replace writing out the
   entire path. The next flag enables compilation for the ``LinearSolver`` example codes,
   some of which require Fortran, and thus require the last flag.

#. Next we can build the executable with,

   .. code-block:: bash

      cmake --build . -j8

   This will tell CMake to use 8 processes to compile the source files.

#. Finally we can run the executable by navigating to the
   ``build/LinearSolvers/ABecLaplacian_C`` folder inside our build directory, and typing
   the name of the executable followed by the inputs file, ``inputs.hypre``.

   .. code-block:: bash

      ./ABecLaplacian_C inputs.hypre

#. To verify the HYPRE solver has been used, look for the line,

   .. code-block:: literal

      HYPRE BoomerAMG: Num. iterations = 4; Relative residual = 3.92236062e-05

   in the output. For this example, the code should quickly reach the end with
   only a handful of iterations at each step. The code successfully completes
   with the AMReX finalized output.

AMReX with HYPRE via Spack
--------------------------


#. Using Spack, install AMReX with HYPRE and Fortran
   options selected.

   .. code-block:: bash

      spack install amrex +hypre +fortran

#. Load the desired version of AMReX.

   .. code-block:: bash

      spack load amrex +hypre +fortran

#. Identify the location of the installed version of AMReX. Because the location is
   usually quite long, we will store the result from Spack as the shell variable,
   ``AMREX_DIR``.

   .. code-block:: bash

      AMREX_DIR=$(spack location -i amrex +hypre +fortran)

#. In this example we will build the ``ABecLaplacian_C`` example code from
   the linear solvers in ``amrex-tutorials``. First navigate to the ``ExampleCodes``
   directory. Then create a build folder to store the compiled files. Inside
   this folder we'll use CMake to compile the code.

   .. code-block:: bash

      cmake .. -DAMReX_DIR=${AMREX_DIR} \
               -DAMReX_HYPRE=ON \
               -DAMReX_FORTRAN=ON \
               -DAMReX_FORTRAN_INTERFACES=ON \
               -DAMREX_LINEAR_SOLVERS=ON

   These configuration commands do the following:

      - AMReX_DIR: Tells CMake where to find the installed version of
        AMReX. If this is not supplied, CMake may be unable to locate
        the AMReX files or it may download the file from the latest release
        from GitHub.

      - AMReX_HYPRE: Enables AMReX to use HYPRE.

      - AMReX_FORTRAN: Enables Fortran for AMReX.

      - AMReX_FORTRAN_INTERFACES: Enables the Fortran API.

      - AMReX_LINEAR_SOLVERS: This command is specific to the
        ``ExamplesCodes`` install configuration, i.e. CMakeLists.txt. It tells CMake
        to compile all the linear solver examples.

#. After setting up the configuration, we build the executables with
   CMake. This command will build the files according to the configuration
   in the current directory using 4 processes (``-j4``).

   .. code-block:: bash

      cmake --build . -j4

#. To run the HYPRE example navigate to the folder,
   ``path_to_base_dir/ExampleCodes/build/LinearSolvers/ABecLaplacian_C``
   and call the executable with the ``inputs.hypre`` file as input.

   .. code-block:: bash

      ./ABecLaplacian_C inputs.hypre

#. To verify the HYPRE solver has been used, look for the line,

   .. code-block:: literal

      HYPRE BoomerAMG: Num. iterations = 4; Relative residual = 3.92236062e-05

   in the output. For this example, the code should quickly reach the end with
   only a handful of iterations at each step. The code successfully completes
   with the AMReX finalized output.



