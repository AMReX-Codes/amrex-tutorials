.. role:: cpp(code)
   :language: cpp


.. _name_multidim:


Namespaces and Multidimensional Code Tutorial
=============================================

.. admonition:: **Time to Complete**: 15 mins
   :class: warning

   **GOALS:**
     - Organize and Simplify Syntax with Namespaces
     - Modify a 3D Code for 2D or 3D Simulation Using Macros and Compiler
       Directives

In this tutorial we will start from the file ``main.cpp`` from
the :ref:`HeatEquation_Simple <guided_heat_simple>`
example and make several modifications to the code. First we will
discuss the ``amrex`` namespace, and use commands to
simplify the syntax of the code. Second, we will cover
AMReX macros and compiler directives commonly used to write
multidimensional code. We will insert these commands into our
example code to enable these features. In the end, we will have
a more capable code with cleaner syntax.


Namespaces
----------

AMReX uses namespaces to organize classes, functions, data, types
and templates and avoid naming conflicts with other libraries.
In the simplified Heat Equation guided tutorial, aspects of the ``amrex``
namespace were accessed by using the prefix ``amrex::``. For example,
to use the AMReX Real type we wrote,

.. code-block:: cpp

   amrex::Real dt;

Moreover, this pattern was repeated 34 times throughout the example
leading to a lot of extra typing! In this case, we can simplify the code
by adding the line,

.. code-block:: cpp

   using namespace amrex;

before the sections of code where we want to access the classes, functions,
data, types or templates from the ``amrex`` namespace. Typically, this
is done after the ``#include`` lines at the top of the file. Once this
is done, we can drop the ``amrex::`` prefix, and write only,

.. code-block:: cpp

   Real dt;

to get the same result.

At this point, we should add the ``using namespace amrex;``
line at the top of our file and remove the ``amrex::`` prefix from all the
commands that use it.

.. note::

   When using math functions such as ``min`` or ``max`` from the C++ standard library,
   it is recommend that users prefix these commands with ``std::``, to have ``std::min``
   or ``std::max`` to avoid any conflicts with the ``amrex`` namespace.

The rest of this tutorial will assume we implemented ``using namespace amrex;``
and removed the ``amrex::`` prefixes.


Writing Multidimensional Code
------------------------------

In this section of the tutorial we focus on using compiler or preprocessor
directives and AMReX macros to write code that will run
in 1-, 2- or 3-Dimensions, depending on what was chosen at compile
time. We continue modifications of the ``main.cpp`` file from
the :ref:`guided_heat_simple`.

What is a compiler directive?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A compiler directive allows a programmer to tell the compiler
to take specific actions at compile time. In C++ they are
often called preprocessor directives, because they are interpreted
before the compilation process. In AMReX they are
commonly used in conjunction with macros for writing 2D/3D code.

What is a Macro?
^^^^^^^^^^^^^^^^^

Macros are a form of text substitution. For example,

.. code-block:: cpp

   #define VALUE_OF_PI 3.14

will substitute 3.14 when the string VALUE_OF_PI is found in
the code. They can
also include more advanced functionality. As we'll
see in the tutorial, AMReX has several helpful macros. They
are named in all caps to avoid confusion with
normal variables.

In this section, we'll start from the ``main.cpp`` code used in
the HeatEquation_Simple example, and modify it for 2D/3D
compilation. We'll begin by adding several macros.


AMREX_D_DECL
~~~~~~~~~~~~

The first line we'll modify is

.. code-block:: cpp

   IntVect dom_lo(0,0,0);

to

.. code-block:: cpp

   IntVect dom_lo(AMREX_D_DECL(0,0,0));

The ``AMREX_D_DECL`` macro expands to a comma-separated list of
1, 2, or 3 of the arguments depending on the dimension selected at
compile time. To be explicit,

if compiled with ``DIM=2`` or ``AMReX_SPACEDIM=2`` the line above will
evaluate to,

.. code-block:: cpp

   IntVect dom_lo (0,0)

if compiled with ``DIM=3`` or ``AMReX_SPACEDIM=3`` it will
evaluate to,

.. code-block:: cpp

   IntVect dom_lo (0,0,0)


Next, modify the definitions of ``dom_hi`` and ``real_box`` to use
the ``AMREX_D_DECL`` macro in a similar manner.


AMREX_SPACEDIM
~~~~~~~~~~~~~~

When we arrive at the line

.. code-block:: cpp

   Array<int,3> is_periodic{1,1,1};

we encounter a slightly different situation. This time we need to
change the dimension of the Array as well as the number of inputs.
For this we change the 3 in ``Array<int,3>``, to ``AMREX_SPACEDIM``.
The inputs to ``is_preiodic`` are treated as above, giving:

.. code-block:: cpp

   Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

The ``AMREX_SPACEDIM`` macro in this statement will evaluate to 1, 2 or 3 depending on the
dimension selected at compile time.


Preprocessor Directives
-----------------------

While macros address many of the dimensional needs of our code
sometimes its necessary to use them in conjunction with
preprocessor directives, such as ``#if``, ``#elif``, and ``#endif``,
to allow for algorithmic differences for
different dimensions.  In our code example, this need arises within calls to ``ParallelFor``.

As a first step to writing a multidimensional version of this code, consider what the algorithm
looks like in 3D dimensions
(This is what we see in the code we're starting with.):

.. code-block:: cpp

           ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
           {
               Real x = (i+0.5) * dx[0];
               Real y = (j+0.5) * dx[1];
               Real z = (k+0.5) * dx[2];
               Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
               phiOld(i,j,k) = 1. + std::exp(-rsquared);
           });

If we wanted a similar initialization in 2D, it would be:

.. code-block:: cpp

           ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
           {
               Real x = (i+0.5) * dx[0];
               Real y = (j+0.5) * dx[1];
               Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01;
               phiOld(i,j,k) = 1. + std::exp(-rsquared);
           });


Notice that much of this code is redundant. The declarations of ``x``, ``y``
and value assigned to ``phiOld`` are all the same. The difference comes
with the addition of ``z`` and definition of the distance ``rsquared`` in 2D and 3D.
We can address this by adding preprocessor directives with AMReX macros
to create different sections of code for different numbers of dimensions.


Splitting the Code by Dimensional Dependence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we will separate the parts of the code by the number of dimensions
they exclusively pertain to. Because the ``x``, ``y`` and ``phiOld`` lines are
included in all cases, they will
remain outside the preprocessor directives. Then we can separate the code like this,

.. code-block:: cpp

   // included in all cases
   Real x = (i+0.5) * dx[0];
   Real y = (j+0.5) * dx[1];

   // dimensional dependent code

   // included in all cases
   phiOld(i,j,k) = 1. + std::exp(-rsquared);


Adding the 2-Dimensional Section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To address the 2-dimensional case, we add preprocessor directives ``#if`` and
``#endif`` with the logic, ``#if (AMREX_SPACEDIM == 2)``. This will
check if the value of the macro ``AMREX_SPACEDIM`` is equal to 2. If true,
it will compile the code inside this section. Therefore we write:

.. code-block:: cpp
   :emphasize-lines: 6-8

   // included in all cases
   Real x = (i+0.5) * dx[0];
   Real y = (j+0.5) * dx[1];

   // dimensional dependent code
   #if (AMREX_SPACEDIM == 2)
      Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01;
   #endif

   // included in all cases
   phiOld(i,j,k) = 1. + std::exp(-rsquared);


Adding the 3-dimensional Section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the above additions, the code will work if we compile with the option ``DIM=2`` or
``AMReX_SPACEDIM=2``. For three dimensions, we include
``#elif (AMREX_SPACEDIM == 3)`` and add the lines for ``z`` and the 3D version
of ``rsquared``:

.. code-block:: cpp
   :emphasize-lines: 7-10

   Real x = (i+0.5) * dx[0];
   Real y = (j+0.5) * dx[1];

   #if (AMREX_SPACEDIM == 2)
      Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01;

   #elif (AMREX_SPACEDIM == 3)
      Real z = (k+0.5) * dx[2];
      Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
   #endif

   phiOld(i,j,k) = 1. + std::exp(-rsquared);

This adds another preprocessor directive
which evaluates the statement ``AMREX_SPACEDIM==3``. If true
(and ``AMREX_SPACEDIM==2`` false), it will
compile this section of code and not the 2-dimensional section.


2D/3D Multidimensional Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Altogether the 2D/3D multidimensional version of the call to ``ParallelFor`` is:

.. code-block:: cpp

   ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
   {
      Real x = (i+0.5) * dx[0];
      Real y = (j+0.5) * dx[1];

   #if (AMREX_SPACEDIM == 2)
      Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01;

   #elif (AMREX_SPACEDIM == 3)
      Real z = (k+0.5) * dx[2];
      Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
   #endif

      phiOld(i,j,k) = 1. + std::exp(-rsquared);
   });

The section of code that will be compiled and executed is now determined
by the value of ``DIM`` or ``AMReX_SPACEDIM`` configured at compile
time. The final addition to this tutorial is to see we need to add
another preprocessor directive to modify the ``ParallelFor`` responsible
for advancing the data by dt. In this case all we need is a conditional
statement to isolate the code that updates the third dimension.


A Note About ParallelFor
------------------------

The ``ParallelFor`` function automatically optimizes code execution for
the hardware according to commands given at compile time. However, when writing
multidimensional code the syntax of the ``ParallelFor`` loop is **always
written as if compiling for three dimensions**. Consider,

.. code-block:: cpp

   ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){ ...  });

and note that we included the iterative variables ``i,j,k`` even though
the 2-dimensional code will not use the ``k`` variable. The ``ParallelFor``
function is aware of the dimension specified at compile time, and
automatically makes this adjustment for the convenience of AMReX users.


Conclusion
----------

Congratulations, you should now be able to compile the code for 2- or
3-dimensional simulations. At this point, our modified ``main.cpp`` code should
look similar to the code used in the :ref:`guided_heat`. The commands to compile
for each number of dimensions with GNU Make and CMake are listed in the
table below.

+---------------------+--------------+-----------------------------+
| Compile Commands    | GNU Make     | CMake                       |
+=====================+==============+=============================+
| 2 Dimensions        | make DIM=2   | cmake .. -DAMReX_SPACEDIM=2 |
+---------------------+--------------+-----------------------------+
| 3 Dimensions        | make DIM=3   | cmake .. -DAMReX_SPACEDIM=3 |
+---------------------+--------------+-----------------------------+


Please be aware that the plotfiles generated in each version of the code
will have different requirements due to the difference in dimensions. For example,
a plotfile from the 3D code will need, ``amrviz3d`` while the 2D code
will need ``amrviz2d``.

