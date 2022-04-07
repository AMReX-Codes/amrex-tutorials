.. role:: cpp(code)
   :language: cpp


.. _multifab_tutorial:


MultiFab Tutorial
=================

.. admonition:: **Time to Complete**: 15 mins
   :class: warning

   **GOALS:**
     - Create a MultiFab
     - Write Data to a MultiFab
     - Plot MultiFab Data

In this tutorial we focus on using the MultiFab data structure.
We will learn the steps to create a MultiFab, load it with data, and plot
it. At the end, we make some comments about writing code that works
for both 2 and 3 dimensions.

.. image:: ./images_tutorial/MultiFabTutorialVideo.png
   :height: 120px
   :align: right
   :alt: Image of Youtube video
   :target: https://youtu.be/498VdW2cNB8

A video companion to this tutorial is available on YouTube. You
can access it by clicking on the thumbnail to the right.

|

What is a MultiFab?
~~~~~~~~~~~~~~~~~~~

A MultiFab is a C++ class in AMReX the stores and operates on
multidimensional arrays in parallel. It contains:

   - Grid information in the form of a BoxArray that contains one or more
     components (scalar values) for a single level of the mesh.
   - A distribution map, that allows for parallel processing of data
     in the MultiFab.
   - Ghost cells that facilitate a variety of mesh refinement, boundary
     conditions, and particle algorithms.

Defining a MultiFab
~~~~~~~~~~~~~~~~~~~

To use the MultiFab type we first need to include the header file,

.. code-block:: cpp

   #include <AMReX_MultiFab.H>


To define a MultiFab,

.. code-block:: cpp

   amrex::MultiFab mf(ba, dm, ncomp, ngrow);

Defining a MultiFab in this way requires 4 inputs:

  - ``ba``, a BoxArray.
  - ``dm``, a DistributionMapping.
  - ``ncomp``, the number of components (scalar values) to store in the MultiFab.
  - ``ngrow``, the number of layers of ghost cells.

The value of ``ncomp`` is straight-forward and can be specified when
the MultiFab is defined. The correct number of ghost cells, ``ngrow``, is important
for many of the algorithms in AMReX such as boundary behavior, mesh refinement, etc.
For this reason it should be carefully selected based on application requirements.
The DistributionMapping, ``dm``,  directs parts of the domain to different processes
for parallel computation. In this example, the default suffices and can be easily
defined once the BoxArray is configured. Defining the BoxArray requires a few steps
that are presented in the next section.

Number of Components and Ghost Cells
------------------------------------

The value of ``ncomp`` determines how many scalar values to store in the MultiFab.
For example, if we want to store the value of phi at every point on the grid,
thats one component for the MultiFab. Then we would set:

.. code-block:: cpp

   int ncomp = 1;

The number of layers of ghost cells around the boundary of the domain is determined
by the value of ``ngrow`` and will depend on your application. In this tutorial
we will not be doing any operations that require ghost cells,
so we will set the value to 0.

.. code-block:: cpp

   int ngrow = 0;

BoxArray: Abstract Domain Setup
-------------------------------

The BoxArray contains a list of boxes that cover the domain. To define a
BoxArray, we will therefore first need to define the domain and some of its properties.
In this example, we chose a 3-dimensional domain and therefore have three inputs.
The steps are:

  #. Define the lower and upper indices of the domain:

     .. code-block:: cpp

        amrex::IntVect dom_lo(0,0,0);
        amrex::IntVect dom_hi(n_cell-1, n_cell-1, n_cell-1);

     We use two IntVects to define the high and low indices of the domain. In the case
     of the high indices, we define it as ``n_cell-1``. ``n_cell`` represents the
     number of cells we want in each dimension. Typically its value is
     read from the inputs file.

  #. Define a single Box with the specified domain:

     .. code-block:: cpp

        amrex::Box domain(dom_lo, dom_hi);

  #. Define the BoxArray object using the domain:

     .. code-block:: cpp

        amrex::BoxArray ba(domain);

  #. Next we define a maximum size for grids ("chunks") in the domain:

     .. code-block:: cpp

        ba.maxSize(max_grid_size);

     ``max_grid_size`` will be read from the inputs file at runtime. Its value will effect how
     AMReX processes the data on the MultiFab in parallel.


Distribution Mapping
--------------------

Once the BoxArray is defined. We can define a DisbributionMapping. The
DistributionMapping will determine how parts of the domain are divided among
processes. For us, the default behavior is sufficient. Therefore we write,

.. code-block:: cpp

   amrex::DistributionMapping dm(ba);

At this point, we have defined the four necessary parts of the MultiFab and
can create it with the line,

.. code-block:: cpp

   amrex::MultiFab mf(ba, dm, ncomp, ngrow);

Adding Data to the MultiFab
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have setup an abstract domain, we want to provide physical
information about the boxes within it so that we can fill it with
useful data. To do that, we will define a Geometry object, and
use a nested loop structure that efficiently iterates over
all the boxes in the domain. AMReX will automatically decide how to
parallelize these operations with MPI ranks or GPU accelerators.


Geometry: Add Physical Properties
---------------------------------

For the geometry, we will
define the size of the box, the coordinate system and the boundary conditions.
It is also convenient to derive the dimensions of each cell at this time.

  #. Define the physical dimensions of the box:

     .. code-block:: cpp

        amrex::RealBox real_box ({ 0., 0., 0.} , {1., 1., 1.});

     In this case, we use ``[0,1]`` for all directions.


  #. Define the Geometry object with the properties specified above:

     .. code-block:: cpp

        amrex::Geometry geom(domain, &real_box);

  #. \*In preparation for future operations on the MultiFab, extract the physical
     dimensions of each cell:

     .. code-block:: cpp

        amrex::GpuArray<amrex::Real,3> dx = geom.CellSizeArray();

     This commands creates a 1-dimensional array of ``amrex::Real``--single or double floating point--
     values that correspond to the physical dimensions of each cell side, i.e. ``dx[0]`` contains the length
     of the cell in the x-direction, ``dx[1]``, contains the length in the y-direction, and so on.


Loop Structure: MFIter and ParallelFor
--------------------------------------

To traverse the elements of the MultiFab we will use two loops. The first loop is a MultiFab iterator,
called ``MFIter``, and the second loop is a ``ParallelFor``. The ``MFIter`` loop efficiently divides
and distributes work on MultiFab data among processors. It does this by telling each processor,
"iterate through all the boxes on the MultiFab, but only do work on the boxes assigned
to you."

Inside the ``MFIter`` for loop, the
``ParallelFor`` function behaves like a triple-nested loop over the ``i,j,k`` coordinates in the box.
Beginning and ending values for each index are derived
from the inputs, and do not need to be explicitly stated. When a GPU backend is enabled, ``ParallelFor``
will launch a thread on the GPU to process each ``i,j,k`` iteration in parallel.
On CPU only, ``ParallelFor`` will enable tiling to process the data in the most efficient way for the hardware available.



Below is an example of typical usage of ``MFIter`` and ``ParallelFor`` to fill the MultiFab with data:

.. code-block:: cpp

   for(amrex::MFIter mfi(mf); mfi.isValid(); ++mfi){

       const amrex::Box& bx = mfi.validbox();
       const amrex::Array4<amrex::Real>& mf_array = mf.array(mfi);

       amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){

           amrex::Real x = (i+0.5) * dx[0];
           amrex::Real y = (j+0.5) * dx[1];
           amrex::Real z = (k+0.5) * dx[2];
           amrex::Real rsquared = ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
           mf_array(i,j,k) = 1.0 + std::exp(-rsquared);

       });
   }

In the first line,

.. code-block:: cpp

   for(amrex::MFIter mfi(mf); mfi.isValid(); ++mfi){

the MFIter object ``mfi`` is defined from the MultiFab ``mf``. ``mfi.isValid()``, tells the for loop
to work on only the parts of the domain designated to that processor by ``mf`` including growth or ghost cells.

The following two lines,

.. code-block:: cpp

   const amrex::Box& bx = mfi.validbox();
   const amrex::Array4<amrex::Real>& mf_array = mf.array(mfi);

define a new box, ``bx`` that represents the current section of the grid being iterated on. To access
and store data in the components
of the MultiFab within this section of the grid, we cast ``mf`` as an Array4 object called ``mf_array``.
We declare ``mf_array`` by reference so that
values changed within the ``ParallelFor`` loop, change values in the components of the MultiFab ``mf``.
Within ``ParallelFor`` each component of ``mf`` can be accessed by
calling ``mf_array(i,j,k,n)`` where ``i,j,k`` represents the location and ``n`` represents the
:math:`(n+1)^{\text th}` component in the MultiFab. Moreover, the first or :math:`0^{\text th}` component, can be accessed
by dropping the 4th index, i.e. ``mf_array(i,j,k)`` is equivalent to ``mf_array(i,j,k,0)``.


The line with the ``ParallelFor`` function,

.. code-block:: cpp

   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){

calls a custom looping function that traverses the three dimensions of the box, ``bx``. The range of values
for ``i,j,k`` are determined by the ``ParallelFor`` from the box information passed to it.
When a GPU backend is enabled, code inside ``ParallelFor`` will be executed in parallel on the GPU. Functions called
within this loop are often called "lambdas". Note that, because ``ParallelFor`` is a function call, and
not a simple for loop, the closing brackets are ``});``.


The remaining lines inside ``ParallelFor`` are,

.. code-block:: cpp

   amrex::Real x = (i+0.5) * dx[0];
   amrex::Real y = (j+0.5) * dx[1];
   amrex::Real z = (k+0.5) * dx[2];
   amrex::Real rsquared = ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
   mf_array(i,j,k) = 1.0 + std::exp(-rsquared);

The first three lines translate integer indices to their ``(x,y,z)`` location in the domain. In this case,
the ``+0.5`` indicates that we want the :math:0^{\text th} index for ``i``, to represent the middle of
the first cell in the x-direction. Together, treating all three variables in this way
indicates cell-center data.

The last line, stores the calculated value in the appropriate location of ``mf_array``
which, in turn, stores it in the corresponding location of the MultiFab ``mf`` because
of the way we declared ``mf_array``.


So far, we've defined and filled a MultiFab with a single component with data. The filling
operation will be done in parallel when compiled with a GPU backend. When compiled for CPU-only computation,
AMReX will use tiling
to increase performance (see `AMReX User's Guide -- Tiling`_). Configuring these
performance optimizations would normally require different lines of code, however, AMReX handles these changes
automatically for portable performance.

.. _`AMReX User's Guide -- Tiling`: https://amrex-codes.github.io/amrex/docs_html/Basics.html#mfiter-and-tiling

Plotting MultiFab Data
~~~~~~~~~~~~~~~~~~~~~~

AMReX can plot MultiFab data with a single function call. To access the plotting
functions we include the header,

.. code-block::

   #include <AMReX_PlotFileUtil.H>

at the top of the file. We can then write the line,

.. code-block:: cpp

   WriteSingleLevelPlotfile("plt001", mf, {"comp0"}, geom, 0., 0);

The first input takes the plotfile name, "plt001".
The second, is the MultiFab that contains th e data we want to plot, ``mf``.
The parameter, ``{"comp0"}``, labels the
first component as "comp0". For multiple components, its necessary to pass additional
variable names, such as ``{"comp0","comp1"}``.  The third parameter, ``geom``, is the
Geometry we defined for the MultiFab above. The
last parameter specifies the level of the MultiFab. In this tutorial, we only had a
single MultiFab to represent the zeroth level of the mesh, therefore we pass ``0`` for the level.

Visualizing the Plotfile
------------------------

The call above to ``WriteSingleLevelPlotfile`` will produce a plotfile in the form
of a directory that contains a Header file, and several subdirectories. The data can be
visualized by passing this directory to one of several visualization software packages.
Information on how to do this is available in the AMReX User's Guide `Visualization`_
section.


.. _`Visualization`: https://amrex-codes.github.io/amrex/docs_html/Visualization_Chapter.html


Conclusion
~~~~~~~~~~

In this tutorial we described the steps to define the domain, and physical properties
of our simulation space. We then demonstrated how to initialize and store scalar values
on a grid across
this domain in the form of a 3-dimensional MultiFab data structure with a single component.
Initializing the data involved the ``MFIter``, ``ParallelFor`` nested-"loop" structure that
required casting the MultiFab as an Array4 for ``i,j,k`` access to the active section of the
grid. Finally, we showed the commands to write out the MultiFab data to a plotfile that
can be visualized with several software packages.

The complete code for this tutorial is available `here`_.

.. _`here`: https:://where?

|

Additional Comments
~~~~~~~~~~~~~~~~~~~

AMReX allows for the selection of 2- or 3- dimensional simulation
at compile time. For simplicity, the example above is presented as 3-dimensional only.
The example below shows how to write code to assign data to a MultiFab
in a way that can be automatically adapted for 2 or 3 dimensions:

.. code-block:: cpp

       for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
       {
           const Box& bx = mfi.validbox();

           const Array4<Real>& phiOld = phi_old.array(mfi);

           // set phi = 1 + e^(-(r-0.5)^2)
           amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
           {
               Real x = (i+0.5) * dx[0];
               Real y = (j+0.5) * dx[1];
   #if (AMREX_SPACEDIM == 2)
               Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01;
   #elif (AMREX_SPACEDIM == 3)
               Real z= (k+0.5) * dx[2];
               Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
   #endif
               phiOld(i,j,k) = 1. + std::exp(-rsquared);
           });
       }

The interesting features to point out are:

  - the preprocessor directives, ``#if``, ``#elif``, ``AMREX_SPACEDIM``, etc. are evaluated
    at compile time, and will only process the code within the appropriate
    ``AMREX_SPACEDIM`` section. Moreover, the value of ``AMREX_SPACEDIM`` can be
    set while compiling the AMReX application.

  - The indices of the ``ParallelFor`` function are ``i,j,k`` despite being applicable
    to both 2- or 3- dimensional code. AMReX has the built-in ability to revert to the two
    dimensional ``i,j`` indices, when compiled as a 2D code.

  - Additional code differences, not shown here, are required in the steps to define
    the BoxArray and Geometry of the MultiFab. An example of the needed modifications
    can be found in ``HeatEuation_EX0_C``.

A ParallelFor Only Approach
---------------------------

The method presented above is the most common at the time of writing the guide. However,
other methods exist to achieve similar results. For example, AMReX has added to the
capabilities of ``ParallelFor`` to include the functionality of both the loops in our
previous example. To use this type of ``ParallelFor`` we need to include the header,

.. code-block:: cpp

   #include <AMReX_MFParallelFor.H>


The newer approach reduces the necessary syntax. However, we still need to cast
the MultiFab as a different type, so that we can access it with ``i,j,k`` indices.
We also need to explicitly pass the number of grow or ghost cells in each direction.
These two things are accomplished in the first two lines.

.. code-block:: cpp


   const& amrex::MultiArray4 mf_arrs = mf.arrays();
   const amrex::IntVect ngs(ngrow);


   amrex::ParallelFor(
      mf, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {

         amrex::Real x = (i+0.5) * dx[0];
         amrex::Real y = (j+0.5) * dx[1];
         amrex::Real z = (k+0.5) * dx[2];
         amrex::Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
         mf_arrs[nbx](i,j,k) = 1. + std::exp(-rsquared);
      }
   });

The other change worth mentioning is ``int nbx``, the first iterative variable in the
inputs to the ``ParallelFor``. In comparison to the first example, iterating through
this variable mimics the functionality of the ``MFIter`` for loop in the first
approach shown in this tutorial.

Finally one last comment, the line,

.. code-block:: cpp

   const& amrex::MultiArray4 mf_arrs = mf.arrays();

is often written to take advantage of the compiler's ability to determine the correct type.
To do this, we replace the above line with

.. code-block:: cpp

   auto const& mf_arrs = mf.arrays();


|

What's Next?
~~~~~~~~~~~~

This tutorial provided an introduction to the MultiFab data structure. The :ref:`guided_heat_simple`
and :ref:`guided_heat` tutorials both have source code that demonstrate MultiFab usage in a loop that evolves
the data over time. At the cost of additional complexity, the :ref:`guided_heat` tutorial makes use
of the :cpp:`using namespace` and preprocessor directives to make writing coding easier and add
additional functionality. On the other hand, the source code in :ref:`guided_heat_simple` more
closely resembles the code used in this example.




