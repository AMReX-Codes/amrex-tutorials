Demo Tutorial
=============

..
   Questions*
   What do people need fingers on keys for. What are the core things to have them do.

 
.. admonition:: **Objectives**
   :class: warning

   - Compile an AMReX Code 
   - Introduce Basic Code Structure
   - Generate and Visualize Output     
     
   **Time to Complete**: 20 mins. 


In this tutorial you will take the steps needed to go from download to
visualized output of an AMReX code. We will demonstrate basic building, 
compiling and output generation. We will also examine key AMReX features
in a C++ code and plot the output with Python in a Jupyter notebook.


Setting Up Your Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorial recommends using Gitpod (Requires a GitHub account).  Gitpod
provides an online terminal that is already preconfigured for our development 
environment.

 Click |GitpodLink| to be taken to the Gitpod workspace. 

.. |GitpodLink| raw:: html

   <a href="https://gitpod.io/#https://github.com/atmyers/ecp-tutorials" target="_blank">here</a>

..
    To download and build AMReX yourself see:
    https://amrex-codes.github.io/amrex/docs_html/GettingStarted.html
    and
    https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX_Chapter.html


Building the Project 
~~~~~~~~~~~~~~~~~~~~

In this example, will use CMake to build the project. Navigate to the directory
:code:`/workspace/exp-tutorials`
and type

.. code-block:: 
   
   mkdir Build
   cd Build
   cmake ..

This will run CMake for all the tutorial directories. During this process
CMake will generate the build files it needs to compile each individual
tutorial.


More information about building options, such as disabling MPI, can be found at
https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX_Chapter.html.

Compiling the Code
~~~~~~~~~~~~~~~~~~

After building the project, Navigate to the directory :code:`03_HeatEquation`. 
At the prompt type :code:`make` and
CMake will compile the code and dependencies. The first time you call :code:`make`, 
you should see a list of all the source AMReX files being compiled:

.. image:: ./images_tutorial/Cmake_make.png

|

When CMake finishes you will be left with an executable named :code:`03_HeatEquation`. 
To run the code type:

.. code-block::

   ./03_HeatEquation inputs

This command will run the :code:`03_HeatEquation` code with the :code:`inputs` file as
the input parameters. Parsing of the information in the :code:`inputs` file is done by
:code:`ParmParse`. More details can be found at
https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse

Code Highlights
~~~~~~~~~~~~~~~

At this point we have built, compiled and ran the :code:`03_HeatEquation` code. Now
we will walk through the code and explain some essential features of AMReX syntax.

Basic Structure
^^^^^^^^^^^^^^^
::

   Main
    |---- Declare Simulation Parameters
    |---- Read Parameter Values From Input File
    |---- Define Simulation Setup & Geometry
    |---- Initialize Data Loop
    |     |---- Set Values For Each Cell
    |---- Write Initial Plot File
    |---- Main Time Progression Loop
          |---- Evolve Values For Each Cell
          |---- Increment
          |---- Write Plot File At Given Interval


AMReX Namespace and Required Commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  

The MultiFab Data Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :code:`MultiFab` is a data structure that can 
be distributed among parallel processes. For 



MFIter and ParallelFor
^^^^^^^^^^^^^^^^^^^^^^

These are the commands we use to iterate through
the cells at each time step. The command MFIter
will go iterate through ... 
The command ParallelFor will automatically utilize
parallel computation methods, such as MPI, OMP, GPUs
or HIP, to iterate through the multidimensional array. 

For more information on the basic components of AMReX, please see
https://amrex-codes.github.io/amrex/docs_html/Basics.html


Visualizing Output
~~~~~~~~~~~~~~~~~~

Data Files
^^^^^^^^^^

In :code:`main.cpp` we called a plot function in two places. The
first time was to plot initial data.

.. code-block::

   129     if (plot_int > 0)
   130     {
   131         int step = 0;
   132         const std::string& pltfile = amrex::Concatenate("plt",step,5);
   133         WriteSingleLevelPlotfile(pltfile, phi_old, {"phi"}, geom, time, 0);
   134     }


The second time plots were generated at given intervals during
the main time progression loop.

.. code-block::

   171         if (plot_int > 0 && step%plot_int == 0)
   172         {
   173             const std::string& pltfile = amrex::Concatenate("plt",step,5);
   174             WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, step);
   175         }

Each time we run the code it will create a series of directories which contain 
data for visualization. Now run :code:`03_HeatEquation` with the :code:`inputs`
file. After it finishes your directory should look like this. 

.. image:: ./images_tutorial/plot_dirs.png


Visualization in Jupyter
^^^^^^^^^^^^^^^^^^^^^^^^

We will use Python and the yt package in a Jupyter notebook to generate plots for the data 
in the directories created in the previous step. First launch the Jupyter notebook
with the command:

.. code-block::

   jupyter notebook

When Jupyter starts, it will generate a token at the command line
and ask for a password in the window it opened. Copy the token
to enter to the notebook.

.. image:: ./images_tutorial/token_hl.png


Once the notebook starts, find :code:`Visualization.ipynb` and open it. 
In this file there are additional notes about the
heat equation example, followed by several cells that use :code:`yt` 
commands to read AMReX output files.  

yt
^^

The following commands import the :code:`yt` package and plot
a 2D slice of the output at from the 1000th time step. 

.. code-block::

   import yt
   from yt.frontends.boxlib.data_structures import AMReXDataset
   ds = AMReXDataset("plt01000")
   sl = yt.SlicePlot(ds, 2, ('boxlib', 'phi'))
   sl

In our example, the commands are already written in the notebook.
To run them, select from the menu: `Kernel -> Restart & Run All`.
Once the run is complete, you will get the following plot.


.. image:: ./images_tutorial/heat_eq_plot.png


Tutorial Features
~~~~~~~~~~~~~~~~~

Useful Features:
  - objectives and time listed at the beginning of the tutorial.
  - less explanations, more actions to follow. Longer explanations linked to. 
  - frequent headings and short text blocks.
