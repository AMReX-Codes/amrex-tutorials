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
compiling and output generation. We will also add key AMReX functionality
to a C++ code and plot the output with Python in a jupyter notebook.


Setting Up Your Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorial recommends using Gitpod (Requires a GitHub account).  Gitpod
provides an online terminal that is already preconfigured for our development 
environment.

 Click here_ to be taken to the Gitpod workspace. 

.. _here: https://gitpod.io/#https://github.com/atmyers/ecp-tutorials


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
CMake will generate the files it needs to compile each individual
tutorial.


More information about building options, such as disabling MPI, can be found at
https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX_Chapter.html.

Compiling the Code
~~~~~~~~~~~~~~~~~~

After building the project, Navigate to the directory :code:`03_HeatEquation`. At the prompt type :code:`make` and
CMake will compile the code and dependencies. The first time you call :code:`make`, 
you should see a list of all the source AMReX files being compiled:

.. image:: ./images_tutorial/Cmake_make.png



When CMake finishes you will be left with an executable named :code:`03_HeatEquation`. 
To run the code type:

.. code-block::

   ./03_HeatEquation inputs

This command will run the :code:`03_HeatEquation` code with the :code:`inputs` file as
the input parameters. Parsing of the information in the :code:`inputs` file is done by
:code:`ParmParse`. More details can be found at
https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse

Code Walkthrough
~~~~~~~~~~~~~~~~

At this point we have built, compiled and ran the :code:`03_HeatEquation` code. Now
we will walkthrough the code and explain some essential features of AMReX syntax.


[Note to self: HeatEquation_EX0: main.cpp could be further simplified.]


For more information on the basic components of AMReX, please see
https://amrex-codes.github.io/amrex/docs_html/Basics.html




Compiling Heat Equation EX0
~~~~~~~~~~~~~~~~~~~~~~~~~~~


You can make the example by 



Running The Executable
~~~~~~~~~~~~~~~~~~~~~~

Running the executable requires specifying the inputs file. 

Inputs
^^^^^^

The input file contains the initial conditions. 
