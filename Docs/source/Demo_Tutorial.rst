Demo Tutorial
=============


*Questions*

What do people need fingers on keys for. What are the core things to have them do.

 



Welcome to AMReX tutorials. This tutorial will take about 20 minutes and in it
you will learn the anatomy of a basic AMReX codeset, from build to visualized
output. 



Setting Up Your Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorial recommends using Gitpod. Gitpod requires a GitHub account, but 
its online terminal is already preconfigured for our development 
environment.

https://gitpod.io/#https://github.com/atmyers/ecp-tutorials


To download and build AMReX yourself see:
https://amrex-codes.github.io/amrex/docs_html/GettingStarted.html
and
https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX_Chapter.html


Compiling the Code 
~~~~~~~~~~~~~~~~~~

The examples here depend on the files in amrex. So first we wil download and
install it. Then we need to tell GNUMake where to find the files. This is
done by either setting the 


Navigate to the directory `03_HeatEquation`. At the prompt type `make` and
GNUMake will compile the code and dependencies. 


.. image:: path/filename.png

When GNUMake finishes you should be left with an executable named 


Code Walkthrough
~~~~~~~~~~~~~~~~



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
