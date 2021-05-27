.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/SENSEI
==========================

SENSEI is a middleware that allows one to send data to various visualization and
analysis back ends through a uniform interface. It's data model and API enable
one to chose the desired visualization and analysis back end for a given task
with out limiting ones options, as the back ends can be inter-changed at run
time via a text based config file.

Overview
--------

The following tutorials are available:

+--------------------+---------------------------------------------------------------------------+
| Name               | Description                                                               |
+====================+===========================================================================+
| Advection_AmrCore  | This tutorial illustrates an explicit SENSEI instrumentation of a code    |
|                    | that makes use of `amrex::AmrMesh`.                                       |
+--------------------+---------------------------------------------------------------------------+
| Advection_AmrLevel | This tutorial illustrates 3 scenarios with a code that makes use of       |
|                    | `amrex::Amr`. The first, `ImplcitAmr`, illustrates using SENSEI with the  |
|                    | built-in instrumentation in `amrex::Amr`. The second, `ExplcitAmr`,       |
|                    | illustrates using SENSEI with an explicit instrumentation. The third,     |
|                    | `ExplicitParticlesAndAmr`, illustrates using SENSEI from a simulation     |
|                    | that generates both particle and meshed based data.                       |
+--------------------+---------------------------------------------------------------------------+

Note that the `Advection_AmrLevel` contains code for 3 different scenarios.
Which of these is active/available depends on how AMReX is compiled. See below
for the details on configuring the build.


Setting up the build
--------------------

Compiling the AMReX SENSEI tutorials requires that SENSEI is previously
installed. The options that SENSEI was built with determine the specific in
situ capabilities available. Additional CMake options must also be passed
when compiling AMReX to activate the SENSEI bridge and adaptors bundled
with AMReX.

Build options
^^^^^^^^^^^^^

The options that AMReX is compiled with determine which SENSEI tutorials are
available. The following table summarizes the various combinations and results.

+---------------------------------+--------------------------------------------------------------+
| CMake Options                   | What gets built                                              |
+=================================+==============================================================+
| -DAMReX_SENSEI=ON               | Enables SENSEI features in AMReX. Required to compile SENSEI |
| -DAMReX_FORTRAN=ON              | tutorials. Enables the AmrCore tutorial and AmrLevel         |
| -DSENSEI_DIR=<path to install>  | implicit tutorial.                                           |
+---------------------------------+--------------------------------------------------------------+
| -DAMReX_SENSEI=ON               | Enables the AmrCore tutorial, AmrLevel explicit tutorial,    |
| -DAMReX_PARTICLES=ON            | and particle based tutorials to be compiled.                 |
| -DAMReX_NO_SENSEI_AMR_INST=TRUE |                                                              |
| -DAMReX_FORTRAN=ON              |                                                              |
| -DSENSEI_DIR=<path to install>  |                                                              |
+---------------------------------+--------------------------------------------------------------+
| -DAMReX_SENSEI=ON               | Enables the AmrCore tutorial, AmrLevel explicit tutorial.    |
| -DAMReX_NO_SENSEI_AMR_INST=TRUE |                                                              |
| -DAMReX_FORTRAN=ON              |                                                              |
| -DSENSEI_DIR=<path to install>  |                                                              |
+---------------------------------+--------------------------------------------------------------+

Running the tutorials
---------------------

Once the tutorials are compiled they can be run from their corresponding
directory.  The executable is passed an AMReX parm-parse `inputs` file
configuring the run. Options inside the `inputs` file configure the SENSEI
instrumentation inside AMReX. Additionally SENSEI needs to configure the
back-end that will process the data generated. This is done with a SENSEI XML
file. Within each tutorial the `sensei` directory contains a number of SENSEI
XML configuration files. The `inputs` file must be modified to point to one of
these. Which one depends on how SENSEI was compiled. For instance the following
snippet from an `inputs` file would configure SENSEI to send data to ParaView
Catalyst,

.. highlight:: shell

::

   sensei.enabled = 1                                 # turn SENSEI in situ on/off
   sensei.config = sensei/render_iso_catalyst_3d.xml  # render simulation data with ParaView Catalyst
   sensei.frequency = 1                               # number of level 0 steps between in situ processing

while the following snippet would configure SENSEI to send data to VisIt Libsim,

.. highlight:: shell

::

   sensei.enabled = 1                                 # turn SENSEI in situ on/off
   sensei.config = sensei/render_iso_libsim_3d.xml    # render simulation data with ParaView Catalyst
   sensei.frequency = 1                               # number of level 0 steps between in situ processing

There are a number of XML files providing the configuration for a number of the
available back-ends. A given SENSEI XML configuration is only valid when the
SENSEI install has been compiled with the requisite back-end enabled.

Note that the `Advection_AmrLevel_ExplicitParticlesAndAmr` uses the file
`inputs.tracers` while the others use the file `inputs`.

The tutorials are run by switching into the tutorial's build directory and
issuing the launching command. For instance the
`Advection_AmrLevel_ImplicitAmr` tutorial is launched by a command similar to:

.. highlight:: shell

::

   mpiexec -np 4 ./Advection_AmrLevel_ImplicitAmr inputs



