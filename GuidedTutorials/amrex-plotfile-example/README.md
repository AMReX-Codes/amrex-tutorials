# AMReX Plotfile Example

This sample code initializes a MultiFab with one component "phi" using
the `(i,j,k)` indices of each cell as:

```
phi(i,j,k) = i + 100.0*j + 10000.0*k
```

The MultiFab is then written to a plotfile.

## Compiling

First, compile AMReX as a library using the directions here:
https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-libamrex

This was tested with g++ 5.4.0 with the default options supplied to `./configure` as:

```
./configure --prefix=[AMReX library prefix]
```

Then compile this example as:

```
make AMREX_LIBRARY_HOME=[AMReX library prefix]
```

## Running

If compilation is successful you should have a `main.exe` executable.

This can be run as-is with defaults corresponding to a 32^3 domain and
max grid size of 16. This will create a domain of size 32^3 cells in
3D and divide it up into 8 boxes each of size 16^3.

It will also write a plotfile named `plt_32_32_32_16`. The naming
convention for the plotfiles is:

```
plt_[# cells x]_[# cells y]_[# cellsz]_[max grid size]
```

To change these defaults, you can use command line options. For
example, this will set up a 64^3 domain and divide it up into 8 boxes
of size 32^3:

```
./main.exe n_cells=64 64 64 max_grid_size=32
```

## Setting C++ flags and library paths manually (optional)

The makefile should automatically find and use the correct C++ and
library flags corresponding to the ones AMReX was built with.

If you need for some reason to change this, or if the flags cannot be
automatically found, then set `CFLAGS` and `LFLAGS` appropriately in
the GNUmakefile.

These flags are automatically extracted from the following file:

```
[AMReX library prefix]/lib/pkgconfig/amrex.pc
```

This file lists entries for `Cflags: ...` and `Libs: ...` something
like the following, but for your system:

```
...
Cflags: -I${includedir}  -Werror=return-type -g -O3 -std=c++14
Libs: -L${libdir} -lamrex -L/usr/lib/gcc/x86_64-linux-gnu/5/ -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichfort -lmpich -lgfortran -lquadmath
...
```
