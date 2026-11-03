# Parser ★★

Evaluate functions given as strings at runtime, e.g., initial conditions
from the inputs file, inside GPU kernels. `amrex::Parser` handles
floating-point expressions and `amrex::IParser` integer ones. The program
compares the results with hard-coded C++ versions of the default expressions.

**Key APIs:** `Parser`, `IParser`, `setConstant`, `registerVariables`,
`compile<N>`, `ParmParse::makeParser`

## Notes

- Build the parser on the host, set constants, register the variables in
  the order they will be passed, then `compile<N>()` with `N` variables. The
  returned executor is a small object that is captured by value in device
  lambdas; keep the `Parser` alive while it is used.
- `ParmParse::makeParser` resolves unknown symbols from the ParmParse
  database, here entries such as `my_constants.A`.
- Expressions support `+ - * / ^ **`, comparisons, `and`, `or`,
  `if(cond, a, b)`, `min`, `max`, and the usual math functions. `IParser`
  also has floor division `//`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_Parser -j 8
cd build/Base/Parser && ./Parser inputs
./Parser inputs ic.function="x*y*exp(-t)"
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
./main3d.gnu.MPI.ex inputs
```

## Expected output

```
ic.function = A * sin(2*pi*x) * cos(2*pi*y) * exp(-t)
ic.tag      = max(i, max(j, k)) < n/2
max |phi| = 1.201407003, tagged cells = 4096
Parser vs. makeParser: max diff = 0
Parser vs. C++: max diff = 0
IParser vs. C++: max diff = 0
PASSED
```
