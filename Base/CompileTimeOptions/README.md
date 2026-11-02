# CompileTimeOptions ★★★

Turn runtime options into compile-time constants inside kernels. A kernel
with two runtime options (the order of a Laplacian stencil and whether to
add an absorption term) is written once with runtime branches and once with
`ParallelFor` taking `CompileTimeOptions`, which compiles a specialized
kernel for every combination and launches the one that matches the runtime
values. The program checks that both give the same result and times them
with `BL_PROFILE` (see [HeatEquation](../HeatEquation#profiling)).

**Key APIs:** `ParallelFor(TypeList<CompileTimeOptions<...>...>, {...}, box, f)`,
`CompileTimeOptions`, `if constexpr`, `AMREX_ENUM`, `BL_PROFILE`

## Notes

- The kernel receives one extra argument per option, of type
  `std::integral_constant<int, value>`, so `opt.value` can be used in
  `if constexpr` and as a template argument.
- Specialized kernels can drop unused code and use fewer registers, which
  matters most for large GPU kernels. The cost is compile time and binary
  size, which grow with the number of combinations.
- nvcc does not allow an `if constexpr` block to be the first use of a
  captured variable. Use a plain `if` there; the compiler still removes the
  dead branch.
- `AnyCTO` does the same for other launch functions such as
  `ParallelForRNG`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_CompileTimeOptions -j 8
cd build/Base/CompileTimeOptions && ./CompileTimeOptions inputs order=second
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
./main3d.gnu.MPI.ex inputs
```

For timings, build with TinyProfiler on (`-DAMReX_TINY_PROFILE=ON` or
`make TINY_PROFILE=TRUE`) and use e.g. `n_cell=256 nrepeat=20` and an
optimized GPU build.

## Expected output

```
order = fourth, absorption = 1
max |result| = 117.2234743, max diff between the two kernels = 0
PASSED
```

With TinyProfiler on, the profiling tables include the timers
`compile_time_options` and `runtime_branches`.
