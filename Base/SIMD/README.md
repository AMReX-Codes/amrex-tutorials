# SIMD ★★

Vectorization on CPUs, in two parts:

1. **Implicit.** On CPU, `ParallelFor` marks its innermost loop with
   `AMREX_PRAGMA_SIMD`, a promise to the compiler that the iterations are
   independent (`#pragma GCC ivdep` with GCC; clang and debug builds get
   none). Whether the loop is then vectorized depends on the kernel.
   The same function is written twice: once in a way that stops the
   vectorizer (a call to a function that is not inlined, with a branch
   inside), for illustration only, and once in the recommended way
   (inlined, branch-free).
2. **Explicit.** With `AMReX_SIMD=ON`, `ParallelForSIMD` hands the kernel
   several consecutive elements at once as a SIMD register. The same code
   runs one element at a time without the option, and on GPUs.

Both parts check their results against the scalar version and time the
kernels with `BL_PROFILE` (see [HeatEquation](../HeatEquation#profiling)).

**Key APIs:** `AMREX_FORCE_INLINE`, `AMREX_NO_INLINE`, `ParallelForSIMD`,
`SIMDindex`, `simd::Vectorized`, `simd::load_1d`, `simd::store_1d`,
`simd::native_simd_size_real`, `BL_PROFILE`

## Writing kernels that vectorize

- Keep the kernel inlinable: small helper functions marked
  `AMREX_FORCE_INLINE`, no calls into other compilation units, no early
  returns in the middle of the loop body.
- Replace branches by selects: compute every quantity unconditionally and
  pick with `?:`, `amrex::min`, `amrex::max` or `std::abs`. Every lane then
  evaluates everything, so the arguments must be safe in every lane, e.g.
  `std::sqrt(std::abs(x))` or `std::sqrt(amrex::max(x, 0.0_rt))` instead of
  `(x > 0) ? std::sqrt(x) : 0.0_rt`. A `sqrt` of a negative number in an
  inactive lane raises a floating-point exception, which the compiler may
  also do on its own when it speculates both arms of a select. The result
  is still correct, but the floating-point traps `amrex.fpe_trap_*` (see
  [Debugging](https://amrex-codes.github.io/amrex/tutorials_html/VisualizationAndDebugging.html))
  would stop the program there. This is why the traps are a debugging aid
  for scalar runs, and why the tests do not enable them.
- GCC's default floating-point model (`-ftrapping-math`) does not evaluate
  arithmetic inside the arms of a `?:` speculatively, so a select like
  `(p > c) ? p - c : c - p` keeps the loop scalar; `-fno-trapping-math`
  (also part of `-ffast-math`) lifts this restriction. Clang speculates by
  default.
- Loops whose iterations are not independent, e.g. scatter-style updates
  into the same cells, must use `amrex::For`, which has the `n`- and
  `Box`-based signatures of `ParallelFor` but no SIMD pragma on CPU, or
  `amrex::Loop` inside an `MFIter` loop.
- Build for the CPU at hand. The default x86-64 code uses 2-wide
  double-precision vectors; `-DCMAKE_CXX_FLAGS="-march=native"` uses the
  full vector width.
- See what the compiler did: GCC `-fopt-info-vec-optimized` and
  `-fopt-info-vec-missed`, clang `-Rpass=loop-vectorize` and
  `-Rpass-missed=loop-vectorize`. The loops are reported at the location of
  `ParallelFor`'s implementation in the AMReX headers.

## Explicit SIMD

`ParallelForSIMD<Functor>(n, f)` calls `f(si)` with a `SIMDindex` covering
`simd_width` consecutive elements, and a scalar remainder.
`simd::load_1d(ptr, si)` returns a register of these elements (a plain
reference to one element without SIMD), and `simd::store_1d` writes a
register back. The functor derives from `simd::Vectorized<>` and is written
generically over the value type, which is `amrex::Real` or a SIMD register
of `Real`s. `AMReX_SIMD=ON` needs the header-only
[vir-simd](https://github.com/mattkretz/vir-simd) library, which the CMake
superbuild downloads; the GNU Make build is scalar only.

Here the scalar push `x += v*dt` is already vectorized by the compiler, so
the explicit version is not faster. Explicit SIMD pays off for kernels the
auto-vectorizer cannot handle, e.g. reductions with many components over
particles; see `ParticleReduceSIMD` in the
[AMReX documentation](https://amrex-codes.github.io/amrex/docs_html/Particle.html).

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base -DAMReX_SIMD=ON -DCMAKE_CXX_FLAGS="-march=native"
cmake --build build --target Base_SIMD -j 8
cd build/Base/SIMD && ./SIMD inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
./main3d.gnu.MPI.ex inputs
```

For timings, build with TinyProfiler on (`-DAMReX_TINY_PROFILE=ON` or
`make TINY_PROFILE=TRUE`) and run with `tiny_profiler.print_threshold=0`,
since the kernels are short.

## Expected output

```
Part 1: max difference between the two kernels = 0
Part 2: SIMD width 4 (AMReX_SIMD=ON)
Part 2: max difference between the SIMD and the scalar push = 0
PASSED
```

The width is 1 without `AMReX_SIMD`, and 2, 4 or 8 doubles depending on
the CPU and `-march`. With TinyProfiler on, the profiling tables include the
timers `branch_not_inlined`, `select_inlined`, `push_simd` and
`push_scalar`; the inlined kernel is typically several times faster than the
one with the call.
