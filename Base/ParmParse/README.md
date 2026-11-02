# ParmParse ★

Read runtime parameters from an inputs file and the command line with
`amrex::ParmParse`: scalars, strings, arrays, math expressions, enums, 2D
tables, TOML-like `[section]` headers, and detection of unused (often
misspelled) inputs. The program checks that the shipped `inputs` (or the
built-in defaults) give the documented values.

**Key APIs:** `ParmParse::query`, `get`, `queryarr`, `queryAdd`,
`queryWithParser`, `querytable`, `query_enum_case_insensitive`, `AMREX_ENUM`,
`ParmParse::prettyPrintUnusedInputs`

## Notes

- `query` returns whether the name was found and leaves the variable alone
  otherwise; `get` aborts if the name is missing. `queryAdd` also stores the
  default in the database.
- Arguments after the inputs file, e.g. `nsteps=200`, override the file.
- `dt = 1.0/64` is evaluated by `query`. Expressions with spaces, such as
  `omega = 2 * pi / L`, need `queryWithParser`. Unknown symbols are looked up
  in the database, here with the `my_constants.` prefix.
- `[physics]` starts a table: `diffusivity = 0.1` below it is
  `physics.diffusivity`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_ParmParse -j 8
cd build/Base/ParmParse && ./ParmParse inputs nsteps=200
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
./main3d.gnu.MPI.ex inputs nsteps=200
```

## Expected output

```
nsteps        = 200 (inputs)
...
time_stepper  = RK2
butcher_a     = { 0 0 } { 1 0 }
physics.diffusivity = 0.1 (inputs)
...
Unused inputs:
physics.difusivity = 0.2

PASSED
```
