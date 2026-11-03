# ElectrostaticPIC ★★★ / Blueprint

An electrostatic particle-in-cell (PIC) code for electrons moving through
a fixed, neutralizing ion background in a periodic box. It is organized so
that it can be copied as the starting point of a PIC code.

| File | Contents |
|---|---|
| `PICParticleContainer.H/.cpp` | Pure SoA particle species: initialization, charge deposition, field gather, push |
| `FieldSolver.H/.cpp` | Periodic Poisson solve with MLMG and `E = -grad(phi)` |
| `main.cpp` | Setup, time loop, diagnostics and plotfiles |

Each step:

1. Deposit the electron charge with cloud-in-cell weights
   (`ParticleToMesh`) and add the ion charge.
2. Solve `lap(phi) = -rho/epsilon0` with `MLPoisson` and `MLMG`, and
   compute `E` at cell centers with centered differences.
3. Interpolate `E` to the particles (`MeshToParticle`).
4. Leapfrog push: `v += (q/m) E dt`, then `x += v dt`, then `Redistribute`.

The problem is a Langmuir (plasma) oscillation. In normalized units
(`epsilon0 = m = n0 = 1`, `q = -1`), the plasma frequency is
`omega_p = sqrt(n0 q^2 / (epsilon0 m)) = 1`. The electrons start at rest,
displaced by `xi0 sin(kx)`. Linear theory gives
`E_x = -q n0 xi0 / epsilon0 sin(kx) cos(omega_p t)`. The code measures the
frequency from the zero crossings of the `sin(kx)` mode of `E_x`, and
checks it, the initial amplitude, and the total (field plus kinetic)
energy against theory to 2%.

**Key APIs:** `amrex::ParticleContainerPureSoA`, `amrex::ParticleToMesh`,
`amrex::MeshToParticle`, `amrex::ParticleInterpolator::Linear`,
`amrex::MLPoisson`, `amrex::MLMG`, `amrex::ParReduce`, `Redistribute`,
`WriteSingleLevelPlotfile`, `WritePlotFile`

## Build and run

Requires AMReX with particles and linear solvers, in 2D or 3D. CMake, from
the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Particles
cmake --build build --target Particles_ElectrostaticPIC -j 8
cd build/Particles/ElectrostaticPIC && mpiexec -n 2 ./ElectrostaticPIC inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

`n_cell` and `nppc` take three values; 2D uses the first two. For a
longer or finer run, increase `nsteps`, `n_cell` and `nppc`.

## Expected output

```
Number of particles: 8192
Step 0, t = 0: E mode 0.009936866131, field energy 0.0003827000414, total energy 0.0003836504961
...
Initial E mode amplitude: 0.009936866131 (linear theory 0.01)
Measured frequency: 0.9956005796 (plasma frequency 1)
Relative errors: frequency 0.004399420394, amplitude 0.006313386911, max energy change 0.001909800398
PASSED
```

Plotfiles `plt00000`, `plt00050` and `plt00100` contain `rho`, `phi`,
`Ex Ey Ez` and the electrons.

## Extending

- More species: create one container per species and sum their charge
  densities before the solve.
- Other boundary conditions: change `setDomainBC` in `FieldSolver` (and the
  particle boundary handling), or use `amrex::FFT::Poisson` for periodic
  problems when AMReX is built with FFT support.
