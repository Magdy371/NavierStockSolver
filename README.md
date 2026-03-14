# NavierStokesSolver

A 3D incompressible Navier-Stokes solver in C++ using the **fractional-step (projection) method** on a uniform Cartesian grid. Designed for simulating viscous fluid flow in a lid-driven cavity.

## How It Works

The solver advances the fluid state each time step through four stages:

1. **CFL Check** — Validates the time step satisfies the Courant–Friedrichs–Lewy stability condition.
2. **Intermediate Velocity** — Computes a provisional velocity field using the convection-diffusion equation (central differences) with gravity.
3. **Pressure Poisson** — Solves the pressure correction equation via Jacobi iteration to enforce incompressibility (∇·**u** = 0).
4. **Velocity Correction** — Projects the intermediate velocity onto a divergence-free field using the pressure gradient.

Boundary conditions are applied after each step via a pluggable interface (currently Dirichlet/no-slip walls with a moving top face).

## Project Structure

```
NavierStockSolver/
├── include/
│   ├── grid/
│   │   └── Grid.h              # 3D velocity & pressure field storage
│   ├── boundary/
│   │   └── BoundaryConditions.h # Abstract BC interface + DirichletBC
│   └── solver/
│       ├── FluidSolver.h        # Abstract solver interface
│       └── NavierStokesSolver.h  # Templated N-S solver implementation
├── src/
│   ├── grid/Grid.cpp
│   ├── boundary/BoundaryConditions.cpp
│   ├── solver/NavierStokesSolver.cpp
│   └── main.cpp                 # Entry point — reads config & runs simulation
├── CMakeLists.txt
└── config.txt                   # Runtime parameters
```

## Configuration

Create a `config.txt` in the project root with space-separated values:

```
nx ny nz dx dy dz nu rho gravity dt num_steps
```

| Parameter   | Description                        | Example |
|-------------|------------------------------------|---------|
| `nx ny nz`  | Grid points in x, y, z (≥ 3 each) | `10 10 10` |
| `dx dy dz`  | Grid spacing                       | `0.1 0.1 0.1` |
| `nu`        | Kinematic viscosity                | `0.01` |
| `rho`       | Fluid density                      | `1.0` |
| `gravity`   | Gravitational acceleration (z-dir) | `9.81` |
| `dt`        | Time step size                     | `0.001` |
| `num_steps` | Number of time steps               | `100` |

**Example `config.txt`:**
```
10 10 10 0.1 0.1 0.1 0.01 1.0 9.81 0.001 100
```

## Build & Run

Requires **CMake ≥ 3.31** and a **C++20** compiler.

```bash
cmake -B build -S .
cmake --build build
cd build && ./NavierStockSolver
```

## Key Classes

| Class | Role |
|-------|------|
| `Grid` | Stores 3D `u`, `v`, `w`, `p` fields with bounds-checked accessors |
| `BoundaryConditions` | Abstract interface for applying boundary conditions |
| `DirichletBC` | No-slip walls + moving lid (top face u=1) |
| `FluidSolver` | Abstract solver interface |
| `NavierStokesSolver<T>` | Templated fractional-step Navier-Stokes solver |

## License

This project is provided as-is for educational and research purposes.
