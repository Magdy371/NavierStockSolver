#ifndef NAVIER_STOKES_SOLVER_H
#define NAVIER_STOKES_SOLVER_H

#include "FluidSolver.h"
#include <vector>

template<typename T = double>
class NavierStokesSolver : public FluidSolver {
public:
    // Constructor: Initialize grid, fluid properties, and gravity
    NavierStokesSolver(size_t nx, size_t ny, size_t nz, T dx, T dy, T dz, T nu, T rho, T gravity);

    // Perform one time step
    void step(double dt) override;

    // Set boundary conditions
    void setBoundaryConditions(std::unique_ptr<BoundaryConditions> bc) override;

    // Access the grid
    Grid& grid() override { return grid_; }

private:
    Grid grid_; // 3D grid
    T nu_; // Viscosity
    T rho_; // Density
    T gravity_; // Gravity in z-direction
    std::unique_ptr<BoundaryConditions> bc_; // Boundary conditions
    std::vector<std::vector<std::vector<T>>> u_star_, v_star_, w_star_; // Intermediate velocities

    // Check CFL condition
    void checkCFL(double dt);

    // Compute intermediate velocity
    void computeIntermediateVelocity(double dt);

    // Solve pressure Poisson equation
    void solvePressurePoisson(double dt);

    // Correct velocity
    void correctVelocity(double dt);
};

#endif