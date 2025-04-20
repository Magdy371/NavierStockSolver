#include "NavierStokesSolver.h"
#include <fstream>
#include <iostream>

int main() {
    try {
        // Read configuration
        std::ifstream config("config.txt");
        if (!config.is_open()) {
            throw std::runtime_error("Could not open config.txt");
        }

        size_t nx, ny, nz, num_steps;
        double dx, dy, dz, nu, rho, gravity, dt;
        config >> nx >> ny >> nz >> dx >> dy >> dz >> nu >> rho >> gravity >> dt >> num_steps;
        config.close();

        // Create solver
        NavierStokesSolver<double> solver(nx, ny, nz, dx, dy, dz, nu, rho, gravity);
        solver.setBoundaryConditions(std::make_unique<DirichletBC>());

        // Set initial condition (top face moves in x-direction)
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                solver.grid().u(i, j, nz - 1) = 1.0;
            }
        }

        // Run simulation
        for (int step = 0; step < num_steps; ++step) {
            solver.step(dt);
            std::cout << "Step " << step + 1 << " completed\n";
        }

        // Output result
        std::cout << "Final u-velocity at center: "
                  << solver.grid().u(nx / 2, ny / 2, nz / 2) << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}