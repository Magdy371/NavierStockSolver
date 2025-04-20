#include "BoundaryConditions.h"

void DirichletBC::apply(Grid& grid) {
    size_t nx = grid.nx(), ny = grid.ny(), nz = grid.nz();

    // X-faces (y-z planes at x=0 and x=nx-1)
    for (size_t j = 0; j < ny; ++j) {
        for (size_t k = 0; k < nz; ++k) {
            grid.u(0, j, k) = 0.0;
            grid.u(nx - 1, j, k) = 0.0;
            grid.v(0, j, k) = 0.0;
            grid.v(nx - 1, j, k) = 0.0;
            grid.w(0, j, k) = 0.0;
            grid.w(nx - 1, j, k) = 0.0;
        }
    }

    // Y-faces (x-z planes at y=0 and y=ny-1)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t k = 0; k < nz; ++k) {
            grid.u(i, 0, k) = 0.0;
            grid.u(i, ny - 1, k) = 0.0;
            grid.v(i, 0, k) = 0.0;
            grid.v(i, ny - 1, k) = 0.0;
            grid.w(i, 0, k) = 0.0;
            grid.w(i, ny - 1, k) = 0.0;
        }
    }

    // Z-faces (x-y planes at z=0 and z=nz-1)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            grid.u(i, j, 0) = 0.0;
            grid.v(i, j, 0) = 0.0;
            grid.w(i, j, 0) = 0.0;
            grid.u(i, j, nz - 1) = 1.0; // Top face moves in x-direction
            grid.v(i, j, nz - 1) = 0.0;
            grid.w(i, j, nz - 1) = 0.0;
        }
    }
}