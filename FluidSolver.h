#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include "Grid.h"
#include "BoundaryConditions.h"
#include <memory>

class FluidSolver {
public:
    virtual ~FluidSolver() = default;
    virtual void step(double dt) = 0;
    virtual void setBoundaryConditions(std::unique_ptr<BoundaryConditions> bc) = 0;
    virtual Grid& grid() = 0;
};

#endif