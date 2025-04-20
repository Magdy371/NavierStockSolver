#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "Grid.h"
#include <memory>

// Abstract base class for boundary conditions
class BoundaryConditions {
public:
    virtual ~BoundaryConditions() = default;
    virtual void apply(Grid& grid) = 0;
};

// Dirichlet boundary conditions (no-slip walls, moving top face)
class DirichletBC : public BoundaryConditions {
public:
    void apply(Grid& grid) override;
};

#endif