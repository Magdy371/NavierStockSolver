#ifndef GRID_H
#define GRID_H

#include <vector>
#include <stdexcept>

// Stores 3D velocity (u, v, w) and pressure (p) fields
class Grid {
public:
    // Constructor: Initialize grid with size (nx, ny, nz) and spacing (dx, dy, dz)
    Grid(size_t nx, size_t ny, size_t nz, double dx, double dy, double dz);

    // Accessors for velocity and pressure
    double& u(size_t i, size_t j, size_t k);
    double& v(size_t i, size_t j, size_t k);
    double& w(size_t i, size_t j, size_t k);
    double& p(size_t i, size_t j, size_t k);

    // Const accessors
    double u(size_t i, size_t j, size_t k) const;
    double v(size_t i, size_t j, size_t k) const;
    double w(size_t i, size_t j, size_t k) const;
    double p(size_t i, size_t j, size_t k) const;

    // Get grid properties
    size_t nx() const { return nx_; }
    size_t ny() const { return ny_; }
    size_t nz() const { return nz_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double dz() const { return dz_; }

private:
    size_t nx_, ny_, nz_; // Grid dimensions
    double dx_, dy_, dz_; // Grid spacing
    std::vector<std::vector<std::vector<double>>> u_, v_, w_, p_; // Fields

    // Check for valid indices
    void checkBounds(size_t i, size_t j, size_t k) const;
};

#endif