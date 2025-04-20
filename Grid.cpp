#include "Grid.h"
#include <vector>
using namespace std;
Grid::Grid(size_t nx, size_t ny, size_t nz, double dx, double dy, double dz)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz),
      u_(nx,vector<vector<double>>(ny, vector<double>(nz, 0.0))),
      v_(nx, vector<vector<double>>(ny, vector<double>(nz, 0.0))),
      w_(nx, vector<vector<double>>(ny, vector<double>(nz, 0.0))),
      p_(nx, vector<vector<double>>(ny, vector<double>(nz, 0.0))) {
    if (nx < 3 || ny < 3 || nz < 3) {
        throw std::invalid_argument("Grid size must be at least 3x3x3");
    }
    if (dx <= 0 || dy <= 0 || dz <= 0) {
        throw std::invalid_argument("Grid spacing must be positive");
    }
}
double& Grid::u(size_t i, size_t j, size_t k) {
    checkBounds(i, j, k);
    return u_[i][j][k];
}

double& Grid::v(size_t i, size_t j, size_t k) {
    checkBounds(i, j, k);
    return v_[i][j][k];
}

double& Grid::w(size_t i, size_t j, size_t k) {
    checkBounds(i, j, k);
    return w_[i][j][k];
}

double& Grid::p(size_t i, size_t j, size_t k) {
    checkBounds(i, j, k);
    return p_[i][j][k];
}
double Grid::u(size_t i, size_t j, size_t k) const {
    checkBounds(i, j, k);
    return u_[i][j][k];
}

double Grid::v(size_t i, size_t j, size_t k) const {
    checkBounds(i, j, k);
    return v_[i][j][k];
}

double Grid::w(size_t i, size_t j, size_t k) const {
    checkBounds(i, j, k);
    return w_[i][j][k];
}

double Grid::p(size_t i, size_t j, size_t k) const {
    checkBounds(i, j, k);
    return p_[i][j][k];
}

void Grid::checkBounds(size_t i, size_t j, size_t k) const {
    if (i >= nx_ || j >= ny_ || k >= nz_) {
        throw std::out_of_range("Grid index out of bounds");
    }
}
