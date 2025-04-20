#include "NavierStokesSolver.h"
#include <cmath>
#include <algorithm>

template<typename T>
NavierStokesSolver<T>::NavierStokesSolver(size_t nx, size_t ny, size_t nz, T dx, T dy, T dz, T nu, T rho, T gravity)
    : grid_(nx, ny, nz, dx, dy, dz), nu_(nu), rho_(rho), gravity_(gravity),
      u_star_(nx, std::vector<std::vector<T>>(ny, std::vector<T>(nz, 0.0))),
      v_star_(nx, std::vector<std::vector<T>>(ny, std::vector<T>(nz, 0.0))),
      w_star_(nx, std::vector<std::vector<T>>(ny, std::vector<T>(nz, 0.0))) {
    if (nu <= 0 || rho <= 0) {
        throw std::invalid_argument("Viscosity and density must be positive");
    }
}

template<typename T>
void NavierStokesSolver<T>::step(double dt) {
    checkCFL(dt);
    computeIntermediateVelocity(dt);
    solvePressurePoisson(dt);
    correctVelocity(dt);
    if (bc_) bc_->apply(grid_);
}

template<typename T>
void NavierStokesSolver<T>::setBoundaryConditions(std::unique_ptr<BoundaryConditions> bc) {
    bc_ = std::move(bc);
}

template<typename T>
void NavierStokesSolver<T>::checkCFL(double dt) {
    T max_vel = 0.0;
    for (size_t i = 0; i < grid_.nx(); ++i) {
        for (size_t j = 0; j < grid_.ny(); ++j) {
            for (size_t k = 0; k < grid_.nz(); ++k) {
                max_vel = std::max(max_vel, std::abs(static_cast<T>(grid_.u(i, j, k))));
                max_vel = std::max(max_vel, std::abs(static_cast<T>(grid_.v(i, j, k))));
                max_vel = std::max(max_vel, std::abs(static_cast<T>(grid_.w(i, j, k))));
            }
        }
    }
    T min_dx = std::min({grid_.dx(), grid_.dy(), grid_.dz()});
    if (max_vel > 0 && dt > min_dx / max_vel) {
        throw std::runtime_error("Time step too large for stability (CFL condition)");
    }
}

template<typename T>
void NavierStokesSolver<T>::computeIntermediateVelocity(double dt) {
    for (size_t i = 1; i < grid_.nx() - 1; ++i) {
        for (size_t j = 1; j < grid_.ny() - 1; ++j) {
            for (size_t k = 1; k < grid_.nz() - 1; ++k) {
                // u derivatives
                T dudx = (grid_.u(i + 1, j, k) - grid_.u(i - 1, j, k)) / (2 * grid_.dx());
                T dudy = (grid_.u(i, j + 1, k) - grid_.u(i, j - 1, k)) / (2 * grid_.dy());
                T dudz = (grid_.u(i, j, k + 1) - grid_.u(i, j, k - 1)) / (2 * grid_.dz());
                T d2udx2 = (grid_.u(i + 1, j, k) - 2 * grid_.u(i, j, k) + grid_.u(i - 1, j, k)) / (grid_.dx() * grid_.dx());
                T d2udy2 = (grid_.u(i, j + 1, k) - 2 * grid_.u(i, j, k) + grid_.u(i, j - 1, k)) / (grid_.dy() * grid_.dy());
                T d2udz2 = (grid_.u(i, j, k + 1) - 2 * grid_.u(i, j, k) + grid_.u(i, j, k - 1)) / (grid_.dz() * grid_.dz());

                u_star_[i][j][k] = grid_.u(i, j, k) + dt * (
                    -grid_.u(i, j, k) * dudx - grid_.v(i, j, k) * dudy - grid_.w(i, j, k) * dudz
                    + nu_ * (d2udx2 + d2udy2 + d2udz2)
                );

                // v derivatives
                T dvdx = (grid_.v(i + 1, j, k) - grid_.v(i - 1, j, k)) / (2 * grid_.dx());
                T dvdy = (grid_.v(i, j + 1, k) - grid_.v(i, j - 1, k)) / (2 * grid_.dy());
                T dvdz = (grid_.v(i, j, k + 1) - grid_.v(i, j, k - 1)) / (2 * grid_.dz());
                T d2vdx2 = (grid_.v(i + 1, j, k) - 2 * grid_.v(i, j, k) + grid_.v(i - 1, j, k)) / (grid_.dx() * grid_.dx());
                T d2vdy2 = (grid_.v(i, j + 1, k) - 2 * grid_.v(i, j, k) + grid_.v(i, j - 1, k)) / (grid_.dy() * grid_.dy());
                T d2vdz2 = (grid_.v(i, j, k + 1) - 2 * grid_.v(i, j, k) + grid_.v(i, j, k - 1)) / (grid_.dz() * grid_.dz());

                v_star_[i][j][k] = grid_.v(i, j, k) + dt * (
                    -grid_.u(i, j, k) * dvdx - grid_.v(i, j, k) * dvdy - grid_.w(i, j, k) * dvdz
                    + nu_ * (d2vdx2 + d2vdy2 + d2vdz2)
                );

                // w derivatives
                T dwdx = (grid_.w(i + 1, j, k) - grid_.w(i - 1, j, k)) / (2 * grid_.dx());
                T dwdy = (grid_.w(i, j + 1, k) - grid_.w(i, j - 1, k)) / (2 * grid_.dy());
                T dwdz = (grid_.w(i, j, k + 1) - grid_.w(i, j, k - 1)) / (2 * grid_.dz());
                T d2wdx2 = (grid_.w(i + 1, j, k) - 2 * grid_.w(i, j, k) + grid_.w(i - 1, j, k)) / (grid_.dx() * grid_.dx());
                T d2wdy2 = (grid_.w(i, j + 1, k) - 2 * grid_.w(i, j, k) + grid_.w(i, j - 1, k)) / (grid_.dy() * grid_.dy());
                T d2wdz2 = (grid_.w(i, j, k + 1) - 2 * grid_.w(i, j, k) + grid_.w(i, j, k - 1)) / (grid_.dz() * grid_.dz());

                w_star_[i][j][k] = grid_.w(i, j, k) + dt * (
                    -grid_.u(i, j, k) * dwdx - grid_.v(i, j, k) * dwdy - grid_.w(i, j, k) * dwdz
                    + nu_ * (d2wdx2 + d2wdy2 + d2wdz2)
                    - gravity_ // Gravity in negative z-direction
                );
            }
        }
    }
}

template<typename T>
void NavierStokesSolver<T>::solvePressurePoisson(double dt) {
    std::vector<std::vector<std::vector<T>>> b(
        grid_.nx(), std::vector<std::vector<T>>(grid_.ny(), std::vector<T>(grid_.nz(), 0.0))
    );
    for (size_t i = 1; i < grid_.nx() - 1; ++i) {
        for (size_t j = 1; j < grid_.ny() - 1; ++j) {
            for (size_t k = 1; k < grid_.nz() - 1; ++k) {
                T dudx = (u_star_[i + 1][j][k] - u_star_[i - 1][j][k]) / (2 * grid_.dx());
                T dvdy = (v_star_[i][j + 1][k] - v_star_[i][j - 1][k]) / (2 * grid_.dy());
                T dwdz = (w_star_[i][j][k + 1] - w_star_[i][j][k - 1]) / (2 * grid_.dz());
                b[i][j][k] = (rho_ / dt) * (dudx + dvdy + dwdz);
            }
        }
    }

    std::vector<std::vector<std::vector<T>>> p_new = grid_.p_;
    const int max_iterations = 50;
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<std::vector<std::vector<T>>> p_temp = p_new;
        for (size_t i = 1; i < grid_.nx() - 1; ++i) {
            for (size_t j = 1; j < grid_.ny() - 1; ++j) {
                for (size_t k = 1; k < grid_.nz() - 1; ++k) {
                    p_new[i][j][k] = (1.0 / 6.0) * (
                        p_temp[i + 1][j][k] + p_temp[i - 1][j][k] +
                        p_temp[i][j + 1][k] + p_temp[i][j - 1][k] +
                        p_temp[i][j][k + 1] + p_temp[i][j][k - 1] -
                        grid_.dx() * grid_.dx() * b[i][j][k]
                    );
                }
            }
        }
    }
    grid_.p_ = p_new;
}

template<typename T>
void NavierStokesSolver<T>::correctVelocity(double dt) {
    for (size_t i = 1; i < grid_.nx() - 1; ++i) {
        for (size_t j = 1; j < grid_.ny() - 1; ++j) {
            for (size_t k = 1; k < grid_.nz() - 1; ++k) {
                T dpdx = (grid_.p(i + 1, j, k) - grid_.p(i - 1, j, k)) / (2 * grid_.dx());
                T dpdy = (grid_.p(i, j + 1, k) - grid_.p(i, j - 1, k)) / (2 * grid_.dy());
                T dpdz = (grid_.p(i, j, k + 1) - grid_.p(i, j, k - 1)) / (2 * grid_.dz());
                grid_.u(i, j, k) = u_star_[i][j][k] - (dt / rho_) * dpdx;
                grid_.v(i, j, k) = v_star_[i][j][k] - (dt / rho_) * dpdy;
                grid_.w(i, j, k) = w_star_[i][j][k] - (dt / rho_) * dpdz;
            }
        }
    }
}

// Explicit template instantiation
template class NavierStokesSolver<double>;