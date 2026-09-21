// heat_solver.hpp
//
// Explicit (FTCS) finite-difference solver for 2D transient conduction
//     rho*cp*dT/dt = k*(d2T/dx2 + d2T/dy2)
// over the Grid2D domain, with:
//   - top boundary (y=0): Neumann flux BC combining the absorbed laser
//     flux and convective + radiative losses to ambient, handled with a
//     second-order ghost-node.
//   - all other boundaries: adiabatic (zero-gradient / mirror ghost),
//     appropriate for a domain large enough that the heat-affected zone
//     does not reach the far edges within the simulated scan.
//
// Each timestep does two disjoint parallel passes over all nodes (a
// Jacobi-style update, required for an explicit scheme to remain
// race-free under `parallel_for`):
//   1. compute the new temperature field from the *old* field,
//   2. commit it and advance each node's transformation kinetics.

#pragma once

#include <cmath>
#include <vector>

#include "grid2d.hpp"
#include "kinetics.hpp"
#include "laser_source.hpp"
#include "material.hpp"
#include "parallel_backend.hpp"

namespace laser_sim {

constexpr double kStefanBoltzmann = 5.670374419e-8; // [W/(m^2 K^4)]

struct SolverConfig {
    double h_conv = 15.0;     // convective coefficient to ambient air [W/(m^2 K)]
    double emissivity = 0.7;  // surface emissivity for radiative loss
};

class HeatSolver {
public:
    HeatSolver(const SteelProperties& steel, const LaserSource& laser,
               const MicrostructureParams& microstructure, SolverConfig cfg = {})
        : steel_(steel), laser_(laser), micro_(microstructure), cfg_(cfg) {}

    // Maximum stable explicit timestep for the given grid (2D FTCS stability).
    double max_stable_dt(const Grid2D& grid) const {
        const double alpha = steel_.alpha();
        const double inv_dx2 = 1.0 / (grid.dx() * grid.dx());
        const double inv_dy2 = 1.0 / (grid.dy() * grid.dy());
        return 0.5 / (alpha * (inv_dx2 + inv_dy2));
    }

    void step(Grid2D& grid, double t, double dt, ParBackend backend) {
        const std::size_t n = grid.size();
        if (T_new_.size() != n) T_new_.assign(n, steel_.T_ambient);

        const double alpha = steel_.alpha();
        const double inv_dx2 = 1.0 / (grid.dx() * grid.dx());
        const double inv_dy2 = 1.0 / (grid.dy() * grid.dy());
        const std::size_t nx = grid.nx();
        const std::size_t ny = grid.ny();

        // Pass 1: compute new temperatures from the current field only.
        //
        // Finite-volume form, cell-centered: at a boundary face the
        // conduction term simply drops the missing neighbor (equivalent
        // to a ghost cell equal to the boundary node itself), which is
        // exactly conservative -- it must NOT mirror the interior
        // neighbor's value, which would double-count flux and slowly
        // inject energy at the boundary (caught by the energy-
        // conservation sanity check during development).
        //
        // The sheet is thin enough that all four exposed faces --not
        // just the irradiated top-- lose heat to ambient air by
        // convection and radiation; only the top face additionally
        // receives the laser flux. Treating the domain as a thin sheet
        // radiating from every face (rather than a semi-infinite,
        // adiabatic-sided block) avoids an unphysical cumulative
        // temperature rise as the beam scans across a finite grid.
        parallel_for(n, backend, [&](std::size_t idx) {
            const std::size_t ix = grid.ix_of(idx);
            const std::size_t iy = grid.iy_of(idx);
            const double T_here = grid.at(ix, iy).T;

            auto surface_loss = [&](double T) {
                const double q_conv = cfg_.h_conv * (T - steel_.T_ambient);
                const double q_rad = cfg_.emissivity * kStefanBoltzmann *
                                      (std::pow(T, 4) - std::pow(steel_.T_ambient, 4));
                return q_conv + q_rad; // [W/m^2], positive = leaving the sheet
            };

            double d2Tdx2, source_x = 0.0;
            if (ix == 0) {
                d2Tdx2 = (grid.at(ix + 1, iy).T - T_here) * inv_dx2;
                source_x = -surface_loss(T_here) / (steel_.rho * steel_.cp * grid.dx());
            } else if (ix == nx - 1) {
                d2Tdx2 = (grid.at(ix - 1, iy).T - T_here) * inv_dx2;
                source_x = -surface_loss(T_here) / (steel_.rho * steel_.cp * grid.dx());
            } else {
                d2Tdx2 = (grid.at(ix + 1, iy).T - 2.0 * T_here + grid.at(ix - 1, iy).T) * inv_dx2;
            }

            double d2Tdy2, source_y = 0.0;
            if (iy == 0) {
                d2Tdy2 = (grid.at(ix, iy + 1).T - T_here) * inv_dy2;
                const double x = grid.x_of(ix);
                const double q_in = laser_.flux(x, t, steel_.absorptivity);
                source_y = (q_in - surface_loss(T_here)) / (steel_.rho * steel_.cp * grid.dy());
            } else if (iy == ny - 1) {
                d2Tdy2 = (grid.at(ix, iy - 1).T - T_here) * inv_dy2;
                source_y = -surface_loss(T_here) / (steel_.rho * steel_.cp * grid.dy());
            } else {
                d2Tdy2 = (grid.at(ix, iy + 1).T - 2.0 * T_here + grid.at(ix, iy - 1).T) * inv_dy2;
            }

            T_new_[idx] = T_here + dt * (alpha * (d2Tdx2 + d2Tdy2) + source_x + source_y);
        }, hw_threads_);

        // Pass 2: commit temperatures and advance transformation kinetics.
        parallel_for(n, backend, [&](std::size_t idx) {
            const std::size_t ix = grid.ix_of(idx);
            const std::size_t iy = grid.iy_of(idx);
            NodeState& node = grid.at(ix, iy);

            node.T = T_new_[idx];
            node.T_peak = std::max(node.T_peak, node.T);
            if (node.T_peak >= steel_.T_melt) node.melted = true;

            advance_austenitization(node, node.T, dt, steel_, micro_);
            node.f_austenite_peak = std::max(node.f_austenite_peak, node.f_austenite);
            advance_martensite(node, node.T, steel_);
        }, hw_threads_);
    }

    void set_thread_count(unsigned n) { hw_threads_ = n; }

private:
    const SteelProperties& steel_;
    const LaserSource& laser_;
    const MicrostructureParams& micro_;
    SolverConfig cfg_;
    std::vector<double> T_new_;
    unsigned hw_threads_ = 0; // 0 => hardware_concurrency()
};

}  // namespace laser_sim
