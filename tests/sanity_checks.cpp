// sanity_checks.cpp
//
// Not a full unit-test framework -- a handful of assert-based physical
// sanity checks, meant to catch the most common ways a solver like this
// silently goes wrong:
//   1. Energy conservation under insulated (zero-flux) boundaries.
//   2. JMAK austenite fraction is monotonically non-decreasing with time
//      at constant super-Ac1 temperature, and saturates towards 1.
//   3. Koistinen-Marburger gives 0 martensite at T == Ms and approaches
//      full transformation of the available austenite well below Ms.

#include <cassert>
#include <cmath>
#include <iostream>

#include "grid2d.hpp"
#include "heat_solver.hpp"
#include "kinetics.hpp"
#include "laser_source.hpp"
#include "material.hpp"

using namespace laser_sim;

void check_energy_conservation() {
    SteelProperties steel{};
    MicrostructureParams micro = preset_fine_equiaxed();
    LaserSource laser;
    laser.peak_power_density = 0.0; // no input flux
    SolverConfig cfg;
    cfg.h_conv = 0.0;      // no convective loss
    cfg.emissivity = 0.0;  // no radiative loss -> fully insulated domain

    const std::size_t nx = 20, ny = 20;
    const double dx = 5e-5, dy = 5e-5;
    Grid2D grid(nx, ny, dx, dy, steel.T_ambient);

    // Seed a hot spot in the interior, away from all boundaries.
    grid.at(nx / 2, ny / 2).T = steel.T_ambient + 500.0;
    grid.at(nx / 2, ny / 2).T_peak = grid.at(nx / 2, ny / 2).T;

    auto total_energy = [&](const Grid2D& g) {
        double e = 0.0;
        for (const auto& n : g.nodes()) e += n.T; // proportional to energy (const rho,cp,cell volume)
        return e;
    };

    const double e0 = total_energy(grid);

    HeatSolver solver(steel, laser, micro, cfg);
    const double dt = 0.8 * solver.max_stable_dt(grid);
    for (int i = 0; i < 500; ++i) {
        solver.step(grid, i * dt, dt, ParBackend::JThread);
    }

    const double e1 = total_energy(grid);
    const double rel_drift = std::fabs(e1 - e0) / e0;

    std::cout << "[energy]   initial=" << e0 << " final=" << e1
              << " rel_drift=" << rel_drift << "\n";
    assert(rel_drift < 1e-6 && "insulated domain must conserve total temperature-sum to FP tolerance");
}

void check_jmak_monotonic_and_saturating() {
    SteelProperties steel{};
    MicrostructureParams micro = preset_fine_equiaxed();
    NodeState node{};
    const double T = steel.Ac1 + 40.0; // comfortably between Ac1 and Ac3
    const double dt = 0.01;

    double prev = -1.0;
    for (int i = 0; i < 2000; ++i) {
        advance_austenitization(node, T, dt, steel, micro);
        assert(node.f_austenite >= prev - 1e-12 && "JMAK fraction must be non-decreasing in time");
        prev = node.f_austenite;
    }
    std::cout << "[jmak]     f_austenite after 20s at Ac1+40K = " << node.f_austenite << "\n";
    assert(node.f_austenite > 0.9 && "JMAK should be close to saturation after many seconds above Ac1");
}

void check_koistinen_marburger_bounds() {
    SteelProperties steel{};
    NodeState node{};
    node.f_austenite_peak = 1.0;

    advance_martensite(node, steel.Ms, steel); // exactly at Ms -> zero transformation
    std::cout << "[km]       f_martensite at T=Ms   = " << node.f_martensite << "\n";
    assert(node.f_martensite < 1e-9);

    advance_martensite(node, steel.Ms - 300.0, steel); // well below Ms
    std::cout << "[km]       f_martensite at Ms-300K = " << node.f_martensite << "\n";
    assert(node.f_martensite > 0.9 && node.f_martensite <= 1.0);
}

int main() {
    check_energy_conservation();
    check_jmak_monotonic_and_saturating();
    check_koistinen_marburger_bounds();
    std::cout << "All sanity checks passed.\n";
    return 0;
}
