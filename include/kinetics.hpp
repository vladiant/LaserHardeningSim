// kinetics.hpp
//
// Two classical solid-state transformation relations, applied per grid
// node using each node's own thermal history:
//
// 1) Johnson-Mehl-Avrami-Kolmogorov (JMAK), used here in incremental
//    (non-isothermal, "additivity rule") form to advance the austenite
//    fraction while a node sits between Ac1 and Ac3:
//
//       f(t) = 1 - exp(-k(T) * t^n)
//       k(T) = k0 * exp(-Q/(R*T))
//
//    At each timestep we find the "equivalent time" t_eq that would have
//    produced the current f at the current T, advance it by dt, and
//    recompute f -- the standard additivity-rule trick for handling a
//    continuously changing temperature with an isothermal kinetics law.
//
// 2) Koistinen-Marburger (K-M), applied once a node cools back through
//    Ms, converting the austenite fraction present at peak temperature
//    into martensite:
//
//       f_M(T) = f_gamma_peak * (1 - exp(-alpha * (Ms - T)))   for T < Ms
//
// Both are simplified, additivity-rule engineering approximations, not
// a full CALPHAD-grade transformation model -- adequate for exploring
// *trends* with initial microstructure, which is the point here.

#pragma once

#include <algorithm>
#include <cmath>

#include "grid2d.hpp"
#include "material.hpp"

namespace laser_sim {

constexpr double kKoistinenMarburgerAlpha = 0.011; // [1/K], standard literature value

// Advance node.f_austenite given the node is currently at temperature T
// (Kelvin) and has been for a timestep dt (seconds), using microstructure
// parameters mp. Only called while Ac1 <= T (below Ac1 nothing happens;
// above Ac3 we just clamp to fully austenitic).
inline void advance_austenitization(NodeState& node, double T, double dt,
                                     const SteelProperties& steel,
                                     const MicrostructureParams& mp) {
    if (T < steel.Ac1) return;

    if (T >= steel.Ac3) {
        node.f_austenite = 1.0;
        return;
    }

    const double k = mp.k0 * std::exp(-mp.Q_over_R / T);
    const double n = mp.n_avrami;
    const double f = std::clamp(node.f_austenite, 0.0, 0.999999);

    // Equivalent time at the current T that would already have produced f.
    double t_eq = 0.0;
    if (f > 0.0 && k > 0.0) {
        t_eq = std::pow(-std::log(1.0 - f) / k, 1.0 / n);
    }

    const double t_new = t_eq + dt;
    const double f_new = 1.0 - std::exp(-k * std::pow(t_new, n));
    node.f_austenite = std::clamp(f_new, 0.0, 1.0);
}

// Apply Koistinen-Marburger martensite formation as the node cools
// through Ms. Call once per timestep during cooling; it is written to be
// monotonic non-decreasing so repeated calls as T keeps dropping are safe.
inline void advance_martensite(NodeState& node, double T, const SteelProperties& steel) {
    if (T >= steel.Ms) return;

    const double f_km = 1.0 - std::exp(-kKoistinenMarburgerAlpha * (steel.Ms - T));
    const double f_m = node.f_austenite_peak * std::clamp(f_km, 0.0, 1.0);
    node.f_martensite = std::max(node.f_martensite, f_m);
}

// Rule-of-mixtures hardness estimate from final phase fractions.
inline double estimate_hardness(const NodeState& node, const SteelProperties& steel) {
    const double f_m = std::clamp(node.f_martensite, 0.0, 1.0);
    return (1.0 - f_m) * steel.HV_base + f_m * steel.HV_martensite;
}

}  // namespace laser_sim
