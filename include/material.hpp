// material.hpp
//
// Thermophysical constants for a generic low-carbon mild steel, plus a
// set of "initial microstructure" presets. The presets are the direct
// software analogue of the PhD thesis variable: for a fixed laser
// recipe (power, spot size, scan speed), the *same* thermal cycle
// produces a different case hardness and case depth depending on the
// starting grain size / pearlite morphology, because that starting
// state sets the nucleation-site density (and thus the JMAK rate
// constant and exponent) for austenite formation on heating.
//
// The numbers below are illustrative, not measured -- they are chosen
// to reproduce the *qualitative* trend reported in laser-hardening
// literature (finer initial grain / more finely divided pearlite ->
// faster, more complete austenitization at a given peak temperature
// and dwell time -> deeper, more uniform hardened case). Replacing
// them with values regressed from real dilatometry / metallography
// data (of the kind collected in the thesis) is the natural next step
// -- see README.md, "Calibrating against thesis data".

#pragma once

#include <string>
#include <string_view>

namespace laser_sim {

struct SteelProperties {
    double k = 45.0;        // thermal conductivity [W/(m*K)]
    double rho = 7850.0;    // density [kg/m^3]
    double cp = 500.0;      // specific heat [J/(kg*K)]
    double T_ambient = 293.15;   // [K]
    double T_melt = 1783.0;      // [K], solidus of low-carbon steel
    double Ac1 = 996.0;          // [K] (~723 C) start of austenite formation
    double Ac3 = 1163.0;         // [K] (~890 C) fully austenitic above this
    double Ms = 723.0;           // [K] (~450 C) martensite start, low-C steel
    double absorptivity = 0.35;  // fraction of incident Nd:Glass beam absorbed

    // Hardness end-members for the mixture rule [HV].
    double HV_base = 160.0;      // as-received ferrite/pearlite
    double HV_martensite = 620.0; // low-carbon martensite (limited by %C)

    double alpha() const { return k / (rho * cp); }  // thermal diffusivity
};

// Parameters controlling isothermal JMAK austenitization kinetics:
//   f(t) = 1 - exp(-k(T) * t^n),   k(T) = k0 * exp(-Q / (R*T))
// k0 and Q are set per microstructure preset; n (Avrami exponent) reflects
// nucleation mode (grain-boundary nucleation for fine grains -> n closer
// to 1; sluggish, site-limited growth in coarse/banded structures -> n
// closer to 2).
struct MicrostructureParams {
    std::string name;
    double grain_size_um;       // prior-austenite/ferrite grain size
    double pearlite_fraction;   // volume fraction pearlite in starting steel
    double pearlite_spacing_um; // interlamellar spacing (finer -> faster dissolution)
    double k0;                  // JMAK pre-exponential [1/s^n]
    double Q_over_R;            // activation temperature for k(T) [K]
    double n_avrami;            // Avrami exponent
};

// Three initial microstructures, spanning the range typically produced by
// hot-rolling / annealing schedules for the same nominal low-carbon
// composition -- exactly the kind of variation the thesis compared.
inline MicrostructureParams preset_fine_equiaxed() {
    return MicrostructureParams{
        "fine_equiaxed_ferrite_pearlite",
        /*grain_size_um=*/6.0,
        /*pearlite_fraction=*/0.18,
        /*pearlite_spacing_um=*/0.25,
        /*k0=*/6.0e6,
        /*Q_over_R=*/9000.0,
        /*n_avrami=*/1.4,
    };
}

inline MicrostructureParams preset_coarse_equiaxed() {
    return MicrostructureParams{
        "coarse_equiaxed_ferrite_pearlite",
        /*grain_size_um=*/28.0,
        /*pearlite_fraction=*/0.16,
        /*pearlite_spacing_um=*/0.55,
        /*k0=*/1.2e6,
        /*Q_over_R=*/9800.0,
        /*n_avrami=*/1.8,
    };
}

inline MicrostructureParams preset_banded_pearlite() {
    return MicrostructureParams{
        "banded_pearlite",
        /*grain_size_um=*/18.0,
        /*pearlite_fraction=*/0.30,
        /*pearlite_spacing_um=*/0.75,
        /*k0=*/0.6e6,
        /*Q_over_R=*/10200.0,
        /*n_avrami=*/2.1,
    };
}

inline MicrostructureParams preset_by_name(std::string_view name) {
    if (name == "fine") return preset_fine_equiaxed();
    if (name == "coarse") return preset_coarse_equiaxed();
    if (name == "banded") return preset_banded_pearlite();
    return preset_fine_equiaxed();
}

}  // namespace laser_sim
