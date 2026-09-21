// laser_source.hpp
//
// Surface heat flux from a Nd:Glass laser scanning across the sheet.
// Nd:Glass systems of the type used in the thesis era were typically
// operated in pulsed (not CW) mode; two modes are offered here:
//
//   Continuous : a Gaussian spot moving at constant scan speed v.
//   Pulsed     : the beam is stepped between discrete positions and
//                fires a fixed-duration pulse at each, reproducing the
//                overlapping-spot-weld pattern typical of pulsed
//                surface hardening trials.

#pragma once

#include <cmath>

namespace laser_sim {

enum class LaserMode { Continuous, Pulsed };

struct LaserSource {
    LaserMode mode = LaserMode::Continuous;

    double peak_power_density = 0.0; // incident, pre-absorption [W/m^2] at beam center
    double spot_radius = 5.0e-4;     // 1/e^2 Gaussian radius [m]
    double scan_speed = 0.02;        // [m/s], used in Continuous mode
    double y0 = 0.0;                 // start x-position of the scan [m]

    // Pulsed-mode parameters.
    double pulse_duration = 2.0e-3;  // [s], laser "on" time per spot
    double pulse_period = 8.0e-3;    // [s], time between pulse starts (includes step)
    double pulse_pitch = 4.0e-4;     // [m], center-to-center spacing between pulses

    // Absorbed surface flux [W/m^2] at position x and time t.
    double flux(double x, double t, double absorptivity) const {
        double x_center;
        double on_fraction = 1.0; // 1 = beam firing, 0 = beam off (between pulses)

        if (mode == LaserMode::Continuous) {
            x_center = y0 + scan_speed * t;
        } else {
            int pulse_index = static_cast<int>(t / pulse_period);
            double t_in_pulse = t - pulse_index * pulse_period;
            on_fraction = (t_in_pulse <= pulse_duration) ? 1.0 : 0.0;
            x_center = y0 + pulse_index * pulse_pitch;
        }

        const double dx = x - x_center;
        const double gaussian = std::exp(-2.0 * dx * dx / (spot_radius * spot_radius));
        return absorptivity * peak_power_density * gaussian * on_fraction;
    }
};

}  // namespace laser_sim
