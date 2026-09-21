// main.cpp
//
// Simulates one pass of a scanning Nd:Glass laser hardening beam over a
// low-carbon steel sheet cross-section, for a chosen initial
// microstructure, and reports the resulting hardness/case-depth profile.
//
// Usage:
//   laser_hardening_sim [options]
//
// Options (all optional, defaults shown):
//   --microstructure fine|coarse|banded   (fine)
//   --mode continuous|pulsed              (continuous)
//   --intensity <W/m^2 peak>              (1.8e9)
//   --spot-radius-um <um>                 (400)
//   --scan-speed-mm-s <mm/s>              (450)
//   --length-mm <mm>          domain length along scan (6.0)
//   --thickness-mm <mm>       domain depth (2.0)
//   --nx <int>                grid cells along scan (240)
//   --ny <int>                grid cells through thickness (80)
//   --backend jthread|stl|serial          (jthread)
//   --threads <int>           0 = hardware_concurrency()  (0)
//   --output-dir <path>                   (./output)
//   --compare-backends        also time Serial vs JThread vs Stl and print a table
//   --martensite-threshold <0..1>  case-depth cutoff (0.5)

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "grid2d.hpp"
#include "heat_solver.hpp"
#include "kinetics.hpp"
#include "laser_source.hpp"
#include "material.hpp"
#include "parallel_backend.hpp"

using namespace laser_sim;
namespace fs = std::filesystem;

struct Args {
    std::string microstructure = "fine";
    std::string mode = "continuous";
    double intensity = 1.8e9;         // W/m^2 peak, pre-absorption
    double spot_radius_um = 400.0;
    double scan_speed_mm_s = 450.0;
    double length_mm = 6.0;
    double thickness_mm = 2.0;
    int nx = 240;
    int ny = 80;
    std::string backend = "jthread";
    unsigned threads = 0;
    std::string output_dir = "./output";
    bool compare_backends = false;
    double martensite_threshold = 0.5;
};

double parse_double(std::string_view s) { return std::strtod(std::string(s).c_str(), nullptr); }
int parse_int(std::string_view s) { return std::atoi(std::string(s).c_str()); }

Args parse_args(int argc, char** argv) {
    Args a;
    std::vector<std::string> tok(argv + 1, argv + argc);
    auto next = [&](std::size_t& i) -> std::string { return tok[++i]; };
    for (std::size_t i = 0; i < tok.size(); ++i) {
        const std::string& s = tok[i];
        if (s == "--microstructure") a.microstructure = next(i);
        else if (s == "--mode") a.mode = next(i);
        else if (s == "--intensity") a.intensity = parse_double(next(i));
        else if (s == "--spot-radius-um") a.spot_radius_um = parse_double(next(i));
        else if (s == "--scan-speed-mm-s") a.scan_speed_mm_s = parse_double(next(i));
        else if (s == "--length-mm") a.length_mm = parse_double(next(i));
        else if (s == "--thickness-mm") a.thickness_mm = parse_double(next(i));
        else if (s == "--nx") a.nx = parse_int(next(i));
        else if (s == "--ny") a.ny = parse_int(next(i));
        else if (s == "--backend") a.backend = next(i);
        else if (s == "--threads") a.threads = static_cast<unsigned>(parse_int(next(i)));
        else if (s == "--output-dir") a.output_dir = next(i);
        else if (s == "--compare-backends") a.compare_backends = true;
        else if (s == "--martensite-threshold") a.martensite_threshold = parse_double(next(i));
        else if (s == "--help") {
            std::cout << "See top of main.cpp for full option list.\n";
            std::exit(0);
        }
    }
    return a;
}

ParBackend backend_from_string(std::string_view s) {
    if (s == "stl") return ParBackend::StlParUnseq;
    if (s == "serial") return ParBackend::Serial;
    return ParBackend::JThread;
}

// Runs the coupled thermal + kinetics loop to completion; returns wall-clock seconds.
// Operates on a *copy* of the initial grid so repeated calls (e.g. for backend
// comparison) all start from the same initial condition.
double run_time_loop(Grid2D grid, const LaserSource& laser, const SteelProperties& steel,
                      const MicrostructureParams& micro, double dt, double t_end,
                      ParBackend backend, unsigned threads, Grid2D* result_out) {
    HeatSolver solver(steel, laser, micro);
    solver.set_thread_count(threads);

    const auto t0 = std::chrono::steady_clock::now();
    for (double t = 0.0; t < t_end; t += dt) {
        solver.step(grid, t, dt, backend);
    }
    const auto t1 = std::chrono::steady_clock::now();

    if (result_out) *result_out = std::move(grid);
    return std::chrono::duration<double>(t1 - t0).count();
}

void write_hardness_map(const fs::path& path, const Grid2D& grid, const SteelProperties& steel) {
    std::ofstream f(path);
    f << "x_mm,y_mm,T_peak_C,f_austenite_peak,f_martensite,hardness_HV,melted\n";
    for (std::size_t iy = 0; iy < grid.ny(); ++iy) {
        for (std::size_t ix = 0; ix < grid.nx(); ++ix) {
            const NodeState& n = grid.at(ix, iy);
            f << (grid.x_of(ix) * 1e3) << ',' << (grid.y_of(iy) * 1e3) << ','
              << (n.T_peak - 273.15) << ',' << n.f_austenite_peak << ',' << n.f_martensite
              << ',' << estimate_hardness(n, steel) << ',' << (n.melted ? 1 : 0) << '\n';
        }
    }
}

void write_case_depth_profile(const fs::path& path, const Grid2D& grid,
                               const SteelProperties& steel, double martensite_threshold) {
    std::ofstream f(path);
    f << "x_mm,case_depth_mm,surface_hardness_HV,surface_T_peak_C\n";
    for (std::size_t ix = 0; ix < grid.nx(); ++ix) {
        double depth_mm = 0.0;
        for (std::size_t iy = 0; iy < grid.ny(); ++iy) {
            if (grid.at(ix, iy).f_martensite >= martensite_threshold) {
                depth_mm = grid.y_of(iy) * 1e3;
            } else {
                break; // fractions are (in practice) monotonically decreasing with depth
            }
        }
        const NodeState& surf = grid.at(ix, 0);
        f << (grid.x_of(ix) * 1e3) << ',' << depth_mm << ',' << estimate_hardness(surf, steel)
          << ',' << (surf.T_peak - 273.15) << '\n';
    }
}

int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);

    SteelProperties steel{};
    MicrostructureParams micro = preset_by_name(args.microstructure);

    const double Lx = args.length_mm * 1e-3;
    const double Ly = args.thickness_mm * 1e-3;
    const auto nx = static_cast<std::size_t>(args.nx);
    const auto ny = static_cast<std::size_t>(args.ny);
    const double dx = Lx / static_cast<double>(nx);
    const double dy = Ly / static_cast<double>(ny);

    LaserSource laser;
    laser.mode = (args.mode == "pulsed") ? LaserMode::Pulsed : LaserMode::Continuous;
    laser.peak_power_density = args.intensity;
    laser.spot_radius = args.spot_radius_um * 1e-6;
    laser.scan_speed = args.scan_speed_mm_s * 1e-3;
    laser.y0 = 4.0 * laser.spot_radius; // start with the spot fully inside the domain

    const double margin = 4.0 * laser.spot_radius;
    const double scan_time = (Lx - 2.0 * margin) / laser.scan_speed;
    const double cooldown = scan_time * 0.5; // let quenching finish after beam exits
    double t_end = scan_time + cooldown;
    if (laser.mode == LaserMode::Pulsed) {
        const int n_pulses = static_cast<int>((Lx - 2.0 * margin) / laser.pulse_pitch);
        t_end = n_pulses * laser.pulse_period + cooldown;
    }

    Grid2D probe_grid(nx, ny, dx, dy, steel.T_ambient);
    HeatSolver probe_solver(steel, laser, micro);
    const double dt = 0.8 * probe_solver.max_stable_dt(probe_grid);

    std::cout << "Grid: " << nx << " x " << ny << " (dx=" << dx * 1e6 << " um, dy=" << dy * 1e6
              << " um)\n";
    std::cout << "Microstructure preset: " << micro.name << "\n";
    std::cout << "Stable dt = " << dt << " s, simulated time = " << t_end << " s ("
              << static_cast<long>(t_end / dt) << " steps)\n";

    if (args.compare_backends) {
        std::cout << "\n-- Backend comparison (full run, no I/O) --\n";
        struct Entry { const char* label; ParBackend b; };
        std::vector<Entry> entries = {
            {"serial ", ParBackend::Serial},
            {"jthread", ParBackend::JThread},
#if defined(HAVE_STL_PARALLEL)
            {"stl par", ParBackend::StlParUnseq},
#endif
        };
        for (auto& e : entries) {
            Grid2D g(nx, ny, dx, dy, steel.T_ambient);
            double secs = run_time_loop(std::move(g), laser, steel, micro, dt, t_end, e.b,
                                         args.threads, nullptr);
            std::cout << "  " << e.label << " : " << secs << " s\n";
        }
        std::cout << "\n";
    }

    Grid2D final_grid(nx, ny, dx, dy, steel.T_ambient);
    const ParBackend backend = backend_from_string(args.backend);
    const double elapsed = run_time_loop(std::move(final_grid), laser, steel, micro, dt, t_end,
                                          backend, args.threads, &final_grid);

    std::cout << "Run (" << args.backend << ", " << (args.threads ? args.threads
                                                                   : std::thread::hardware_concurrency())
              << " threads): " << elapsed << " s\n";

    fs::create_directories(args.output_dir);
    write_hardness_map(fs::path(args.output_dir) / "hardness_map.csv", final_grid, steel);
    write_case_depth_profile(fs::path(args.output_dir) / "case_depth_profile.csv", final_grid,
                              steel, args.martensite_threshold);

    double max_case_depth = 0.0, max_surface_hv = 0.0, max_T_peak_C = -273.15;
    bool any_melted = false;
    for (std::size_t ix = 0; ix < final_grid.nx(); ++ix) {
        double depth_mm = 0.0;
        for (std::size_t iy = 0; iy < final_grid.ny(); ++iy) {
            const NodeState& n = final_grid.at(ix, iy);
            if (n.melted) any_melted = true;
            max_T_peak_C = std::max(max_T_peak_C, n.T_peak - 273.15);
            if (n.f_martensite >= args.martensite_threshold) depth_mm = final_grid.y_of(iy) * 1e3;
        }
        max_case_depth = std::max(max_case_depth, depth_mm);
        max_surface_hv = std::max(max_surface_hv, estimate_hardness(final_grid.at(ix, 0), steel));
    }

    std::cout << "\n-- Summary --\n"
              << "Peak temperature reached : " << max_T_peak_C << " C\n"
              << "Surface melted           : " << (any_melted ? "yes (reduce intensity or speed up scan)" : "no") << "\n"
              << "Max case depth (>= " << args.martensite_threshold << " martensite) : "
              << max_case_depth << " mm\n"
              << "Max surface hardness     : " << max_surface_hv << " HV\n"
              << "Outputs written to       : " << fs::absolute(args.output_dir).string() << "\n";

    return 0;
}
