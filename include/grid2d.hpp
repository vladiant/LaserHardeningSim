// grid2d.hpp
//
// A flat 2D grid over a vertical cross-section of the steel sheet:
//   x -> along the surface, in the laser scan direction   [0, Lx]
//   y -> depth into the sheet, y = 0 is the irradiated surface [0, Ly]
//
// Each node carries the state needed to reproduce (in simplified form)
// the microstructural evolution studied in the PhD work: temperature
// history, progress of austenitization, and the martensite fraction
// that results on quenching.

#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace laser_sim {

struct NodeState {
    double T = 293.15;         // current temperature [K]
    double T_peak = 293.15;    // maximum temperature ever reached [K]
    double f_austenite = 0.0;  // current austenitized fraction [0,1] (JMAK)
    double f_austenite_peak = 0.0; // austenite fraction "locked in" at peak
    double f_martensite = 0.0; // martensite fraction after quench [0,1]
    bool   melted = false;     // true if T_peak exceeded the melting point
};

class Grid2D {
public:
    Grid2D(std::size_t nx, std::size_t ny, double dx, double dy, double T0)
        : nx_(nx), ny_(ny), dx_(dx), dy_(dy), nodes_(nx * ny, NodeState{T0, T0, 0, 0, 0, false}) {}

    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    std::size_t size() const { return nodes_.size(); }

    std::size_t index(std::size_t ix, std::size_t iy) const { return iy * nx_ + ix; }
    std::size_t ix_of(std::size_t idx) const { return idx % nx_; }
    std::size_t iy_of(std::size_t idx) const { return idx / nx_; }

    NodeState& at(std::size_t ix, std::size_t iy) { return nodes_[index(ix, iy)]; }
    const NodeState& at(std::size_t ix, std::size_t iy) const { return nodes_[index(ix, iy)]; }

    std::span<NodeState> nodes() { return nodes_; }
    std::span<const NodeState> nodes() const { return nodes_; }

    // Physical coordinates of a node, cell-centered.
    double x_of(std::size_t ix) const { return (static_cast<double>(ix) + 0.5) * dx_; }
    double y_of(std::size_t iy) const { return (static_cast<double>(iy) + 0.5) * dy_; }

private:
    std::size_t nx_, ny_;
    double dx_, dy_;
    std::vector<NodeState> nodes_;
};

}  // namespace laser_sim
