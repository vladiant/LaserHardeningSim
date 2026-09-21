# Laser Surface Hardening Simulator (C++20)

A small, self-contained simulation of Nd:Glass laser transformation hardening of a low-carbon mild steel sheet -- built as a direct software analogue of the question at the center of the thesis *"Influence of the initial microstructure on the level of modification of properties of low carbon mild sheet steel during laser treatment with Nd:Glass laser."*

It solves the transient 2D heat-conduction problem for a scanning laser spot, tracks phase transformation kinetics at every grid point, and reports the resulting hardness / case-depth profile -- for a choice of **initial microstructure**, which is the thesis's central independent variable.

## Why this is the right computational bridge

The thesis measured how a given Nd:Glass laser treatment produces different hardening outcomes depending on the steel's starting microstructure (grain size, pearlite fraction and morphology). That is fundamentally a statement about **transformation kinetics**: a finer initial grain structure offers more ferrite/pearlite grain-boundary area, i.e. more nucleation sites for austenite, so at a given temperature and dwell time it austenitizes faster and more completely than a coarse or banded structure -- which then governs how much martensite (and therefore how much hardness) the subsequent rapid self-quench can produce.

The simulation reproduces exactly that causal chain, numerically:

```
initial microstructure  ->  JMAK rate constant / exponent
                          -> austenitization completeness at each node's local peak-temperature dwell
                          -> martensite fraction on quench (Koistinen-Marburger)
                          -> hardness (rule of mixtures)
```

with the laser's power, spot size, and scan speed held fixed across runs, so any difference in the output is attributable to the microstructure parameter alone -- the same controlled-variable logic as the original experimental thesis work, just run as a numerical experiment instead of a metallography one.

## Physical model

- **Heat conduction**: explicit (FTCS) finite-difference solution of `rho*cp*dT/dt = k*(d2T/dx2 + d2T/dy2)` on a 2D cross-section (x = scan direction, y = depth). See `include/heat_solver.hpp`.
- **Laser source**: a moving Gaussian surface flux (`include/laser_source.hpp`), with a `--mode pulsed` option that reproduces the overlapping discrete-spot pattern typical of historical pulsed Nd:Glass hardening rigs (as opposed to CW).
- **Boundary conditions**: the sheet loses heat to ambient air by convection + radiation on *all four* exposed faces (it's a thin sheet, not a semi-infinite block); the top face additionally receives the absorbed laser flux. This was a deliberate, verified choice -- see "A numerical pitfall worth knowing about" below.
- **Austenitization**: incremental Johnson-Mehl-Avrami-Kolmogorov (JMAK) kinetics, `f(t) = 1 - exp(-k(T)*t^n)`, applied with the standard non-isothermal additivity rule so a changing temperature history is handled correctly (`include/kinetics.hpp`).
- **Martensite formation**: Koistinen-Marburger relation on cooling through Ms, applied to the fraction of austenite each node reached at its own peak temperature.
- **Hardness**: linear rule-of-mixtures between a base ferrite/pearlite hardness and a martensite hardness end-member.

Three **initial microstructure presets** are provided in `include/material.hpp` (`fine`, `coarse`, `banded`), each with its own JMAK rate constant, activation temperature, and Avrami exponent, loosely reflecting the grain-boundary nucleation density and pearlite interlamellar spacing you'd characterize by metallography.

These parameter values are illustrative defaults, not measurements. The one place where a numerical model like this earns its keep is if its rate constants are regressed against real time-temperature and hardness data -- see "Calibrating against thesis data" below.

## A numerical pitfall worth knowing about (and why the code avoids it)

Early in development, the solver used a "mirror ghost node" for insulated boundaries -- reflecting the interior neighbor's temperature across the boundary. That looks reasonable but silently **double-counts flux at every boundary face**, injecting energy from nowhere. It only showed up as a slow drift once heat actually reached a domain edge.

The `tests/sanity_checks.cpp` energy-conservation check exists specifically to catch this class of bug: it runs the solver on a fully insulated grid and asserts the total temperature sum stays constant to floating-point tolerance. The fix was to use the correct finite-volume closure (drop the missing neighbor term entirely, rather than mirroring it) -- the kind of bug that is easy to introduce and easy to miss without an explicit conservation check, which is exactly why the check is there and why the accompanying CMake target runs it as `ctest`.

A related, more physical pitfall: a first version used adiabatic (zero-loss) side/bottom boundaries, which is the right approximation for a *bulk, semi-infinite* solid but wrong for a thin, finite *sheet* being scanned end-to-end -- the whole sheet's background temperature crept upward over the course of the scan with no way to shed the input energy. Letting every exposed face lose heat by convection and radiation (as a real sheet in air does) fixed it.

## Building

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build          # runs the sanity checks
```

CMake looks for Intel TBB (`find_package(TBB)`); if found, the `--backend stl` option becomes a real `std::execution::par_unseq` parallel run. If not found, `--backend stl` transparently falls back to the jthread backend, and the build still succeeds. Note that on libstdc++, `std::execution::par_unseq` requires TBB to actually run in parallel; it silently degrades to serial without it, which is exactly why this project treats the "stl" backend as opt-in and pairs it with a hand-rolled `std::jthread` fan-out/fan-in backend that has no external dependency (`include/parallel_backend.hpp`). Both operate on the same stencil kernel, so `--compare-backends` is an apples-to-apples comparison of "roll your own" vs "let the standard library do it" parallelism.

## Running

```bash
# A representative continuous-wave scan, no melting, ~0.06 mm case depth:
./build/laser_hardening_sim --microstructure fine

# Same laser recipe, coarse starting grain structure -- compare the
# hardness_map.csv / case_depth_profile.csv against the run above:
./build/laser_hardening_sim --microstructure coarse

# Pulsed-mode Nd:Glass operation (closer to the historical rig):
./build/laser_hardening_sim --mode pulsed --intensity 9e8

# Compare the jthread / stl / serial backends on identical physics:
./build/laser_hardening_sim --compare-backends
```

Full option list is documented at the top of `src/main.cpp`
(`--microstructure`, `--mode`, `--intensity`, `--spot-radius-um`,
`--scan-speed-mm-s`, `--length-mm`, `--thickness-mm`, `--nx`, `--ny`,
`--backend`, `--threads`, `--output-dir`, `--compare-backends`,
`--martensite-threshold`).

Each run writes two CSVs to `--output-dir` (default `./output`):

- `hardness_map.csv` -- per-node peak temperature, austenite fraction, martensite fraction, hardness, and melt flag.
- `case_depth_profile.csv` -- case depth and surface hardness as a function of scan position.

`scripts/plot_hardness.py output_dir` (needs numpy + matplotlib, not required to build or run the simulator itself) renders both as PNGs.

### Where the microstructure effect is visible

At the default laser settings, most of the near-surface material heats well above Ac3 and austenitizes essentially instantly in all three presets, so the *peak* hardness and case depth look similar. The microstructure signal is concentrated exactly where the physics says it should be: in the **intercritical band** (Ac1 < T < Ac3) at greater depth or off-center in x, where the JMAK kinetics haven't saturated and the rate constant actually controls the outcome. Compare the `f_austenite_peak` and `hardness_HV` columns across presets at those grid points (or lower `--intensity` toward ~1.0e9 W/m^2 with `--scan-speed-mm-s 400` to push most of the hardened zone into that kinetically-limited regime, making the preset-to-preset difference visible directly in `case_depth_profile.csv`).

## Calibrating against thesis data

To turn this from a qualitative demonstrator into something that reproduces the thesis's actual numbers, the natural next step is to replace the illustrative constants in `include/material.hpp` (`k0`, `Q_over_R`, `n_avrami` per preset, and the `SteelProperties` thermal constants) with values regressed from real dilatometry/metallography data of the kind collected in the thesis work -- e.g. fitting `k0` and `n` per starting microstructure against measured case-depth vs. laser-parameter curves, and cross-checking the resulting hardness predictions against Vickers traverses across the heat-affected zone.

## Limitations

- 2D cross-section, not a full 3D volume (reasonable for a scan long compared to the spot size, less so near the very start/end of a pass or for a stationary spot).
- Constant thermal properties (no temperature dependence of k/rho/cp, no latent heat of melting/solidification modeled beyond flagging melted nodes).
- JMAK/Koistinen-Marburger with additivity-rule handling of a varying temperature history is a standard engineering approximation, not a CALPHAD-grade transformation model; retained austenite and tempering effects during multi-pass overlap are not modeled.
- Explicit FTCS timestepping is simple and easy to reason about, but restricts the timestep via the CFL-like stability bound the code computes automatically; an implicit scheme would allow larger steps at the cost of more code.
