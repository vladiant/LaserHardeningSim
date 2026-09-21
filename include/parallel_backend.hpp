// parallel_backend.hpp
//
// Two interchangeable parallel-for backends, selected at compile time:
//
//   - JTHREAD:   manual std::jthread fan-out/fan-in, no external deps.
//   - STL_PAR:   std::for_each(std::execution::par_unseq, ...), requires
//                a parallel STL backend (Intel TBB on libstdc++). Enabled
//                by defining HAVE_STL_PARALLEL (CMake does this when TBB
//                is found).
//
// The point of keeping both is pedagogical: the CMake build can toggle
// PAR_BACKEND=jthread|stl and the two give an apples-to-apples comparison
// of "roll your own" vs "let the standard library do it" parallelism on
// the same stencil kernel.

#pragma once

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#if defined(HAVE_STL_PARALLEL)
#include <execution>
#endif

namespace laser_sim {

enum class ParBackend { Serial, JThread, StlParUnseq };

// Runs f(i) for i in [0, n) using the requested backend.
// f must be safe to call concurrently for distinct i (no shared mutable
// state other than disjoint writes indexed by i).
template <typename F>
void parallel_for(std::size_t n, ParBackend backend, F&& f, unsigned hw_threads = 0) {
    if (n == 0) return;

    switch (backend) {
        case ParBackend::Serial: {
            for (std::size_t i = 0; i < n; ++i) f(i);
            return;
        }

        case ParBackend::StlParUnseq: {
#if defined(HAVE_STL_PARALLEL)
            std::vector<std::size_t> idx(n);
            for (std::size_t i = 0; i < n; ++i) idx[i] = i;
            std::for_each(std::execution::par_unseq, idx.begin(), idx.end(),
                          [&](std::size_t i) { f(i); });
            return;
#else
            // Fall through to JThread if not compiled with parallel STL support.
            [[fallthrough]];
#endif
        }

        case ParBackend::JThread: {
            const unsigned nthreads = hw_threads ? hw_threads
                                     : std::max(1u, std::thread::hardware_concurrency());
            const std::size_t chunk = (n + nthreads - 1) / nthreads;

            {
                std::vector<std::jthread> workers;
                workers.reserve(nthreads);
                for (unsigned t = 0; t < nthreads; ++t) {
                    const std::size_t begin = t * chunk;
                    const std::size_t end = std::min(n, begin + chunk);
                    if (begin >= end) break;
                    workers.emplace_back([begin, end, &f] {
                        for (std::size_t i = begin; i < end; ++i) f(i);
                    });
                }
                // std::jthread destructors join automatically here.
            }
            return;
        }
    }
}

}  // namespace laser_sim
