cmake_minimum_required(VERSION 3.16)
project(laser_hardening_sim CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Threads REQUIRED)

# Optional: enable std::execution::par_unseq via Intel TBB, if present.
# On (most) libstdc++ builds, parallel <algorithm> policies require TBB
# to actually run in parallel; without it they silently fall back to
# serial execution. Give an explicit HAVE_STL_PARALLEL define so the
# code only offers the "stl" backend when it will really parallelize.
find_package(TBB QUIET)

add_executable(laser_hardening_sim
  src/main.cpp
)
target_include_directories(laser_hardening_sim PRIVATE include)
target_link_libraries(laser_hardening_sim PRIVATE Threads::Threads)

if(TBB_FOUND)
  message(STATUS "TBB found: enabling std::execution::par_unseq backend (--backend stl)")
  target_compile_definitions(laser_hardening_sim PRIVATE HAVE_STL_PARALLEL)
  target_link_libraries(laser_hardening_sim PRIVATE TBB::tbb)
else()
  message(STATUS "TBB not found: '--backend stl' will fall back to the jthread backend")
endif()

# Small standalone sanity-check binary (not a full test framework -- just
# assert-based checks on the kinetics laws and energy behavior).
add_executable(sanity_checks
  tests/sanity_checks.cpp
)
target_include_directories(sanity_checks PRIVATE include)
target_link_libraries(sanity_checks PRIVATE Threads::Threads)
if(TBB_FOUND)
  target_compile_definitions(sanity_checks PRIVATE HAVE_STL_PARALLEL)
  target_link_libraries(sanity_checks PRIVATE TBB::tbb)
endif()

enable_testing()
add_test(NAME sanity_checks COMMAND sanity_checks)
