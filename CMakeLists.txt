cmake_minimum_required(VERSION 3.31)
project(saigon C)

set(CMAKE_C_STANDARD 17)

add_executable(saigon saigon.c)
add_executable(saigon_util saigon_util.c)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)
target_link_libraries(saigon PRIVATE PkgConfig::GMP)
target_link_libraries(saigon_util PRIVATE PkgConfig::GMP)

# Add compiler flags for Release configuration (e.g. for production builds)
# Use: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -DNDEBUG" CACHE STRING "Flags used by the C compiler for Release builds.")

# Enable Link-Time Optimization (LTO) for Release builds if supported
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON CACHE BOOL "Enable LTO for Release builds")
endif()
