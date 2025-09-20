# Saigon Project

This project contains a C program for analyzing electrical circuits. 
It uses the GMP library for arbitrary-precision arithmetic.
This is a solver for "Saigon Towers" https://oeis.org/A384183 and for 
https://oeis.org/A339548 and https://oeis.org/A339808 (resistor networks 
with resistances approaching 1)

## Prerequisites

Before building, you need the following tools installed:

*   A C compiler (e.g., GCC or Clang)
*   CMake (version 3.31 or higher)
*   pkg-config
*   GNU Multiple Precision Arithmetic Library (`gmp`)

### Prerequisite Installation

**On Debian/Ubuntu:**
```sh
sudo apt-get update
sudo apt-get install build-essential cmake pkg-config libgmp-dev
```

**On macOS (using Homebrew):**
```sh
brew install cmake pkg-config gmp
```

## Building

This project now builds two executables:
- saigon
- saigon_util (a utility copy of saigon.c)

1.  Configure the project using CMake.

    For a debug build, run this from the project root:
    ```sh
    cmake -S . -B build
    ```

    For an optimized release build:
    ```sh
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    ```

2.  Compile all targets:
    ```sh
    cmake --build build
    ```

    Or build a specific target, for example only saigon_util:
    ```sh
    cmake --build build --target saigon_util
    ```

## Running

The executables will be created in the build directory. Examples:

```sh
./build/saigon
./build/saigon_util
```

## Clean, build and run in one step:

```sh
rm -r build ; cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build  && ./build/saigon
```