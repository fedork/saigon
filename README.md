# Saigon Project

This project contains a C program for analyzing electrical circuits, likely related to finding equivalent resistance in complex networks. It uses the GMP library for arbitrary-precision arithmetic.

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

1.  **Configure the project using CMake.**

    For a debug build, run this from the project root:
    ```sh
    cmake -S . -B build
    ```

    For an optimized release build:
    ```sh
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    ```

2.  **Compile the code:**
    ```sh
    cmake --build build
    ```

## Running

The executable `saigon` will be created in the `build` directory. You can run it from the project root like this:

```sh
./build/saigon
```
