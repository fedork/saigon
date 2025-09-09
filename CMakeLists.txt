cmake_minimum_required(VERSION 3.31)
project(saigon C)

set(CMAKE_C_STANDARD 17)

add_executable(saigon main.c)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)
target_link_libraries(saigon PRIVATE PkgConfig::GMP)
