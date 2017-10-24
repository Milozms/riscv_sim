cmake_minimum_required(VERSION 3.8)
project(riscv_sim)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES Simulation.cpp Simulation.h Reg_def.h)
add_executable(riscv_sim ${SOURCE_FILES})