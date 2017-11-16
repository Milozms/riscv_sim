cmake_minimum_required(VERSION 3.8)
project(riscv_sim)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES Simulation-pipe.cpp Simulation-pipe.h Reg_def.h Read_Elf.cpp Read_Elf.h)
add_executable(sim ${SOURCE_FILES})