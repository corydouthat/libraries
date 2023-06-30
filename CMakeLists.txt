cmake_minimum_required(VERSION 3.13)

project(libraries C CXX ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

add_library(libraries INTERFACE)

target_include_directories(libraries INTERFACE .)
