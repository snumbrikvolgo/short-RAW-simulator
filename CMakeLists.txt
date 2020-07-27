cmake_minimum_required(VERSION 3.12)
project(IITP)

set(CMAKE_CXX_STANDARD 17)

add_executable(IITP main.cpp station.cpp station.h frame.h counter.cpp counter.h)