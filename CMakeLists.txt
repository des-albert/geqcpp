cmake_minimum_required(VERSION 3.19)
project(geq)

set(CMAKE_CXX_STANDARD 20)
include_directories("/usr/local/dislin")
link_libraries(discpp)
add_executable(geq geq.cpp geq.h bndmat.cpp eqsil.cpp flux.cpp splnco.cpp intpol.cpp startt.cpp compar.cpp xcur.cpp gelg.cpp saddle.cpp curf.cpp topol.cpp plotit.cpp)