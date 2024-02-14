cmake_minimum_required(VERSION 3.20)
project(GEO1004_ass_1)

set(CMAKE_CXX_STANDARD 17)

find_package( CGAL REQUIRED )

add_executable(GEO1004_ass_1 main.cpp)

target_link_libraries(GEO1004_ass_1 CGAL::CGAL)
