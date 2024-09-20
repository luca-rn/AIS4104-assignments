cmake_minimum_required(VERSION 3.29)
project(AIS4104-assignments)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(assignment-1)
add_subdirectory(assignment-2)

add_subdirectory(math)
add_subdirectory(session-1)

add_dependencies(assignment2 math)
add_dependencies(session-1 math)