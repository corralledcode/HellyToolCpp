cmake_minimum_required(VERSION 3.28)
project(HellyToolCpp)

set(CMAKE_CXX_STANDARD 23)

add_executable(HellyToolCpp main.cpp
        Batchprocesseddata.cpp
        Graph.cpp
        Formatdata.cpp
        Formatgraph.cpp
        EdgesforHelly.cpp
)

