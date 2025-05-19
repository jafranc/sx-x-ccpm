cmake_minimum_required (VERSION 2.8.12)
project(test)

set(exe "test")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build)

find_package(Boost REQUIRED)

add_executable(${exe} main.cpp)
target_link_libraries(${exe} ${Boost_GRAPH_LIBRARY})
