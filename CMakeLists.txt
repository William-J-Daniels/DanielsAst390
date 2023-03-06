cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(Ast390 VERSION 1.0.0)

option(build_docs "Build homework writeups" ON)

add_subdirectory(HomeworkOne)
add_subdirectory(HomeworkTwo)
add_subdirectory(HomeworkThree)
