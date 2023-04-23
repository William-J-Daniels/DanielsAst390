cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(Ast390 VERSION 1.0.0)

find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(OpenMP REQUIRED)

option(build_docs "Build homework writeups" ON)
option(debugging "Complier settings for debugging" OFF)

if(debugging)
  add_compile_options(-Wall)
  add_compile_options(-fsanitize=address)
  link_libraries(-fsanitize=address)
endif()

find_package(Python COMPONENTS Interpreter REQUIRED)
# add checks for reqiured packages later

add_subdirectory(HomeworkOne)
add_subdirectory(HomeworkTwo)
add_subdirectory(HomeworkThree)
add_subdirectory(HomeworkFour)
add_subdirectory(HomeworkFive)
add_subdirectory(HomeworkSix)
