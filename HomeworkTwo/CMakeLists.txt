cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(HomeworkTwo VERSION 1.0.0)

add_library(${PROJECT_NAME} SHARED
  include/simpsons.h
  src/simpsons.cpp
)

target_include_directories(${PROJECT_NAME}
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set(public_headers
  include/simpsons.h
)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/")

add_executable(ProblemOne examples/problem_one.cpp)
target_link_libraries(ProblemOne PUBLIC ${PROJECT_NAME})

add_executable(ProblemTwo examples/problem_two.cpp)
target_link_libraries(ProblemTwo PUBLIC ${PROJECT_NAME})
