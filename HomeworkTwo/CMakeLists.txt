cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(HomeworkTwo VERSION 1.0.0)

add_library(${PROJECT_NAME} SHARED
  include/simpsons.h
  src/simpsons.cpp
  include/indefinite_integrals.h
  src/indefinite_integrals.cpp
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
  include/indefinite_integrals.h
)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")

add_executable(Hw2P1 examples/problem_one.cpp)
target_link_libraries(Hw2P1 PUBLIC ${PROJECT_NAME})

add_executable(Hw2P2 examples/problem_two.cpp)
target_link_libraries(Hw2P2 PUBLIC ${PROJECT_NAME})

add_subdirectory(docs/)
