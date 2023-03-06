cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(HomeworkOne VERSION 1.0.0)

add_library(${PROJECT_NAME} SHARED
  include/differentiator.h
  src/differentiator.cpp
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set(public_headers
  include/differentiator.h
)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")

add_executable(Hw1P1 examples/question_one.cpp)
target_link_libraries(Hw1P1 PUBLIC ${PROJECT_NAME})

add_executable(Hw1P2 examples/question_two.cpp)

add_subdirectory(docs/)
