cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(HomeworkThree VERSION 1.0.0)

add_library(${PROJECT_NAME} SHARED
  include/ode_integrator.h
  src/ode_integrator.cpp
  include/eulermethod.h
  src/eulermethod.cpp
)

target_include_directories(${PROJECT_NAME}
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set(public_headers
  include/ode_integrator.h
  include/eulermethod.h
)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")

add_executable(Hw3P1 examples/problem_one.cpp)
target_link_libraries(Hw3P1 PUBLIC ${PROJECT_NAME})
#
# add_executable(ProblemTwo examples/problem_two.cpp)
# target_link_libraries(ProblemTwo PUBLIC ${PROJECT_NAME})

if(build_docs)
  add_subdirectory(docs/)
endif()
