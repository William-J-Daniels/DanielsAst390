cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(HomeworkFive VERSION 1.0.0)

add_library(${PROJECT_NAME} SHARED
  include/twodiminterp.h
  src/twodiminterp.cpp
  include/logspace.h
)

target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(${PROJECT_NAME} PUBLIC TBB::tbb)

target_include_directories(${PROJECT_NAME}
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set(public_headers
  include/twodiminterp.h
  include/logspace.h
)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

add_executable(Hw5P1 examples/problem.cpp)
target_link_libraries(Hw5P1 PUBLIC ${PROJECT_NAME})
