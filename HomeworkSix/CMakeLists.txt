cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(HomeworkSix VERSION 1.0.0)

if(build_docs)
  add_subdirectory(docs)
endif()
