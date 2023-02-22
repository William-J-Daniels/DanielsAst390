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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/")

add_executable(QuestionOne examples/question_one.cpp)
target_link_libraries(QuestionOne PUBLIC ${PROJECT_NAME})

add_executable(QuestionTwo examples/question_two.cpp)

FIND_PACKAGE(LATEX)
IF(LATEX_FOUND AND LATEX_PDFLATEX_FOUND)
    add_custom_target(latex-PRE_BUILD
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR}/docs/ ${CMAKE_CURRENT_SOURCE_DIR}/docs/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs/
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/docs/writeup.tex
    )
    add_custom_target(latex-pdf
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR}/docs/ ${CMAKE_CURRENT_SOURCE_DIR}/docs/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs/
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/docs/writeup.tex
    )

    add_custom_target(all-formats ALL)
    add_dependencies(all-formats latex-pdf)
ELSE()
    MESSAGE(ERROR "Latex might not be installed, not building docs.")
ENDIF()
