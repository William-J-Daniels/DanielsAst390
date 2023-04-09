FIND_PACKAGE(LATEX)

IF(LATEX_FOUND AND LATEX_PDFLATEX_FOUND)
    add_custom_target(Hw5PreBuild
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    add_custom_target(Hw5Latex
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_custom_target(Hw5Writeup ALL)
    add_dependencies(Hw5Writeup Hw5Latex)
ELSE()
    MESSAGE(ERROR "Latex might not be installed, not building docs.")
ENDIF()
