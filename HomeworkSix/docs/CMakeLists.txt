FIND_PACKAGE(LATEX)

IF(LATEX_FOUND AND LATEX_PDFLATEX_FOUND)
    add_custom_target(Hw6PreBuild
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    add_custom_target(Hw6Latex
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_custom_target(Hw6Writeup ALL)
    add_dependencies(Hw6Writeup Hw6Latex)
ELSE()
    MESSAGE(ERROR "Latex might not be installed, not building docs.")
ENDIF()
