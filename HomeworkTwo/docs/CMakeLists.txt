FIND_PACKAGE(LATEX)

IF(LATEX_FOUND AND LATEX_PDFLATEX_FOUND)
    add_custom_target(Hw2PreBuild
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
    )
    add_custom_target(Hw2Latex
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
    )

    add_custom_target(Hw2Writeup ALL)
    add_dependencies(Hw2Writeup Hw2Latex)
ELSE()
    MESSAGE(ERROR "Latex might not be installed, not building docs.")
ENDIF()
