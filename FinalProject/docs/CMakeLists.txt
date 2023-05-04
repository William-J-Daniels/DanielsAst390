FIND_PACKAGE(LATEX)

IF(LATEX_FOUND AND LATEX_PDFLATEX_FOUND)
    add_custom_target(FpPreBuild
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    add_custom_target(FpLatex
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/writeup.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_custom_target(FpWriteup ALL)
    add_dependencies(FpWriteup FpLatex)
ELSE()
    MESSAGE(ERROR "Latex might not be installed, not building docs.")
ENDIF()
