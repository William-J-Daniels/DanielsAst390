# make the part 1 csv
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/TenDeg.csv
         ${CMAKE_CURRENT_SOURCE_DIR}/HundredDeg.csv
  DEPENDS Hw3P1
  COMMAND Hw3P1
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tools/
)

add_custom_target(Hw3P1_csv ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/TenDeg.csv
          ${CMAKE_CURRENT_SOURCE_DIR}/HundredDeg.csv
)

# make the part 1 plot
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/P1Plots.jpeg
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../tools/P1_plots.py
          ${CMAKE_CURRENT_SOURCE_DIR}/TenDeg.csv
          ${CMAKE_CURRENT_SOURCE_DIR}/HundredDeg.csv
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tools/P1_plots.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tools/
)

add_custom_target(Hw3P1Plots ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/P1Plots.jpeg
)
