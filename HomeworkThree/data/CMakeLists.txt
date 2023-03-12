# make the csv
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/P1.csv
  DEPENDS Hw3P1
  COMMAND Hw3P1
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tools/
)

add_custom_target(Hw3P1_csv ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/P1.csv
)

# make the plot
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/P1Plots.jpeg
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../tools/P1_plots.py
          ${CMAKE_CURRENT_SOURCE_DIR}/P1.csv
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tools/P1_plots.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tools/
)

add_custom_target(Hw3P1Plots ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/P1Plots.jpeg
)
