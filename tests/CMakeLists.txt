add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

add_executable(
    unit_tests
    example_add.cpp
    example_subtract.cpp
    )

  #add_library(example STATIC
  #  example_add.cpp
  #  example_subtract.cpp
  #)

target_link_libraries(
    unit_tests
    gtest_main
    sleipnir
    )

add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/unit_tests
  )
