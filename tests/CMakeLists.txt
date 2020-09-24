add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

add_executable(
    unit_tests
    # example_add.cpp
    # example_subtract.cpp
    seekplatform_test.cpp
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

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_link_libraries( unit_tests stdc++fs)
endif()

add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/unit_tests
  )
