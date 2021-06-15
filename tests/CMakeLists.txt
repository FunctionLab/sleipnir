add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

add_executable(
    unit_tests
    # example_add.cpp
    # example_subtract.cpp
    seekplatform_test.cpp
    seekhelper_test.cpp
    seekerror_test.cpp
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    find_library(LIBOMP libomp.a)
    list(APPEND EXTRA_LIBRARIES ${LIBOMP})
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    list(APPEND EXTRA_LIBRARIES stdc++fs gomp)
endif()


target_link_libraries(
    unit_tests
    gtest_main
    sleipnir
    ${EXTRA_LIBRARIES}
)

add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/unit_tests
)
