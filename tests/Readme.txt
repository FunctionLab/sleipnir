# To run the tests
# Build the code and in the binaries from the build run test/unit_tests
cd Debug/
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
./tests/unit_tests

# To add a new unit test 
# Add the new test file to tests/ directory
# Edit tests/CMakeLists.txt and add the test file to add_executable section

