# To run the tests
# Build the code and in the binaries from the build run test/unit_tests
cd Debug/
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
./tests/unit_tests

# To add a new unit test 
# Add the new test file to tests/ directory
# Edit tests/CMakeLists.txt and add the test file to add_executable section

# To run seekRPC and seekMiner system Tests (pytest based):
1) Install the conda environment (first time only)
conda env create --file conda_environment.yml

2) Activate the conda environment
conda activate genomics

3) Run the system tests
# From the sleipnir/ directory
python -m pytest -s -v tests