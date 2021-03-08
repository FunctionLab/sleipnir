#include <iostream>
#include <fstream>
#include <filesystem>
#include <ctime>
#include <chrono>
#include "gtest/gtest.h"
#include "seekhelper.h"
#include <vector>

using namespace std;
using namespace Sleipnir;
namespace fs = std::filesystem;


// Setup a fixture class of data for the tests
class SeekHelperTest : public ::testing::Test
{
public:

protected:
};

TEST_F(SeekHelperTest, loadOneColumnTextFile)
{
    fs::path inputFilePath = fs::path(__BASE_FILE__).remove_filename();
    inputFilePath /= "test_inputs/oneColTextFile.txt";
    vector<string> vals;
    loadOneColumnTextFile(inputFilePath, vals);
    ASSERT_EQ(vals.size(), 7);
}

TEST_F(SeekHelperTest, loadTwoColumnTextFile)
{
    fs::path inputFilePath = fs::path(__BASE_FILE__).remove_filename();
    inputFilePath /= "test_inputs/twoColTextFile.txt";
    map<string, string> vals;
    loadTwoColumnTextFile(inputFilePath, vals);
    ASSERT_EQ(vals.size(), 4);
}
