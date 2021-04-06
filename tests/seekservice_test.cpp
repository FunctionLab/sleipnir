#include <iostream>
#include <fstream>
#include <filesystem>
#include "gtest/gtest.h"
#include "seekhelper.h"
#include "seekdataset.h"
#include "sampleConfigFiles.h"

using namespace std;
using namespace Sleipnir;
using namespace SeekRPC;

// Setup a fixture class of data for the tests
class SeekRPCTest : public ::testing::Test
{
public:
    string testDir = "/tmp/test/";
    string humanTomlFile = testDir + humanTomlFilename;
    string yeastTomlFile = testDir + yeastTomlFilename;

protected:
    void SetUp() override
    {
        createConfigFiles(testDir);
    }
};

// TEST - load toml config files
TEST_F(SeekHelperTest, loadTomlConfigs)
{
}