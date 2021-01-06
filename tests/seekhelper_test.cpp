#include <iostream>
#include <fstream>
#include <filesystem>
#include <ctime>
#include <chrono>
#include "gtest/gtest.h"
#include "seekhelper.h"
#include "seekdataset.h"
#include "sampleConfigFiles.h"

using namespace std;
using namespace Sleipnir;

// Forward declarations
bool compareSeekDBSettings(CSeekDBSetting* cc1, CSeekDBSetting* cc2);
bool compareSeekSettings(SeekSettings &s1, SeekSettings &s2);
void printingAndTimingGetConfig(vector<string> &configFiles);


// Setup a fixture class of data for the tests
class SeekHelperTest : public ::testing::Test
{
public:
    string testDir = "/tmp/test/";
    string humanTomlFile = testDir + humanTomlFilename;
    string yeastTomlFile = testDir + yeastTomlFilename;
    string testLegacyConfigFile = testDir + testLegacyConfigFilename;

protected:
    void SetUp() override
    {
        createConfigFiles(testDir);
    }
};

// TEST saveLoadData - save out the seekPlatform to files and reload it
TEST_F(SeekHelperTest, loadTomlConfigs)
{
    SeekSettings tomlSettings;
    bool res;

    // load the toml file
    res = parseTomlConfig(humanTomlFile,  tomlSettings);
    ASSERT_TRUE(res);
    ASSERT_EQ(tomlSettings.species, "human");
    ASSERT_EQ(tomlSettings.port, 1234);
    ASSERT_EQ(tomlSettings.numThreads, 12);
    ASSERT_EQ(tomlSettings.numBufferedDBs, 23);
    ASSERT_EQ(tomlSettings.scoreCutoff, -3.1);
    ASSERT_TRUE(tomlSettings.squareZ);
    ASSERT_TRUE(tomlSettings.isNibble);
    ASSERT_TRUE(tomlSettings.scoreCutoff);

    // load the from the legacy config format
    vector <CSeekDBSetting*> legacyDBSettings;
    res = legacyReadDBConfigFile(testLegacyConfigFile, legacyDBSettings);
    ASSERT_TRUE(res);
    ASSERT_EQ(tomlSettings.dbs.size(), 2);
    ASSERT_EQ(legacyDBSettings.size(), 2);
    ASSERT_TRUE(compareSeekDBSettings(tomlSettings.dbs[0], legacyDBSettings[0]));
    ASSERT_TRUE(compareSeekDBSettings(tomlSettings.dbs[1], legacyDBSettings[1]));
}

TEST_F(SeekHelperTest, loadSpeciesConfigs)
{
    vector<string> configFiles;
    map<string, SeekSettings> speciesConfigs;

    configFiles.push_back(humanTomlFile);
    configFiles.push_back(yeastTomlFile);
    getConfigs(configFiles, speciesConfigs);

    // load the from the legacy config format
    vector <CSeekDBSetting*> legacyDBSettings;
    bool res = legacyReadDBConfigFile(testLegacyConfigFile, legacyDBSettings);
    ASSERT_TRUE(res);
    ASSERT_EQ(speciesConfigs["human"].dbs.size(), 2);
    ASSERT_EQ(legacyDBSettings.size(), 2);
    ASSERT_TRUE(compareSeekDBSettings(speciesConfigs["human"].dbs[0], legacyDBSettings[0]));
    ASSERT_TRUE(compareSeekDBSettings(speciesConfigs["human"].dbs[1], legacyDBSettings[1]));

    // load yeast separately and compare
    SeekSettings yeastTomlSettings;
    res = parseTomlConfig(yeastTomlFile,  yeastTomlSettings);
    ASSERT_TRUE(compareSeekSettings(yeastTomlSettings, speciesConfigs["yeast"]));

    printingAndTimingGetConfig(configFiles);
}

void printingAndTimingGetConfig(vector<string> &configFiles) {
    using namespace chrono;
    int iters = 10;
    map<string, SeekSettings> speciesConfigs;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (int i=0; i<iters; i++) {
        getConfigs(configFiles, speciesConfigs);
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "Time load config ("<< iters << " iters): " << time_span.count() << " seconds" << endl;

    for( auto const& [key, val] : speciesConfigs ) {
        cout << "CONFIG: " << key << endl;
        cout << val << endl;
    }
}

bool compareSeekDBSettings(CSeekDBSetting* cc1, CSeekDBSetting* cc2) {
    // Note could also add an == operator to CSeekDBSetting class
    // bool operator==(CSeekDBSetting &rhs) {
    //     if (m_geneMapFile == rhs.geneMapFile && ...) {
    //         return true;
    //     }
    //     return false;
    // }
    if (cc1->geneMapFile == cc2->geneMapFile &&
        // cc1->geneSymbolFile == cc2->geneSymbolFile &&
        cc1->datasetFile == cc2->datasetFile &&
        cc1->quantFile == cc2->quantFile &&
        cc1->gvarDir == cc2->gvarDir &&
        cc1->sinfoDir == cc2->sinfoDir &&
        cc1->dbDir == cc2->dbDir &&
        cc1->prepDir == cc2->prepDir &&
        cc1->platDir == cc2->platDir &&
        cc1->dsetSizeFile == cc2->dsetSizeFile &&
        cc1->GetNumDB() == cc2->GetNumDB() ) {
           return true;
       }
       return false;
}

bool compareSeekSettings(SeekSettings &s1, SeekSettings &s2) {
    if (s1.species == s2.species &&
        s1.port == s2.port &&
        s1.numThreads == s2.numThreads &&
        s1.numBufferedDBs == s2.numBufferedDBs &&
        s1.squareZ == s2.squareZ &&
        s1.scoreCutoff == s2.scoreCutoff &&
        s1.outputAsText == s2.outputAsText) {
            if (s1.dbs.size() != s2.dbs.size()) {
                return false;
            }
            for (int i=0; i<s1.dbs.size(); i++) {
                if (compareSeekDBSettings(s1.dbs[i], s2.dbs[i]) == false) {
                    return false;
                } 
            }
            return true;
        }
    return false;
}