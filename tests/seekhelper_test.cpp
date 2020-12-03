#include <iostream>
#include <fstream>
#include <filesystem>
#include "gtest/gtest.h"
#include "seekhelper.h"
#include "seekdataset.h"

using namespace std;
using namespace Sleipnir;

// Forward declarations
bool compareSeekDBSettings(CSeekDBSetting* cc1, CSeekDBSetting* cc2);

string tomlConfigData =
R"""(port = 1234
numThreads = 12
numBufferedDBs = 23  # Number of Databaselets to store in memory
scoreCutoff = -3.1  # The gene-gene score cutoff to add, default: no cutoff
squareZ = true  # If using z-score, square-transform z-scores. Usually used in conjunction with --score-cutoff
isNibble = true  # If CDatabase is stored in nibble or byte format
outputAsText = true  # Output results (gene list and dataset weights) as text

[Database]
    [Database.DB1]
    DB_DIR        = "/path/db1/db"
    PREP_DIR      = "/path/db1/prep"
    PLATFORM_DIR  = "/path/db1/plat"
    SINFO_DIR     = "/path/db1/sinfo"
    GVAR_DIR      = "/path/db1/gvar"
    QUANT_FILE    = "/path/db1/quant2"
    DSET_MAP_FILE = "/path/db1/dataset.map"
    GENE_MAP_FILE = "/path/db1/gene_map.txt"
    DSET_SIZE_FILE = "/path/db1/dataset_size"
    NUMBER_OF_DB  = 1000

    [Database.DB2]
    PREP_DIR      = "/path/db2/prep"
    PLATFORM_DIR  = "/path/db2/plat"
    DB_DIR        = "/path/db2/db"
    DSET_MAP_FILE = "/path/db2/dataset.map"
    GENE_MAP_FILE = "/path/db2/gene_map.txt"
    QUANT_FILE    = "/path/db2/quant2"
    SINFO_DIR     = "/path/db2/sinfo"
    GVAR_DIR      = "/path/db2/gvar"
    DSET_SIZE_FILE = "/path/db2/dataset_size"
    NUMBER_OF_DB  = 2000
)""";

// Note: legacy file must use tab separations
string legacyDBConfigData = 
R"""(START
DB_DIR	/path/db1/db
PREP_DIR	/path/db1/prep
PLATFORM_DIR	/path/db1/plat
SINFO_DIR	/path/db1/sinfo
GVAR_DIR	/path/db1/gvar
QUANT_FILE	/path/db1/quant2
DSET_MAP_FILE	/path/db1/dataset.map
GENE_MAP_FILE	/path/db1/gene_map.txt
DSET_SIZE_FILE	/path/db1/dataset_size
NUMBER_OF_DB	1000
END
START
DB_DIR	/path/db2/db
PREP_DIR	/path/db2/prep
PLATFORM_DIR	/path/db2/plat
SINFO_DIR	/path/db2/sinfo
GVAR_DIR	/path/db2/gvar
QUANT_FILE	/path/db2/quant2
DSET_MAP_FILE	/path/db2/dataset.map
GENE_MAP_FILE	/path/db2/gene_map.txt
DSET_SIZE_FILE	/path/db2/dataset_size
NUMBER_OF_DB	2000
END
)""";

// Setup a fixture class of data for the tests
class SeekHelperTest : public ::testing::Test
{
public:
    string testDir = "/tmp/test/";
    string testTomlFile = testDir + "seekConfig.toml";
    string testLegacyConfigFile = testDir + "seekAdditionalDB";

protected:
    void SetUp() override
    {
        if (!filesystem::exists(testDir))
        {
            filesystem::create_directories(testDir);
        }
        ofstream tmpFile;
        tmpFile.open(testTomlFile);
        tmpFile << tomlConfigData;
        tmpFile.close();

        tmpFile.open(testLegacyConfigFile);
        tmpFile << legacyDBConfigData;
        tmpFile.close();
    }
};

// TEST saveLoadData - save out the seekPlatform to files and reload it
TEST_F(SeekHelperTest, loadTomlConfigs)
{
    SeekSettings tomlSettings;
    bool res;

    // load the toml file
    res = parseTomlConfig(testTomlFile,  tomlSettings);
    ASSERT_TRUE(res);
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
    ASSERT_TRUE(compareSeekDBSettings(tomlSettings.dbs[0], legacyDBSettings[0]));
}

bool compareSeekDBSettings(CSeekDBSetting* cc1, CSeekDBSetting* cc2) {
    // Note could also add an == operator to CSeekDBSetting class
    // bool operator==(CSeekDBSetting &rhs) {
    //     if (m_geneMapFile == rhs.GetValue("gene") && ...) {
    //         return true;
    //     }
    //     return false;
    // }
    if (cc1->GetValue("gene") == cc2->GetValue("gene") &&
        cc1->GetValue("dset") == cc2->GetValue("dset") &&
        cc1->GetValue("quant") == cc2->GetValue("quant") &&
        cc1->GetValue("gvar") == cc2->GetValue("gvar") &&
        cc1->GetValue("sinfo") == cc2->GetValue("sinfo") &&
        cc1->GetValue("db") == cc2->GetValue("db") &&
        cc1->GetValue("prep") == cc2->GetValue("prep") &&
        cc1->GetValue("platform") == cc2->GetValue("platform") &&
        cc1->GetValue("dset_size") == cc2->GetValue("dset_size") &&
        cc1->GetNumDB() == cc2->GetNumDB() ) {
           return true;
       }
       return false;
}
