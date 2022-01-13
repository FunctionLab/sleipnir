#include <iostream>
#include <fstream>
#include <filesystem>

using namespace std;

string humanTomlFilename = "humanConfig.toml";
string yeastTomlFilename = "yeastConfig.toml";
string testLegacyConfigFilename = "seekAdditionalDB.txt";

string tomlConfigHuman =
R"""(species = "human"
port = 1234
numThreads = 12
numBufferedDBs = 23  # Max num genes per query and max num of Databaselets to store in memory
pclCacheSize = 32  # Max CPCLs to keep open in cached memory (for PCL Server queries)
scoreCutoff = -3.1  # The gene-gene score cutoff to add, default: no cutoff
squareZ = true  # If using z-score, square-transform z-scores. Usually used in conjunction with --score-cutoff
isNibble = true  # If CDatabase is stored in nibble or byte format
outputAsText = true  # Output results (gene list and dataset weights) as text

[Database]
    [Database.DB1]
    DB_DIR        = "/path/human/db1/db"
    PREP_DIR      = "/path/human/db1/prep"
    PLATFORM_DIR  = "/path/human/db1/plat"
    SINFO_DIR     = "/path/human/db1/sinfo"
    GVAR_DIR      = "/path/human/db1/gvar"
    PCL_DIR      = "/path/human/db1/pcl"
    PVAL_DIR      = "/path/human/db1/random"
    QUANT_FILE    = "/path/human/db1/quant2"
    DSET_MAP_FILE = "/path/human/db1/dataset.map"
    GENE_MAP_FILE = "/path/human/db1/gene_map.txt"
    GENE_SYMBOL_FILE = "/path/human/db1/gene_symbols.txt"
    DSET_SIZE_FILE = "/path/human/db1/dataset_size"
    NUMBER_OF_DB  = 1000

    [Database.DB2]
    PREP_DIR      = "/path/human/db2/prep"
    PLATFORM_DIR  = "/path/human/db2/plat"
    DB_DIR        = "/path/human/db2/db"
    PCL_DIR        = "/path/human/db2/pcl"
    PVAL_DIR        = "/path/human/db2/random"
    DSET_MAP_FILE = "/path/human/db2/dataset.map"
    GENE_MAP_FILE = "/path/human/db2/gene_map.txt"
    GENE_SYMBOL_FILE = "/path/human/db2/gene_symbols.txt"
    QUANT_FILE    = "/path/human/db2/quant2"
    SINFO_DIR     = "/path/human/db2/sinfo"
    GVAR_DIR      = "/path/human/db2/gvar"
    DSET_SIZE_FILE = "/path/human/db2/dataset_size"
    NUMBER_OF_DB  = 2000
)""";

string tomlConfigYeast =
R"""(species = "yeast"
port = 5678
numThreads = 34
numBufferedDBs = 45  #  Max num genes per query and max num of Databaselets to store in memory
pclCacheSize = 16  # Max CPCLs to keep open in cached memory (for PCL Server queries)
scoreCutoff = 8.9  # The gene-gene score cutoff to add, default: no cutoff
squareZ = true  # If using z-score, square-transform z-scores. Usually used in conjunction with --score-cutoff
isNibble = true  # If CDatabase is stored in nibble or byte format
outputAsText = true  # Output results (gene list and dataset weights) as text

[Database]
    [Database.DB1]
    DB_DIR        = "/path2/yeast/db"
    PREP_DIR      = "/path2/yeast/prep"
    PLATFORM_DIR  = "/path2/yeast/plat"
    SINFO_DIR     = "/path2/yeast/sinfo"
    GVAR_DIR      = "/path2/yeast/gvar"
    PCL_DIR      = "/path2/yeast/pcl"
    PVAL_DIR      = "/path2/yeast/random"
    QUANT_FILE    = "/path2/yeast/quant2"
    DSET_MAP_FILE = "/path2/yeast/dataset.map"
    GENE_MAP_FILE = "/path2/yeast/gene_map.txt"
    GENE_SYMBOL_FILE = "/path2/yeast/gene_symbols.txt"
    DSET_SIZE_FILE = "/path2/yeast/dataset_size"
    NUMBER_OF_DB  = 1000

)""";

// Note: legacy file must use tab separations
string legacyDBConfigData = 
R"""(START
DB_DIR	/path/human/db1/db
PREP_DIR	/path/human/db1/prep
PLATFORM_DIR	/path/human/db1/plat
SINFO_DIR	/path/human/db1/sinfo
GVAR_DIR	/path/human/db1/gvar
PCL_DIR	/path/human/db1/pcl
PVAL_DIR	/path/human/db1/random
QUANT_FILE	/path/human/db1/quant2
DSET_MAP_FILE	/path/human/db1/dataset.map
GENE_MAP_FILE	/path/human/db1/gene_map.txt
GENE_SYMBOL_FILE	/path/human/db1/gene_symbols.txt
DSET_SIZE_FILE	/path/human/db1/dataset_size
NUMBER_OF_DB	1000
END
START
DB_DIR	/path/human/db2/db
PREP_DIR	/path/human/db2/prep
PLATFORM_DIR	/path/human/db2/plat
SINFO_DIR	/path/human/db2/sinfo
PCL_DIR	/path/human/db2/pcl
PVAL_DIR	/path/human/db2/random
GVAR_DIR	/path/human/db2/gvar
QUANT_FILE	/path/human/db2/quant2
DSET_MAP_FILE	/path/human/db2/dataset.map
GENE_MAP_FILE	/path/human/db2/gene_map.txt
GENE_SYBMOL_FILE	/path/human/db2/gene_symbols.txt
DSET_SIZE_FILE	/path/human/db2/dataset_size
NUMBER_OF_DB	2000
END
)""";

bool createConfigFiles(string testDir) {
    string humanTomlPath = testDir + humanTomlFilename;
    string yeastTomlPath = testDir + yeastTomlFilename;
    string testLegacyConfigPath = testDir + testLegacyConfigFilename;

    if (!filesystem::exists(testDir))
    {
        filesystem::create_directories(testDir);
    }
    ofstream tmpFile;
    tmpFile.open(humanTomlPath);
    tmpFile << tomlConfigHuman;
    tmpFile.close();

    tmpFile.open(yeastTomlPath);
    tmpFile << tomlConfigYeast;
    tmpFile.close();

    tmpFile.open(testLegacyConfigPath);
    // if (tmpFile.fail()) {
    //     cerr << "Error: file open " << testLegacyConfigFilename << ": " << strerror(errno) << endl;
    //     return false;
    // }
    tmpFile << legacyDBConfigData;
    tmpFile.close();
    return true;
}