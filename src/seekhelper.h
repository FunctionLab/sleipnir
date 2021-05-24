#ifndef SEEKHELPER_H
#define SEEKHELPER_H

#include <mutex>
#include <condition_variable>
#include "seekdataset.h"

using namespace Sleipnir;

uint32_t omp_enabled_test();

// Settings that specify a Seek Database configurations. There
//  is typically a different database per species.
struct SeekSettings {
    vector <CSeekDBSetting*> dbs;
    string species;
    int64_t port = 9000;
    int64_t numThreads = 8; 
    int64_t numBufferedDBs = 20;
    double scoreCutoff = -9999;
    bool squareZ = false;
    bool isNibble = false;
    bool outputAsText = false;
    friend ostream& operator<<(ostream& os, const SeekSettings& settings) {
        os << "Species: " << settings.species << endl;
        os << "Port: " << settings.port << endl;
        os << "NumThreads: " << settings.numThreads << endl;
        os << "NumBufferedDBs: " << settings.numBufferedDBs << endl;
        os << "SquareZ: " << settings.squareZ << endl;
        for (const CSeekDBSetting *db : settings.dbs) {
            os << *db;
        }
        return os;
    }
};

/*
Example Toml Config Format:

species = "human"
port = 9000
numThreads = 8
numBufferedDBs = 20  # Max num genes per query and max num of Databaselets to store in memory
scoreCutoff = -9999.0  # The gene-gene score cutoff to add, default: no cutoff
squareZ = false  # If using z-score, square-transform z-scores. Usually used in conjunction with --score-cutoff
isNibble = false  # If CDatabase is stored in nibble or byte format
outputAsText = false  # Output results (gene list and dataset weights) as text

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
    GENE_SYMBOL_FILE = "/path/db1/gene_symbols.txt"
    DSET_SIZE_FILE = "/path/db1/dataset_size"
    NUMBER_OF_DB  = 1000

    [Database.DB2]
    PREP_DIR      = "/path/db2/prep"
    PLATFORM_DIR  = "/path/db2/plat"
    DB_DIR        = "/path/db2/db"
    DSET_MAP_FILE = "/path/db2/dataset.map"
    GENE_MAP_FILE = "/path/db2/gene_map.txt"
    GENE_SYMBOL_FILE = "/path/db2/gene_symbols.txt"
    QUANT_FILE    = "/path/db2/quant2"
    SINFO_DIR     = "/path/db2/sinfo"
    GVAR_DIR      = "/path/db2/gvar"
    DSET_SIZE_FILE = "/path/db2/dataset_size"
    NUMBER_OF_DB  = 1000
*/
// Read database config files
bool parseTomlConfig(string tomlConfigFile, SeekSettings &settings);

// Read in a set of config files for different species
// map [speciesName -> Configs]
void getConfigs(vector<string> &configFiles, map<string, SeekSettings> &configs);
void getConfigs_old(vector<string> &configFiles, map<string, SeekSettings> &configs);

// Parse config file for invocation parameters
// New db instances are pushed onto CSeekDBSetting cc
// Make optional parameters distanceMeasure and check_dset_size_flag and set 
// them to the most strict case in terms of checking meta-data is available.
bool legacyReadDBConfigFile(string dbConfigFile,
                            vector<CSeekDBSetting*> &cc, 
                            CSeekDataset::DistanceMeasure eDistMeasure = CSeekDataset::CORRELATION,
                            bool check_dset_size_flag = true);

void loadOneColumnTextFile(string filename, vector<string> &vals);

void loadTwoColumnTextFile(string filename, map<string, string> &vals);


// Implementation of a basic semaphore class for resource management
class Semaphore
{
public:
    Semaphore(int count, int maxCount) : 
        _count(count), _maxCount(maxCount) {}
    void notify();
    void wait();
    bool try_wait();
    void lock() { this->wait(); }
    void unlock() { this->notify(); }
private:
    mutex _mutex;
    condition_variable _condition;
    unsigned long _count = 0;
    unsigned long _maxCount = 0;

};


// A class that implements BasicLockable and can be used with lock_guard primative
//  to initalize a boolean false and automatically set it when the context exits.
class BoolFlag {
public:
    BoolFlag(bool &flag) : _flag(flag) {}
    void lock() { this->_flag = false; }
    void unlock() { this->_flag = true; }
private:
    bool &_flag;
};

#endif  // SEEKHELPER_H
