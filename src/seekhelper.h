#ifndef SEEKHELPER_H
#define SEEKHELPER_H

#include <list>
#include <mutex>
#include <queue>
#include <shared_mutex>
#include <unordered_map>
#include <condition_variable>
#include "seekdataset.h"
#include "seekerror.h"

using namespace Sleipnir;

uint32_t omp_enabled_test();

// Settings that specify a Seek Database configurations. There
//  is typically a different database per species.
struct SeekSettings {
    vector <CSeekDBSetting*> dbs;
    string species;
    int32_t port = 9000;
    int32_t numThreads = 8;
    int32_t numBufferedDBs = 20;
    int32_t pclCacheSize = 100;
    double scoreCutoff = -9999;
    bool squareZ = false;
    bool isNibble = false;
    bool outputAsText = false;
    friend ostream& operator<<(ostream& os, const SeekSettings& settings) {
        os << "Species: " << settings.species << endl;
        os << "Port: " << settings.port << endl;
        os << "NumThreads: " << settings.numThreads << endl;
        os << "NumBufferedDBs: " << settings.numBufferedDBs << endl;
        os << "PclCacheSize: " << settings.pclCacheSize << endl;
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
    PCL_DIR      = "/path/db1/pcl"
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
    PCL_DIR     = "/path/db2/pcl"
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


template <typename T>
class ThreadSafeQueue {
public:
    void enqueue(T element) {
        lock_guard tlock(m_mutex);
        m_queue.push(element);
    }
    T dequeue() {
        lock_guard tlock(m_mutex);
        if (m_queue.empty()) {
            throw state_error("ThreadSafeQueue: Dequeue called on an empty queue");
        }
        T element = m_queue.front();
        m_queue.pop();
        return element;
    }
    uint32_t size() {
        lock_guard tlock(m_mutex);
        return m_queue.size();
    }
    bool empty() {
        lock_guard tlock(m_mutex);
        return (m_queue.size() == 0);
    }
private:
    queue<T> m_queue;
    mutex m_mutex;
};


template <typename K, typename V = K>
class LRUCache
{
public:
    LRUCache(uint32_t s) :csize(s) {};
    void set(const K key, const V value) {
        // Take unique lock. It is automatically released
        //  on scope exit
        unique_lock ulock(this->cacheMutex);
        auto pos = keyValuesMap.find(key);
        if (pos == keyValuesMap.end()) {
            items.push_front(key);
            keyValuesMap[key] = { value, items.begin() };
            if (keyValuesMap.size() > csize) {
                keyValuesMap.erase(items.back());
                items.pop_back();
            }
        } else {
            items.erase(pos->second.second);
            items.push_front(key);
            keyValuesMap[key] = { value, items.begin() };
        }
    }
    bool get(const K key, V &value) {
        // Take shared lock. It is automatically released
        //  on scope exit
        shared_lock slock(this->cacheMutex);
        auto pos = keyValuesMap.find(key);
        if (pos == keyValuesMap.end()) {
            return false;
        }
        items.erase(pos->second.second);
        items.push_front(key);
        keyValuesMap[key] = { pos->second.first, items.begin() };
        value = pos->second.first;
        return true;
    }
    uint32_t cacheSize() { return csize; }

private:
    list<K>items;
    unordered_map <K, pair<V, typename list<K>::iterator>> keyValuesMap;
    uint32_t csize;
    shared_mutex cacheMutex;
};

#endif  // SEEKHELPER_H
