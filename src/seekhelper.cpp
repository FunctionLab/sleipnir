#include <map>
#include <regex>
#include <vector>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "toml.hpp"
#include "seekhelper.h"

using namespace std;
using namespace Sleipnir;


// Helper function used by parseTomlConfig to get a named setting
template <typename T>
void tomlGetValue(toml::table tbl, string key, T &retVal) {
    if (tbl[key]) {
        optional<T> opval = tbl[key].value<T>();
        if (opval.has_value()) {
            retVal = opval.value();
            return;
        }
        cerr << "Toml Config: unable to parse key: " << key << endl;
    }
    cerr << "Toml: No entry found for \'" << key << "\'" << ", default to: " << retVal << endl;
}


/*
 * ParseTomlConfig:
 * Seek can input it's invocation parameters from a TOML configuration
 *   file as opposed to using command line switches.
 * Parse the TOML config file and return the settings in SeekSettings
 */
bool parseTomlConfig(string tomlConfigFile, SeekSettings &settings) {
    toml::table tbl;
    try
    {
        tbl = toml::parse_file(tomlConfigFile);
    }
    catch (const toml::parse_error& err)
    {
        cerr << "readTomlConf parsing failed:\n" << err << "\n";
        // return false;
        throw err;
    }
    // populate the top level settings
    tomlGetValue<string>(tbl, "species", settings.species);
    tomlGetValue<int32_t>(tbl, "port", settings.port);
    tomlGetValue<int32_t>(tbl, "numThreads", settings.numThreads);
    tomlGetValue<int32_t>(tbl, "numBufferedDBs", settings.numBufferedDBs);
    tomlGetValue<int32_t>(tbl, "pclCacheSize", settings.pclCacheSize);
    tomlGetValue<double>(tbl, "scoreCutoff", settings.scoreCutoff);
    tomlGetValue<bool>(tbl, "squareZ", settings.squareZ);
    tomlGetValue<bool>(tbl, "isNibble", settings.isNibble);
    tomlGetValue<bool>(tbl, "outputAsText", settings.outputAsText);

    // Database settings will be under the 'database' section
    // get the database settings
    if (tbl["Database"].is_table()) {
        auto databaseTbl = tbl["Database"];
        string key;
        int idx = 0;
        // Multiple databases can be described, each in a 'DB' subsection
        while(idx++, key="DB"+to_string(idx), databaseTbl[key].is_table()) {
            // auto dbTbl = databaseTbl[key].value<toml::v2::table>().value();
            toml::table *dbTbl = databaseTbl[key].as_table();
            string db_dir = "NA";
            string prep_dir = "NA";
            string platform_dir = "NA";
            string sinfo_dir = "NA";
            string gvar_dir = "NA";
            string pcl_dir = "NA";
            string quant_file = "NA";
            string gene_map_file = "NA";
            string gene_symbol_file = "NA";
            string dset_map_file = "NA";
            string dset_size_file = "NA";
            int32_t num_db = -1;

            tomlGetValue<string>(*dbTbl, "DB_DIR", db_dir);
            tomlGetValue<string>(*dbTbl, "PREP_DIR", prep_dir);
            tomlGetValue<string>(*dbTbl, "PLATFORM_DIR",platform_dir);
            tomlGetValue<string>(*dbTbl, "SINFO_DIR", sinfo_dir);
            tomlGetValue<string>(*dbTbl, "GVAR_DIR", gvar_dir);
            tomlGetValue<string>(*dbTbl, "PCL_DIR", pcl_dir);
            tomlGetValue<string>(*dbTbl, "QUANT_FILE", quant_file);
            tomlGetValue<string>(*dbTbl, "GENE_MAP_FILE", gene_map_file);
            tomlGetValue<string>(*dbTbl, "GENE_SYMBOL_FILE", gene_symbol_file);
            tomlGetValue<string>(*dbTbl, "DSET_MAP_FILE", dset_map_file);
            tomlGetValue<string>(*dbTbl, "DSET_SIZE_FILE", dset_size_file);
            tomlGetValue<int32_t>(*dbTbl, "NUMBER_OF_DB", num_db);

            CSeekDBSetting *dbSetting2 = 
                new CSeekDBSetting(gvar_dir, sinfo_dir, platform_dir, 
                                prep_dir, db_dir, gene_map_file, gene_symbol_file,
                                quant_file, dset_map_file,
                                dset_size_file, num_db);
            dbSetting2->setPclDir(pcl_dir);
            settings.dbs.push_back(dbSetting2);
        }
    } else {
        cerr << "Toml: expecting a section \'Database\' to have subsections for each DB" << endl;
        return false;
    }
    if (settings.dbs.size() == 0) {
        cerr << "Toml: No DB sections found, check spelling and capatilization, i.e. [Database.DB1]" << endl;
        return false;
    }
    return true;
}

// Loop through a set of config files and create map: species name --> config settings
void getConfigs(vector<string> &configFiles, map<string, SeekSettings> &configs) {
    for (int i=0; i<configFiles.size(); i++) {
        SeekSettings settings;
        parseTomlConfig(configFiles[i], settings);
        if (!settings.species.empty()) {
            // move/copy the local settings variable to the map
            configs[settings.species] = settings;
        }
    }
}

// Parse a legacy style seek db config file and load the CSeekDBSettings
// An example of the legacy db config is in tests/sampleConfigFiles.h:legacyDBConfigData
bool legacyReadDBConfigFile(string dbConfigFile,
                            vector<CSeekDBSetting*> &cc, 
                            CSeekDataset::DistanceMeasure eDistMeasure,
                            bool check_dset_size_flag) {
    const int lineSize = 1024;
    ifstream ifsm;
    ifsm.open(dbConfigFile.c_str());
    if (!ifsm.is_open()) {
        fprintf(stderr, "Error opening file %s\n", dbConfigFile.c_str());
        return false;
    }
    char acBuffer[lineSize];
    uint32_t c_iBuffer = lineSize;
    vector <map<string, string>> parameters; //an array of CDatabase's
    while (!ifsm.eof()) {
        ifsm.getline(acBuffer, c_iBuffer - 1);
        if (acBuffer[0] == 0) break;
        acBuffer[c_iBuffer - 1] = 0;
        string strB = acBuffer;
        if (strB == "START") {
            map <string, string> p;
            while (!ifsm.eof()) {
                ifsm.getline(acBuffer, c_iBuffer - 1);
                if (acBuffer[0] == 0) {
                    fprintf(stderr, "Invalid line (empty)\n");
                    return false;
                }
                strB = acBuffer;
                if (strB == "END") break;
                vector <string> tok;
                CMeta::Tokenize(acBuffer, tok); //separator is tab
                p[tok[0]] = tok[1];
            }
            parameters.push_back(p);
        }
    }
    ifsm.close();
    if (parameters.size() == 0) {
        fprintf(stderr, "Error, extra_db setting file must begin with START and end with END lines\n");
        return false;
    }

    for (int i = 0; i < parameters.size(); i++) {
        string sinfo_dir = "NA";
        string gvar_dir = "NA";
        string platform_dir = "NA";
        string prep_dir = "NA";
        string pcl_dir = "NA";
        string db_dir = "NA";
        string dset_map_file = "NA";
        string gene_map_file = "NA";
        string gene_symbol_file = "NA";
        string quant_file = "NA";
        string dset_size_file = "NA";
        int num_db = -1;

        if (eDistMeasure == CSeekDataset::CORRELATION) {
            if (parameters[i].find("SINFO_DIR") == parameters[i].end() ||
                parameters[i].find("SINFO_DIR")->second == "NA") {
                fprintf(stderr, "Please specify an sinfo directory for the extra db\n");
                return false;
            }
            sinfo_dir = parameters[i].find("SINFO_DIR")->second;
        }
        if (parameters[i].find("GVAR_DIR") != parameters[i].end())
            gvar_dir = parameters[i].find("GVAR_DIR")->second;

        if (parameters[i].find("PCL_DIR") != parameters[i].end())
            pcl_dir = parameters[i].find("PCL_DIR")->second;

        if (check_dset_size_flag == true) {
            if (parameters[i].find("DSET_SIZE_FILE") == parameters[i].end() ||
                parameters[i].find("DSET_SIZE_FILE")->second == "NA") {
                fprintf(stderr, "Please specify the dataset size file for the extra db\n");
                return false;
            }
            dset_size_file = parameters[i].find("DSET_SIZE_FILE")->second;
        }

        if (parameters[i].find("PREP_DIR") == parameters[i].end() ||
            parameters[i].find("PLATFORM_DIR") == parameters[i].end() ||
            parameters[i].find("DB_DIR") == parameters[i].end() ||
            parameters[i].find("DSET_MAP_FILE") == parameters[i].end() ||
            parameters[i].find("GENE_MAP_FILE") == parameters[i].end() ||
            parameters[i].find("QUANT_FILE") == parameters[i].end() ||
            parameters[i].find("NUMBER_OF_DB") == parameters[i].end()) {
            fprintf(stderr, "Some arguments are missing. Please make sure the following are provided:\n");
            fprintf(stderr, "PREP_DIR, DB_DIR, DSET_MAP_FILE, GENE_MAP_FILE, QUANT_FILE, NUMBER_OF_DB\n");
        }

        platform_dir = parameters[i].find("PLATFORM_DIR")->second;
        db_dir = parameters[i].find("DB_DIR")->second;
        prep_dir = parameters[i].find("PREP_DIR")->second;
        dset_map_file = parameters[i].find("DSET_MAP_FILE")->second;
        gene_map_file = parameters[i].find("GENE_MAP_FILE")->second;
        quant_file = parameters[i].find("QUANT_FILE")->second;
        num_db = atoi(parameters[i].find("NUMBER_OF_DB")->second.c_str());

        CSeekDBSetting *dbSetting2 = new CSeekDBSetting(gvar_dir, sinfo_dir,
                                                        platform_dir, prep_dir, db_dir, 
                                                        gene_map_file, gene_symbol_file,
                                                        quant_file, dset_map_file,
                                                        dset_size_file, num_db);
        dbSetting2->setPclDir(pcl_dir);
        cc.push_back(dbSetting2);
    }
    return true;
}


// Loads a vector from a text file containing a single column of entries
void loadOneColumnTextFile(string filename, vector<string> &vals) {
    ifstream fileHandle(filename);
    string line;
    vals.clear();
    while (getline(fileHandle, line)) {
        boost::trim(line);
        if (line.length() > 0) {
            vals.push_back(line);
        }
    }
}

// Loads a map from a text file containing two colunms of entries
// Creates map: col1_vals --> col2_vals
void loadTwoColumnTextFile(string filename, map<string, string> &vals) {
    ifstream fileHandle(filename);
    string line;
    int lineNum = 0;
    vals.clear();
    regex whitespace_regex("\\s+");
    while (getline(fileHandle, line)) {
        boost::trim(line);
        if (line.length() > 0) {
            vector<string> items {
                sregex_token_iterator(line.begin(), line.end(), whitespace_regex, -1), {}
            };
            // vector<string> items;
            // boost::split(items, line, boost::is_any_of("\t "));
            if (items.size() != 2) {
                throw runtime_error("Expecting 2 cols: " + filename + ":" + to_string(lineNum));
            }
            vals[items[0]] = items[1];
        }
        lineNum++;
    }
}


void Semaphore::notify() {
    lock_guard<mutex> lock(_mutex);
    ++_count;
    if (_count > _maxCount) { _count = _maxCount; }
    _condition.notify_one();
}

void Semaphore::wait() {
    unique_lock<mutex> lock(_mutex);
    while(!_count) // Handle spurious wake-ups.
        _condition.wait(lock);
    --_count;
}

bool Semaphore::try_wait() {
    lock_guard<mutex> lock(_mutex);
    if(_count) {
        --_count;
        return true;
    }
    return false;
}


#include <omp.h>
uint32_t omp_enabled_test() {
    uint32_t thread_count = 0;
    #pragma omp parallel num_threads(4)
    {
        #pragma omp critical
        {
            thread_count++;
            std::cout << "tid = " << omp_get_thread_num() << std::endl;
        }
    }
    return thread_count;
}


// Note: If the function implementation is here instead of 
//  in the header, then you must explicitly instantiate
//  each template type you will use, such as:
// template class LRUCache<string, int>;

// template <typename K, typename V>
// void LRUCache<K, V>::set(const K key, const V value) {
//     auto pos = keyValuesMap.find(key);
//     if (pos == keyValuesMap.end()) {
//         items.push_front(key);
//         keyValuesMap[key] = { value, items.begin() };
//         if (keyValuesMap.size() > csize) {
//             keyValuesMap.erase(items.back());
//             items.pop_back();
//         }
//     }
//     else {
//         items.erase(pos->second.second);
//         items.push_front(key);
//         keyValuesMap[key] = { value, items.begin() };
//     }
// }

// template <typename K, typename V>
// bool LRUCache<K, V>::get(const K key, V &value) {
//     auto pos = keyValuesMap.find(key);
//     if (pos == keyValuesMap.end())
//         return false;
//     items.erase(pos->second.second);
//     items.push_front(key);
//     keyValuesMap[key] = { pos->second.first, items.begin() };
//     value = pos->second.first;
//     return true;
// }