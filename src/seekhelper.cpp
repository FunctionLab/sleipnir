#include <fstream>
#include <map>
#include <vector>
#include "seekdataset.h"

using namespace std;
using namespace Sleipnir;

void readTomlConfig() {

}

bool legacyReadDBConfigFile(string dbConfigFile,
                            vector<CSeekDBSetting*> &cc, 
                            enum CSeekDataset::DistanceMeasure eDistMeasure = CSeekDataset::CORRELATION,
                            bool check_dset_size_flag = true) {
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
        string db_dir = "NA";
        string dset_map_file = "NA";
        string gene_map_file = "NA";
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
                                                        platform_dir, prep_dir, db_dir, gene_map_file, quant_file,
                                                        dset_map_file,
                                                        dset_size_file, num_db);
        cc.push_back(dbSetting2);
    }
    return true;
}