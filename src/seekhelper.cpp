#include <map>
#include <vector>
#include <fstream>
#include <regex>
#include <boost/algorithm/string.hpp>
#include "seekhelper.h"

using namespace std;


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

