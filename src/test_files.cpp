#include <fstream>
#include <sstream> //std::stringstream
#include <map>
#include <vector>

# For testing the speed of creating and apending files

class testfile {

  std::fstream m_fstm;
}

int numFiles = 2
char* baseDir = "/tmp"

bool create_file(const std::string &strFilename, const std::string &text) {
  std::fstream fp;
  fp.open(strFilename.c_str(), ios_base::out);
  if (!fp.is_open()) {
    return false;
  }
  fp.write((char *) text.c_str(), text.size);
  fp.close()
  return true;
}

bool append_file(const std::string &strFilename, const std::string &text) {
  std::fstream fp;
  fp.open(strFilename.c_str(), ios_base::out | ios_base::app)
  if (!fp.is_open()) {
    return false;
  }
  fp.write((char *) text.c_str(), text.size);
  fp.close()
  return true;
}

std::string &resStr read_file(const std::string &strFilename) {
  std::ifstream inFile;
  std::stringstream strStream;
  inFile.open(strFilename.c_str());
  if (!fp.is_open()) {
    return false;
  }
  strStream << inFile.rdbuf();
  return strStream.str();
}

int main(void) {
  std::string fileName;
  for (i=0; i<numFiles; i++) {
    fileName = baseDir + "/test" + i;
    cout << fileName;
  }

}
