#include "gtest/gtest.h"
#include "seekplatform.h"
#include <filesystem>

// forward declarations
void calcPlatformStats(Sleipnir::SeekPlatforms &seekPlatforms, vector<vector<vector<float>>> &rawData, 
                       vector<string> &platformNames, map<string, uint32_t> &platformMap); 
            
void printPlatformStats(string name, Sleipnir::SeekPlatforms &seekPlatforms);


// Setup a fixture class of data for the tests
class SeekPlatformsTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    char* username = getenv("USER");
    assert(username != NULL);
    testDir = "/tmp/" + string(username) + "/test/platform";

    if (!filesystem::exists(testDir)) {
      filesystem::create_directories(testDir);
    }
    calcPlatformStats(seekPlatforms1, rawData1, platformNames1, platformMap1);
    calcPlatformStats(seekPlatforms2, rawData2, platformNames2, platformMap2);
    calcPlatformStats(seekPlatforms3, rawData3, platformNames3, platformMap3);
  }

  // void TearDown() override {}

  // Sleipnir::SeekPlatforms seekPlatforms;
  string testDir;

  // SeekPlat1: data for 3 platforms, 3 genes, 4 samples per gene
  Sleipnir::SeekPlatforms seekPlatforms1;
  vector<string> platformNames1 = {"zero", "one", "two"};
  map<string, uint32_t> platformMap1 = {{"zero", 0} , {"one", 1}, {"two", 2}};
  vector<vector<vector<float>>> rawData1{
      {{0.1, 0.2, 0.3, 0.4}, {0.11, 0.12, 0.13, 0.14}, {0.21, 0.22, 0.23, 0.24}},
      {{1, 2, 3, 4}, {11, 12, 13, 14}, {21, 22, 23, 24}},
      {{20, 22, 23, 24}, {31, 32, 33, 34}, {41, 42, 43, 44}}};

  // SeekPlat2:  data for 3 platforms, two of which overlap with rawData1 platforms
  Sleipnir::SeekPlatforms seekPlatforms2;
  vector<string> platformNames2 = {"one", "two", "three"};
  map<string, uint32_t> platformMap2 = {{"one", 0}, {"two", 1}, {"three", 2}};
  vector<vector<vector<float>>> rawData2{
      {{11, 12, 13, 14}, {21, 22, 23, 24}, {31, 32, 33, 34}},
      {{24, 31, 26, 30}, {38, 42, 31, 44}, {32, 31, 30, 40}},
      {{301, 302, 303, 304}, {311, 312, 313, 314}, {321, 322, 323, 324}}};

  // combined data for platforms 1 and 2 above
  Sleipnir::SeekPlatforms seekPlatforms3;
  vector<string> platformNames3 = {"zero", "one", "two", "three"};
  map<string, uint32_t> platformMap3 = {{"zero", 0} , {"one", 1}, {"two", 2}, {"three", 3}};
  vector<vector<vector<float>>> rawData3{
      {{0.1, 0.2, 0.3, 0.4}, {0.11, 0.12, 0.13, 0.14}, {0.21, 0.22, 0.23, 0.24}},
      {{1, 2, 3, 4, 11, 12, 13, 14}, {11, 12, 13, 14, 21, 22, 23, 24}, {21, 22, 23, 24, 31, 32, 33, 34}},
      {{20, 22, 23, 24, 24, 31, 26, 30}, {31, 32, 33, 34, 38, 42, 31, 44}, {41, 42, 43, 44, 32, 31, 30, 40}},
      {{301, 302, 303, 304}, {311, 312, 313, 314}, {321, 322, 323, 324}}};

};


// Helper function to initialize a SeekPlatforms object from raw sample data
void calcPlatformStats(Sleipnir::SeekPlatforms &seekPlatforms,
                       vector<vector<vector<float>>> &rawData, 
                       vector<string> &platformNames,
                       map<string, uint32_t> &platformMap)
{
  uint32_t numPlatforms = rawData.size();
  uint32_t numGenes = rawData[0].size();

  // *** Set the SeekPlatform with the data ***
  seekPlatforms.initialize(numPlatforms, numGenes);
  seekPlatforms.setPlatformNames(platformNames);
  seekPlatforms.setPlatformNameMap(platformMap);

  assert(seekPlatforms.getNumPlatforms() == numPlatforms);
  assert(platformNames.size() == numPlatforms);
  assert(platformMap.size() == numPlatforms);

  // calculate the avg and stdev of the gene data per platform
  for (int plat = 0; plat < numPlatforms; plat++)
  {
    for (int gene = 0; gene < numGenes; gene++)
    {
      float sum = 0;
      float sum_sq = 0;
      uint32_t samplesPerGene = rawData[plat][gene].size();

      for (int samp = 0; samp < samplesPerGene; samp++)
      {
        sum += rawData[plat][gene][samp];
        sum_sq += pow(rawData[plat][gene][samp], 2);
      }
      double mean = sum / samplesPerGene;
      double var = sum_sq / samplesPerGene - pow(mean, 2);
      double stdev = sqrt(var);
      seekPlatforms.platformAvgMatrix.Set(plat, gene, mean);
      seekPlatforms.platformStdevMatrix.Set(plat, gene, stdev);
      seekPlatforms.platformCountMatrix.Set(plat, gene, samplesPerGene);
    }
  }
  seekPlatforms.setCSeekPlatformData();
}

void printPlatformStats(string name, Sleipnir::SeekPlatforms &seekPlatforms)
{
  uint32_t numPlatforms = seekPlatforms.getNumPlatforms();
  uint32_t numGenes = seekPlatforms.getNumGenes();
  vector<string> &platNames = seekPlatforms.getPlatformNames();
  vector<Sleipnir::CSeekPlatform> &platVec = seekPlatforms.getCSeekPlatforms();
  cout << "PlATFORM STATS: " << name << endl;
  for (int platIdx = 0; platIdx < numPlatforms; platIdx++)
  {
    cout << "STATS: platform " << platIdx << ", " << platNames[platIdx] << endl;
    for (int geneIdx = 0; geneIdx < numGenes; geneIdx++)
    {
      printf("STATS: Gene %d, pA avg: %.02f\n", geneIdx, platVec[platIdx].GetPlatformAvg(geneIdx));
    }
  }
}

// Helper function to compare two seekPlatforms
void compareSeekPlatforms(Sleipnir::SeekPlatforms &seekPlatformsA, Sleipnir::SeekPlatforms &seekPlatformsB, bool includesCounts=true) 
{
  uint32_t numPlatforms = seekPlatformsA.getNumPlatforms();
  uint32_t numGenes = seekPlatformsA.getNumGenes();
  // bool includesCounts = seekPlatformsA.bIncludesCounts();

  // loop through and compare the platform names
  vector<string> &platNamesA = seekPlatformsA.getPlatformNames();
  vector<string> &platNamesB = seekPlatformsB.getPlatformNames();
  for (int platIdx = 0; platIdx < numPlatforms; platIdx++)
  {
    ASSERT_EQ(platNamesA[platIdx], platNamesB[platIdx]);
  }

  // loop through and compare the platform order map
  map<string, utype> &platMapA = seekPlatformsA.getPlatformMap();
  map<string, utype> &platMapB = seekPlatformsB.getPlatformMap();
  for (int platIdx = 0; platIdx < numPlatforms; platIdx++)
  {
    string name = platNamesA[platIdx];
    ASSERT_EQ(platMapA[name], platMapB[name]);
  }

  // loop through and compare values between platforms 1 and 2
  vector<Sleipnir::CSeekPlatform> &platVecA = seekPlatformsA.getCSeekPlatforms();
  vector<Sleipnir::CSeekPlatform> &platVecB = seekPlatformsB.getCSeekPlatforms();
  for (int platIdx = 0; platIdx < numPlatforms; platIdx++)
  {
    // cout << "COMPARE platform " << platIdx << ", " << platNamesA[platIdx] << endl;
    for (int geneIdx = 0; geneIdx < numGenes; geneIdx++)
    {
      // printf("Gene %d, pA avg: %.02f, pB avg: %.02f\n", geneIdx, platVecA[platIdx].GetPlatformAvg(geneIdx), platVecB[platIdx].GetPlatformAvg(geneIdx));
      // printf("Gene %d, pA stdev: %.02f, pB stdev: %.02f\n", geneIdx, platVecA[platIdx].GetPlatformStdev(geneIdx), platVecB[platIdx].GetPlatformStdev(geneIdx));
      ASSERT_EQ(platVecA[platIdx].GetPlatformAvg(geneIdx), platVecB[platIdx].GetPlatformAvg(geneIdx));
      ASSERT_EQ(platVecA[platIdx].GetPlatformStdev(geneIdx), platVecB[platIdx].GetPlatformStdev(geneIdx));
      if (includesCounts) {
        // printf("Gene %d, pA count: %d, pB count: %d\n", geneIdx, platVecA[platIdx].GetPlatformCount(geneIdx), platVecB[platIdx].GetPlatformCount(geneIdx));
        ASSERT_EQ(platVecA[platIdx].GetPlatformCount(geneIdx), platVecB[platIdx].GetPlatformCount(geneIdx));
      }
    }
  }

}


// TEST saveLoadData - save out the seekPlatform to files and reload it
TEST_F(SeekPlatformsTest, saveLoadData)
{
  // Save current platform data to files
  seekPlatforms1.savePlatformDataToFiles(testDir);

  // Load platform data from files to a new platform instance
  Sleipnir::SeekPlatforms seekPlatformsReload;
  seekPlatformsReload.loadPlatformDataFromFiles(testDir);

  // Make sure the loaded platform matches the original
  compareSeekPlatforms(seekPlatforms1, seekPlatformsReload);
}


// TEST loadNoCounts - save out seekPlatform but without a counts file and reload it
TEST_F(SeekPlatformsTest, loadNoCounts)
{
  // Save current platform data to files
  seekPlatforms1.savePlatformDataToFiles(testDir);

  // Remove the counts file
  filesystem::remove(testDir + "/all_platforms.gplatcount");

  // Load platform data from files to a new platform instance
  Sleipnir::SeekPlatforms seekPlatformsReload;
  seekPlatformsReload.loadPlatformDataFromFiles(testDir);
  ASSERT_FALSE(seekPlatformsReload.bIncludesCounts());

  // Make sure the loaded platform matches the original
  compareSeekPlatforms(seekPlatforms1, seekPlatformsReload, false);
}


// TEST SeekPlatform Copy
TEST_F(SeekPlatformsTest, copyPlatform)
{
  auto *srcPlatform = new Sleipnir::SeekPlatforms;
  Sleipnir::SeekPlatforms platDst1;
  Sleipnir::SeekPlatforms platDst2;

  calcPlatformStats(*srcPlatform, rawData1, platformNames1, platformMap1);

  platDst1.copy(*srcPlatform);
  platDst2.copy(platDst1);

  // delete original data to make sure copy didn't just copy references
  delete srcPlatform;

  // Compare the two copies
  compareSeekPlatforms(platDst1, platDst2);
  compareSeekPlatforms(platDst2, seekPlatforms1);
}


// TEST combine function
TEST_F(SeekPlatformsTest, combinePlatforms)
{
  // printPlatformStats("Platform1", seekPlatforms1);
  // printPlatformStats("Platform2", seekPlatforms2);
  // printPlatformStats("Platform3", seekPlatforms3);

  seekPlatforms1.combineWithPlatform(seekPlatforms2);
  // cout << "COMPARE PLATFORMS 1, 3" << endl;
  compareSeekPlatforms(seekPlatforms1, seekPlatforms3);
}


// TEST combine an empty seekPlatforms with a new seekPlatforms
TEST_F(SeekPlatformsTest, combineEmpty)
{
  // Test combining a filled platform into an empty platform
  Sleipnir::SeekPlatforms emptyPlatform;
  emptyPlatform.combineWithPlatform(seekPlatforms1);
  compareSeekPlatforms(emptyPlatform, seekPlatforms1);

  // Test combining an empty platform into a filled platform
  Sleipnir::SeekPlatforms emptyPlatform2;
  Sleipnir::SeekPlatforms p2test;
  p2test.copy(seekPlatforms2);
  p2test.combineWithPlatform(emptyPlatform2);
  compareSeekPlatforms(p2test, seekPlatforms2);
}


// TEST combine function in SeekPrep
TEST_F(SeekPlatformsTest, combineUtilityFunction)
{
  string p1_dir = testDir + "/p1";
  string p2_dir = testDir + "/p2";
  string combined_dir = testDir + "/combine12";
  filesystem::create_directories(p1_dir);
  filesystem::create_directories(p2_dir);
  filesystem::create_directories(combined_dir);

  // Save current platform data to files
  seekPlatforms1.savePlatformDataToFiles(p1_dir);
  seekPlatforms2.savePlatformDataToFiles(p2_dir);

  // Call the utility function that combines the statistic files
  combinePlatformStatistics(p1_dir, p2_dir, combined_dir);

  // Load the combined stats files and compare to platform3.
  Sleipnir::SeekPlatforms seekPlatformsCombined;
  seekPlatformsCombined.loadPlatformDataFromFiles(combined_dir);
  compareSeekPlatforms(seekPlatforms3, seekPlatformsCombined);

  // If SeekPrep is used outside of these tests and the results put in
  // testDir/SeekPrepCombine then read is those results and compare to platform3
  string seekPrep_dir = testDir + "/SeekPrepCombine";
  if (filesystem::exists(seekPrep_dir)) {
    cout << "Verify SeekPrepCombine resutls" << endl;
    Sleipnir::SeekPlatforms seekPrepCombine;
    seekPrepCombine.loadPlatformDataFromFiles(seekPrep_dir);
    compareSeekPlatforms(seekPlatforms3, seekPrepCombine);
  }
}

