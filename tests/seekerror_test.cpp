#include <iostream>
#include "gtest/gtest.h"
#include "seekerror.h"

using namespace std;

void func3();

// Setup a fixture class of data for the tests
class SeekErrorTest : public ::testing::Test
{
public:

protected:
    void SetUp() override
    {
        return;
    }
};

// TEST - load toml config files
TEST_F(SeekErrorTest, errorStackTrace)
{
    try {
        func3();
    } catch (const exception &err) {        
        print_exception_stack(err);
    } catch (...) {
        FAIL();
    }
}


void func0() {
  throw init_error(FILELINE + "No arg provided");
}

void func1() {
  try {
    func0();
  } catch(...) {
    throw_with_nested(config_error(FILELINE + "func1 failed with call func0"));
  }
}

void func2() {
  try {
    func1();
  } catch(...) {
    throw_with_nested(logic_error("func2 failed with call func1"));
  }
}

void func3() {
  try {
    func2();
  } catch(...) {
    throw_with_nested(query_error(FILELINE + "func3 failed with call func2"));
  }
}