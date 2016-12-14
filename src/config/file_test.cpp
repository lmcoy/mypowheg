#include "gtest/gtest.h"

#include "config/file.h"

#include <sstream>
#include <iostream>

class FileTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        strstr << "[cuts]\n\teta = (-2.5, 2.5)\n\tpT = (25, inf)\n\tmll "
                  "= (50, inf)\n\n[recombination]\n\tdR = 0.1\n";
    }

    std::stringstream strstr;
};

using namespace Config;

TEST_F(FileTest, Read) {
    File file;
    File::ReadError err = file.Read(strstr);   
    ASSERT_EQ(err, File::ReadError::NoError);

    double dR = -1.0;
    File::Error e = file.GetDouble("recombination", "dR", &dR);
    ASSERT_EQ(e, File::Error::NoError);
    EXPECT_DOUBLE_EQ(dR, 0.1);


    double min = -1.0;
    double max = -1.0;
    e = file.GetDoubleInterval("cuts", "eta", &min, &max);
    ASSERT_EQ(e, File::Error::NoError);
    EXPECT_DOUBLE_EQ(min, -2.5);
    EXPECT_DOUBLE_EQ(max, 2.5);

    e = file.GetDoubleInterval("cuts", "pT", &min, &max);
    ASSERT_EQ(e, File::Error::NoError);
    EXPECT_DOUBLE_EQ(min, 25);
    EXPECT_TRUE(1e9 < max);
}
