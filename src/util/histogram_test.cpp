#include "gtest/gtest.h"

#include "util/histogram.h"

TEST(HistogramGetBin, TwoBinsRight) {
    Util::Histogram hist(2, 0.0, 10.0);

    int index = hist.GetBin(7.5);
    EXPECT_EQ(index, 1);
}

TEST(HistogramGetBin, TwoBinsLeft) {
    Util::Histogram hist(2, 0.0, 10.0);

    int index = hist.GetBin(1.5);
    EXPECT_EQ(index, 0);
}

TEST(HistogramGetBin, TwoBinsMid) {
    Util::Histogram hist(2, 0.0, 10.0);

    int index = hist.GetBin(5.0);
    EXPECT_EQ(index, 0);
}

TEST(HistogramGetBin, Bin) {
    Util::Histogram hist(8, 0.1, 20.0);

    int index = hist.GetBin(1.0);
    EXPECT_EQ(index, 0);

    index = hist.GetBin(3.0);
    EXPECT_EQ(index, 1);

    index = hist.GetBin(5.5);
    EXPECT_EQ(index, 2);

    index = hist.GetBin(8.5);
    EXPECT_EQ(index, 3);

    index = hist.GetBin(10.5);
    EXPECT_EQ(index, 4);

    index = hist.GetBin(13.5);
    EXPECT_EQ(index, 5);

    index = hist.GetBin(15.5);
    EXPECT_EQ(index, 6);

    index = hist.GetBin(18.5);
    EXPECT_EQ(index, 7);
}

TEST(HistogramGetBin, OutOfBound) {
    Util::Histogram hist(8, 0.1, 20.0);
    int index = hist.GetBin(0.0);
    EXPECT_EQ(index, -1);
    index = hist.GetBin(30.0);
    EXPECT_EQ(index, -2);
}

TEST(Histogram, GetValue) {
    Util::Histogram hist(8, 0.1, 20.0);
    double v = -1.0;
    bool r = hist.Get(2.5, &v);
    EXPECT_EQ(r, true);
    EXPECT_DOUBLE_EQ(v, 0.0);
}

TEST(Histogram, SetValue) {
    Util::Histogram hist(8, 0.0, 20.0);
    double v = 12.3;
    bool r = hist.SetValue(3.5, v);
    EXPECT_EQ(r, true);
    double v2 = -1.0;
    r = hist.Get(3.5, &v2);
    EXPECT_DOUBLE_EQ(v, v2);
}

TEST(Histogram, GetXLower) {
    Util::Histogram hist(8, 0.0, 20.0);

    EXPECT_DOUBLE_EQ(hist.GetXLower(0), 0.0);
    EXPECT_DOUBLE_EQ(hist.GetXLower(1), 2.5);
    EXPECT_DOUBLE_EQ(hist.GetXLower(2), 5.0);
    EXPECT_DOUBLE_EQ(hist.GetXLower(3), 7.5);
    EXPECT_DOUBLE_EQ(hist.GetXLower(4), 10.0);
    EXPECT_DOUBLE_EQ(hist.GetXLower(5), 12.5);
    EXPECT_DOUBLE_EQ(hist.GetXLower(6), 15.0);
    EXPECT_DOUBLE_EQ(hist.GetXLower(7), 17.5);
}

TEST(Histogram, GetXUpper) {
    Util::Histogram hist(8, 0.0, 20.0);

    EXPECT_DOUBLE_EQ(hist.GetXUpper(0), 2.5);
    EXPECT_DOUBLE_EQ(hist.GetXUpper(1), 5.0);
    EXPECT_DOUBLE_EQ(hist.GetXUpper(2), 7.5);
    EXPECT_DOUBLE_EQ(hist.GetXUpper(3), 10.0);
    EXPECT_DOUBLE_EQ(hist.GetXUpper(4), 12.5);
    EXPECT_DOUBLE_EQ(hist.GetXUpper(5), 15.0);
    EXPECT_DOUBLE_EQ(hist.GetXUpper(6), 17.5);
    EXPECT_DOUBLE_EQ(hist.GetXUpper(7), 20.0);
}
        
TEST(HistogramGetBin, Between2) {
    Util::Histogram hist(5, 0.0, 250.0);

    EXPECT_EQ(hist.GetBin(105.466), 2);
}

TEST(Histogram, AddHistogram) {
    Util::Histogram hist1(5, 0.0, 250.0);
    Util::Histogram hist2(5, 0.0, 250.0);

    hist1.SetValue(0.5, 1.3);
    hist2.SetValue(0.5, 1.4);

    hist1.SetValue(125.0, 1.0);
    hist2.SetValue(125.0, 2.0);

    hist1.AddHistogram( &hist2);
    double v;
    hist1.Get(0.5, &v);
    EXPECT_DOUBLE_EQ( v, 2.7);

    hist1.Get(125.0, &v);
    EXPECT_DOUBLE_EQ(v, 3.0);
}
