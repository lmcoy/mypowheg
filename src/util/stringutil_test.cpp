#include "gtest/gtest.h"

#include <cmath>

#include "util/stringutil.h"

TEST(Format, NoFormatString) {
    std::string s1 = Strings::Format("abcdef");
    EXPECT_EQ(s1, "abcdef");

    s1 = Strings::Format(
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");
}

TEST(Format, SimpleFormatString) {
    EXPECT_EQ(Strings::Format("%d", 4), "4");
    EXPECT_EQ(Strings::Format("%.6f", 0.1), "0.100000");
    EXPECT_EQ(Strings::Format("%d,%d,%d", 4, 5, 6), "4,5,6");
}

TEST(StringFormat, Error) {
    EXPECT_EQ(Strings::Format("%"), "error in Strings::Format (1)");
}

// test whether StringFormat works when the internal buffer is too small and a
// larger one has to be used.
TEST(Format, LongSubstitute) {
    std::string s("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567"
                  "89abcdefghijklmnopqrstuvwxyz");
    EXPECT_EQ(Strings::Format("%s", s.c_str()), s);
}

TEST(TrimLeft, StringWithWhitespaces) {
    std::string s(" \t  \t\t  hallo  ");
    EXPECT_EQ(Strings::TrimLeft(s), std::string("hallo  "));
}

TEST(TrimLeft, StringWithoutWhitespaces) {
    std::string s("hallo welt!");
    EXPECT_EQ(Strings::TrimLeft(s), s);
}

TEST(TrimLeft, EmptyString) {
    std::string s("");
    EXPECT_EQ(Strings::TrimLeft(s), s);
}

TEST(TrimLeft, OnlyWhitespaces) {
    std::string s("  \t   \t\t\n");
    EXPECT_EQ(Strings::TrimLeft(s), std::string(""));
}


TEST(TrimRight, StringWithWhitespaces) {
    std::string s("  hallo \n \t   \r\n");
    EXPECT_EQ(Strings::TrimRight(s), std::string("  hallo"));
}

TEST(TrimRight, StringWithoutWhitespaces) {
    std::string s("hallo welt!");
    EXPECT_EQ(Strings::TrimRight(s), s);
}

TEST(TrimRight, EmptyString) {
    std::string s("");
    EXPECT_EQ(Strings::TrimRight(s), s);
}

TEST(TrimRight, OnlyWhitespaces) {
    std::string s("  \t   \t\t\n");
    EXPECT_EQ(Strings::TrimRight(s), std::string(""));
}


TEST(Trim, StringWithWhitespaces) {
    std::string s("  hallo \n \t   \r\n");
    EXPECT_EQ(Strings::Trim(s), std::string("hallo"));
}

TEST(Trim, StringWithoutWhitespaces) {
    std::string s("hallo welt!");
    EXPECT_EQ(Strings::Trim(s), s);
}

TEST(Trim, EmptyString) {
    std::string s("");
    EXPECT_EQ(Strings::Trim(s), s);
}

TEST(Trim, OnlyWhitespaces) {
    std::string s("  \t   \t\t\n");
    EXPECT_EQ(Strings::Trim(s), std::string(""));
}

TEST(StripComment, EmptyString) {
    std::string s("");
    EXPECT_EQ( Strings::StripComment(s, "#"), std::string("") );
}

TEST(StripComment, NoComment) {
    std::string s("hallo welt!");
    EXPECT_EQ(Strings::StripComment(s, "#"), s);
}

TEST(StripComment, WithComment) {
    std::string s("hallo welt! # comment");
    EXPECT_EQ(Strings::StripComment(s, "#"), std::string("hallo welt! "));
}

TEST(StripComment, WithTwoComments) {
    std::string s("hallo welt! # comment 1 # comment 2");
    EXPECT_EQ(Strings::StripComment(s, "#"), std::string("hallo welt! ") );
}

TEST(StripComment, WithCPPLikeComment) {
    std::string s("hallo welt! // comment");
    EXPECT_EQ(Strings::StripComment( s, "//"), "hallo welt! ");
}

TEST(Split, EmptyString) {
    std::string s("");
    std::vector<std::string> res = Strings::Split(s, "#");
    ASSERT_EQ(res.size(), 1u);
    EXPECT_EQ(res[0], std::string(""));
}

TEST(Split, MulitpleSep) {
    std::string s("a,b,c");
    std::vector<std::string> res = Strings::Split(s, ",");
    ASSERT_EQ(res.size(), 3u);
    EXPECT_EQ(res[0], "a");
    EXPECT_EQ(res[1], "b");
    EXPECT_EQ(res[2], "c");
}

TEST(Split, NoSep) {
    std::string s = "hallo welt";
    std::vector<std::string> res = Strings::Split(s, ",");
    ASSERT_EQ(res.size(), 1u);
    EXPECT_EQ(res[0], s);
}

TEST(Split, EmptySep) {
    std::string s = "abc";
    std::vector<std::string> res = Strings::Split(s, "");
    ASSERT_EQ(res.size(), 1u);
    EXPECT_EQ(res[0], s);
}

TEST(SplitN, EmptyString) {
    std::string s("");
    std::vector<std::string> res = Strings::SplitN(s, "#", 2);
    ASSERT_EQ(res.size(), 1u);
    EXPECT_EQ(res[0], std::string(""));
}

TEST(SplitN, MulitpleSep) {
    std::string s("a,b,c");
    std::vector<std::string> res = Strings::SplitN(s, ",", 2);
    ASSERT_EQ(res.size(), 2u);
    EXPECT_EQ(res[0], "a");
    EXPECT_EQ(res[1], "b,c");
}

TEST(SplitN, NoSep) {
    std::string s = "hallo welt";
    std::vector<std::string> res = Strings::SplitN(s, ",", 2);
    ASSERT_EQ(res.size(), 1u);
    EXPECT_EQ(res[0], s);
}

TEST(SplitN, EmptySep) {
    std::string s = "abc";
    std::vector<std::string> res = Strings::SplitN(s, "", 2);
    ASSERT_EQ(res.size(), 1u);
    EXPECT_EQ(res[0], s);
}

TEST(HasPrefix, PrefixTooLong) {
    std::string s = "hallo";
    std::string prefix = "hallo welt";
    EXPECT_FALSE( Strings::HasPrefix(s, prefix) );
}

TEST(HasPrefix, StringWithPrefix) {
    std::string s1 = "prefix string";
    std::string prefix1 = "prefix";
    EXPECT_TRUE(Strings::HasPrefix(s1, prefix1));

    std::string s2 = "string";
    std::string prefix2 = "str";
    EXPECT_TRUE(Strings::HasPrefix(s2, prefix2));
}

TEST(HasPrefix, StringWithoutPrefix) { 
    std::string s1 = "string without prefix";
    std::string prefix1 = "prefix";
    EXPECT_FALSE(Strings::HasPrefix(s1, prefix1));

    std::string s2 = "astring";
    std::string prefix2 = "str";
    EXPECT_FALSE(Strings::HasPrefix(s2, prefix2));
}

TEST(HasPrefix, EmptyPrefix) {
    std::string s = "hallo welt";
    std::string prefix = "";
    EXPECT_TRUE(Strings::HasPrefix(s, prefix));
}

TEST(HasPrefix, EmptyStrings) {
    EXPECT_TRUE(Strings::HasPrefix("", ""));
}

TEST(HasPrefix, StringEqualsPrefix) {
    EXPECT_TRUE(Strings::HasPrefix("string", "string"));
}

TEST(ParseDouble, Numbers) {
    double d = 0.0;
    EXPECT_EQ(Strings::ParseDouble("13.4", &d), 0);
    EXPECT_DOUBLE_EQ(d, 13.4);
    EXPECT_EQ(Strings::ParseDouble("1.0e-5", &d), 0);
    EXPECT_DOUBLE_EQ(d, 1.0e-5);
}

TEST(ParseDouble, NotANumber) {
    double d = -1.0;
    EXPECT_EQ(Strings::ParseDouble("hallo", &d), 0);
    EXPECT_DOUBLE_EQ(d, 0.0);
}

TEST(ParseDouble, NULLPTR) {
    EXPECT_EQ(Strings::ParseDouble("1.0", NULL), -1);
}

TEST(ParseDouble, StartsWithNumber) {
    double d = -1.0;
    EXPECT_EQ(Strings::ParseDouble("1.0hallo", &d), 3);
    EXPECT_DOUBLE_EQ(d, 1.0);
}

TEST(ParseDouble, Infinity) {
    double d = -1.0;
    EXPECT_EQ(Strings::ParseDouble("Infinity", &d), 0);
    EXPECT_TRUE( isinf(d) );
}

TEST(ParseDouble, EmptyString) {
    double d = -1.0;
    EXPECT_EQ(Strings::ParseDouble("", &d), -2);
    EXPECT_DOUBLE_EQ(d, 0.0);
}

TEST(ParseInt, Numbers) {
    long int d = 0.0;
    EXPECT_EQ(Strings::ParseInt("13", &d), 0);
    EXPECT_EQ(d, 13);
    EXPECT_EQ(Strings::ParseInt("-1", &d), 0);
    EXPECT_EQ(d, -1);
    EXPECT_EQ(Strings::ParseInt("0xfa", &d), 0);
    EXPECT_EQ(d, 0xfa);
}

TEST(ParseInt, NotANumber) {
    long int d = -1;
    EXPECT_EQ(Strings::ParseInt("hallo", &d), 0);
    EXPECT_EQ(d, 0);
}

TEST(ParseInt, NULLPTR) {
    EXPECT_EQ(Strings::ParseInt("1.0", NULL), -1);
}

TEST(ParseInt, StartsWithNumber) {
    long int d = -1;
    EXPECT_EQ(Strings::ParseInt("1hallo", &d), 1);
    EXPECT_EQ(d, 1);
}

TEST(ParseInt, EmptyString) {
    long int d = 1;
    EXPECT_EQ(Strings::ParseInt("", &d), -2);
    EXPECT_EQ(d, 0);
}


TEST(ParseIntegerList, SimpleList) {
    size_t end;
    std::vector<int> list;
    auto e = Strings::ParseIntegerList("1,2,3", &list, &end);
    EXPECT_EQ(e, Strings::IntegerListError::NoError);
    ASSERT_EQ(list.size(), 3u);
    EXPECT_EQ(list[0], 1);
    EXPECT_EQ(list[1], 2);
    EXPECT_EQ(list[2], 3);
}

TEST(ParseIntegerList, Range) {
    size_t end;
    std::vector<int> list;
    auto e = Strings::ParseIntegerList("1-3", &list, &end);
    EXPECT_EQ(e, Strings::IntegerListError::NoError);
    ASSERT_EQ(list.size(), 3u);
    EXPECT_EQ(list[0], 1);
    EXPECT_EQ(list[1], 2);
    EXPECT_EQ(list[2], 3);
}

TEST(ParseIntegerList, EmptyField) {
    size_t end;
    std::vector<int> list;
    std::string s = "1,   ,2";
    auto e = Strings::ParseIntegerList(s, &list, &end);
    EXPECT_EQ(e, Strings::IntegerListError::EmptyField);
    EXPECT_EQ(end, 2u);
    ASSERT_EQ(list.size(), 1u);
} 

TEST(ParseIntegerList, List) {
    size_t end;
    std::vector<int> list;
    std::string s = "1, -2 --1  ,2";
    auto e = Strings::ParseIntegerList(s, &list, &end);
    EXPECT_EQ(e, Strings::IntegerListError::NoError);
    EXPECT_EQ(end, 0u);
    ASSERT_EQ(list.size(), 4u);
    EXPECT_EQ(list[0], 1);
    EXPECT_EQ(list[1], -2);
    EXPECT_EQ(list[2], -1);
    EXPECT_EQ(list[3], 2);
}

TEST(ParseIntegerList, UnknownChar) {
    size_t end;
    std::vector<int> list;
    std::string s = "1 aa";
    auto e = Strings::ParseIntegerList(s, &list, &end);
    EXPECT_EQ(e, Strings::IntegerListError::UnknownChar);
    EXPECT_EQ(end, 1u);
    ASSERT_EQ(list.size(), 0u);
}

TEST(ParseIntegerList, UnknownCharInRange) {
    size_t end;
    std::vector<int> list;
    std::string s = "1  - 3 aa";
    auto e = Strings::ParseIntegerList(s, &list, &end);
    EXPECT_EQ(e, Strings::IntegerListError::UnknownChar);
    EXPECT_EQ(end, 6u);
    ASSERT_EQ(list.size(), 0u);
}

TEST(IsPrint, NormalString) {
    EXPECT_TRUE(Strings::IsPrint("hallo welt!") );
    EXPECT_TRUE(Strings::IsPrint("aab"));
}

TEST(IsPrint, NotPrintable) {
    EXPECT_FALSE(Strings::IsPrint("\taad"));
}
