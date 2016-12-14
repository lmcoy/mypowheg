#ifndef _STRINGUTIL_H_
#define _STRINGUTIL_H_

#include <string>
#include <vector>

/**
\brief functions for strings.

\ingroup Util
 */
namespace Strings {

typedef std::vector<std::string> StringList;
typedef std::vector<int> IntList;
typedef std::vector<double> DoubleList;

std::string Format(const char *, ...);

std::string TrimLeft(const std::string &s); 
std::string TrimRight(const std::string &s); 
std::string Trim(const std::string &s);

std::string StripComment(const std::string &s, const std::string &del);

StringList Split(const std::string &s, const std::string del);

StringList SplitN(const std::string &s, const std::string sep,
                                int n);

bool HasPrefix(const std::string &s, const std::string &prefix);

int ParseDouble(const std::string &s, double * d);
int ParseInt(const std::string &s, long int *i);

enum class IntegerListError {
    NoError,
    EmptyField,
    UnknownChar,
    MissingMinus,
    InvalidRange
};

IntegerListError ParseIntegerList(const std::string &s, IntList *dest,
                                  size_t *endindex);

enum class DoubleListError {
    NoError,
    EmptyField,
    UnknownChar
};

DoubleListError ParseDoubleList(const std::string &s, DoubleList *dest,
                                  size_t *endindex);

bool IsPrint( const std::string& s );

} // namespace

#endif
