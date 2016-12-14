#include "util/stringutil.h"

#include <cstring>
#include <cstdarg>
#include <cstdio>
#include <cassert>
#include <cctype>

namespace Strings {

/**
StringFormat returns a std::string which contains the formatted text.

StringFormat uses the vsnprintf from C to format the string. Therefore, the same
formatting string syntax as for printf is used.

### Note
The code relies on the fact that vsnprintf returns the string length of the
untruncated string if the buffer was not larger enough. Therefore, glibc > 2.0.6
is needed!
*/
std::string Format(const char *fmt, ...) {
    // guess the buffer size
    size_t buffersize = 16 + 2 * strlen(fmt);

    char buffer[buffersize];

    va_list args;

    va_start(args, fmt);
    int printed = vsnprintf(&buffer[0], buffersize, fmt, args);
    va_end(args);

    if (printed < 0) {
        return std::string("error in Strings::Format (1)");
    }

    if ((size_t)printed >= buffersize) {
        buffersize = printed + 1;
    } else {
        return std::string(buffer);
    }

    char buffer2[buffersize];
    va_start(args, fmt);
    printed = vsnprintf(&buffer2[0], buffersize, fmt, args);
    va_end(args);

    if (printed < 0) {
        return std::string("error in Strings::Format (2)");
    }
    if ((size_t)printed >= buffersize) {
        return std::string("error in Strings::Format (3)");
    }

    return std::string(buffer2);
}

/**
 * @brief Remove leading white spaces
 *
 * TrimLeft returns a string with all leading white spaces removed.
 * White spaces are " \t\f\n\r\v".
 */
std::string TrimLeft(const std::string &s) {
    size_t pos = s.find_first_not_of(" \t\f\n\r\v", 0);
    if (pos == std::string::npos) {
        return "";
    };
    size_t len = s.size() - pos;
    return s.substr(pos, len);
}

/**
 * @brief Remove trailing white spaces
 *
 * TrimRight returns a string with all leading white spaces removed.
 * White spaces are " \t\f\n\r\v"
 */
std::string TrimRight(const std::string &s) {
    size_t pos = s.find_last_not_of(" \t\f\n\r\v");
    if (pos == std::string::npos) {
        return "";
    }
    return s.substr(0, pos + 1);
}

/**
 * @brief Remove leading and trailing white spaces
 *
 * TrimRight returns a string with all leading and trailing white spaces
 * removed. White spaces are " \t\f\n\r\v"
 */
std::string Trim(const std::string &s) {
    size_t pos1 = s.find_first_not_of(" \t\f\n\r\v", 0);
    size_t pos2 = s.find_last_not_of(" \t\f\n\r\v");
    if (pos1 == std::string::npos) {
        return "";
    }
    size_t len = pos2 - pos1 + 1;
    return s.substr(pos1, len);
}

/**
 * @brief Remove comment from string
 *
 * StripComment returns a string without any comments starting with sep.
 */
std::string StripComment(const std::string &s, const std::string &sep) {
    size_t pos = s.find(sep, 0);
    if (pos == std::string::npos) {
        return s;
    }
    return s.substr(0, pos);
}

/**
 * @brief split a string
 *
 * Split returns a list of substrings of s which are separated by sep.
 * If sep or s is empty, the returned list contains only s.
 */
StringList Split(const std::string &s, const std::string sep) {
    if (sep.empty()) {
        std::vector<std::string> r;
        r.push_back(s);
        return r;
    }
    size_t start = 0;
    std::vector<std::string> res;
    while (true) {
        size_t pos = s.find(sep, start);
        if (pos == std::string::npos) {
            res.push_back(s.substr(start, s.size() - start));
            break;
        }
        res.push_back(s.substr(start, pos - start));
        start = pos + sep.size();
    }
    return res;
}

/**
 * @brief split a string in N substrings
 *
 * SplitN returns a list of substrings of s which are separated by sep.
 * The returned list has at maximum n elements. The last substring is the
 * unsplit remainder.
 * If sep or s is empty, the returned list contains only s.
 */
StringList SplitN(const std::string &s, const std::string sep, int n) {
    assert(n > 0);
    if (sep.empty() || n == 1) {
        std::vector<std::string> r;
        r.push_back(s);
        return r;
    }

    size_t start = 0;
    std::vector<std::string> res;
    res.reserve(n);
    while (true) {
        size_t pos = s.find(sep, start);
        if (pos == std::string::npos) {
            res.push_back(s.substr(start, s.size() - start));
            break;
        }
        res.push_back(s.substr(start, pos - start));
        start = pos + sep.size();
        if (res.size() == (size_t)n - 1) {
            res.push_back(s.substr(start, std::string::npos));
            break;
        }
    }
    return res;
}

/**
 * @brief Test if string has prefix
 *
 * HasPrefix tests if the string begins with prefix.
 */
bool HasPrefix(const std::string &s, const std::string &prefix) {
    return prefix.size() <= s.size() &&
           s.compare(0, prefix.size(), prefix) == 0;
}

/**
 * Convert string to double
 *
 * ParseDouble uses strtod to parse the string. The parsed
 * value is stored in d.
 *
 * @return
 * It returns 0 if the string was parsed without any errors.
 * If the string contains a double + string, it returns the
 * index of the string, e.g. "12.0hello" sets d to 12.0 and returns
 * 4. It returns -1 if i is NULL. If s is empty, it returns -2 and
 * sets d to 0.0. If s does not start with a number, it returns -3 and sets d
 * =0.0.
 */
int ParseDouble(const std::string &s, double * d) {
    if (d == NULL || d == 0) {
        return -1;
    }
    if (s.empty()) {
        *d = 0.0;
        return -2;
    }
    char *endptr;

    const char *cs = s.c_str();
    *d = strtod(cs, &endptr);
    if (endptr == cs) {
        // no conversion at all
        return -3;
    }
    if (*endptr != '\0') {
        return endptr - cs;
    }
    return 0;
}

/**
 * Convert string to int
 *
 * ParseInt uses strtol to parse the string. The parsed
 * value is stored in i.
 *
 * @return
 * It returns 0 if the string was parsed without any errors.
 * If the string contains an integer + string, it returns the 
 * index of the string, e.g. "12hello" sets i to 12 and returns
 * 2. It returns -1 if i is NULL. If s is empty, it returns -2
 * and sets i to 0.
 */
int ParseInt(const std::string &s, long int *i) {
    if (i == NULL || i == 0) {
        return -1;
    }
    if (s.empty()) {
        *i = 0;
        return -2;
    }
    char *endptr;

    const char *cs = s.c_str();
    *i = strtol(cs, &endptr, 0);
    if (*endptr != '\0') {
        return endptr - cs;
    }
    return 0;
}

/**
 * Parse a comma separated list of integers
 *
 * This functions allows to use ranges, i.e. you can use 1-3 or -1-1, -5--2.
 * 
 * @param s string which is parsed
 * @param dest output lost
 * @param endindex if an error occurrs, the position of it is stored to endindex.
 *
 *
 * @return
 * ParseIntegerList returns 0 on success and endindex ist set to 0.
 *
 */
IntegerListError ParseIntegerList(const std::string &s, IntList *dest,
                                  size_t *endindex) {
    std::vector<std::string> fields = Split(s, ",");
    size_t parsed_len = 0;
    for (size_t n = 0; n < fields.size(); n++) {
        std::string s1 = TrimRight(fields[n]);
        if (s1.empty()) {
            *endindex = parsed_len;
            return IntegerListError::EmptyField;
        }
        long int number = -1;
        int e = ParseInt(s1, &number);
        if (e == 0) {
            dest->push_back(number);
            parsed_len += fields[n].size() + 1;
        } else {
            size_t index_sep = s1.find("-", e);
            if (!Strings::Trim(s1.substr(e, index_sep-e)).empty()) {
               *endindex = parsed_len + e;
              return IntegerListError::UnknownChar;
            } 
            if (index_sep == std::string::npos) {
                *endindex = parsed_len + e;
                return IntegerListError::MissingMinus;
            }
            long int number2 = -1;
            int e2 = ParseInt(s1.substr(index_sep + 1, std::string::npos), &number2);
            if (e2 != 0) {
                *endindex = parsed_len + index_sep + e2 + 1;
                return IntegerListError::UnknownChar;
            }
            if (number > number2) {
                *endindex = parsed_len + e;
                return IntegerListError::InvalidRange;
            }
            for( int i = 0; i < number2 - number + 1; i++) {
                dest->push_back(number + i);
                parsed_len += fields[n].size() + 1;
            }
        }
    }
    *endindex = 0;
    return IntegerListError::NoError;
}

bool IsPrint( const std::string& s ) {
    for( char c : s ) {
        if ( isprint(c) == 0 ) {
            return false;
        }
    }
    return true;
}

DoubleListError ParseDoubleList(const std::string &s, DoubleList *dest,
                                  size_t *endindex) {
  std::vector<std::string> fields = Split(s, ",");
  size_t parsed_len = 0;
  for (size_t n = 0; n < fields.size(); n++) {
    std::string s1 = TrimRight(fields[n]);
    if (s1.empty()) {
            *endindex = parsed_len;
            return DoubleListError::EmptyField;
        }
        double number = -1;
        int e = ParseDouble(s1, &number);
        if (e == 0) {
            dest->push_back(number);
            parsed_len += fields[n].size() + 1;
	} else {
	    *endindex = parsed_len + e;
	    return DoubleListError::UnknownChar;
	}
  }
  *endindex = 0;
  return DoubleListError::NoError;
}

} // namespace
