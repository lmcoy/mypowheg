#ifndef CONFIG_FILE_H
#define CONFIG_FILE_H

#include <map>
#include <string>
#include <vector>
#include <istream>

#include "util/stringutil.h"

namespace Config {

class File {
  public:
    enum class Error {
        NoError,
        CategoryNotFound,
        KeyNotFound,
        DoubleParser,
        IntParser,
        IntListParser,
        DoubleIntervalParser,
    };

    enum class ReadError {
        NoError,
        ParserError,
        FileError,
    };

    ReadError ReadFromFile(const char *filename);

    ReadError Read(std::istream &);

    Error GetDouble(const std::string &category, const std::string &name,
                    double *dest);

    Error GetInt(const std::string &category, const std::string &name,
                 long int *dest);

    Error GetString(const std::string &category, const std::string &name,
                    std::string *dest);

    Error GetStringList(const std::string &category, const std::string &key,
                        Strings::StringList *dest);

    Error GetIntList(const std::string &category, const std::string &key,
                     Strings::IntList *dest);

    Error GetDoubleInterval(const std::string &category, const std::string &key,
                            double *min, double *max);

    const std::string &ErrorMsg() const;

  private:
    typedef std::map<std::string, std::string> FieldList;
    std::map<std::string, FieldList> fields_;
    std::string errormsg_;
};

} // namespace Config

#endif

