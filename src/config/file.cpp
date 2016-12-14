#include "config/file.h"

#include <cfloat>
#include <fstream>

namespace Config {

File::ReadError File::ReadFromFile(const char *filename) {
    std::ifstream istr(filename);
    if (!istr.is_open()) {
        errormsg_ = Strings::Format("Couldn't open file \"%s\"", filename);
        return File::ReadError::FileError;
    }
    return Read(istr);
}

File::ReadError File::Read(std::istream &is) {
    std::string line = "";
    std::string category = "";
    bool incat = false;
    int linenr = 0;
    while (!is.eof()) {
        std::getline(is, line);
        linenr += 1;
        line = Strings::Trim(line);
        if (line.empty()) {
            continue;
        }
        if (line[0] == '[') {
            if (line[line.size() - 1] != ']') {
                errormsg_ = Strings::Format("l %d, no ] found", linenr);
                return ReadError::ParserError;
            }
            std::string tmp = line.substr(1, line.size() - 2);
            if (!Strings::IsPrint(tmp)) {
                errormsg_ = Strings::Format(
                    "l %d, category name is not printable", linenr);
                return ReadError::ParserError;
            }
            category = tmp;
            incat = true;
            continue;
        }
        auto fields = Strings::SplitN(line, "=", 2);
        if (fields.size() != 2) {
            errormsg_ = Strings::Format("l %d, no = found", linenr);
            return ReadError::ParserError;
        }
        if (!incat) {
            errormsg_ = Strings::Format("l %d, no category specified", linenr);
            return ReadError::ParserError;
        }
        std::string key = Strings::Trim(fields[0]);
        std::string value = Strings::Trim(fields[1]);
        if (!Strings::IsPrint(key)) {
            errormsg_ =
                Strings::Format("l %d, key name is not printable", linenr);
            return ReadError::ParserError;
        }
        if (!Strings::IsPrint(value)) {
            errormsg_ = Strings::Format("l %d, value is not printable", linenr);
            return ReadError::ParserError;
        }

        auto ret = fields_[category].insert(make_pair(key, value));
        if (!ret.second) {
            errormsg_ = Strings::Format("l %d, key %s exists already", linenr,
                                        key.c_str());
            return ReadError::ParserError;
        }
    }

    return ReadError::NoError;
}

const std::string &File::ErrorMsg() const { return errormsg_; }

File::Error File::GetString(const std::string &category, const std::string &key,
                            std::string *dest) {
    auto category_it = fields_.find(category);
    if (category_it == fields_.end()) {
        *dest = "";
        errormsg_ = "unknown category: \"" + category + "\"";
        return Error::CategoryNotFound;
    }

    FieldList *fields = &category_it->second;
    auto key_it = fields->find(key);
    if (key_it == fields->end()) {
        *dest = "";
        errormsg_ =
            "key \"" + key + "\" not found in category \"" + category + "\"";
        return Error::KeyNotFound;
    }

    *dest = key_it->second;

    errormsg_ = "";
    return Error::NoError;
}

File::Error File::GetDouble(const std::string &category, const std::string &key,
                            double *dest) {
    std::string str;
    Error error = GetString(category, key, &str);
    if (error != Error::NoError) {
        *dest = 0.0;
        return error;
    }

    int ret = Strings::ParseDouble(str, dest);
    if (ret > 0) {
        errormsg_ = Strings::Format(
            "error in %s.%s: \"%s\" of \"%s\" is not a valid double.",
            category.c_str(), key.c_str(),
            str.substr(ret, std::string::npos).c_str(), str.c_str());
        return Error::DoubleParser;
    }
    if (ret == -2) {
        errormsg_ = Strings::Format("error in %s.%s: value is empty.",
                                    category.c_str(), key.c_str());
        return Error::DoubleParser;
    }

    errormsg_ = "";
    return Error::NoError;
}

File::Error File::GetInt(const std::string &category, const std::string &key,
                         long int *dest) {
    std::string str;
    Error error = GetString(category, key, &str);
    if (error != Error::NoError) {
        *dest = 0;
        return error;
    }

    int ret = Strings::ParseInt(str, dest);
    if (ret > 0) {
        errormsg_ = Strings::Format(
            "error in %s.%s: \"%s\" of \"%s\" is not a valid int.",
            category.c_str(), key.c_str(),
            str.substr(ret, std::string::npos).c_str(), str.c_str());
        return Error::DoubleParser;
    }
    if (ret == -2) {
        errormsg_ = Strings::Format("error in %s.%s: value is empty.",
                                    category.c_str(), key.c_str());
        return Error::DoubleParser;
    }

    errormsg_ = "";
    return Error::NoError;
}

File::Error File::GetStringList(const std::string &category,
                                const std::string &key, Strings::StringList *dest) {
    std::string str;
    Error error = GetString(category, key, &str);
    if (error != Error::NoError) {
        *dest = Strings::StringList();
        return error;
    }

    *dest = Strings::Split(str, ",");
    for (auto &s : *dest) {
        s = Strings::Trim(s);
    }

    errormsg_ = "";
    return Error::NoError;
}

File::Error File::GetIntList(const std::string &category,
                             const std::string &key, Strings::IntList *dest) {
    std::string str;
    Error error = GetString(category, key, &str);
    if (error != Error::NoError) {
        *dest = Strings::IntList();
        return error;
    }

    size_t endindex = 0;
    Strings::IntegerListError err =
        Strings::ParseIntegerList(str, dest, &endindex);
    using namespace Strings;
    switch (err) {
    case IntegerListError::InvalidRange:
        errormsg_ = Format("invalid range at %d", endindex);
        return Error::IntListParser;
    case IntegerListError::MissingMinus:
        errormsg_ = Format("missing minus or comma at %d", endindex);
        return Error::IntListParser;
    case IntegerListError::UnknownChar:
        errormsg_ =
            Format("unknown char \"%c\" at %d", str[endindex], endindex);
        return Error::IntListParser;
    case IntegerListError::EmptyField:
        errormsg_ = Format("empty field at %d", endindex);
        return Error::IntListParser;
    case IntegerListError::NoError:
        break;
    }
    errormsg_ = "";
    return Error::NoError;
}

File::Error File::GetDoubleInterval(const std::string &category,
                                    const std::string &key, double *min,
                                    double *max) {
    std::string str;
    Error error = GetString(category, key, &str);
    if (error != Error::NoError) {
        *min = 0.0;
        *max = 0.0;
        return error;
    }
    str = Strings::Trim(str);
    if (str[0] != '(') {
        errormsg_ = Strings::Format("\"%s\"->\"%s\" doesn't start with (",
                                    category.c_str(), key.c_str());
        return Error::DoubleIntervalParser;
    }
    if (str[str.size() - 1] != ')') {
        errormsg_ = Strings::Format("\"%s\"->\"%s\" doesn't end with )",
                                    category.c_str(), key.c_str());
        return Error::DoubleIntervalParser;
    }
    auto fields = Strings::Split(str.substr(1, str.size() - 2), ",");

    if (fields.size() < 2) {
        errormsg_ = Strings::Format("missing , in \"%s\"->\"%s\"",
                                    category.c_str(), key.c_str());
        return Error::DoubleIntervalParser;
    }
    if (fields.size() > 2) {
        errormsg_ = Strings::Format("too many , in \"%s\"->\"%s\"",
                                    category.c_str(), key.c_str());
        return Error::DoubleIntervalParser;
    }
    double tmin = -1.0;
    int ret = 0;
    if (fields[0] == "-inf") {
        tmin = DBL_MIN;
    } else {
        ret = Strings::ParseDouble(Strings::Trim(fields[0]), &tmin);
        if (ret > 0) {
            errormsg_ = Strings::Format(
                "error in %s.%s: \"%s\" of \"%s\" is not a valid double.",
                category.c_str(), key.c_str(),
                fields[0].substr(ret, std::string::npos).c_str(),
                fields[0].c_str());
            return Error::DoubleIntervalParser;
        }
        if (ret == -2) {
            errormsg_ = Strings::Format("error in %s.%s: value is empty.",
                                        category.c_str(), key.c_str());
            return Error::DoubleIntervalParser;
        }
    }
    double tmax = -1.0;
    if (fields[1] == "inf") {
        tmax = DBL_MAX;
    } else {
        ret = Strings::ParseDouble(Strings::Trim(fields[1]), &tmax);
        if (ret > 0) {
            errormsg_ = Strings::Format(
                "error in %s.%s: \"%s\" of \"%s\" is not a valid double.",
                category.c_str(), key.c_str(),
                fields[1].substr(ret, std::string::npos).c_str(),
                fields[1].c_str());
            return Error::DoubleIntervalParser;
        }
        if (ret == -2) {
            errormsg_ = Strings::Format("error in %s.%s: value is empty.",
                                        category.c_str(), key.c_str());
            return Error::DoubleIntervalParser;
        }
    }

    if (tmax <= tmin) {
        errormsg_ = Strings::Format("min = %f >= %f = max in %s.%s", tmin, tmax,
                                    category.c_str(), key.c_str());
        return Error::DoubleIntervalParser;
    }
    *min = tmin;
    *max = tmax;

    errormsg_ = "";
    return Error::NoError;
}

} // end namespace config
