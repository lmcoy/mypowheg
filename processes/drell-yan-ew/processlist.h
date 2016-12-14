#ifndef PROCESSLIST_H_
#define PROCESSLIST_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "config/file.h"

struct Process {
    Process(int process, const std::string &type) : Id(process), Type(type) {}
    int Id;
    std::string Type;
};

typedef std::vector<Process> ProcessList;

int ReadProcessFile(ProcessList *list, const char *filename) {
    Config::File cfile;
    auto cerror = cfile.ReadFromFile(filename);
    if (cerror != Config::File::ReadError::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    Strings::StringList slist;
    std::map<std::string, int> process_list;
    process_list["u~u>mu+mu-"] = 0;
    process_list["uu~>mu+mu-"] = 1;
    process_list["d~d>mu+mu-"] = 2;
    process_list["dd~>mu+mu-"] = 3;

    for (auto &p : process_list) {
        if (cfile.GetStringList("processes", p.first, &slist) !=
            Config::File::Error::NoError) {
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (slist.size() == 1 && slist[0].empty() ) {
            continue;
        }
        for (auto &s : slist) {
            if (s != "b" && s != "r" && s != "v") {
                std::cerr << "unknown type \"" << s
                          << "\": {b, v, r} expected\n";
                return 1;
            }
            list->push_back(Process(p.second, s));
        }
    }
    return 0;
}

#endif

