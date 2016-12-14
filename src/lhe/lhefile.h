#ifndef LHEFILE_H_LRQ7JVWO
#define LHEFILE_H_LRQ7JVWO

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

namespace LHE {

class Process {
  public:
    int ID = 0;
    double XSec = 0.0;
    double XSecErr = 0.0;
    double XMax = 0.0;
};

class File {
  public:
    ~File() {
        if (stage_ > 0) {
            Close();
        }
    }
    int Init(const char *filename);

    int WriteHeader(std::istream &input, bool xml);

    int WriteProcessInfo(double SqrtS, int lhaid,
                         const std::vector<LHE::Process> &plist);


    int WriteEvents(size_t len, const char* buffer);
    int Close();
    std::string Filename() const { return filename_; }

  private:
    std::ofstream output_;
    int stage_ = 0;
    std::string filename_;
};

} // end namespace LHE

#endif /* end of include guard: LHEFILE_H_LRQ7JVWO */
