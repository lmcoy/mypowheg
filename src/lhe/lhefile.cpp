#include "lhe/lhefile.h"
#include "libconfig.h"

namespace LHE {

int File::Init(const char *filename) {
    if (stage_ != 0) {
        return 2;
    }
    output_.open(filename);
    if (!output_) {
        return 1;
    }

    output_ << "<LesHouchesEvents version=\"1.0\">\n";
    stage_ = 1;
    filename_ = std::string(filename);
    return 0;
}

int File::WriteHeader(std::istream &input, bool xml) {
    if (stage_ != 1) {
        return 2;
    }
    std::string start_tag;
    std::string end_tag;
    std::string git_start;
    std::string git_end;
    if (xml) {
        start_tag = "<header>";
        end_tag = "</header>";
        git_start = "<git>";
        git_end = "</git>";
    } else {
        start_tag = "<!--";
        end_tag = "-->";
        git_start = "Git:";
        git_end = "Git End";
    }
    output_ << start_tag << "\n";
    for (std::string line; std::getline(input, line);) {
        output_ << line << "\n";
    }
    output_ << git_start << "\n";
    output_ << GITVERSION << "\n";
    output_ << GITDIFF << "\n";
    output_ << git_end << "\n";
    output_ << end_tag << "\n";
    return 0;
}

int File::WriteProcessInfo(double SqrtS, int lhaid,
                           const std::vector<LHE::Process> &plist) {
    if (stage_ != 1) {
        return 2;
    }
    double beamE = SqrtS / 2.0;
    int N = plist.size();
    output_ << "<init>\n";
    output_ << " 2212 2212 "; // proton proton collider
    output_ << beamE << " " << beamE << " ";
    output_ << "0 0 "; // pdf (lhapdf)
    output_ << lhaid << " " << lhaid << " ";
    output_ << "3 "; // weighting: use every event
    output_ << N << "\n";

    for (const auto &process : plist) {
        output_ << " " << process.XSec << " " << process.XSecErr << " "
                << process.XMax << " " << process.ID << "\n";
    }
    output_ << "</init>\n";

    stage_ = 2;
    return 0;
}

int File::WriteEvents(size_t len, const char *buffer) {
    if (stage_ != 2) {
        return 2;
    }
    output_.write(buffer, len);
    return 0;
}

int File::Close() {
    if (stage_ < 1) {
        return 2;
    }
    output_ << "</LesHouchesEvents>\n";
    output_.close();
    stage_ = 0;
    return 0;
}

} /* LHE */
