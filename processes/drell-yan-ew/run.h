#ifndef RUN_H_Z8FLX5GK
#define RUN_H_Z8FLX5GK

#include <array>
#include <vector>
#include <cassert>

struct Run {
    int Number = 0;
    std::vector<std::array<double, FKS::RadiationRegion::NPDF> > Norms;
    double Maximum = 0.0;
    double MaximumVoverB = 0.0;
    double XSec = 0.0;
    double XSecErr = 0.0;
    bool Write(const char *filename) {
        std::ofstream out(filename);
        if (!out) {
            std::cerr << "could not open " << filename << "\n";
            return false;
        }
        out << Number << "\n";
        out << Norms.size() << "\n";
        for (const auto &n : Norms) {
            for (size_t i = 0; i < n.size(); i++) {
                out << " " << n[i];
                if (i != n.size() - 1) {
                    out << ",";
                }
            }
            out << "\n";
        }
        out << Maximum << "\n";
        out << MaximumVoverB << "\n";
        out << XSec << "\n";
        out << XSecErr << "\n";

        return true;
    }

    bool Read(const char *filename) {
        std::ifstream in(filename);
        if (!in) {
            std::cerr << "could not open " << filename << "\n";
            return false;
        }
        std::string line;
        std::getline(in, line);
        long N = 0;
        Strings::ParseInt(line, &N);
        Number = N;
        std::getline(in, line);
        Strings::ParseInt(line, &N);
        std::vector<double> out;
        for (long i = 0; i < N; i++) {
            out.clear();
            std::array<double, FKS::RadiationRegion::NPDF> array;
            std::getline(in, line);
            size_t end = 0;
            Strings::ParseDoubleList(line, &out, &end);
            assert(array.size() == out.size());
            for (size_t j = 0; j < array.size(); j++) {
                array[j] = out[j];
            }
            Norms.push_back(array);
        }
        std::getline(in, line);
        Strings::ParseDouble(line, &Maximum);
        std::getline(in, line);
        Strings::ParseDouble(line, &MaximumVoverB);
        std::getline(in, line);
        Strings::ParseDouble(line, &XSec);
        std::getline(in, line);
        Strings::ParseDouble(line, &XSecErr);
        return true;
    }
};


#endif /* end of include guard: RUN_H_Z8FLX5GK */ 
