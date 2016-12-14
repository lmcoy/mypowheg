#ifndef PROCESS_GENBORNPS_H_
#define PROCESS_GENBORNPS_H_

namespace Phasespace { class Phasespace; }
namespace UserProcess { struct Data; }

void GenBornPhasespace(Phasespace::Phasespace *ps_out, const int ndim,
                       const double x[], const UserProcess::Data *userdata);

#endif

