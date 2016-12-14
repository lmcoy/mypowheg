#ifndef WDATA_H_SF30FW39
#define WDATA_H_SF30FW39

#include "fks/process.h"
#include "process/data.h"

class WData : public UserProcess::Data {
  public:
    virtual ~WData();

  protected:
    virtual int ProcessInit(Config::File &);
    virtual void ProcessPrint() const;
};

#endif
