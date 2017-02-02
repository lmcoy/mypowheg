#ifndef WJDATA_H_AHSJC7SM
#define WJDATA_H_AHSJC7SM

#include "fks/process.h"
#include "process/data.h"

class WJData : public UserProcess::Data {
  public:
    virtual ~WJData();

  protected:
    virtual int ProcessInit(Config::File &);
    virtual void ProcessPrint() const;
};

#endif