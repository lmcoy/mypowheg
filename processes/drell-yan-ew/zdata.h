#ifndef ZDATA_H_G8GSJ8CQ
#define ZDATA_H_G8GSJ8CQ

#include "fks/process.h"
#include "process/data.h"

class ZData : public UserProcess::Data {
public:
  virtual ~ZData();
protected:
  virtual int ProcessInit(Config::File &);
  virtual void ProcessPrint() const;
};

#endif
