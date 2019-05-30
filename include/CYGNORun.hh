#ifndef CYGNORun_h
#define CYGNORun_h 1

#include "G4Run.hh"

class G4Event;

class CYGNORun : public G4Run 
{
public:
  CYGNORun();
  virtual ~CYGNORun() {};
  
  virtual void RecordEvent(const G4Event*);
  virtual void Merge(const G4Run*);
private:
  G4int CYGNOID;
};
#endif
