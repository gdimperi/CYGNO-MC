#include "CYGNOHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4ThreadLocal G4Allocator<CYGNOHit>* CYGNOHitAllocator=0;

CYGNOHit::CYGNOHit() {}

CYGNOHit::~CYGNOHit() {}

CYGNOHit::CYGNOHit(const CYGNOHit& right)
  : G4VHit()
{
  parentID   = right.parentID;
  particleID   = right.particleID;
  trackID   = right.trackID;
  globalTime = right.globalTime;
  kinEne = right.kinEne;
  processIni = right.processIni;
  processFin = right.processFin;
  edep      = right.edep;
  pos       = right.pos;
  trackLen       = right.trackLen;
}

const CYGNOHit& CYGNOHit::operator=(const CYGNOHit& right)
{
  parentID   = right.parentID;
  particleID  = right.particleID;
  trackID   = right.trackID;
  globalTime = right.globalTime;
  kinEne = right.kinEne;
  processIni = right.processIni;
  processFin = right.processFin;
  edep      = right.edep;
  pos       = right.pos;
  trackLen       = right.trackLen;

  return *this;
}

G4int CYGNOHit::operator==(const CYGNOHit& right) const
{
  return (this==&right) ? 1 : 0;
}

void CYGNOHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void CYGNOHit::Print()
{
//  G4cout << "  parent id : " << parentID
//	 << "  particleID : " << particleID
//	 << "  track id : " << trackID
// 	 << "  global time : " << G4BestUnit(globalTime,"Time")
//	 << "  kinetic energy : " << G4BestUnit(kinEne,"Energy")
//	 << "  process ini : " << processIni 
//	 << "  process fin : " << processFin 
//     << "  detector number : " << detN
//	 << "  energy deposit : " << G4BestUnit(edep,"Energy")
//	 << "  position : " << G4BestUnit(pos,"Length") << G4endl;
}

