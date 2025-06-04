#include "CYGNOHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "TMath.h"  

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


static Double_t f_helium(Double_t E);
static Double_t fprime_helium(Double_t E);
static Double_t QF_helium(Double_t E);
static Double_t QFprime_helium(Double_t E);
static Double_t dRdE_helium(Double_t E);

static Double_t f_carbon(Double_t E);
static Double_t fprime_carbon(Double_t E);
static Double_t QF_carbon(Double_t E);
static Double_t QFprime_carbon(Double_t E);
static Double_t dRdE_carbon(Double_t E);

static Double_t f_flor(Double_t E);
static Double_t fprime_flor(Double_t E);
static Double_t QF_flor(Double_t E);
static Double_t QFprime_flor(Double_t E);
static Double_t dRdE_flor(Double_t E);

// New function to apply the quenching factor (derivative of E*QF)
// when the hit corresponds to a nucleus.
void CYGNOHit::ApplyQuenching()
{
  // Only apply quenching if the hit belongs to a nucleus.
  // (Here we check the particleID for helium, carbon, or fluorine.)
  double dRdE = 1.0;
  if (particleID == 1000020040) {         // Helium ion (example PDG code)
    dRdE = dRdE_helium(kinEne);
  } else if (particleID == 1000060120) {  // Carbon ion (example PDG code)
    dRdE = dRdE_carbon(kinEne);
  } else if (particleID == 1000090190) {  // Fluorine ion (example PDG code)
    dRdE = dRdE_flor(kinEne);
  }
  edep *= dRdE;
}
void CYGNOHit::ApplyQuenchingAvg()
{
  // Only apply quenching if the hit belongs to a nucleus.
  // (Here we check the particleID for helium, carbon, or fluorine.)
  double QFavg = 1.0;
  if (particleID == 1000020040) {         // Helium ion (example PDG code)
    QFavg = QF_helium(kinEne);
  } else if (particleID == 1000060120) {  // Carbon ion (example PDG code)
    QFavg = QF_carbon(kinEne);
  } else if (particleID == 1000090190) {  // Fluorine ion (example PDG code)
    QFavg = QF_flor(kinEne);
  }
  edep *= QFavg;
}

// --- QF functions for different ions ---

// Helium from F. Di Giambattista PhD thesis
static const double k_helium = 0.117;
static const double a_helium = 3.9;
static const double b_helium = 0.44;

Double_t f_helium(Double_t E) {
    return k_helium * (E + a_helium * pow(E, b_helium));
}

Double_t fprime_helium(Double_t E) {
    return k_helium * (1.0 + a_helium * b_helium * pow(E, b_helium - 1.0));
}

Double_t QF_helium(Double_t E) {
    Double_t val = f_helium(E);
    return val / (1.0 + val);
}

Double_t QFprime_helium(Double_t E) {
    Double_t fp = fprime_helium(E);
    Double_t fE = f_helium(E);
    return fp / pow(1.0 + fE, 2.0);
}

Double_t dRdE_helium(Double_t E) {
    return QF_helium(E) + E * QFprime_helium(E);
}

// Carbon from F. Di Giambattista PhD thesis
static const double k_carbon = 0.0195;
static const double a_carbon = 14.7;
static const double b_carbon = 0.33;

Double_t f_carbon(Double_t E) {
    return k_carbon * (E + a_carbon * pow(E, b_carbon));
}

Double_t fprime_carbon(Double_t E) {
    return k_carbon * (1.0 + a_carbon * b_carbon * pow(E, b_carbon - 1.0));
}

Double_t QF_carbon(Double_t E) {
    Double_t val = f_carbon(E);
    return val / (1.0 + val);
}

Double_t QFprime_carbon(Double_t E) {
    Double_t fp = fprime_carbon(E);
    Double_t fE = f_carbon(E);
    return fp / pow(1.0 + fE, 2.0);
}

Double_t dRdE_carbon(Double_t E) {
    return QF_carbon(E) + E * QFprime_carbon(E);
}

// Fluorine from F. Di Giambattista PhD thesis
static const double k_flor   = 0.0083;
static const double a_flor   = 27.4;
static const double b_flor   = 0.303;

Double_t f_flor(Double_t E) {
    return k_flor * (E + a_flor * pow(E, b_flor));
}

Double_t fprime_flor(Double_t E) {
    return k_flor * (1.0 + a_flor * b_flor * pow(E, b_flor - 1.0));
}

Double_t QF_flor(Double_t E) {
    Double_t val = f_flor(E);
    return val / (1.0 + val);
}

Double_t QFprime_flor(Double_t E) {
    Double_t fp = fprime_flor(E);
    Double_t fE = f_flor(E);
    return fp / pow(1.0 + fE, 2.0);
}

Double_t dRdE_flor(Double_t E) {
    return QF_flor(E) + E * QFprime_flor(E);
}

