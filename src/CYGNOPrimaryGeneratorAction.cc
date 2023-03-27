#include "CYGNOPrimaryGeneratorAction.hh"
#include "CYGNODetectorConstruction.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <math.h>
#include "G4GeneralParticleSource.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"
#include "G4Gamma.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#define PI 3.14159265
using namespace std;

CYGNOPrimaryGeneratorAction::CYGNOPrimaryGeneratorAction(CYGNODetectorConstruction* myDC):myDetector(myDC),fGenerator(""),fMessenger(nullptr)
{
  n_particle = 1;

  MuonParticleGun = new G4ParticleGun(n_particle);
  ER_Gun = new G4ParticleGun(n_particle);
  particleGun = new G4GeneralParticleSource();

  //fMessenger = new CYGNOPrimaryGeneratorActionMessenger(this);

  // default particle
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  particleGun->SetParticleDefinition(particle);
  particle_energy = 1.0*MeV;

  DefineCommands();
  fFileName = " ";

  }

CYGNOPrimaryGeneratorAction::~CYGNOPrimaryGeneratorAction()
{
  //delete fMessenger;
  delete particleGun;
  delete MuonParticleGun;
  delete ER_Gun;
  delete fMessenger;
  if(fInputFile.is_open())
    fInputFile.close();
}

void CYGNOPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 

//----------MUSUN---------------
/*  if(fGenerator == "Musun")
  {
    G4int    nEvent = 0;
    G4double time   = 0.0;
    G4double energy = 0.0 * MeV;
    G4double px, py, pz;
    G4double theta, phi;
    G4double x = 0, y = 0, z = 0;
    G4int    particleID = 0;

    //fInputFile >> nEvent >> particleID >> energy >> x >> y >> z >> theta >> phi;
    fInputFile >> nEvent >> particleID >> energy >> x >> y >> z >> px >> py >> pz;

    // G4cout  << nEvent << " " << x << " " << y << " " << z << G4endl;
    if(fInputFile.eof())
    {
      fInputFile.close();
      G4cerr << "File over: not enough events! Debugoutput" << G4endl;
      G4Exception("CYGNOPrimaryGeneratorAction::GeneratePrimaryVertex()", "err001",
                  FatalException, "Exit");
      return;
    }

    G4double particle_time = time * s;
    energy                 = energy * GeV;
    theta                  = theta * rad;
    phi                    = phi * rad;
    x                      = x * cm;
    y                      = y * cm;
    z                      = z * cm;

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();

    G4String particleName = " ";

    if(particleID == 10)
      particleName = "mu+";
    else
      particleName = "mu-";

    G4ParticleDefinition* muon = theParticleTable->FindParticle(particleName);
    G4double theMass     = muon->GetPDGMass();
    G4double totMomentum = std::sqrt(energy * energy + 2 * theMass * energy);
    //pz                   = -1 * std::cos(theta);
    //px                   = std::sin(theta) * cos(phi);
    //py                   = std::sin(theta) * sin(phi);
    G4ThreeVector momentumDir(px, py, pz);

    MuonParticleGun->SetParticleDefinition(muon);

    MuonParticleGun->SetParticleMomentumDirection(momentumDir);

    MuonParticleGun->SetParticleEnergy(energy);

    MuonParticleGun->SetParticlePosition(G4ThreeVector(y, z, x));

    MuonParticleGun->GeneratePrimaryVertex(anEvent);
    
       G4cout << "Primary coordinates: " << position/m << " m" << G4endl;
       G4cout << "Primary coordinates: " << x/cm << " " <<  y/cm << " " << z/cm << " "
       << G4endl; G4cout << "Primary energy: " << energy/GeV << " GeV" << G4endl; G4cout
       << "Theta: " << theta/deg << " deg; Phi: " << phi/deg << " deg" << G4endl;
  }*/
  /*else if (fGenerator=='txt_list')
  {
    G4int    nEvent = 0;
    G4double time   = 0.0;
    G4double energy = 0.0 * MeV;
    G4double px, py, pz;
    G4double theta, phi;
    G4double x = 0, y = 0, z = 0;
    G4int    particleID = 0;

    fInputFile >> nEvent >> particleID >> energy >> x >> y >> z >> theta >> phi;
    //fInputFile >> nEvent >> particleID >> energy >> x >> y >> z >> px >> py >> pz;

     G4cout  << nEvent << " " << energy << " " << x << " " << y << " " << z << G4endl;
    if(fInputFile.eof())
    {
      fInputFile.close();
      G4cerr << "File over: not enough events! Debugoutput" << G4endl;
      G4Exception("CYGNOPrimaryGeneratorAction::GeneratePrimaryVertex()", "err001",
                  FatalException, "Exit");
      return;
    }


    if (particleID>1e9) return;
    G4double particle_time = time * s;
    energy                 = energy * keV;
    theta                  = theta * rad;
    phi                    = phi * rad;
    x                      = x * mm;
    y                      = y * mm;
    z                      = z * mm;

    //   G4cout << "Primary coordinates: " << position/m << " m" << G4endl;
    //   G4cout << "Primary coordinates: " << x/cm << " " <<  y/cm << " " << z/cm << " "
    //   << G4endl; G4cout << "Primary energy: " << energy/GeV << " GeV" << G4endl; G4cout
    //   << "Theta: " << theta/deg << " deg; Phi: " << phi/deg << " deg" << G4endl;

    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();

    G4ParticleDefinition* recoil = theParticleTable->FindParticle(particleID);
    G4cout<<recoil->GetParticleName()<<G4endl;
    G4double theMass     = recoil->GetPDGMass();
    //G4double totMomentum = std::sqrt(energy * energy + 2 * theMass * energy);
    pz                   = -1 * std::cos(theta);
    px                   = std::sin(theta) * cos(phi);
    py                   = std::sin(theta) * sin(phi);
    G4ThreeVector momentumDir(px, py, pz);

    ER_Gun->SetParticleDefinition(recoil);

    ER_Gun->SetParticleMomentumDirection(momentumDir);

    ER_Gun->SetParticleEnergy(energy);

    ER_Gun->SetParticlePosition(G4ThreeVector(x, y, z));

    ER_Gun->GeneratePrimaryVertex(anEvent);
	  
  }*/
  //else{
  G4int numParticles=1;
  particleGun->SetNumberOfParticles(numParticles);
  particleGun->GeneratePrimaryVertex(anEvent);
  //}
  
}

void CYGNOPrimaryGeneratorAction::DefineCommands()
{
  // Define /CYGNO/generator command directory using generic messenger class
  fMessenger =
    new G4GenericMessenger(this, "/CYGNO/generator/", "Primary generator control");

  // musun file command
  auto& musunfileCmd =
    fMessenger
      ->DeclareMethod("setInputFile", &CYGNOPrimaryGeneratorAction::ChangeFileName)
      .SetGuidance("Set input file name")
      .SetParameterName("filename", false)
      .SetDefaultValue("./musun_gs_100M.dat");
  // generator command
  // switch command
  fMessenger->DeclareMethod("setGenerator", &CYGNOPrimaryGeneratorAction::SetGenerator)
    .SetGuidance("Set generator of primary muons from MUSUN or from a previous GEANT4 simulation")
    .SetDefaultValue(" ")
    .SetCandidates(
      "Musun txt_list");
}

void CYGNOPrimaryGeneratorAction::ChangeFileName(G4String newFile)
{
  if(fFileName != newFile)  // check if the new file is equal to the other
  {
    if(fInputFile.is_open())
      fInputFile.close();  // close the old file
    G4cout << "opening file: " << newFile << G4endl;
    fFileName = newFile;
    // open the new one
    fInputFile.open(fFileName, std::ifstream::in);
    if(!(fInputFile.is_open()))
    {  
       G4cerr << "file not valid! Name: " << fFileName << G4endl;
    }
  }
}

void CYGNOPrimaryGeneratorAction::SetGenerator(const G4String& name)
{
  std::set<G4String> knownGenerators = {
    "Musun", "txt_list"
  };
  if(knownGenerators.count(name) == 0)
  {
    G4Exception("CYGNOPrimaryGeneratorAction::SetGenerator", "", JustWarning,
                ("Invalid generator name '" + name + "'").c_str());
    return;
  }
  fGenerator = name;
}

