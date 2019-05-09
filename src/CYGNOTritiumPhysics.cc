//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "CYGNOTritiumPhysics.hh"

#include "G4Triton.hh"
//#include "G4DecayPhysics.hh"//Already included in any physics list
//#include "G4RadioactiveDecayPhysics.hh"//Already included in Shielding
#include "G4RadioactiveDecay.hh"

#include "G4ProcessManager.hh"
//#include "G4LossTableManager.hh"
//#include "G4ParticleDefinition.hh"

// factory
//#include "G4PhysicsConstructorFactory.hh"
//G4_DECLARE_PHYSCONSTR_FACTORY(CYGNOTritiumPhysics);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOTritiumPhysics::CYGNOTritiumPhysics() 
  : G4VPhysicsConstructor()
	//  : G4VUserPhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOTritiumPhysics::~CYGNOTritiumPhysics()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void CYGNOTritiumPhysics::ConstructParticle()
{
  G4Triton::TritonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CYGNOTritiumPhysics::ConstructProcess()
{
  // Define transportation process
  //AddTransportation();
  // Decay process
  ConstructDecay();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CYGNOTritiumPhysics::ConstructDecay()
{
  // Make tritium unstable, so it can be decayed
  G4ParticleDefinition* tr = G4Triton::TritonDefinition();
  tr->SetPDGStable(false);
  // Remove G4Decay process, which requires a registered decay table
  G4VProcess* decay = 0;
  G4ProcessManager* pman = tr->GetProcessManager();
  G4ProcessVector* pvec = pman->GetAtRestProcessVector();
  for (G4int i=0; i<pvec->size() && decay==0; i++) {
    if ((*pvec)[i]->GetProcessName() == "Decay") decay = (*pvec)[i];
  }
  if (decay) pman->RemoveProcess(decay);
  // Attach RDM, which is a rest-discrete process
  tr->GetProcessManager()->AddProcess(new G4RadioactiveDecay(), 1000, -1, 1000);//the -1 means that nothing is done along the step, something is done at rest and after step instead
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void CYGNOTritiumPhysics::SetCuts()
{
  //   the G4VUserPhysicsList::SetCutsWithDefault() method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
