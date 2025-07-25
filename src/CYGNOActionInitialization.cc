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

#include "CYGNOActionInitialization.hh"
#include "CYGNOPrimaryGeneratorAction.hh"
#include "CYGNORunAction.hh"
#include "CYGNOEventAction.hh"
//#include "CYGNOSteppingAction.hh"
//#include "CYGNOStackingAction.hh"
//#include "CYGNOTrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOActionInitialization::CYGNOActionInitialization(CYGNODetectorConstruction* det)
  : G4VUserActionInitialization(),fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CYGNOActionInitialization::~CYGNOActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CYGNOActionInitialization::BuildForMaster() const
{
//  SetUserAction(new CYGNORunAction(fDetector));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CYGNOActionInitialization::Build() const
{
  CYGNOPrimaryGeneratorAction* gen_action = new CYGNOPrimaryGeneratorAction(fDetector);
  SetUserAction(gen_action);

  CYGNORunAction* run_action = new CYGNORunAction(fDetector);
  SetUserAction(run_action);
  
  CYGNOEventAction* event_action = new CYGNOEventAction(run_action,fDetector);
  SetUserAction(event_action);
  
  //G4UserSteppingAction* stepping_action = new CYGNOSteppingAction(fDetector,event_action);
  //SetUserAction(stepping_action);

  //SetUserAction(new CYGNOStackingAction);
  //SetUserAction(new CYGNOTrackingAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
