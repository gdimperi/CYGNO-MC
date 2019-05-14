/* ************************************************
 * GEANT4 CAD INTERFACE - template
 *
 * File:      DetectorConstruction.hh
 *
 * Author:    Christopher M Poole,
 * Email:     mail@christopherpoole.net
 *
 * Date:      17th August, 2017
 **************************************************/

#pragma once

// STL //
#include <string>

// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <vector>


class CYGNODetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    CYGNODetectorConstruction();
    ~CYGNODetectorConstruction();

    G4VPhysicalVolume* Construct();

    void SetExternalRockThickness(G4double rockthick) {rockThicknessOuter = rockthick;}
    void SetProductionRockThickness(G4double rockthick) {productionLayerThickness = rockthick;}
    void SetInternalRockThickness(G4double rockthick) {rockThicknessInner = rockthick;}

    void SetCYGNOLab(G4String lab) {CYGNOLab = lab;}
    G4String GetCYGNOLab() {return CYGNOLab;}
  
  private:
    
    G4String CYGNOLab;

    //G4VSolid * world_solid;
    //G4LogicalVolume* world_logical;
    //G4VPhysicalVolume* world_physical;
    

    G4double InsideVolume_OR;
    G4double InsideVolume_Z;
    G4double rockThicknessOuter;
    G4double rockThicknessInner;
    G4double productionLayerThickness;
  
    //Building blocks: logic volumes, sizes and positions
    G4ThreeVector  tr_Tot;
    G4LogicalVolume* Rock_log;
    G4ThreeVector size_Rock;
    G4ThreeVector tr_Rock;
    G4RotationMatrix rot_Rock;
    G4RotationMatrix absrot_Rock;
    G4LogicalVolume* Laboratory_log;
    G4ThreeVector size_Laboratory;
    G4ThreeVector tr_Laboratory;
    G4RotationMatrix rot_Laboratory;
    G4LogicalVolume* Shielding_log;
    G4ThreeVector size_Shielding;
    G4ThreeVector tr_Shielding;
    G4RotationMatrix rot_Shielding;
    G4RotationMatrix absrot_Shielding;
    G4LogicalVolume* InsideVolume_log;
    G4ThreeVector size_InsideVolume;
    G4ThreeVector tr_InsideVolume;
    G4RotationMatrix rot_InsideVolume;
    G4ThreeVector size_cad;
    G4ThreeVector tr_cad;
    G4RotationMatrix rot_cad;
    G4RotationMatrix absrot_cad;
    G4ThreeVector tr_CYGNO_gas;
    G4RotationMatrix rot_CYGNO_gas;
    G4RotationMatrix absrot_CYGNO_gas;
    G4ThreeVector tr_cad_internal;
    
    
    //Solids and meshes
    G4VSolid * cad_shell_solid;
    G4VSolid * cad_camera_carter_solid;
    G4VSolid * cad_cameras_all_solid;
    G4VSolid * cad_window_solid;
    G4VSolid * cad_internal_structure_solid;
    G4VSolid * cad_gem_support_solid;
    G4VSolid * cad_gem_solid;
    G4VSolid * cad_cathode_frame_solid;
    G4VSolid * cad_square_turn_solid;
   
    
    // Logical volumes
    G4LogicalVolume* WorldVolume_log;
    G4LogicalVolume * cad_shell_logical;
    G4LogicalVolume * cad_camera_carter_logical;
    G4LogicalVolume * cad_cameras_all_logical;
    G4LogicalVolume * cad_window_logical;
    G4LogicalVolume * CYGNO_log;
    G4LogicalVolume * cad_internal_structure_logical;
    G4LogicalVolume * cad_gem_support_logical;
    G4LogicalVolume * cad_gem_logical;
    G4LogicalVolume * cad_cathode_frame_logical;
    G4LogicalVolume * cad_square_turn_logical;
  
    // Physical volumes
    G4VPhysicalVolume* WorldVolume_phys;
    G4VPhysicalVolume* productionRockThinTube_phys;
    G4VPhysicalVolume * cad_shell_physical;
    G4VPhysicalVolume * cad_camera_carter_physical;
    G4VPhysicalVolume * cad_cameras_all_physical;
    G4VPhysicalVolume * cad_window_physical;
    G4VPhysicalVolume * CYGNO_phys;
    G4VPhysicalVolume * cad_internal_structure_physical;
    G4VPhysicalVolume * cad_gem_support_physical;
    G4VPhysicalVolume * cad_gem_physical;
    G4VPhysicalVolume * cad_cathode_frame_physical;
    G4VPhysicalVolume * cad_square_turn_physical;


};

