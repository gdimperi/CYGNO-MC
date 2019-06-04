#pragma once

// STL //
#include <string>

class CYGNODetectorConstructionMessenger;
// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

// USER //
class CYGNOSensitiveDetector;

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
    void SaveMassAndDensity();
    void UpdateGeometry();

    void SetExternalRockThickness(G4double rockthick) {rockThicknessOuter = rockthick;}
    void SetProductionRockThickness(G4double rockthick) {productionLayerThickness = rockthick;}
    void SetInternalRockThickness(G4double rockthick) {rockThicknessInner = rockthick;}

    void SetCYGNOLab(G4String lab) {CYGNOLab = lab;}
    G4String GetCYGNOLab() {return CYGNOLab;}

    void SetCYGNOShielding(G4String shield) {CYGNOShielding = shield;}
    G4String GetCYGNOShielding() {return CYGNOShielding;}
    
    void SetShieldThick0(G4double thick) {thick0 = thick;}
    void SetShieldThick1(G4double thick) {thick1 = thick;}
    void SetShieldThick2(G4double thick) {thick2 = thick;}
    void SetShieldThick3(G4double thick) {thick3 = thick;}    
    void SetShield0Material(G4String shm) {Mat0 = shm;}
    void SetShield1Material(G4String shm) {Mat1 = shm;}
    void SetShield2Material(G4String shm) {Mat2 = shm;}
    void SetShield3Material(G4String shm) {Mat3 = shm;}
    void SetInsideVolumeRadius(G4double r) {InsideVolume_OR = r;}
    void SetInsideVolumeHeight(G4double h) {InsideVolume_Z = h;}

  private:
    
    CYGNODetectorConstructionMessenger* fMessenger;

    //G4VSolid * world_solid;
    //G4LogicalVolume* world_logical;
    //G4VPhysicalVolume* world_physical;
    

    G4double InsideVolume_OR;
    G4double InsideVolume_Z;
    G4double rockThicknessOuter;
    G4double rockThicknessInner;
    G4double productionLayerThickness;
    
    G4double gasHeFrac;
    G4double gasCF4Frac;

    G4String CYGNOLab;
    G4String CYGNOShielding;
    G4double thick0;
    G4double thick1;
    G4double thick2;
    G4double thick3;
    G4String Mat0;
    G4String Mat1;
    G4String Mat2;
    G4String Mat3;

  
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
    //G4VSolid * cad_window_solid;
    G4VSolid * cad_fc_support_solid;
    G4VSolid * cad_turns_support_solid;
    G4VSolid * cad_gem_support_solid;
    G4VSolid * cad_gem_solid;
    G4VSolid * cad_cathode_frame_solid;
    G4VSolid * cad_cathode_solid;
    G4VSolid * cad_field_cage_solid;
   
    
    // Logical volumes
    G4LogicalVolume* WorldVolume_log;
    G4LogicalVolume* Shield0_log; 
    G4LogicalVolume* Shield1_log; 
    G4LogicalVolume* Shield2_log; 
    G4LogicalVolume* Shield3_log; 
    G4LogicalVolume* AirBox_log;

    G4LogicalVolume * cad_shell_logical;
    G4LogicalVolume * cad_camera_carter_logical;
    G4LogicalVolume * cad_cameras_all_logical;
    //G4LogicalVolume * cad_window_logical;
    G4LogicalVolume * CYGNO_log;
    G4LogicalVolume * cad_fc_support_logical;
    G4LogicalVolume * cad_turns_support_logical;
    G4LogicalVolume * cad_gem_support_logical;
    G4LogicalVolume * cad_gem_logical;
    G4LogicalVolume * cad_cathode_frame_logical;
    G4LogicalVolume * cad_cathode_logical;
    G4LogicalVolume * cad_field_cage_logical;
  
    // Physical volumes
    G4VPhysicalVolume* WorldVolume_phys;
    G4VPhysicalVolume* productionRockThinTube_phys;
    G4VPhysicalVolume* externalRock_phys;
    G4VPhysicalVolume* InnerAirSphere_phys;
    G4VPhysicalVolume* Shield0_phys;
    G4VPhysicalVolume* Shield1_phys;
    G4VPhysicalVolume* Shield2_phys;
    G4VPhysicalVolume* Shield3_phys;
    G4VPhysicalVolume* AirBox_phys;

    G4VPhysicalVolume * cad_shell_physical;
    G4VPhysicalVolume * cad_camera_carter_physical;
    G4VPhysicalVolume * cad_cameras_all_physical;
    //G4VPhysicalVolume * cad_window_physical;
    G4VPhysicalVolume * CYGNO_phys;
    G4VPhysicalVolume * cad_fc_support_physical;
    G4VPhysicalVolume * cad_turns_support_physical;
    G4VPhysicalVolume * cad_gem_support_physical;
    G4VPhysicalVolume * cad_gem_physical;
    G4VPhysicalVolume * cad_cathode_frame_physical;
    G4VPhysicalVolume * cad_cathode_physical;
    G4VPhysicalVolume * cad_field_cage_physical;

    //CYGNO sensitive detector
    CYGNOSensitiveDetector * CYGNOSD;



};

