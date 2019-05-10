//SF6 DETECTOR HEADER FILE

#ifndef CYGNODetectorConstruction_H
#define CYGNODetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class AllSD;

class CYGNODetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "G4SDManager.hh" 

class CYGNODetectorConstruction : public G4VUserDetectorConstruction {

 public:
  
   CYGNODetectorConstruction();
  ~CYGNODetectorConstruction();

  G4VPhysicalVolume* Construct();
   
private:

  G4LogicalVolume*   world_log;
  G4VPhysicalVolume* world_phys;
  
  G4LogicalVolume* vessel_top_log;
  G4VPhysicalVolume* vessel_top_phys;
 
  G4LogicalVolume* vessel_bottom_log;
  G4VPhysicalVolume* vessel_bottom_phys;
  
  G4LogicalVolume* vessel_front_log;
  G4VPhysicalVolume* vessel_front_phys;
  
  G4LogicalVolume* vessel_back_log;
  G4VPhysicalVolume* vessel_back_phys;
  
  G4LogicalVolume* glass_left_log;
  G4VPhysicalVolume* glass_left_phys;

  G4LogicalVolume* glass_right_log;
  G4VPhysicalVolume* glass_right_phys;

  G4LogicalVolume* GAS_log;
  G4VPhysicalVolume* GAS_phys;  
  
  G4LogicalVolume* camera_shield_bottom_log;
  G4VPhysicalVolume* camera_shield_bottom_phys;

  G4LogicalVolume* camera_shield_left_log;
  G4VPhysicalVolume* camera_shield_left_phys;

  G4LogicalVolume* camera_shield_right_log;
  G4VPhysicalVolume* camera_shield_right_phys;

  G4LogicalVolume* camera_shield_front_log;
  G4VPhysicalVolume* camera_shield_front_phys;

  G4LogicalVolume* camera_shield_back_log;
  G4VPhysicalVolume* camera_shield_back_phys;
  
  G4LogicalVolume* camera_shield_log;
  G4VPhysicalVolume* camera_shield_phys;

  G4LogicalVolume* camera_box_log;
  G4VPhysicalVolume* camera_box_phys;

  G4LogicalVolume* GEM_Kapton_log;
  G4VPhysicalVolume* GEM_Kapton_phys;

  G4LogicalVolume* GEM_Copper_log;
  G4VPhysicalVolume* GEM_Copper_phys;

  G4LogicalVolume* shield_left_log;
  G4VPhysicalVolume* shield_left_phys;

  G4LogicalVolume* camera_hole_log;
  G4VPhysicalVolume* camera_hole_phys;

 

  // G4SubtractionSolid* GEM_Copper;
  //G4UnionSolid* GEM;
 
#define world_length (3*m) // Boulby : add 25 cm on each side
#define world_width (3*m)
#define world_height (3*m)
   
#define vessel_length (GAS_length+4*cm) // 2cm either side
#define vessel_width (GAS_width+4*cm)
#define vessel_height (GAS_height+4*cm)

#define GAS_length (1*m) 
#define GAS_width (1*m)
#define GAS_height (1*m)  

#define camera_box_length (100*mm) 
#define camera_box_width (100*mm) 
#define camera_box_height (150*mm) 

#define camera_shield_length (140*mm)
#define camera_shield_width (140*mm)
#define camera_shield_height (170*mm)
 
  // Sensitive Detector
  
  G4SDManager* SDman; 

  AllSD* allSD;
 
  public  :
  G4int NbOfDrift;
  G4int GetNbOfDrift() {return NbOfDrift;};
  G4VPhysicalVolume* GetDriftPhys(G4int ){return GAS_phys;}
  
};

#endif
