#include "CYGNODetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4VVisManager.hh"
#include "G4PVParameterised.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "G4UserSpecialCuts.hh"
#include "G4ios.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PVReplica.hh"

// Includes Physical Constants and System of Units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

CYGNODetectorConstruction::CYGNODetectorConstruction() :
   world_log(0),     world_phys(0)

{
  SDman = G4SDManager::GetSDMpointer();
  allSD = NULL;  
}

CYGNODetectorConstruction::~CYGNODetectorConstruction() {}

G4VPhysicalVolume* CYGNODetectorConstruction::Construct() {

 
////************************************************************************
//// MATERIALS: Data taken from http://www.webelements.com/

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density,shieldDensity, solid_shieldDensity, cardboardDensity, woodDensity, pressure, temperature, fractionmass;
  G4double He_frac, CF4_frac;
  G4String name, symbol;
  G4int ncomponents, natoms;

//// Define all elements used

  a = 10.811*g/mole;
  G4Element* elB = new G4Element(name="Boron", symbol="B", z=5, a);
  
  a = 1.00794*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1, a);

  a = 4.002602*g/mole;
  G4Element* elHe = new G4Element(name="Helium", symbol="He", z=2, a);
  
  a = 12.0107*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6, a);

  a = 14.0*g/mole;
  G4Element* elC14 = new G4Element(name="Carbon14", symbol="C14", z=6, a);

  a = 14.0067*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7, a);

  a = 15.9994*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8, a);

  a = 22.9897*g/mole;
  G4Element* elNa = new G4Element(name="Sodium", symbol="Na", z=11, a);

  a = 28.0855*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14, a);

  a = 30.9738*g/mole;
  G4Element* elP = new G4Element(name="Phosphorus", symbol="P", z=15, a);

  a = 32.066*g/mole;
  G4Element* elS = new G4Element(name="Sulphur", symbol="S", z=16, a);

  a = 35.453*g/mole;
  G4Element* elCl = new G4Element(name="Chlorine", symbol="Cl", z=17, a);

  a = 51.9961*g/mole;
  G4Element* elCr = new G4Element(name="Chromium", symbol="Cr", z=24, a);

  a = 54.938*g/mole;
  G4Element* elMn = new G4Element(name="Manganese", symbol="Mn", z=25, a);

  a = 55.845*g/mole;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", z=26, a);

  a = 58.6934*g/mole;
  G4Element* elNi = new G4Element(name="Nickel", symbol="Ni", z=28, a);

  a = 92.9064*g/mole;
  G4Element* elNb = new G4Element(name="Niobium", symbol="Nb", z=41, a);

  a = 207.2*g/mole;
  G4Element* elPb = new G4Element(name="Lead", symbol="Pb", z=82, a);

  a = 63.546*g/mole;
  G4Element* elCu = new G4Element(name="Copper", symbol="Cu", z=29, a);

  a=40.078*g/mole;
  G4Element* elCa = new G4Element(name="Calcium",symbol="Ca" , z= 20., a);
  
  a=26.98*g/mole;
  G4Element* elAl = new G4Element(name="Aluminium"  ,symbol="Al" , z= 13., a);  
	
  a=18.998*g/mole;
  G4Element* elF = new G4Element(name="Flourine"  ,symbol="F" , z= 9., a);

  a= 47.867*g/mole;
  G4Element*   elTi = new G4Element(name="Titanium", symbol="Ti", z=22.,a);

  //// Temperature of DRIFT lab
  temperature = 300*kelvin;

  //Boubly Rock

  // Create rock salt from sodium and chlorine
  density = 2170.*kg/m3;
  G4Material* RockSalt = new G4Material(name="RockSalt", density, ncomponents=2);
  RockSalt->AddElement(elNa, natoms=1);
  RockSalt->AddElement(elCl, natoms=1);
  
  // Create air from 80% nitrogen and 20% oxygen
  density = 1.290*mg/cm3;         // density of air
  pressure = 1.0*atmosphere;      // pressure of air in JIF
  G4Material* Air = new G4Material(name="Air", density, ncomponents=2,
				   kStateGas, temperature, pressure);
  Air->AddElement(elN, fractionmass=80*perCent);
  Air->AddElement(elO, fractionmass=20*perCent);

  //Borated HDPE. 95% HDPE (High density Polyethylene ...), 5% Boron ->Used in medical physics to shield against neutrons and gammas.
  density = 0.95*g/cm3; //Density from http://www.emcoplastics.com/borated-hdpe/
  G4Material* HDPE = new G4Material(name="HDPE", density, ncomponents=4);
  HDPE->AddElement(elC, fractionmass=81.4*perCent);
  HDPE->AddElement(elH, fractionmass=13.6*perCent);
  HDPE->AddElement(elB, fractionmass=5.0*perCent);
    

  // Gran sasso limestone CaCO3, limestone 2.7 g/cm3
  // https://arxiv.org/pdf/1509.00770.pdf


  density = 2.7*g/cm3;
  G4Material* RockLime = new G4Material(name="RockLime", density, ncomponents=3);
  RockLime->AddElement(elCa, natoms=1);
  RockLime->AddElement(elC, natoms=1);
  RockLime->AddElement(elO, natoms=3); 

 
  // Create stainless steel
  density = 8.027*g/cm3;
  G4Material* SS = new G4Material(name="Stainless_steel", density,
				  ncomponents=9);
  SS->AddElement(elFe, fractionmass=67.445*perCent);
  SS->AddElement(elCr, fractionmass=18*perCent);
  SS->AddElement(elNi, fractionmass=11*perCent);
  SS->AddElement(elMn, fractionmass=2*perCent);
  SS->AddElement(elSi, fractionmass=1*perCent);
  SS->AddElement(elNb, fractionmass=0.4*perCent);
  SS->AddElement(elS, fractionmass=0.08*perCent);
  SS->AddElement(elP, fractionmass=0.045*perCent);
  SS->AddElement(elS, fractionmass= 0.03*perCent);


  // Create Copper

  density = 8.96*g/cm3; //density updated from wikipedia
  G4Material* Copper = new G4Material(name="Copper", density, ncomponents=1);
  Copper->AddElement(elCu, natoms=1);  

  density = 2700*kg/m3;
  G4Material* Aluminium = new G4Material(name="Aluminium", density, ncomponents=1);
  Aluminium->AddElement(elAl, natoms=1);
 
  // Create lead for layer beneath detector
  density = 11340*kg/m3;
  G4Material* Lead = new G4Material(name="Lead", density, ncomponents=1);
  Lead->AddElement(elPb, natoms=1);

  // Create wax for layer beneath detector
  density = 0.8088*g/cm3;
  G4Material* Wax = new G4Material(name="Wax", density, ncomponents=2);
  Wax->AddElement(elH, natoms=54);
  Wax->AddElement(elC, natoms=26);

 //water
  G4Material* water = new G4Material(name="water", density=1.00*g/cm3, ncomponents=2);
  water->AddElement(elH , 2);
  water->AddElement(elO , 1);

 // Copper
 // density = 8.92*g/cm3;
 // G4Material* copper = new G4Material(name="copper", density, ncomponents=1);
 // copper->AddElement(elCu, natoms=1);

  // making quartz
  G4Material* quartz = new G4Material
  (name="quartz", density=2.200*g/cm3, ncomponents=2);
  quartz->AddElement(elSi, 1);
  quartz->AddElement(elO , 2);

 // Perspex (Acrylic)
  density = 1.180*g/cm3; //1180kg/m^3
  G4Material* Perspex = new G4Material(name="Perspex", density, ncomponents=3);
  Perspex->AddElement(elH, natoms=8);
  Perspex->AddElement(elC, natoms=5);
  Perspex->AddElement(elO, natoms=2);

  //PVC
  density = 1.3*g/cm3;
  G4Material* PVC = new G4Material(name="PVC",density, ncomponents=3);
  PVC->AddElement(elH, natoms=3);
  PVC->AddElement(elC, natoms=2);
  PVC->AddElement(elCl, natoms=1);

  //Polyethylene (plastic) LDPE on UKDMC
  density = 0.93*g/cm3;
  G4Material* Polyethylene = new G4Material(name="Polyethylene", density, ncomponents =2);
  Polyethylene -> AddElement(elH, natoms=4);
  Polyethylene -> AddElement(elC, natoms=2);

 // CH2 Shielding
  //shieldDensity = 1000.*kg/m3;
  shieldDensity = 600.*kg/m3;  
  G4Material* CH2 = new G4Material(name="CH2", shieldDensity, ncomponents=2);
  CH2->AddElement(elC, natoms=1);
  CH2->AddElement(elH, natoms=2);

// solid CH2 slab Sgielding
  solid_shieldDensity = 1000.*kg/m3;
  G4Material* CH2_solid = new G4Material(name="CH2", solid_shieldDensity, ncomponents=2);
  CH2_solid->AddElement(elC, natoms=1);
  CH2_solid->AddElement(elH, natoms=2); //

  // wood
  woodDensity = 0.9*g/cm3;
  G4Material* wood = new G4Material(name="wood", woodDensity, ncomponents=3);
  wood->AddElement(elH, natoms=4);
  wood->AddElement(elO, natoms=1);
  wood->AddElement(elC, natoms=2);


  // cardboard
  cardboardDensity = 689.0*kg/m3;
  G4Material* cardboard = new G4Material(name="cardboard", cardboardDensity, ncomponents=3);
  cardboard->AddElement(elH, natoms=4);
  cardboard->AddElement(elO, natoms=1);
  cardboard->AddElement(elC, natoms=2);
 
  //concrete
 
  G4Material* concrete = new G4Material(name="Concrete", density=2.3*g/cm3, ncomponents=6);
  concrete->AddElement(elSi, 0.227915);
  concrete->AddElement(elO, 0.60541);
  concrete->AddElement(elH, 0.09972);
  concrete->AddElement(elCa, 0.04986);
  concrete->AddElement(elAl, 0.014245);
  concrete->AddElement(elFe, 0.00285);

  // plaster board  , gypsum
  G4Material* gypsum = new G4Material("gypsum", density= 708.25*kg/m3, ncomponents =4);
  gypsum->AddElement(elCa, natoms = 1);
  gypsum->AddElement(elS, natoms = 1);
  gypsum->AddElement(elO, natoms = 6);
  gypsum->AddElement(elH, natoms = 4);


  // Alumina for ceramics
  G4Material* Alumina = new G4Material("Alumina", density=3.95*g/cm3, ncomponents=2);
  Alumina->AddElement(elAl, natoms=2);
  Alumina->AddElement(elO, natoms=3);

  // Bakelite as glass filled polymid for mupic sim
  G4Material* Bakelite = new G4Material("Bakelite", density=1.4*g/cm3,ncomponents=3);
  Bakelite->AddElement(elC,natoms=7);
  Bakelite->AddElement(elH, natoms=8);
  Bakelite->AddElement(elO, natoms=2);

  // Silica glass

  G4Material* Silica = new G4Material("Silica", density=2.203*g/cm3,ncomponents=1);
  Silica->AddElement(elSi, natoms=1);

  // SF6_gas 

  density = 156.14*g/m3; //sf6 20 torr
  pressure = 0.0263158*atmosphere; //sf6 20 torr
  G4Material* SF6_gas = new G4Material(name="SF6_gas", density, ncomponents=2, kStateGas, temperature,pressure);
  SF6_gas->AddElement(elS, natoms =1);
  SF6_gas->AddElement(elF, natoms =6);

  G4cout << "Enter CF4 fraction :" << G4endl;
  G4cin >> CF4_frac;
  G4cout << G4endl;
  G4cout << "Enter He fraction :" << G4endl;
  G4cin >> He_frac;
  G4cout << G4endl;


  // He_frac = 0.70;
  //CF4_frac = 0.30;

  //He_gas
  density = 162.488*He_frac*g/m3;
  pressure = 1*He_frac*atmosphere;
  G4Material* He_gas = new G4Material(name="He_gas", density, ncomponents=1, kStateGas, temperature,pressure);
  He_gas->AddElement(elHe, natoms=1);

  //CF4_gas
  density = 3574.736*CF4_frac*g/m3;
  pressure = 1*CF4_frac*atmosphere;
  G4Material* CF4_gas = new G4Material(name="CF4_gas", density, ncomponents=2, kStateGas, temperature, pressure);
  CF4_gas->AddElement(elC, natoms=1);
  CF4_gas->AddElement(elF, natoms=4);

  //CYGNO_gas
  density = He_gas->GetDensity()+CF4_gas->GetDensity();
  pressure = He_gas->GetPressure()+CF4_gas->GetPressure();
  G4Material* CYGNO_gas = new G4Material(name="CYGNO_gas", density, ncomponents=2, kStateGas, temperature, pressure);
  CYGNO_gas->AddMaterial(He_gas, fractionmass = He_gas->GetDensity()/density*100*perCent);
  CYGNO_gas->AddMaterial(CF4_gas, fractionmass = CF4_gas->GetDensity()/density*100*perCent);

  G4cout << "\n" << "GAS INFO:" << G4endl;
  G4cout << "CF4" << "\t" <<
    "Pressure: " << CYGNO_gas->GetMaterial("CF4_gas")->GetPressure()/(atmosphere) << " atm " << "\t" <<
    "Density: "  << CYGNO_gas->GetMaterial("CF4_gas")->GetDensity()/(g/m3) << " g/m3" <<  G4endl;
  G4cout << "He " << "\t" <<  
    "Pressure: " << CYGNO_gas->GetMaterial("He_gas")->GetPressure()/(atmosphere) << " atm " << "\t" <<
    "Density: "  << CYGNO_gas->GetMaterial("He_gas")->GetDensity()/(g/m3) << " g/m3" <<  G4endl; 
  G4cout << "Mix" << "\t" << 
    "Pressure: " << CYGNO_gas->GetPressure()/(atmosphere) << " atm" << "\t" << "\t"
    "Density: "  << CYGNO_gas->GetDensity()/(g/m3) << " g/m3" << "\n" << G4endl;

  char temp;
  G4cout << "Check gas info and press y then enter to continue..." << G4endl;
  G4cin >> temp;
  

  //kapton

  density = 1.42*g/cm3;
  G4Material* Kapton = new G4Material("Kapton", density, ncomponents=4);
  Kapton->AddElement(elC, natoms=22);
  Kapton->AddElement(elH, natoms=10);
  Kapton->AddElement(elN, natoms=2);
  Kapton->AddElement(elO, natoms=5);

  //Titanium

  density = 4.506*g/cm3;
  G4Material* TitaniumMetal = new G4Material(name="TitaniumMetal", density, ncomponents=1);
  TitaniumMetal->AddElement(elTi,natoms=1);

  //CS2
  //density = 166.9*g/m3; // ideal gas law for 41torr 300K, molar mass = 166.9g/m3
  density = 4172.5*g/m3; // x25 for more stats : 25*41torr
  G4Material* CS2 = new G4Material(name="CS2_gas", density, ncomponents=2,
        			   kStateGas, temperature);
  CS2->AddElement(elC, natoms=1);
  CS2->AddElement(elS, natoms=2);


//*************************************************************************
// VOLUMES
 
  NbOfDrift=1; //What is this?

  G4bool* CheckOverlap=0;

  //----------------------------
  // Creation of the volumes
  //----------------------------
 
  G4Box* world = new G4Box("world", world_length/2, world_width/2, world_height/2);
  world_log = new G4LogicalVolume(world,CYGNO_gas,"world_log",0,0,0); // Boulby Version
  // world_log = new G4LogicalVolume(world,RockLime,"world_log",0,0,0); // LGNS Version
  world_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0), "world_phys",world_log,0,CheckOverlap,0);
  
  G4Box* vessel_top = new G4Box("vessel_top", (vessel_length/2)+2*cm, (vessel_width/2), 1*cm);
  vessel_top_log = new G4LogicalVolume(vessel_top,Aluminium,"vessel_top_log",0,0,0);
  vessel_top_phys = new G4PVPlacement(0,G4ThreeVector(0,0,(vessel_height/2)-1*cm),"vessel_top_phys",vessel_top_log,world_phys,CheckOverlap,0);
				      
  G4Box* vessel_bottom = new G4Box("vessel_bottom", (vessel_length/2)+2*cm, (vessel_width/2), 1*cm);
  vessel_bottom_log = new G4LogicalVolume(vessel_bottom,Aluminium,"vessel_bottom_log",0,0,0);
  vessel_bottom_phys = new G4PVPlacement(0,G4ThreeVector(0,0,(-vessel_height/2)+1*cm),"vessel_bottom_phys",vessel_bottom_log,world_phys,CheckOverlap,0);
				       
  G4Box* vessel_front = new G4Box("vessel_front", (vessel_length/2)+2*cm, 1*cm, GAS_height/2);
  vessel_front_log = new G4LogicalVolume(vessel_front,Aluminium,"vessel_front_log",0,0,0);
  vessel_front_phys = new G4PVPlacement(0,G4ThreeVector(0,(vessel_width/2)-1*cm,0),"vessel_front_phys",vessel_front_log,world_phys,CheckOverlap,0);

  G4Box* vessel_back = new G4Box("vessel_back", (vessel_length/2)+2*cm, 1*cm, GAS_height/2);
  vessel_back_log = new G4LogicalVolume(vessel_back,Aluminium,"vessel_back_log",0,0,0);
  vessel_back_phys = new G4PVPlacement(0,G4ThreeVector(0,(-vessel_width/2)+1*cm,0),"vessel_back_phys",vessel_back_log,world_phys,CheckOverlap,0);
   
  G4Box* glass_left = new G4Box("glass_left", 2*cm, GAS_width/2, GAS_height/2);
  glass_left_log = new G4LogicalVolume(glass_left,Silica,"glass_left_log",0,0,0);
  glass_left_phys = new G4PVPlacement(0,G4ThreeVector(vessel_length/2,0,0),"glass_left_phys",glass_left_log,world_phys,CheckOverlap,0);

  G4Box* glass_right = new G4Box("glass_right", 2*cm, GAS_width/2, GAS_height/2);
  glass_right_log = new G4LogicalVolume(glass_right,Silica,"glass_right_log",0,0,0);
  glass_right_phys = new G4PVPlacement(0,G4ThreeVector(-vessel_length/2,0,0),"glass_right_phys",glass_right_log,world_phys,CheckOverlap,0);
  
  G4Box* GAS = new G4Box("GAS",GAS_length/2,GAS_width/2,GAS_height/2);
  GAS_log = new G4LogicalVolume(GAS,CYGNO_gas,"GAS_log",0,0,0);
  GAS_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),"GAS_phys",GAS_log,world_phys,CheckOverlap,0);

  G4Box* shield_left = new G4Box("shield_left", 2*cm, GAS_width/2, GAS_height/2);
  shield_left_log = new G4LogicalVolume(shield_left,Copper,"shield_left_log",0,0,0);
  shield_left_phys = new G4PVPlacement(0,G4ThreeVector(1215,0,0),"shield_left_phys",shield_left_log,world_phys,CheckOverlap,0);

  G4RotationMatrix* yrot = new G4RotationMatrix();;
  yrot->rotateY(90*deg);

  G4Tubs* camera_hole = new G4Tubs("camera_hole",0,2*cm,1*cm,0,2*M_PI);
  camera_hole_log = new G4LogicalVolume(camera_hole,Silica,"camera_hole_log",0,0,0);
  camera_hole_phys = new G4PVPlacement(yrot,G4ThreeVector(0,0,0),"camera_hole_phys",camera_hole_log,shield_left_phys,CheckOverlap,0);

  /*
  G4Box* camera_shield = new G4Box("camera_shield",camera_shield_length/2,camera_shield_width/2,camera_shield_height/2);
  camera_shield_log = new G4LogicalVolume(camera_shield,Copper,"camera_shield_log",0,0,0);
  camera_shield_phys = new G4PVPlacement(0,G4ThreeVector(1217*mm,0,-10*mm),"camera_shield_phys",camera_shield_log,world_phys,CheckOverlap,0);
 
  G4Box* camera_box = new G4Box("camera_box",camera_box_length/2,camera_box_width/2,camera_box_height/2);
  camera_box_log = new G4LogicalVolume(camera_box,CYGNO_gas,"camera_box_log",0,0,0);
  camera_box_phys = new G4PVPlacement(0,G4ThreeVector(1217*mm,0,0),"camera_box_phys",camera_box_log,world_phys,CheckOverlap,0);
  */
  /*

  G4Box* GEM_Kapton = new G4Box("GEM_Kapton",10*cm,12*cm,0.0025*cm);
  G4Box* GEM_copper = new G4Box("GEM_Copper",9.5*cm,11.5*cm,0.00275*cm); 
  G4SubtractionSolid* GEM_Copper = new G4SubtractionSolid("GEM_Copper", GEM_Kapton, GEM_copper);
  GEM_Kapton_log = new G4LogicalVolume(GEM_Kapton,Kapton,"GEM_Kapton_log",0,0,0);
  GEM_Kapton_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),"GEM_Kapton_phys",GEM_Kapton_log,GAS_phys,CheckOverlap,0);
  GEM_Copper_log = new G4LogicalVolume(GEM_Copper,Copper,"GEM_Copper_log",0,0,0);
  GEM_Copper_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),"GEM_Copper_phys",GEM_Copper_log,GAS_phys,CheckOverlap,0);
  G4UnionSolid* GEM = new G4UnionSolid("GEM", GEM_Kapton, GEM_Copper);
  
  */
  G4RunManager* theRunManager = G4RunManager::GetRunManager();
  theRunManager->DefineWorldVolume(world_phys);
  theRunManager->GeometryHasBeenModified();
  
//  G4SDManager* SDman = G4SDManager::GetSDMpointer();
//  allSD = new AllSD("allSD",this);
//  SDman->AddNewDetector( allSD );
//
//  GAS_log -> SetSensitiveDetector (allSD); //CHANGE WHAT PART OF THE DETECTOR IS SENSITIVE HERE!!!

  
  return world_phys;
}
