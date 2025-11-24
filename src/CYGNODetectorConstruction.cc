#include <iostream>
#include <string.h>
#include <fstream> 

// GEANT4 //
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
//#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"
#include "G4EllipticalTube.hh"
#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"

#include "G4UserLimits.hh"

#include "G4RotationMatrix.hh"

// USER //
#include "CYGNODetectorConstruction.hh"
#include "CYGNODetectorConstructionMessenger.hh"
#include "CYGNODetectorLNGS.hh"
#include "CYGNODetectorMaterial.hh"
#include "CYGNODetectorProperty.hh"
#include "CYGNOSensitiveDetector.hh"
#include "CYGNOVolumes.hh"
#include "CADMesh.hh"
//#include "CYGNOBiasMultiParticleChangeCrossSection.hh"

CYGNODetectorConstruction::CYGNODetectorConstruction() :
   CYGNOGeomPath("../geometry/cygno_04_v3_ASCII/"),
   rockThicknessOuter(-999*m),
   rockThicknessInner(-999*m),
   //rockThicknessInner(4.*m),
   //productionLayerThickness(2.5*m),
   productionLayerThickness(-999*m),
   CYGNOLab("NoCave"),
   //CYGNOLab("LNGS"),
   //CYGNOLab("MuonLNGS"),
   CYGNOShielding("FullShield"),
  //  CYGNOShielding("NoShield"),
   thick0(0.90*m), thick1(0.40*m), thick2(0.20*m), thick3(0.05*m), 
   Mat0("Water"), Mat1("PE"), Mat2("Pb"), Mat3("Cu")
{
	
     fMessenger = new CYGNODetectorConstructionMessenger(this);
}

CYGNODetectorConstruction::~CYGNODetectorConstruction()
{
	delete fMessenger;

}

G4VPhysicalVolume* CYGNODetectorConstruction::Construct()
{
   
    //UpdateGeometryPath(CYGNOGeomPath); 
	
    G4NistManager * nist_manager = G4NistManager::Instance();
   
    //-----------------------------
    // construction of materials
    //-----------------------------
    G4cout << "Constructing materials..." << G4endl;

    //G4Material * air = nist_manager->FindOrBuildMaterial("G4_AIR");
    //G4Material * water = nist_manager->FindOrBuildMaterial("G4_WATER");
    CYGNODetectorMaterial* CYGNOMaterials = CYGNODetectorMaterial::GetInstance();
    G4cout << "... done" << G4endl;

    ////-----------------------------
    //// construction of general properties
    ////-----------------------------
    //G4cout << "Constructing general properties...";
    //CYGNODetectorProperty* CYGNOProperties = CYGNODetectorProperty::GetInstance();
    //G4cout << "... done" << G4endl;




    //**********************************************************************
    //   DEFINITION OF THE GEOMETRY
    //**********************************************************************
      
    //INITIALIZING TRANSLATION VECTORS TO 0:
    Rock_log = 0;
    size_Rock = G4ThreeVector();
    tr_Rock = G4ThreeVector();
    rot_Rock = G4RotationMatrix();
    absrot_Rock = G4RotationMatrix();
    Laboratory_log = 0;
    size_Laboratory = G4ThreeVector();
    tr_Laboratory = G4ThreeVector();
    rot_Laboratory = G4RotationMatrix();
    Shielding_log = 0;
    size_Shielding = G4ThreeVector();
    tr_Shielding = G4ThreeVector();
    rot_Shielding = G4RotationMatrix();
    absrot_Shielding = G4RotationMatrix();
    InsideVolume_log = 0;
    size_InsideVolume = G4ThreeVector();
    tr_InsideVolume = G4ThreeVector();
    rot_InsideVolume = G4RotationMatrix();

    //Name of the volumes
    G4String name_solid="";
    G4String name_log="";
    G4String name_phys="";

    G4ThreeVector tr;
    G4RotationMatrix rot;

    //**********************************************************************
    // WORLD ***************************************
    //**********************************************************************
    G4double world_x = 100.0*m;
    G4double world_y = 100.0*m;
    G4double world_z = 400.0*m;
      
    name_phys="WorldVolume";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* WorldVolume_box = new G4Box(name_solid,0.5*world_x,0.5*world_y,0.5*world_z);
    WorldVolume_log = new G4LogicalVolume(WorldVolume_box,CYGNOMaterials->Material("Vacuum"),name_log,0,0,0);
    //THE WORLD CANNOT BE TRANSLATED, therefore the first volume inside it (Rock_log) is the volume that must be translated in order to adjust the origin. Rock_log cannot be as large as the world to avoid that this volume is moved out of the world
    WorldVolume_phys = new G4PVPlacement(0,G4ThreeVector(),WorldVolume_log,name_phys,0,false,0,true);//The world volume cannot be translated
          
    //These variables are used to set the thin tube for depth studies in the externalRock_log
    G4double rockdist_z;
    G4double rockdepth_z;
        
    //**********************************************************************
    // LABORATORY ***************************************
    //**********************************************************************
    G4cout << "Constructing laboratory..." << G4endl;
    G4bool isThinTubeCompatible=false;
    // ---------------------------------- LNGS
    if (CYGNOLab == "LNGS"){
    
        /////////////////////////////////////////////////////////
	// In this configuration the translations are wrong, WIP
	// /////////////////////////////////////////////////////
	G4cout << "========== WARNING =============" << G4endl;
	G4cout << "Translation vectors in this configuration ar wrong, WIP, use 'NoCave'" << G4endl;
	G4cout << "================================" << G4endl;
	CYGNODetectorLNGS* LNGS = CYGNODetectorLNGS::GetInstance();
        if (rockThicknessOuter != -999*m)
              LNGS->SetExternalRockThickness(rockThicknessOuter);
        if (productionLayerThickness != -999*m)
              LNGS->SetProductionRockThickness(productionLayerThickness);
        if (rockThicknessInner != -999*m)
              LNGS->SetInternalRockThickness(rockThicknessInner);	  
        
	LNGS->SetDetectorMaterial(CYGNOMaterials);
	LNGS->ConstructRock();
        Rock_log=LNGS->GetRock();
        size_Rock=LNGS->GetRockSizeXYZ();
        absrot_Rock=LNGS->GetRockAbsRotation();
        Laboratory_log=LNGS->GetLaboratory();
        size_Laboratory=LNGS->GetLaboratorySizeXYZ();
        tr_Laboratory=LNGS->GetLaboratoryTranslation();
        rot_Laboratory=LNGS->GetLaboratoryRotation();
        
        ////for thin tube
        //if(LNGS->GetProductionRockThickness()==0.*cm && LNGS->GetInternalRockThickness()==0.*cm)
        //  isThinTubeCompatible=true;
        //rockdist_z=LNGS->rockdist_z;
        //rockdepth_z=LNGS->rockdepth_z;
    }
    // ---------------------------------- NoCave
    else if (CYGNOLab == "NoCave") 
	{
	  //**********************************************************************
	  // Double Air sphere surrounding the whole detector and shielding
	  //**********************************************************************        
	  G4double airInnerRadius = 25.0*m;
	  G4double airThickness = 20.0*m;
	  if (productionLayerThickness != -999*m)
		airThickness = productionLayerThickness;
	  if (rockThicknessInner != -999*m)
		airInnerRadius=rockThicknessInner;

	  //Air permeates the are around the detector
	  name_phys="OuterAirSphere";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  G4Sphere* OuterAirSphere = new G4Sphere(name_solid, 0., airInnerRadius+airThickness, 0*degree, 360*degree, 0*degree,180*degree);
	  G4LogicalVolume* OuterAirSphere_log = new G4LogicalVolume(OuterAirSphere,CYGNOMaterials->Material("Air"),name_log);
	  Rock_log=OuterAirSphere_log;
	  size_Rock=G4ThreeVector(airInnerRadius+airThickness,airInnerRadius+airThickness,airInnerRadius+airThickness);
	  absrot_Rock = G4RotationMatrix();
		
	  name_phys="InnerAirSphere";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  G4Sphere* InnerAirSphere = new G4Sphere(name_solid,  0., airInnerRadius,  0*degree, 360*degree, 0*degree,180*degree);
	  G4LogicalVolume* InnerAirSphere_log = new G4LogicalVolume(InnerAirSphere, CYGNOMaterials->Material("Air"), name_log);
	  Laboratory_log=InnerAirSphere_log;
	  size_Laboratory=G4ThreeVector(airInnerRadius,airInnerRadius,airInnerRadius);
	  tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
	  tr_Laboratory+=(rot_Laboratory*tr);
	  rot = G4RotationMatrix();// rotation of daughter volume
	  rot_Laboratory*=rot; //equivalent to rot_Laboratory=rot_Laboratory*rot
	  InnerAirSphere_phys = new G4PVPlacement(G4Transform3D(rot,tr),InnerAirSphere_log,name_phys,OuterAirSphere_log,false,0,true);        
    }
    else if (CYGNOLab == "MuonLNGS")
      {
	  //**********************************************************************
	  // Air parallelepiped inside a parallelepiped made of rock
	  //**********************************************************************        
	  G4double airInnerRadius = 3.0*m; //3.0*m for Muon simulation
	  G4double airThickness = 5*cm; //5.0*m for Muon simulation
	  if (productionLayerThickness != -999*m)
		airThickness = productionLayerThickness;
	  if (rockThicknessInner != -999*m)
		airInnerRadius=rockThicknessInner;

	  //Air permeates the are around the detector
	  name_phys="OuterAirSphere";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  //G4Sphere* OuterAirSphere = new G4Sphere(name_solid, 0., airInnerRadius+airThickness, 0*degree, 360*degree, 0*degree,180*degree); 
	  
	  //LNGS rock:
	  G4Box* OuterAirSphere = new G4Box(name_solid,(airInnerRadius+airThickness)/2.,(airInnerRadius+airThickness)/2.,(airInnerRadius+airThickness)/2.);
	  G4LogicalVolume* OuterAirSphere_log = new G4LogicalVolume(OuterAirSphere,CYGNOMaterials->Material("LNGSRock"),name_log); //LNGS rock
	  Rock_log=OuterAirSphere_log;
	  size_Rock=G4ThreeVector(airInnerRadius+airThickness,airInnerRadius+airThickness,airInnerRadius+airThickness); //seems it's not used
	  absrot_Rock = G4RotationMatrix();
		
	  name_phys="InnerAirSphere";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  G4Box* InnerAirSphere = new G4Box(name_solid,airInnerRadius/2.,airInnerRadius/2.,airInnerRadius/2.);
	  G4LogicalVolume* InnerAirSphere_log = new G4LogicalVolume(InnerAirSphere, CYGNOMaterials->Material("Air"), name_log);
	  Laboratory_log=InnerAirSphere_log;
	  size_Laboratory=G4ThreeVector(airInnerRadius,airInnerRadius,airInnerRadius);
	  tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
	  tr_Laboratory+=(rot_Laboratory*tr);
	  rot = G4RotationMatrix();// rotation of daughter volume
	  rot_Laboratory*=rot; //equivalent to rot_Laboratory=rot_Laboratory*rot
	  InnerAirSphere_phys = new G4PVPlacement(G4Transform3D(rot,tr),InnerAirSphere_log,name_phys,OuterAirSphere_log,false,0,true);

    }
    else
    {
      G4cout << "ERROR: Something went wrong with the definition of the variable CYGNOLab" << G4endl;
      //throw std::exception();
      exit(1);		
    }


    G4cout << "Laboratory done." << G4endl;
 
    //**********************************************************************
    // SHIELDING **********************************************************
    
          
    if (CYGNOShielding == "FullShield") 
    {
        G4cout<<"Shielding Construction Started"<<G4endl;
        
        if (thick3==-999*m || thick2==-999*m || thick1==-999*m || thick0==-999*m || Mat0=="" || Mat1=="" || Mat2=="" || Mat3=="") 
        {
            G4cout << "ERROR: Please specify the thickennses and the materials of the shielding layers if you want to use the shielding design: FullShield" << G4endl;
            //throw std::exception();
            exit(1);
        }

        // ----------------------------------- Inner room dimensions
        G4double AirBox_x;
        G4double AirBox_y;
        G4double AirBox_z;
        G4Box* AirBox;
        AirBox_x = 1.2*m;  //1.2*m; //cygno 2.65*m; lime 2.0*m; inner cu shield 1.2*m; inner water shield 1.8*m;
        AirBox_y = 0.7*m;  //0.7*m; //cygno 1.45*m; lime 0.8*m; inner cu shield 0.7*m; inner water shield 1.0*m;
        AirBox_z = 0.6*m;  //0.6*m; //cygno 1.45*m; lime 0.8*m; inner cu shield 0.6*m; inner water shield 1.0*m;
        tr_InsideVolume = G4ThreeVector(0.,0.,0.);
        rot_InsideVolume = G4RotationMatrix();		
        size_InsideVolume = G4ThreeVector(AirBox_x/2.,
              								AirBox_y/2.,
              								AirBox_z/2.);		
        size_Shielding = G4ThreeVector(AirBox_x/2. + thick3 + thick2 + thick1 + thick0,
              							 AirBox_y/2. + thick3 + thick2 + thick1 + thick0,
              							 AirBox_z/2. + thick3 + thick2 + thick1 + thick0);
        absrot_Shielding = G4RotationMatrix();

        // ----------------------------------- Shield 0
        G4double Shield0_x = AirBox_x + 2.*thick3 + 2.*thick2 + 2.*thick1 + 2.*thick0 ;
        G4double Shield0_y = AirBox_y + 2.*thick3 + 2.*thick2 + 2.*thick1 + 2.*thick0 ;
        G4double Shield0_z = AirBox_z + 2.*thick3 + 2.*thick2 + 2.*thick1 + 2.*thick0 ;        
        G4Material* Shield0Mat = CYGNOMaterials->Material(Mat0);
        name_phys="Shield0";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield0 = new G4Box(name_solid,0.5*Shield0_x,0.5*Shield0_y,0.5*Shield0_z);
        Shield0_log = new G4LogicalVolume(Shield0,Shield0Mat,name_log);
        Shielding_log = Shield0_log;
        Shield0_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat0));

        // ----------------------------------- Shield 1
        G4double Shield1_x = AirBox_x + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
        G4double Shield1_y = AirBox_y + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
        G4double Shield1_z = AirBox_z + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
        G4Material* Shield1Mat = CYGNOMaterials->Material(Mat1);
        name_phys="Shield1";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield1 = new G4Box(name_solid,0.5*Shield1_x,0.5*Shield1_y,0.5*Shield1_z);
        Shield1_log = new G4LogicalVolume(Shield1,Shield1Mat,name_log);
        Shield1_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat1));
        
        // ----------------------------------- Shield 2        
        G4double Shield2_x = AirBox_x + 2.*thick3 + 2.*thick2 ;
        G4double Shield2_y = AirBox_y + 2.*thick3 + 2.*thick2 ;
        G4double Shield2_z = AirBox_z + 2.*thick3 + 2.*thick2 ;
        G4Material* Shield2Mat = CYGNOMaterials->Material(Mat2);
        name_phys="Shield2";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield2 = new G4Box(name_solid,0.5*Shield2_x,0.5*Shield2_y,0.5*Shield2_z);
        Shield2_log = new G4LogicalVolume(Shield2,Shield2Mat,name_log);
        Shield2_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat2));
        
        // ----------------------------------- Shield 3        
        G4double Shield3_x = AirBox_x + 2.*thick3 ;
        G4double Shield3_y = AirBox_y + 2.*thick3 ;
        G4double Shield3_z = AirBox_z + 2.*thick3 ;
        G4Material* Shield3Mat = CYGNOMaterials->Material(Mat3);
        name_phys="Shield3";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        G4Box* Shield3 = new G4Box(name_solid,0.5*Shield3_x,0.5*Shield3_y,0.5*Shield3_z);
        Shield3_log = new G4LogicalVolume(Shield3,Shield3Mat,name_log);
        Shield3_log->SetVisAttributes(CYGNOMaterials->VisAttributes(Mat3));
        
        // ----------------------------------- Airbox
        name_phys="AirBox";
        name_log=name_phys+"_log";
        name_solid=name_phys+"_solid";
        AirBox = new G4Box(name_solid,0.5*AirBox_x,0.5*AirBox_y,0.5*AirBox_z);
        AirBox_log = new G4LogicalVolume(AirBox,CYGNOMaterials->Material("Air"),name_log);
        AirBox_log->SetVisAttributes(CYGNOMaterials->VisAttributes("Air"));
	InsideVolume_log = AirBox_log;
    
    }
    // ---------------------------------- No shielding
    else if (CYGNOShielding == "NoShield") 
    {
	  G4double AirBox_x;
	  G4double AirBox_y;
	  G4double AirBox_z;
	  G4Box* AirBox;
          AirBox_x = 5.*m;
          AirBox_y = 3.*m;
          AirBox_z = 3.*m;       
	  name_phys="AirBox";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  AirBox = new G4Box(name_solid,0.5*AirBox_x,0.5*AirBox_y,0.5*AirBox_z);
	  AirBox_log = new G4LogicalVolume(AirBox,CYGNOMaterials->Material("Air"),name_log,0,0,0);

	  Shielding_log=AirBox_log;
	  size_Shielding=G4ThreeVector(AirBox_x/2.,AirBox_y/2.,AirBox_z/2.);
	  absrot_Shielding = G4RotationMatrix();
	  size_InsideVolume=G4ThreeVector(AirBox_x/2.,AirBox_y/2.,AirBox_z/2.);
	  tr_InsideVolume=G4ThreeVector(0.,0.,0.);
	  rot_InsideVolume=G4RotationMatrix();
    }
    else if (CYGNOShielding == "CYGNO_04") 
    {
	  G4double AirBox_x;
	  G4double AirBox_y;
	  G4double AirBox_z;
	  G4Box* AirBox;
          AirBox_x = 7.*m;
          AirBox_y = 7.*m;
          AirBox_z = 7.*m;        
	  name_phys="AirBox";
	  name_log=name_phys+"_log";
	  name_solid=name_phys+"_solid";
	  AirBox = new G4Box(name_solid,0.5*AirBox_x,0.5*AirBox_y,0.5*AirBox_z);
	  AirBox_log = new G4LogicalVolume(AirBox,CYGNOMaterials->Material("Air"),name_log,0,0,0);

	  Shielding_log=AirBox_log;
	  size_Shielding=G4ThreeVector(AirBox_x/2.,AirBox_y/2.,AirBox_z/2.);
	  absrot_Shielding = G4RotationMatrix();
	  size_InsideVolume=G4ThreeVector(AirBox_x/2.,AirBox_y/2.,AirBox_z/2.);
	  tr_InsideVolume=G4ThreeVector(0.,0.,0.);
	  rot_InsideVolume=G4RotationMatrix();
    }
    else
    {
      G4cout << "ERROR: Something went wrong with the definition of the variable CYGNOShielding" << G4endl;
      //throw std::exception();
      exit(1);		
    }
    
    G4double tolerance = 1*mm;
    

    //**********************************************************************
    // ********* CYGNO volumes form CADMesh *****************************
    //**********************************************************************
    
    char namestl[70];
    
    snprintf(namestl, sizeof(namestl), "%s/Cathode.stl", CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;   
    ifstream infile(CYGNOGeomPath.c_str());
    if (infile.good())
      mesh_Cathode = CADMesh::TessellatedMesh::FromSTL(namestl);     
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;

    snprintf(namestl, sizeof(namestl), "%s/FCSupport.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_FCSupport = CADMesh::TessellatedMesh::FromSTL(namestl);
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;
    
    snprintf(namestl, sizeof(namestl), "%s/FieldCage.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_FieldCage = CADMesh::TessellatedMesh::FromSTL(namestl);  
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;
    
    snprintf(namestl, sizeof(namestl), "%s/GEM.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_GEM = CADMesh::TessellatedMesh::FromSTL(namestl);
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;
    
    snprintf(namestl, sizeof(namestl), "%s/GemFrame.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_GemFrame = CADMesh::TessellatedMesh::FromSTL(namestl); 
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;
    
    snprintf(namestl, sizeof(namestl), "%s/InnerShieldCu.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;     
    if (infile.good())
      mesh_InnerShieldCu = CADMesh::TessellatedMesh::FromSTL(namestl);
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;

    snprintf(namestl, sizeof(namestl), "%s/OuterShieldCu.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_OuterShieldCu = CADMesh::TessellatedMesh::FromSTL(namestl); 
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;
 
    snprintf(namestl, sizeof(namestl), "%s/PEBase.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl;
    if (infile.good())
      mesh_PEBase = CADMesh::TessellatedMesh::FromSTL(namestl); 
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;
    
    snprintf(namestl, sizeof(namestl), "%s/PMMABox.stl",CYGNOGeomPath.c_str());
    G4cout << namestl << G4endl; 
    if (infile.good())
      mesh_PMMABox = CADMesh::TessellatedMesh::FromSTL(namestl);
    else G4cout << "########## WARNING: CAD geometry from "<< namestl  <<" not loaded! ##############" << G4endl;
    
        

    if (infile.good()){
    
      //Cathode
      cad_Cathode_solid = mesh_Cathode->GetSolid();
      cad_Cathode_logical = new G4LogicalVolume(cad_Cathode_solid, CYGNOMaterials->Material("Cu"), "cad_Cathode_logical");
      cad_Cathode_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));


      //FCSupport
      cad_FCSupport_solid = mesh_FCSupport->GetSolid();
      cad_FCSupport_logical = new G4LogicalVolume(cad_FCSupport_solid, CYGNOMaterials->Material("Perspex"), "cad_FCSupport_logical", 0, 0, 0);
      cad_FCSupport_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
      
      
      //FIeldCage
      cad_FieldCage_solid = mesh_FieldCage->GetSolid();
      cad_FieldCage_logical = new G4LogicalVolume(cad_FieldCage_solid, CYGNOMaterials->Material("FieldCage"), "cad_FieldCage_logical", 0, 0, 0);
      cad_FieldCage_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
 

      //GEM  
      cad_GEM_solid = mesh_GEM->GetSolid();
      cad_GEM_logical = new G4LogicalVolume(cad_GEM_solid, CYGNOMaterials->Material("GEM"), "cad_GEM_logical", 0, 0, 0);
      cad_GEM_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));

      //GEM frame
      cad_GemFrame_solid = mesh_GemFrame->GetSolid();
      cad_GemFrame_logical = new G4LogicalVolume(cad_GemFrame_solid, CYGNOMaterials->Material("Perspex"), "cad_GemFrame_logical", 0, 0, 0);
      cad_GemFrame_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Perspex"));
      
      //InnerShield Cu
      cad_InnerShieldCu_solid = mesh_InnerShieldCu->GetSolid();
      cad_InnerShieldCu_logical = new G4LogicalVolume(cad_InnerShieldCu_solid, CYGNOMaterials->Material("Cu"), "cad_InnerShieldCu_logical", 0, 0, 0);
      cad_InnerShieldCu_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));


      //OuterShield Cu
      cad_OuterShieldCu_solid = mesh_OuterShieldCu->GetSolid();
      cad_OuterShieldCu_logical = new G4LogicalVolume(cad_OuterShieldCu_solid, CYGNOMaterials->Material("Cu"), "cad_OuterShieldCu_logical", 0, 0, 0); 
      cad_OuterShieldCu_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));
      
        
      //PEbase
      cad_PEBase_solid = mesh_PEBase->GetSolid();
      cad_PEBase_logical = new G4LogicalVolume(cad_PEBase_solid, CYGNOMaterials->Material("PE"), "cad_PEBase_logical", 0, 0, 0);
      cad_PEBase_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("PE"));
        
      //PMMABox
      cad_PMMABox_solid = mesh_PMMABox->GetSolid();
      cad_PMMABox_logical = new G4LogicalVolume(cad_PMMABox_solid, CYGNOMaterials->Material("Perspex"), "cad_PMMABox_logical", 0, 0, 0);
      cad_PMMABox_logical->SetVisAttributes(CYGNOMaterials->VisAttributes("Cu"));

    }
  
    //TPC gas
    G4double TPC_x = 1000.*mm;
    G4double TPC_y = 800.*mm;
    G4double TPC_z = 500.*mm;
      
    name_phys="TPC";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* TPC_box = new G4Box(name_solid,0.5*TPC_x,0.5*TPC_y,0.5*TPC_z);
    TPC_log = new G4LogicalVolume(TPC_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
    
    //CYGNO fiducial gas
    G4double CYGNO_x = 480.*mm;
    G4double CYGNO_y = 750.*mm;
    G4double CYGNO_z = 450.*mm;
      
    name_phys="CYGNO";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* CYGNO_box = new G4Box(name_solid,0.5*CYGNO_x,0.5*CYGNO_y,0.5*CYGNO_z);
    CYGNO_log = new G4LogicalVolume(CYGNO_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
  
////FIXME remove step limiter    
//    G4double maxStep = 10*um;
//    fStepLimit = new G4UserLimits(maxStep);
//    CYGNO_log->SetUserLimits(fStepLimit); 

    CYGNO_log->SetVisAttributes(CYGNOMaterials->VisAttributes("CYGNO_gas"));
    
    //FIXME

    //tr_cad=G4ThreeVector(-3554*mm,-3845*mm,230.*mm);
    tr_cad=G4ThreeVector(-35.5*cm,-25.8*cm,-26*cm);
   
    G4double ztr_cam = 785.*mm ; //893.4mm
    G4ThreeVector trcam0(ztr_cam,-20.*mm,0.);
    G4RotationMatrix* rotcam0 = new G4RotationMatrix();
   
    G4ThreeVector trlens0(ztr_cam+1.*cm,-20.*mm,0.);
    G4RotationMatrix* rotlens0 = new G4RotationMatrix();
    rotlens0->rotateY(90*deg);

    
    if (CYGNOLab == "LNGS"){
	tr+=G4ThreeVector(0.,-1*size_Laboratory.y()+size_Shielding.y(),size_Laboratory.z()-10*m);
	tr_cad+=G4ThreeVector(0.,1.0*m-1*size_Laboratory.y()+size_Shielding.y(),size_Laboratory.z()-10*m);	  
    
	rot = G4RotationMatrix();// rotation of daughter volume
	tr_Shielding+=(rot_Shielding*tr);
    }
    else if (CYGNOLab == "NoCave" || CYGNOLab == "MuonLNGS") {
	tr=G4ThreeVector(0.,0.,0.);
	tr_cad+=G4ThreeVector(0.,0.,0.);
	rot = G4RotationMatrix();
	tr_Shielding+=(rot_Shielding*tr);

    }

    if (CYGNOShielding == "NoShield")
    {
        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        AirBox_log = Shielding_log;
    }
    
    else if (CYGNOShielding == "FullShield") 
    {   
        // ----------------------------------- Volume placements

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        Shield1_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield1_log,"Shield1",Shielding_log,false,0,true);

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        Shield2_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield2_log,"Shield2",Shield1_log,false,0,true);

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        Shield3_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield3_log,"Shield3",Shield2_log,false,0,true);

        tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
        tr_InsideVolume+=(rot_InsideVolume*tr);
        rot = G4RotationMatrix();// rotation of daughter volume
        rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
        AirBox_phys = new G4PVPlacement(G4Transform3D(rot,tr), AirBox_log, "AirBox", Shield3_log, false, 0,true); 
    
    }
    
    G4ThreeVector  size;
    

    if (infile.good()){
      cad_PEBase_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_PEBase_logical,"cad_PEBase_physical", AirBox_log, false, 0, true); 
      cad_InnerShieldCu_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_InnerShieldCu_logical,"cad_InnerShieldCu_physical", AirBox_log, false, 0, true);
      cad_OuterShieldCu_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_OuterShieldCu_logical,"cad_OuterShieldCu_physical", AirBox_log, false, 0, true);
      cad_PMMABox_physical = new G4PVPlacement(G4Transform3D(rot_cad_shield,tr_cad), 
        	    cad_PMMABox_logical,"cad_PMMABox_physical", AirBox_log, false, 0, true); 

      cad_FieldCage_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_FieldCage_logical,"cad_FieldCage_physical", AirBox_log, false, 0, true);
      cad_FCSupport_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_FCSupport_logical,"cad_FCSupport_physical", TPC_log, false, 0, true);
      cad_GEM_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_GEM_logical,"cad_GEM_physical", TPC_log, false, 0, true);
      cad_GemFrame_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
        	    cad_GemFrame_logical,"cad_GemFrame_physical", TPC_log, false, 0, true);
    }
    tr=G4ThreeVector(0.,0.,0.);
    //FIXME
    TPC_phys = new G4PVPlacement(G4Transform3D(rot, tr),
      	    TPC_log,"TPC_gas", AirBox_log, false, 0, true);
    
    tr_CYGNO_gas_1= tr;
    //tr_CYGNO_gas_1=G4ThreeVector(TPC_x/2.-CYGNO_x/2.-50.*mm,-20.*mm,0.);
    CYGNO_phys = new G4PVPlacement(G4Transform3D(rot,tr_CYGNO_gas_1),
      	    CYGNO_log,"CYGNO_gas", TPC_log, false, 0, true);
          
    tr=G4ThreeVector(0.,0.,0.);
    rot = G4RotationMatrix();
    //if (infile.good()){
    //  cad_Cathode_physical = new G4PVPlacement(G4Transform3D(rot,tr_cad), 
    //    	    cad_Cathode_logical,"cad_Cathode_physical", TPC_log, false, 0, true);
    //}  
    
	camera_phys = new G4PVPlacement(rotcam0,trcam0,camera_log,"camera",AirBox_log, false, 0, true);
        camera_lens_phys = new G4PVPlacement(rotlens0,trlens0,camera_lens_log,"camera_lens",AirBox_log, false, 0, true);

    //
    //**********************************************************************
    // GLOBAL TRANSLATIONS ***************************************
    G4cout<<"Placement of Laboratory in the World started"<<G4endl;
    
    tr_Rock=-1*(tr_Laboratory+rot_Laboratory*(tr+rot_cad*(tr_InsideVolume+rot_InsideVolume*(tr_Shielding))));//The shift of Rock_log in the world volume to make the origin be the center of the detector

    G4RotationMatrix rot_check = absrot_Rock*(rot_Laboratory*(rot_cad*(rot_InsideVolume)));  
    name_phys="externalRock";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    externalRock_phys = new G4PVPlacement(G4Transform3D(absrot_Rock,tr_Rock),Rock_log,name_phys,WorldVolume_log,false,0,true);

    G4cout << "The main CYGNO volume is translated w.r.t the center of the rock volume of:\t x="<< -1*tr_Rock.x()/cm << " cm\t y=" << -1*tr_Rock.y()/cm << " cm\t z=" << -1*tr_Rock.z()/cm << " cm"<< G4endl;

    G4cout << "The Rock volume has been translated to put the main CYGNO volume in the center of the coordinate system"<< G4endl;
    G4cout<<"Placement of Laboratory in the World ended"<<G4endl;
 
    //======= Save volumes mass and density ======
    
    G4cout<<"Saving masses and densities of the volumes"<<G4endl;
    SaveMassAndDensity();


    //===========
    return WorldVolume_phys;


}

void CYGNODetectorConstruction::SaveMassAndDensity()
{
  CYGNODetectorProperty* CYGNOProperties = CYGNODetectorProperty::GetInstance();

  G4cout << "Saving masses and densities of the volumes of the CYGNODetectorConstruction class"<< G4endl;
  //CYGNOProperties->AddVolumeNameMassAndDensity(Rock_log);
  //CYGNOProperties->AddVolumeNameMassAndDensity(Laboratory_log);
  if (CYGNOShielding=="FullShield"){
      CYGNOProperties->AddVolumeNameMassAndDensity(Shielding_log);
      CYGNOProperties->AddVolumeNameMassAndDensity(Shield1_log);
      CYGNOProperties->AddVolumeNameMassAndDensity(Shield2_log);
      CYGNOProperties->AddVolumeNameMassAndDensity(Shield3_log);
  }
  CYGNOProperties->AddVolumeNameMassAndDensity(AirBox_log);
  ifstream infile(CYGNOGeomPath.c_str());
  if (infile.good()) {
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_Cathode_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_FCSupport_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(camera_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(camera_lens_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(TPC_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(CYGNO_log);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_FieldCage_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_GemFrame_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_InnerShieldCu_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_OuterShieldCu_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_Cathode_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_GEM_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_PEBase_logical);
    CYGNOProperties->AddVolumeNameMassAndDensity(cad_PMMABox_logical);
  //  CYGNOProperties->AddVolumeNameMassAndDensity(cad_WaterShielding_logical);
  }

//  if ( productionRockThinTube_phys )
//	{
//	  CYGNOProperties->AddPhysVolumeNameMassAndDensity(productionRockThinTube_phys);
//	}
  G4cout << "All volume masses and densities saved"<< G4endl;
}

void CYGNODetectorConstruction::UpdateGeometry()
{

  G4cout << "Updating the Geometry"<< G4endl;
  CYGNODetectorProperty* CYGNOProperties = CYGNODetectorProperty::GetInstance();
  CYGNOProperties->Refresh();

  //Removing sensitive detectors
  CYGNO_log->SetSensitiveDetector(0);
  CYGNO_log=0;

  //Deleting all the solids, logical and physical objects
  //G4RunManager::GetRunManager()->ReinitializeGeometry(true);
  //Equivalent to
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore* PhysicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
  PhysicalVolumeStore->Clean();
  G4LogicalVolumeStore* LogicalVolumeStore = G4LogicalVolumeStore::GetInstance();
  LogicalVolumeStore->Clean();
  G4SolidStore* SolidStore = G4SolidStore::GetInstance();
  SolidStore->Clean();
  
  CYGNODetectorMaterial* CYGNOMaterials = CYGNODetectorMaterial::GetInstance();
  CYGNOMaterials->Refresh();


  //The memory for these pointers has been freed by the above Clean() methods
  CYGNODetectorLNGS* CYGNOLNGS = CYGNODetectorLNGS::GetInstance();
  CYGNOLNGS->Refresh();

  Shield0_log=0; 
  Shield1_log=0; 
  Shield2_log=0; 
  Shield3_log=0; 
  AirBox_log=0;

  cad_Cathode_logical=0;
  cad_FCSupport_logical=0;
  TPC_log=0;
  CYGNO_log=0;
  cad_FieldCage_logical=0;
  cad_InnerShieldCu_logical=0;
  cad_GemFrame_logical=0;
  cad_OuterShieldCu_logical=0;
  cad_Cathode_logical=0;
  cad_GEM_logical=0;
  cad_PEBase_logical=0;
  cad_PMMABox_logical=0;
  camera_log=0;
  camera_lens_log=0;

  InsideVolume_log=0;
  Shielding_log=0;
  Laboratory_log=0;
  Rock_log=0;
  WorldVolume_log=0;
  productionRockThinTube_phys=0;
  //log->ClearDaughters();

  // Delete all the geometry you had defined and build everything from scratch
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
void CYGNODetectorConstruction::UpdateGeometryPath(G4String newpath)
{
  G4cout << "Updating the Geometry path to "<< newpath << G4endl;
  SetGeomPath(newpath); 

}

void CYGNODetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //
  auto cygnoSD
    = new CYGNOSensitiveDetector("cygnoSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(cygnoSD);
  SetSensitiveDetector("CYGNO_log",cygnoSD);

}
