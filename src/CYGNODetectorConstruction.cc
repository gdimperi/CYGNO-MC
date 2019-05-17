/* ************************************************
 * GEANT4 CAD INTERFACE - template
 *
 * File:      DetectorConstruction.cc
 *
 * Author:    Christopher M Poole,
 * Email:     mail@christopherpoole.net
 *
 * Date:      13th August, 2017
 **************************************************/

// CADMESH //
#include "CADMesh.hh"

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

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

// USER //
#include "CYGNODetectorConstruction.hh"
#include "CYGNODetectorLNGS.hh"
#include "CYGNODetectorMaterial.hh"



CYGNODetectorConstruction::CYGNODetectorConstruction() :
   rockThicknessOuter(-999*m),
   rockThicknessInner(-999*m),
   CYGNOLab("NoCave")

{
}

CYGNODetectorConstruction::~CYGNODetectorConstruction()
{
}

G4VPhysicalVolume* CYGNODetectorConstruction::Construct()
{
    G4NistManager * nist_manager = G4NistManager::Instance();

    //-----------------------------
    // construction of materials
    //-----------------------------
    G4cout << "Constructing materials...";

    //G4Material * air = nist_manager->FindOrBuildMaterial("G4_AIR");
    //G4Material * water = nist_manager->FindOrBuildMaterial("G4_WATER");
    CYGNODetectorMaterial* CYGNOMaterials = CYGNODetectorMaterial::GetInstance();
    G4cout << "... done" << G4endl;

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
	  G4double airInnerRadius = 10.0*m;
	  G4double airThickness = 5.0*m;
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
	  G4LogicalVolume* InnerAirSphere_log = new G4LogicalVolume(InnerAirSphere,
																CYGNOMaterials->Material("Air"),
																name_log);
	  Laboratory_log=InnerAirSphere_log;
	  size_Laboratory=G4ThreeVector(airInnerRadius,airInnerRadius,airInnerRadius);
	  tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
	  tr_Laboratory+=(rot_Laboratory*tr);
	  rot = G4RotationMatrix();// rotation of daughter volume
	  rot_Laboratory*=rot; //equivalent to rot_Laboratory=rot_Laboratory*rot
	  G4PVPlacement* InnerAirSphere_phys = new G4PVPlacement(G4Transform3D(rot,tr),InnerAirSphere_log,name_phys,OuterAirSphere_log,false,0,true);        
    }

    G4cout << "Laboratory done." << G4endl;
 
    //**********************************************************************
    // SHIELDING **********************************************************
          
    //Choose shielding, default is UPoP
    //
    // ---------------------------------- Full
    G4cout<<"Shielding Construction Started"<<G4endl;
    
    if (thick3==-999*m || thick2==-999*m || thick1==-999*m || thick0==-999*m || Mat0=="" || Mat1=="" || Mat2=="" || Mat3=="") 
    {
        G4cout << "ERROR: Please specify the thickennses and the materials of the shielding layers if you want to use the shielding design: Full" << G4endl;
	//throw std::exception();
	exit(1);
    }

    // ----------------------------------- Inner room dimensions
    G4double AirBox_x;
    G4double AirBox_y;
    G4double AirBox_z;
    G4Box* AirBox;
    AirBox_x = 1.6*m;
    AirBox_y = 1.8*m;
    AirBox_z = 2.2*m;        
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
    G4Material* Shield0Mat = SABREMaterials->FindOrBuildMaterial(Mat0);
    name_phys="Shield0";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* Shield0 = new G4Box(name_solid,0.5*Shield0_x,0.5*Shield0_y,0.5*Shield0_z);
    Shield0_log = new G4LogicalVolume(Shield0,Shield0Mat,name_log,0,0,0);
    Shielding_log = Shield0_log;
    
    // ----------------------------------- Shield 1
    G4double Shield1_x = AirBox_x + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
    G4double Shield1_y = AirBox_y + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
    G4double Shield1_z = AirBox_z + 2.*thick3 + 2.*thick2 + 2.*thick1 ;
    G4Material* Shield1Mat = SABREMaterials->FindOrBuildMaterial(Mat1);
    name_phys="Shield1";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* Shield1 = new G4Box(name_solid,0.5*Shield1_x,0.5*Shield1_y,0.5*Shield1_z);
    Shield1_log = new G4LogicalVolume(Shield1,Shield1Mat,name_log,0,0,0);
    
    // ----------------------------------- Shield 2        
    G4double Shield2_x = AirBox_x + 2.*thick3 + 2.*thick2 ;
    G4double Shield2_y = AirBox_y + 2.*thick3 + 2.*thick2 ;
    G4double Shield2_z = AirBox_z + 2.*thick3 + 2.*thick2 ;
    G4Material* Shield2Mat = SABREMaterials->FindOrBuildMaterial(Mat2);
    name_phys="Shield2";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* Shield2 = new G4Box(name_solid,0.5*Shield2_x,0.5*Shield2_y,0.5*Shield2_z);
    Shield2_log = new G4LogicalVolume(Shield2,Shield2Mat,name_log,0,0,0);
    
    // ----------------------------------- Shield 3        
    G4double Shield3_x = AirBox_x + 2.*thick3 ;
    G4double Shield3_y = AirBox_y + 2.*thick3 ;
    G4double Shield3_z = AirBox_z + 2.*thick3 ;
    G4Material* Shield3Mat = SABREMaterials->FindOrBuildMaterial(Mat3);
    name_phys="Shield3";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* Shield3 = new G4Box(name_solid,0.5*Shield3_x,0.5*Shield3_y,0.5*Shield3_z);
    Shield3_log = new G4LogicalVolume(Shield3,Shield3Mat,name_log,0,0,0);
    
    // ----------------------------------- Airbox
    name_phys="AirBox";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    AirBox = new G4Box(name_solid,0.5*AirBox_x,0.5*AirBox_y,0.5*AirBox_z);
    G4LogicalVolume* AirBox_log = new G4LogicalVolume(AirBox,SABREMaterials->Material("Air"),name_log,0,0,0);
    AirBox_log->SetVisAttributes(G4Color(1.,0.,0.));
    InsideVolume_log = AirBox_log;
    
    // ----------------------------------- Volume placements
    tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
    tr_InsideVolume+=(rot_InsideVolume*tr);
    rot = G4RotationMatrix();// rotation of daughter volume
    rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
    G4PVPlacement* Shield1_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield1_log,"Shield1",Shield0_log,false,0,true);
    tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
    tr_InsideVolume+=(rot_InsideVolume*tr);
    rot = G4RotationMatrix();// rotation of daughter volume
    rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
    G4PVPlacement* Shield2_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield2_log,"Shield2",Shield1_log,false,0,true);
    tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
    tr_InsideVolume+=(rot_InsideVolume*tr);
    rot = G4RotationMatrix();// rotation of daughter volume
    rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
    G4PVPlacement* Shield3_phys = new G4PVPlacement(G4Transform3D(rot,tr),Shield3_log,"Shield3",Shield2_log,false,0,true);
    tr = G4ThreeVector(0.,0.,0.);//translation in mother frame
    tr_InsideVolume+=(rot_InsideVolume*tr);
    rot = G4RotationMatrix();// rotation of daughter volume
    rot_InsideVolume*=rot; //equivalent to rot_InsideVolume=rot_InsideVolume*rot
    G4PVPlacement* AirBox_phys = new G4PVPlacement(G4Transform3D(rot,tr), AirBox_log, "AirBox", Shield3_log, false, 0,true);
    

    //**********************************************************************
    // ********* CYGNO volumes form CADMesh *****************************
    //**********************************************************************
    CADMesh * mesh_shell = new CADMesh("../geometry/v2/shell.stl");    
    CADMesh * mesh_camera_carter = new CADMesh("../geometry/v2/camera_carter.stl");    
    CADMesh * mesh_camera = new CADMesh("../geometry/v2/cameras.stl");    
    //CADMesh * mesh_window = new CADMesh("../geometry/v1/glass_windows.stl");    
    CADMesh * mesh_fc_support = new CADMesh("../geometry/v2/FC_support.stl");    
    CADMesh * mesh_turns_support = new CADMesh("../geometry/v2/turns_support.stl");    
    CADMesh * mesh_field_cage = new CADMesh("../geometry/v2/field_cage.stl");    
    CADMesh * mesh_gem_support = new CADMesh("../geometry/v2/gem_frame.stl");    
    CADMesh * mesh_gem = new CADMesh("../geometry/v2/gem.stl");    
    CADMesh * mesh_cathode_frame = new CADMesh("../geometry/v2/cathode_frame.stl");   

    mesh_shell->SetScale(mm);
    mesh_camera_carter->SetScale(mm);
    mesh_camera->SetScale(mm);
    //mesh_window->SetScale(mm);
    mesh_fc_support->SetScale(mm);
    mesh_field_cage->SetScale(mm);
    mesh_gem_support->SetScale(mm);
    mesh_gem->SetScale(mm);
    mesh_cathode_frame->SetScale(mm);


    //shell
    cad_shell_solid = mesh_shell->TessellatedMesh();
    cad_shell_logical = new G4LogicalVolume(cad_shell_solid, CYGNOMaterials->Material("Perspex"), "cad_shell_logical", 0, 0, 0);
    
    //camera carter
    cad_camera_carter_solid = mesh_camera_carter->TessellatedMesh();
    cad_camera_carter_logical = new G4LogicalVolume(cad_camera_carter_solid, CYGNOMaterials->Material("Perspex"), "cad_camera_carter_logical", 0, 0, 0);
    
    //cameras
    cad_cameras_all_solid = mesh_camera->TessellatedMesh();
    cad_cameras_all_logical = new G4LogicalVolume(cad_cameras_all_solid, CYGNOMaterials->Material("Perspex"), "cad_cameras_all_logical", 0, 0, 0);
    
    //window
    //cad_window_solid = mesh_window->TessellatedMesh();
    //cad_window_logical = new G4LogicalVolume(cad_window_solid, CYGNOMaterials->Material("Silica"), "cad_window_logical", 0, 0, 0);
    
    //CYGNO gas
    G4double CYGNO_x = 1250.*mm;
    G4double CYGNO_y = 1252.*mm;
    G4double CYGNO_z = 1310.*mm;
      
    name_phys="CYGNO";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* CYGNO_box = new G4Box(name_solid,0.5*CYGNO_x,0.5*CYGNO_y,0.5*CYGNO_z);
    CYGNO_log = new G4LogicalVolume(CYGNO_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
   
    //fc support
    cad_fc_support_solid = mesh_fc_support->TessellatedMesh();
    cad_fc_support_logical = new G4LogicalVolume(cad_fc_support_solid, CYGNOMaterials->Material("Cu"), "cad_fc_support_logical", 0, 0, 0);

    //turns support
    cad_turns_support_solid = mesh_turns_support->TessellatedMesh();
    cad_turns_support_logical = new G4LogicalVolume(cad_turns_support_solid, CYGNOMaterials->Material("Cu"), "cad_turns_support_logical", 0, 0, 0);
    
    //field cage  
    cad_field_cage_solid = mesh_field_cage->TessellatedMesh();
    cad_field_cage_logical = new G4LogicalVolume(cad_field_cage_solid, CYGNOMaterials->Material("Cu"), "cad_field_cage_logical", 0, 0, 0);
    
    //GEM support
    cad_gem_support_solid = mesh_gem_support->TessellatedMesh();
    cad_gem_support_logical = new G4LogicalVolume(cad_gem_support_solid, CYGNOMaterials->Material("Cu"), "cad_gem_support_logical", 0, 0, 0);
    
    //GEM
    cad_gem_solid = mesh_gem->TessellatedMesh();
    cad_gem_logical = new G4LogicalVolume(cad_gem_solid, CYGNOMaterials->Material("Kapton"), "cad_gem_logical", 0, 0, 0);
    
    //cathode
    cad_cathode_frame_solid = mesh_cathode_frame->TessellatedMesh();
    cad_cathode_frame_logical = new G4LogicalVolume(cad_cathode_frame_solid, CYGNOMaterials->Material("Cu"), "cad_cathode_frame_logical", 0, 0, 0);
    


    tr_CYGNO_gas+=G4ThreeVector(0.,0., 0.);	  

    if (CYGNOLab == "LNGS"){
    	tr_cad+=G4ThreeVector(0.,-1*size_Laboratory.y()+1*m,size_Laboratory.z()-10*m);	  
    }
    else if (CYGNOLab == "NoCave") {
    	tr_cad_internal= -1*tr_CYGNO_gas;	   
    }
    G4ThreeVector  size;
    rot_cad=(rot_Laboratory.inverse()*absrot_Rock.inverse())*absrot_cad;//Rotation of the CYGNO outer volume
    
    cad_shell_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_shell_logical,"cad_shell_physical", Laboratory_log, false, 0, true);
    cad_camera_carter_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_camera_carter_logical,"cad_camera_carter_physical", Laboratory_log, false, 0, true);
    cad_cameras_all_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_cameras_all_logical,"cad_cameras_all_physical", Laboratory_log, false, 0, true);
//    //cad_window_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
////		    cad_window_logical,"cad_window_physical", Laboratory_log, false, 0, true);
    CYGNO_phys = new G4PVPlacement(G4Transform3D(rot_CYGNO_gas,tr_CYGNO_gas),
		    CYGNO_log,"CYGNO_gas", Laboratory_log, false, 0, true);
    cad_fc_support_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_fc_support_logical,"cad_fc_support_physical", CYGNO_log, false, 0, true);
    cad_turns_support_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_turns_support_logical,"cad_turns_support_physical", CYGNO_log, false, 0, true);
    cad_field_cage_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad_internal), 
		    cad_field_cage_logical,"cad_field_cage_physical", CYGNO_log, false, 0, true);
    cad_gem_support_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad_internal), 
		    cad_gem_support_logical,"cad_gem_support_physical", CYGNO_log, false, 0, true);
    cad_gem_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad_internal), 
		    cad_gem_logical,"cad_gem_physical", CYGNO_log, false, 0, true);
    cad_cathode_frame_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad_internal), 
		    cad_cathode_frame_logical,"cad_cathode_frame_physical", CYGNO_log, false, 0, true);

    //
    //**********************************************************************
    // GLOBAL TRANSLATIONS ***************************************
    G4cout<<"Placement of Laboratory in the World started"<<G4endl;
    
    tr_Rock=-1*(tr_Laboratory+rot_Laboratory*(tr_cad+rot_cad*(tr_InsideVolume+rot_InsideVolume*(tr_cad))));//The shift of Rock_log in the world volume to make the origin be the center of the detector

    G4RotationMatrix rot_check = absrot_Rock*(rot_Laboratory*(rot_cad*(rot_InsideVolume)));  
    name_phys="externalRock";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4PVPlacement* externalRock_phys = new G4PVPlacement(G4Transform3D(absrot_Rock,tr_Rock),Rock_log,name_phys,WorldVolume_log,false,0,true);

    G4cout << "The main CYGNO volume is translated w.r.t the center of the rock volume of:\t x="<< -1*tr_Rock.x()/cm << " cm\t y=" << -1*tr_Rock.y()/cm << " cm\t z=" << -1*tr_Rock.z()/cm << " cm"<< G4endl;
    G4cout << "The Rock volume has been translated to put the main CYGNO volume in the center of the coordinate system"<< G4endl;
    G4cout<<"Placement of Laboratory in the World ended"<<G4endl;
    
   

   
    return WorldVolume_phys;
}

