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
    // ********* CYGNO volumes form CADMesh *****************************
    //**********************************************************************
    CADMesh * mesh_shell = new CADMesh("/storage/giulia/CYGNO/geometry/v1/shell.stl");    
    CADMesh * mesh_camera_carter = new CADMesh("/storage/giulia/CYGNO/geometry/v1/carters.stl");    
    CADMesh * mesh_camera = new CADMesh("/storage/giulia/CYGNO/geometry/v1/cameras_all.stl");    
    CADMesh * mesh_window = new CADMesh("/storage/giulia/CYGNO/geometry/v1/glass_windows.stl");    
    CADMesh * mesh_internal_structure = new CADMesh("/storage/giulia/CYGNO/geometry/v1/internal_structure.stl");    
    CADMesh * mesh_gem_support = new CADMesh("/storage/giulia/CYGNO/geometry/v1/gem_support.stl");    
    CADMesh * mesh_cathode_frame = new CADMesh("/storage/giulia/CYGNO/geometry/v1/cathode_frame.stl");   
    CADMesh * mesh_square_turn = new CADMesh("/storage/giulia/CYGNO/geometry/v1/square_turn.stl");    

    mesh_shell->SetScale(mm);
    mesh_camera_carter->SetScale(mm);
    mesh_camera->SetScale(mm);
    mesh_window->SetScale(mm);
    mesh_internal_structure->SetScale(mm);
    mesh_gem_support->SetScale(mm);
    mesh_cathode_frame->SetScale(mm);
    mesh_square_turn->SetScale(mm);


    //shell
    cad_shell_solid = mesh_shell->TessellatedMesh();
    cad_shell_logical = new G4LogicalVolume(cad_shell_solid, CYGNOMaterials->Material("Perspex"), "cad_shell_logical", 0, 0, 0);
    
    //camera carter
    cad_camera_carter_solid = mesh_camera_carter->TessellatedMesh();
    cad_camera_carter_logical = new G4LogicalVolume(cad_camera_carter_solid, CYGNOMaterials->Material("Perspex"), "cad_camera_carter_logical", 0, 0, 0);
    
    //cameras
    cad_cameras_all_solid = mesh_camera->TessellatedMesh();
    cad_cameras_all_logical = new G4LogicalVolume(cad_cameras_all_solid, CYGNOMaterials->Material("Perspex"), "cad_cameras_all_logical", 0, 0, 0);
    
    //glass window
    cad_window_solid = mesh_window->TessellatedMesh();
    cad_window_logical = new G4LogicalVolume(cad_window_solid, CYGNOMaterials->Material("Silica"), "cad_window_logical", 0, 0, 0);
    
    //CYGNO gas
    G4double CYGNO_x = 1250.*mm;
    G4double CYGNO_y = 1250.*mm;
    G4double CYGNO_z = 1300.*mm;
      
    name_phys="CYGNO";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4Box* CYGNO_box = new G4Box(name_solid,0.5*CYGNO_x,0.5*CYGNO_y,0.5*CYGNO_z);
    CYGNO_log = new G4LogicalVolume(CYGNO_box,CYGNOMaterials->Material("CYGNO_gas"),name_log,0,0,0);
   
    //internal structure
    cad_internal_structure_solid = mesh_internal_structure->TessellatedMesh();
    cad_internal_structure_logical = new G4LogicalVolume(cad_internal_structure_solid, CYGNOMaterials->Material("Cu"), "cad_internal_structure_logical", 0, 0, 0);

    //GEM support
    cad_gem_support_solid = mesh_gem_support->TessellatedMesh();
    cad_gem_support_logical = new G4LogicalVolume(cad_gem_support_solid, CYGNOMaterials->Material("Cu"), "cad_gem_support_logical", 0, 0, 0);
    
    //cathode
    cad_cathode_frame_solid = mesh_cathode_frame->TessellatedMesh();
    cad_cathode_frame_logical = new G4LogicalVolume(cad_cathode_frame_solid, CYGNOMaterials->Material("Cu"), "cad_cathode_frame_logical", 0, 0, 0);
    
    //fili 
    cad_square_turn_solid = mesh_square_turn->TessellatedMesh();
    cad_square_turn_logical = new G4LogicalVolume(cad_square_turn_solid, CYGNOMaterials->Material("Cu"), "cad_square_turn_logical", 0, 0, 0);
    
    //GEM



    if (CYGNOLab == "LNGS"){
    	tr_cad+=G4ThreeVector(0.,-1*size_Laboratory.y()+1*m,size_Laboratory.z()-10*m);	  
    }
    else if (CYGNOLab == "NoCave") {}
    G4ThreeVector  size;
    rot_cad=(rot_Laboratory.inverse()*absrot_Rock.inverse())*absrot_cad;//Rotation of the CYGNO outer volume
    //cad_physical = new G4PVPlacement(0, G4ThreeVector(), cad_logical,
    cad_shell_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_shell_logical,"cad_shell_physical", Laboratory_log, false, 0, true);
    cad_camera_carter_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_camera_carter_logical,"cad_camera_carter_physical", Laboratory_log, false, 0, true);
    cad_cameras_all_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_cameras_all_logical,"cad_cameras_all_physical", Laboratory_log, false, 0, true);
    cad_window_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_window_logical,"cad_window_physical", Laboratory_log, false, 0, true);
    CYGNO_phys = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad),
		    CYGNO_log,"CYGNO_gas", Laboratory_log, false, 0, true);
    cad_internal_structure_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_internal_structure_logical,"cad_internal_structure_physical", CYGNO_log, false, 0, true);
    cad_gem_support_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_gem_support_logical,"cad_gem_support_physical", CYGNO_log, false, 0, true);
    cad_cathode_frame_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_cathode_frame_logical,"cad_cathode_frame_physical", CYGNO_log, false, 0, true);
    cad_square_turn_physical = new G4PVPlacement(G4Transform3D(rot_cad,tr_cad), 
		    cad_square_turn_logical,"cad_square_turn_physical", CYGNO_log, false, 0, true);

    //
    //**********************************************************************
    // GLOBAL TRANSLATIONS ***************************************
    G4cout<<"Placement of Laboratory in the World started"<<G4endl;
    
//    tr_Rock=-1*(tr_Laboratory+rot_Laboratory*(tr_Shielding+rot_Shielding*(tr_InsideVolume+rot_InsideVolume*(tr_Vessel+rot_Vessel*(tr_Scintillator+rot_Scintillator*(tr_Assembly+rot_Assembly*(rot_Enclosure*(tr_Crystal))))))));//The shift of Rock_log in the world volume to make the origin be the center of the crystal in the middle of the Assembly volume
    tr_Rock=-1*(tr_Laboratory+rot_Laboratory*(tr_cad+rot_cad*(tr_InsideVolume+rot_InsideVolume*(tr_cad))));//The shift of Rock_log in the world volume to make the origin be the center of the detector

    //G4RotationMatrix rot_check = absrot_Rock*(rot_Laboratory*(rot_Shielding*(rot_InsideVolume*(rot_Vessel*(rot_Scintillator*(rot_Assembly*rot_Enclosure))))));  
    G4RotationMatrix rot_check = absrot_Rock*(rot_Laboratory*(rot_cad*(rot_InsideVolume)));  
    //if (rot_check != absrot_Enclosure)
    //      {
    //        G4cout << "ERROR: the product of all the rotations does not coincide with the desired absolute rotation of the enclosures" << G4endl;
    //        //throw std::exception();
    //        exit(1);		
    //      }

    //Translating the rock in order to have the origin of the coordinate system in the middle of the crystal
    name_phys="externalRock";
    name_log=name_phys+"_log";
    name_solid=name_phys+"_solid";
    G4PVPlacement* externalRock_phys = new G4PVPlacement(G4Transform3D(absrot_Rock,tr_Rock),Rock_log,name_phys,WorldVolume_log,false,0,true);

    G4cout << "The main CYGNO volume is translated w.r.t the center of the rock volume of:\t x="<< -1*tr_Rock.x()/cm << " cm\t y=" << -1*tr_Rock.y()/cm << " cm\t z=" << -1*tr_Rock.z()/cm << " cm"<< G4endl;
    G4cout << "The Rock volume has been translated to put the main CYGNO volume in the center of the coordinate system"<< G4endl;
    G4cout<<"Placement of Laboratory in the World ended"<<G4endl;
    
   

   
    return WorldVolume_phys;
}

