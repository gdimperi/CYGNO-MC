Introduction
------------

This repository hosts the Monte Carlo simulation for CYGNO experiment.
The prerequisite to run this simulations is to have ROOT (v6.X.X), GEANT4 (v10.05.X) and CADMesh (v1.1) software installed.
For more informations about the software see:
* ROOT:  https://root.cern.ch/ 
* GEANT4: http://geant4.cern.ch/
* CADMesh: https://github.com/christopherpoole/CADMesh


Setup ROOT, GEANT4 and CADMesh
----------------------

Setup all the environment variables of ROOT and GEANT4.
In `farm-login.roma1.infn.it` you can do
```
## root and other libraries
source /chirone/soft/brew/set.sh
### geant 4
source /ua9/soft/geant4.10.05.p01-install/bin/geant4.sh
alias g4cmake='cmake -DGeant4_DIR=/ua9/soft/geant4.10.05.p01-install/lib64/Geant4-10.5.1/'
## cadmesh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ua9/soft/CADMesh-install/lib

```

In general:

```
source path-to-root-install/bin/thisroot.sh
source path-to-geant-install/bin/geant.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path-to-CADmesh-install/lib
```


You can put these commands in your `.bashrc` to execute them automatically every time you open a bash shell.

Download CYGNO-MC repository
---------------------------

Run 
```
git clone git@github.com:CYGNUS-RD/CYGNO-MC.git
```
or, if you don't want to configure a ssh key in baltig use https protocol
```
git clone https://github.com/CYGNUS-RD/CYGNO-MC.git
```

Now you have downloaded the code in `CYGNO-MC/` directory

Setup CYGNO-MC code
------------------

Create a build directory 
```
mkdir CYGNO-MC-build
```
and compile CYGNO-MC code
```
cd CYGNO-MC-build
g4cmake -Dcadmesh_DIR=/ua9/soft/CADMesh-install/lib/cmake/cadmesh-1.1.0/ ../CYGNOMC
make -j`nproc`
```

Now you have the CYGNO executable in the build directory.
You can run it in graphic mode:
```
./CYGNO
```
or specify a macro to run instructions
```
./CYGNO macro.mac
```
Some example macros are available in the `macro` directory.


A guide of CYGNO commands is in  CYGNOCommandsREADME file.

