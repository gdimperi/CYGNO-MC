import os, sys, random
import BackgroundReader
from optparse import OptionParser

## python3 scripts/create_submit_condor.py -m <macroname> -n <tot_events> -j <events_per_job> -c <code_directory> -b <executable_directory> -g <geometry_path> -s <external_or_radioactivity> -f <isotopes_file> 

#example: python3 scripts/create_submit_condor.py -m macro/LIMEtest_condor -n 10000 -j 1000 -c /jupyter-workspace/private/Simulation/CYGNO-MC/ -b /jupyter-workspace/private/Simulation/CYGNO-MC-build/ -g /jupyter-workspace/private/Simulation/geometry/lime_new -s radioactivity -f test_bkg.txt

#notice: the submit files are written in a "submit_macros" directory, parallel to the code directory

def parseInputArgs():
    parser = OptionParser()

    parser.add_option('-m', '--macro', default=None, help='Macro that will be split (omit .mac)')
    parser.add_option('-n', '--nevents', default=10, help='Specify number of events (ignored if -f/--bkgfile is specified)')
    parser.add_option('-j', '--nperjob', default=10, help='Maximum number of events per job (the macro will be split in nevents/nperjob macros)')
    parser.add_option('-c', '--codedir', default='/jupyter-workspace/private/CYGNO-MC', help='Path to CYGNO-MC code')
    parser.add_option('-b', '--builddir', default='/jupyter-workspace/private/CYGNO-MC-build', help='Path to CYGNO-MC build directory')
    parser.add_option('-g', '--geopath', default='/jupyter-workspace/private/geometry/lime', help='Path to geometry files')
    parser.add_option('-s', '--typesim', default='external', help='Type of simulation: external or radioactivity')
    parser.add_option('-f', '--bkgfile', default=None, help='File with list of isotopes to be simulated')
    
    (options, args) = parser.parse_args()
    return options

def query_yes_no():
    yes = set(['yes','y'])
    no = set(['no','n'])

    choice = input()
    if choice in yes:
        return True
    elif choice in no:
        return False
    else:
        print ('Please respond with yes or no')
        query_yes_no()

def ModifyMacro(newpath,i_split,isotope='0',Z='0',A='0'):

    options = parseInputArgs()

    oldmac = options.macro.split("/")[-1]
    #For external background simulation:
    if options.typesim=='external':
        newname = oldmac+'_part'+str(i_split)
    elif options.typesim == 'radioactivity' and options.bkgfile:
        newname = isotope+'_'+oldmac+'_part'+str(i_split)

    new_file = open(newpath+'/'+newname+'.mac','w')
    old_file = open(options.macro+'.mac')

    for line in old_file:
        line.strip()
        elements=line.split()
        if len(elements)==0:
            continue
        if line.split()[0]=='/gps/ion' and options.typesim=='radioactivity' and options.bkgfile:
            line='/gps/ion %s %s\n'%(Z,A)
        #setting outname equal to macro name
        elif line.split()[0]=='/CYGNO/outfile':
            line=line.replace(elements[-1],newname)
        elif line.split()[0]=='/run/beamOn':
            #line=line.replace(elements[-1],str(int(int(elements[-1])/int(n_split))))
            line=line.replace(elements[-1],options.nperjob)
        elif line.split()[0]=='/random/setSeeds':
            seed1=int(random.random()*1000000000)
            seed2=int(random.random()*1000000000)
            line=line.replace(elements[-2],str(seed1))
            line=line.replace(elements[-1],str(seed2))
        new_file.write(line)
    new_file.close()
    old_file.close()


def CreateSubmitFile(newpath,njob,isotope='0'):
    
    options = parseInputArgs()

    oldmac = options.macro.split("/")[-1]
    if options.typesim=='external':
        print('Splitting macro for external background simulation')
        newname = 'submit_'+oldmac
        macrobase = oldmac
    elif options.typesim=='radioactivity' and options.bkgfile:
        macrobase = isotope+'_'+oldmac
        newname = 'submit_'+isotope+'_'+oldmac

    new_file = open(newpath+'/'+newname,'w')
    
    new_file.write('universe = vanilla\nexecutable = {build}/CYGNO\narguments = {base}_part$(ProcId).mac\nerror = {base}_part$(ProcId).error\noutput = {base}_part$(ProcId).out\nlog = {base}_part$(ProcId).log\ninitialdir = {code}\ngetenv = True\nshould_transfer_files = yes\nwhen_to_transfer_output = ON_EXIT\ntransfer_input_files = {build}/CYGNO, {macropath}/, {geopath}, /usr/local/lib/libcadmesh.so\ntransfer_output_files = {base}_part$(ProcId).root\n+OWNER = "condor"\nqueue {N}' .format(base = macrobase, build = options.builddir, code = options.codedir, macropath = newpath, geopath = options.geopath, N = njob))

def SplitMacros():
    #initializing random seed with the time
    random.seed()
    options = parseInputArgs()
    MACRO=options.macro
    CODEDIR=options.codedir
    nevents=int(options.nevents)
    nperjob=int(options.nperjob)
    typesim=options.typesim
    if options.bkgfile: bkg_file=options.bkgfile

    #set path and name of new macros
    oldmac = MACRO.split("/")[-1]
    newpath = CODEDIR+'../submit_macros/'+oldmac
    if os.path.exists(newpath):
        print('%s already exists. Do you want to write there? (y/n)'%(newpath))
        if not query_yes_no():
            sys.exit()
    else: os.system("mkdir -p %s"%(newpath))

    #split the macro!
    if typesim=='external':
        CreateSubmitFile(newpath,nevents//nperjob)
        for i in range(nevents//nperjob):
            ModifyMacro(newpath,i)
    elif typesim=='radioactivity' and options.bkgfile:
        print('simulating radioactivity from file %s'%(bkg_file))
        sources=BackgroundReader.ReadFile(bkg_file, CODEDIR+'/backgrounds/')
        print('we have %d isotopes'%(len(sources['Name'])))
        for j in range(0,len(sources['Name'])):
            CreateSubmitFile(newpath,int(sources['NEvents'][j])//nperjob,sources['Name'][j])
            for i in range(int(sources['NEvents'][j])//nperjob):
                ModifyMacro(newpath,i,sources['Name'][j],sources['Z'][j],sources['A'][j])
    else:
        print("Specify type of simulation (external or radioactivity) with -s/--typesim, and (only for radioactivity) the background file with -f/--bkgfile")
        sys.exit()

    return

SplitMacros()