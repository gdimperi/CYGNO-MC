import os, subprocess, sys, optparse

def parseInputArgs():
    parser = optparse.OptionParser(description='Simplified job submission script.')
    parser.add_option('-m', '--macro', default=None, help='Macro configuration to run')
    parser.add_option('-n', '--nevts', default=None, help='Number of events')
    parser.add_option('-e', '--nevtsperjob', default=None, help='Max events per job')
    parser.add_option('-f', '--file', default=None, help='File with radioactive processes')
    parser.add_option('--tag', default=None, help='Simulation version identifier')
    parser.add_option('--codedir', default=None, help='Simulation code directory')
    parser.add_option('--builddir', default=None, help='Compiled code directory')
    parser.add_option('--tmpdir', default=None, help='Temporary working directory')
    parser.add_option('--outdir', default=None, help='Directory to store output files and logs')
    (options, args) = parser.parse_args()
    return options

def SubmitJob(MACRO, CODEDIR, BUILDDIR, TMPDIR, TAG, OUTDIR):
    WORKDIR = os.path.join(TMPDIR, TAG, MACRO)
    os.makedirs(WORKDIR, exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, 'pbs_logs', TAG), exist_ok=True)

    script_template = f"""#!/bin/bash
export TAG={TAG}
export CODEDIR={CODEDIR}
export BUILDDIR={BUILDDIR}
export TMPDIR={TMPDIR}
export MACRO={MACRO}
export OUTDIR={OUTDIR}

mkdir -p $TMPDIR
cd $TMPDIR
cp $CODEDIR/macro/$MACRO.mac .
$BUILDDIR/CYGNO $MACRO.mac > $OUTDIR/pbs_logs/{TAG}/${MACRO}.log
"""

    script_filename = os.path.join(WORKDIR, f"{MACRO}.sh")
    with open(script_filename, 'w') as f:
        f.write(script_template)

    subprocess.call(['chmod', '+x', script_filename])
    subprocess.call(['qsub', '-q', 'cygno', '-d', TMPDIR, '-o', os.path.join(OUTDIR, 'pbs_logs', TAG, f"{MACRO}.log"), script_filename])
    print(f"Submitted job for macro: {MACRO}")

def main():
    options = parseInputArgs()

    required = [options.codedir, options.builddir, options.tmpdir, options.outdir, options.tag, options.macro, options.nevts]
    if not all(required):
        print("ERROR: Missing required arguments.")
        return

    CODEDIR = options.codedir.rstrip('/')
    BUILDDIR = options.builddir.rstrip('/')
    TMPDIR = options.tmpdir.rstrip('/') + '/'
    OUTDIR = options.outdir.rstrip('/') + '/'
    TAG = options.tag
    MacrosList = [m.strip() for m in options.macro.split(',')]
    NEvts = [n.strip() for n in options.nevts.split(',')]

    if len(MacrosList) != len(NEvts):
        print("ERROR: Macro list and event list sizes do not match.")
        return

    for i in range(len(MacrosList)):
        MACRO = MacrosList[i]
        SubmitJob(MACRO, CODEDIR, BUILDDIR, TMPDIR, TAG, OUTDIR)

if __name__ == "__main__":
    main()

