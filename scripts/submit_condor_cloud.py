import os
from optparse import OptionParser

parser=OptionParser()

parser.add_option('-p', '--path', default=None, help='Path where macros and condor submit files are')

(options, args) = parser.parse_args()

#first run create_submit_condor.py to split the macros and create the submit files

for filename in os.listdir(options.path):
    if filename.startswith('submit'):
        os.system('condor_submit -spool {path}/{submit_file} >> jobID_{submit_file}.txt'.format(path=options.path,submit_file=filename))
        
#after the job finishes, in order to get the output files you must do manually: condor_transfer_data <JOB_ID>