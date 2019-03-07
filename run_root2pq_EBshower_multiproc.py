import os, glob, re
from multiprocessing import Pool

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [int(c) if c.isdigit() else c for c in re.split('([0-9]+)',s)]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

def run_process(process):
    os.system('python %s'%process)

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-d', '--genDR', default=10, type=int, help='gen-level dR.')
parser.add_argument('-p', '--p_drop', default=0.80, type=float, help='p(drop) scale.')
args = parser.parse_args()

genDR = args.genDR
p_drop_scale = args.p_drop

xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'

#pu = 'noPU'
pu = '2016_25ns_Moriond17MC_PoissonOOTPU'
decay = 'DoublePi0Pt15To100_m0To1600_pythia8_%s_mlog_ptexp'%pu
#decay = '%s_genDR%d_recoDR16_IMG'%(decay, genDR)
decay = '%s_genDR%d_recoDR16_seedPos_phoVars_IMG'%(decay, genDR)
date_str = '190301_012801' # DR10

# Paths to input files 
rhFileLists = '%s/%s/%s/*/output_*.root'%(eosDir, decay, date_str)
print(" >> Input file list: %s"%rhFileList)
rhFileLists = glob.glob(rhFileLists)
assert len(rhFileLists) > 0
print(" >> %d files found"%len(rhFileLists))
rhFileLists = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileLists]
print(' >> Input File[0]: %s'%rhFileList[0])
sort_nicely(rhFileLists)

# Weights file
wgt_file = '%s_mvpt_weights_pdrop%.2f.npz'%(decay, p_drop_scale)
assert os.path.isfile(wgt_file) 

# Output path
outDir='/uscms/physics_grp/lpcml/nobackup/mandrews' # NOTE: Space here is limited, transfer files to EOS after processing
outDir='%s/%s'%(outDir, decay)
if not os.path.isdir(outDir):
    os.makedirs(outDir)
print(' >> Output directory: %s'%outDir)

proc_file = 'convert_root2pq_EBshower.py'
processes = ['%s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFile, outDir, decay, i+1, wgt_file) for i,rhFile in enumerate(rhFileLists)]
#processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhFile, outDir, decay, i+1) for i,rhFile in enumerate(rhFileLists)]
print(' >> Process[0]: %s'%processes[0])

#os.system('python %s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFileLists[0], outDir, decay, 1, wgt_file))
pool = Pool(processes=len(processes))
pool.map(run_process, processes)
