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

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-d', '--genDR', default=10, type=int, help='gen-level dR.')
parser.add_argument('-p', '--p_drop', default=0.95, type=float, help='p(drop) scale.')
args = parser.parse_args()

eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'
#eosDir='/store/user/lpcml/mandrews/IMG'
#outDir='~lpcml/nobackup/mandrews' # NOTE: Space here is limited, transfer files to EOS after processing
xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN

genDR = args.genDR
p_drop_scale = args.p_drop

decay = 'DoublePi0Pt15To100_m0To1600_pythia8_noPU_mlog_ptexp'
#decay = '%s_genDR%d_recoDR16_IMG'%(decay, genDR)
decay = '%s_genDR%d_recoDR16_seedPos_IMG'%(decay, genDR)
if genDR == 5:
    date_str = '190219_012120' #DR5
else:
    #date_str = '190219_012301' #DR10
    date_str = '190221_191131' #DR10
#decay = 'DoublePi0Pt15To100_m0To1600_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_IMG'
#date_str = '190214_155902'
#decay = 'DoublePi0Pt15To100_m000_pythia8_noPU_genDR10_recoDR16_IMG'
#date_str = '190220_160854'
#decay = 'DoublePhotonPt15To100_pythia8_noPU_IMG'
#date_str = '190220_155357'
#decay = 'DoublePhotonPt15To100_pythia8_noPU_seedPos_IMG'
#date_str = '190221_195631'

# Load input TTrees into TChain
rhTreeStrs = '%s/%s/%s/*/output_*.root'%(eosDir, decay, date_str)
print(" >> Input string: %s"%rhTreeStrs)
rhTreeStrs = glob.glob(rhTreeStrs)
assert len(rhTreeStrs) > 0
print(" >> %d files found"%len(rhTreeStrs))
rhTreeStrs = [('%s/%s'%(xrootd, rhtree)).replace('/eos/uscms','') for rhtree in rhTreeStrs]
print(' >> infile[0]: %s'%rhTreeStrs[0])
sort_nicely(rhTreeStrs)

wgt_file = '%s_mvpt_weights_pdrop%.2f.npz'%(decay, p_drop_scale)
assert os.path.isfile(wgt_file) 

outDir='/uscms/physics_grp/lpcml/nobackup/mandrews/%s'%(decay)
if not os.path.isdir(outDir):
    os.makedirs(outDir)
print(' >> outdir: %s'%outDir)

proc_file = 'convert_root2pq_EBshower.py'
#proc_file = 'convert_root2pq_EBshower_lhood.py'
processes = ['%s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhtree, outDir, decay, i+1, wgt_file) for i,rhtree in enumerate(rhTreeStrs)]
#processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhtree, outDir, decay, i+1) for i,rhtree in enumerate(rhTreeStrs)]
print(' >> process[0]: %s'%processes[0])

def run_process(process):
    os.system('python %s'%process)

pool = Pool(processes=len(processes))
pool.map(run_process, processes)
