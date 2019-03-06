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

xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'

decay = 'QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU_IMG'

# Paths to input files
rhFileList = '%s/%s/*/*/output_*.root'%(eosDir, decay)
print(" >> Input file list: %s"%rhFileList)
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print(" >> %d files found"%len(rhFileList))
rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
print(' >> Input File[0]: %s'%rhFileList[0])
sort_nicely(rhFileList)

# Output path
outDir='/uscms/physics_grp/lpcml/nobackup/mandrews' # NOTE: Space here is limited, transfer files to EOS after processing
outDir='%s/%s'%(outDir, decay)
if not os.path.isdir(outDir):
    os.makedirs(outDir)
print(' >> Output directory: %s'%outDir)

proc_file = 'convert_root2pq_jet.py'
processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhtree, outDir, decay, i+1) for i,rhtree in enumerate(rhFileList)]
print(' >> Process[0]: %s'%processes[0])

pool = Pool(processes=len(processes))
pool.map(run_process, processes)
