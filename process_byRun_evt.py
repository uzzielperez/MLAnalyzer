import h5py
import numpy as np
from dask.delayed import delayed
import dask.array as da
import glob

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-l', '--label', required=True, type=int, help='Decay label.')
parser.add_argument('-f', '--list_idx', required=True, type=int, help='List index, 0: train, 1:test')
parser.add_argument('-r', '--run', default=0, type=int, help='First run')
args = parser.parse_args()

#eosDir='/eos/uscms/store/user/mba2012/IMG'
eosDir='/eos/cms/store/user/mandrews/IMG'
#eosDir='~/work/IMG'
#eosDir='~lpcml/nobackup/IMG'
list_idx = '0000%d'%args.list_idx
label = args.label
decay = 'QCD_Pt_80_170_%s_IMG'%list_idx
#ver = ''
ver = 'v2'
print " >> Doing", decay
files = glob.glob('%s/%s%s/%s_*_label%d_list%s_*.hdf5'%(eosDir,decay,ver,decay,label,list_idx))
print " >> N files:",len(files)
dsets = [h5py.File(f) for f in files]
if list_idx == '00000':
    tgtNjets = [107778, 148990, 140182] #00000
elif list_idx == '00001':
    tgtNjets = [18136, 23770, 27747] #00001
else:
    quit() 

# runId
jetEventId_ = da.concatenate([da.from_array(dset['jetEventId'], chunks=dset['jetEventId'].chunks) for dset in dsets])
jetRunId_ = da.concatenate([da.from_array(dset['jetRunId'], chunks=dset['jetRunId'].chunks) for dset in dsets])
print " >> %s: %s"%('jetRunId', jetRunId_.shape)
X_ECAL_stacked_ = da.concatenate([da.from_array(dset['X_ECAL_stacked'], chunks=dset['X_ECAL_stacked'].chunks) for dset in dsets])
y_jets_ = da.concatenate([da.from_array(dset['y_jets'], chunks=dset['y_jets'].chunks) for dset in dsets])
print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked_.shape)

runs = [194533, 200519, 206859]
for i,r in enumerate(runs):

    if i < args.run:
        continue

    print " >> Run[%d]: %d"%(i,r)
    runMask = (jetRunId_==r).compute()
    nJets = len(runMask[runMask==True])
    print " >> nJets: %d -> %d"%(nJets, tgtNjets[i])
    nJets = tgtNjets[i]
    jetEventId = jetEventId_[runMask][:nJets]
    print " >> %s: %s"%('jetEventId', jetEventId.shape)
    jetRunId = jetRunId_[runMask][:nJets]
    print " >> %s: %s"%('jetRunId', jetRunId.shape)
    X_ECAL_stacked = X_ECAL_stacked_[runMask][:nJets]
    print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)
    y_jets = y_jets_[runMask][:nJets]
    print " >> %s: %s"%('y_jets', y_jets.shape)

    #file_out_str = "test_jets.hdf5"
    file_out_str = "%s/%s%s/%s_n%d_label%d_run%d.hdf5"%(eosDir, decay, ver, decay, nJets, label, i)
    #file_out_str = "%s/%s/%s_n%d_label%d_jet%d_run%d.hdf5"%(eosDir, decay, decay, nJets, label, ijet, i)
    print " >> Writing to:", file_out_str
    da.to_hdf5(file_out_str, {
                            #'runId': runId,
                            #'lumiId': lumiId,
                            #'eventId': eventId,
                            'X_ECAL_stacked': X_ECAL_stacked,
                            #'y': y,
                            'jetRunId': jetRunId,
                            'jetEventId': jetEventId,
                            #'jetSeed_iphi': jetSeed_iphi,
                            #'jetSeed_ieta': jetSeed_ieta,
                            #'jetM': jetM,
                            #'jetPt': jetPt,
                            #'X_jets': X_jets,
                            'y_jets': y_jets 
                            }, compression='lzf')

    print " >> Done.\n"
