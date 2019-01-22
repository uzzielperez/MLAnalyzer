import h5py
import numpy as np
from dask.delayed import delayed
import dask.array as da
import glob

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-l', '--label', required=True, type=int, help='Decay label.')
parser.add_argument('-f', '--list_idx', required=True, type=int, help='List index, 0: train, 1:test')
args = parser.parse_args()

#eosDir='/eos/uscms/store/user/mba2012/IMG'
eosDir='/eos/cms/store/user/mandrews/IMG'
#eosDir='~/work/IMG'
#eosDir='~lpcml/nobackup/IMG'
list_idx = '0000%d'%args.list_idx
label = args.label
decay = 'QCD_Pt_80_170_%s_IMGjet'%list_idx
if list_idx == '00000':
    tgtNjets = [107778, 148990, 140182] #00000
elif list_idx == '00001':
    tgtNjets = [18136, 23770, 27747] #00001
else:
    quit() 

ijet = 0 
files = glob.glob('%s/%s/%s_*_label%d_jet%d_list%s_*.hdf5'%(eosDir,decay,decay,label,ijet,list_idx))
print " >> N files:",len(files)
dsets0 = [h5py.File(f) for f in files]

# runId
jetEventId0_ = da.concatenate([da.from_array(dset['jetEventId'], chunks=dset['jetEventId'].chunks) for dset in dsets0])
jetRunId0_ = da.concatenate([da.from_array(dset['jetRunId'], chunks=dset['jetRunId'].chunks) for dset in dsets0])
print " >> %s: %s"%('jetRunId', jetRunId0_.shape)
#jetM_ = da.concatenate([da.from_array(dset['jetM'], chunks=dset['jetM'].chunks) for dset in dsets0])
#jetPt_ = da.concatenate([da.from_array(dset['jetPt'], chunks=dset['jetPt'].chunks) for dset in dsets0])
X_jets0_ = da.concatenate([da.from_array(dset['X_jets'], chunks=dset['X_jets'].chunks) for dset in dsets0])
X_FC_ = da.concatenate([da.from_array(dset['X_FC'], chunks=dset['X_FC'].chunks) for dset in dsets0])
#X_ECAL_stacked_ = da.concatenate([da.from_array(dset['X_ECAL_stacked'], chunks=dset['X_ECAL_stacked'].chunks) for dset in dsets0])
y_jets0_ = da.concatenate([da.from_array(dset['y_jets'], chunks=dset['y_jets'].chunks) for dset in dsets0])
print " >> %s: %s"%('X_jets0', X_jets0_.shape)
#print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked_.shape)

ijet = 1
files = glob.glob('%s/%s/%s_*_label%d_jet%d_list%s_*.hdf5'%(eosDir,decay,decay,label,ijet,list_idx))
print " >> N files:",len(files)
#dsets1 = [h5py.File(f) for f in files]
dsets1 = [h5py.File(f) for f in reversed(files)]
jetEventId1_ = da.concatenate([da.from_array(dset['jetEventId'], chunks=dset['jetEventId'].chunks) for dset in dsets1])
jetRunId1_ = da.concatenate([da.from_array(dset['jetRunId'], chunks=dset['jetRunId'].chunks) for dset in dsets1])
X_jets1_ = da.concatenate([da.from_array(dset['X_jets'], chunks=dset['X_jets'].chunks) for dset in dsets1])
y_jets1_ = da.concatenate([da.from_array(dset['y_jets'], chunks=dset['y_jets'].chunks) for dset in dsets1])
print " >> %s: %s"%('X_jets1', X_jets1_.shape)

runs = [194533, 200519, 206859]
for i,r in enumerate(runs):

    print " >> Run[%d]: %d"%(i,r)
    runMask0 = (jetRunId0_==r).compute()
    runMask1 = (jetRunId1_==r).compute()
    nJets = len(runMask0[runMask0==True])
    print " >> nJets: %d -> %d"%(nJets, tgtNjets[i])
    nJets = tgtNjets[i]
    jetEventId0 = jetEventId0_[runMask0][:nJets]
    print " >> %s: %s"%('jetEventId0', jetEventId0.shape)
    jetEventId1 = jetEventId1_[runMask1][:nJets]
    print " >> %s: %s"%('jetEventId1', jetEventId1.shape)
    jetRunId0 = jetRunId0_[runMask0][:nJets]
    print " >> %s: %s"%('jetRunId0', jetRunId0.shape)
    jetRunId1 = jetRunId1_[runMask1][:nJets]
    print " >> %s: %s"%('jetRunId1', jetRunId1.shape)
    #jetM = jetM_[runMask0][:nJets]
    #print " >> %s: %s"%('jetM', jetM.shape)
    #jetPt = jetPt_[runMask0][:nJets]
    #print " >> %s: %s"%('jetPt', jetPt.shape)
    X_jets0 = X_jets0_[runMask0][:nJets]
    print " >> %s: %s"%('X_jets0', X_jets0.shape)
    X_jets1 = X_jets1_[runMask1][:nJets]
    print " >> %s: %s"%('X_jets1', X_jets1.shape)
    #X_FC = X_FC_[runMask0][:nJets]
    #print " >> %s: %s"%('X_FC', X_FC.shape)
    #X_ECAL_stacked = X_ECAL_stacked_[runMask0][:nJets]
    #print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)
    y_jets0 = y_jets0_[runMask0][:nJets]
    print " >> %s: %s"%('y_jets0', y_jets0.shape)
    y_jets1 = y_jets1_[runMask1][:nJets]
    print " >> %s: %s"%('y_jets1', y_jets1.shape)

    assert da.all(jetEventId0 != jetEventId1)
    assert da.all(jetRunId0 == jetRunId1)
    assert da.all(y_jets0 == y_jets1)

    #file_out_str = "test_jets.hdf5"
    file_out_str = "%s/%s/%s_n%d_label%d_jetcombo-rev_run%d.hdf5"%(eosDir, decay, decay, nJets, label, i)
    print " >> Writing to:", file_out_str
    da.to_hdf5(file_out_str, {
                            #'runId': runId,
                            #'lumiId': lumiId,
                            #'eventId': eventId,
                            #'X_ECAL_stacked': X_ECAL_stacked,
                            #'y': y,
                            'jetRunId0': jetRunId0,
                            'jetRunId1': jetRunId1,
                            'jetEventId0': jetEventId0,
                            'jetEventId1': jetEventId1,
                            #'jetSeed_iphi': jetSeed_iphi,
                            #'jetSeed_ieta': jetSeed_ieta,
                            #'jetM': jetM,
                            #'jetPt': jetPt,
                            #'X_FC': X_FC,
                            'X_jets0': X_jets0,
                            'X_jets1': X_jets1,
                            'y_jets': y_jets0 
                            }, compression='lzf')

    print " >> Done.\n"
