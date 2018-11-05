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
decay = 'QCD_Pt_80_170_%s'%list_idx
print " >> Doing", decay
files = glob.glob('%s/%s_FC/%s_label%d_ijet0_*.hdf5'%(eosDir,decay,decay,label))
print " >> N files:",len(files)
dsets = [h5py.File(f) for f in files]
if list_idx == '00000':
    tgtNjets = [107778, 148990, 140182] #00000
elif list_idx == '00001':
    tgtNjets = [18136, 23770, 27747] #00001
else:
    quit()

# runId
jetEventId_ = np.concatenate([dset['jetEventId'][:] for dset in dsets])
jetRunId_ = np.concatenate([dset['jetRunId'][:] for dset in dsets])
X_ = np.concatenate([dset['events'][:] for dset in dsets])
y_ = np.concatenate([dset['y'][:] for dset in dsets])
print " >> %s: %s"%('X', X_.shape)

runs = [194533, 200519, 206859]
for i,r in enumerate(runs):

    print " >> Run[%d]: %d"%(i,r)
    runMask = (jetRunId_==r)
    nJets = len(runMask[runMask==True])
    print " >> nJets: %d -> %d"%(nJets, tgtNjets[i])

    nJets = tgtNjets[i]
    jetEventId = jetEventId_[runMask][:nJets]
    print " >> %s: %s"%('jetEventId', jetEventId.shape)
    jetRunId = jetRunId_[runMask][:nJets]
    print " >> %s: %s"%('jetRunId', jetRunId.shape)
    #jetM = jetM_[runMask][:nJets]
    #print " >> %s: %s"%('jetM', jetM.shape)
    #jetPt = jetPt_[runMask][:nJets]
    #print " >> %s: %s"%('jetPt', jetPt.shape)
    #X_jets = X_jets_[runMask][:nJets]
    #print " >> %s: %s"%('X_jets', X_jets.shape)
    X = X_[runMask][:nJets]
    print " >> %s: %s"%('X', X.shape)
    y = y_[runMask][:nJets]
    print " >> %s: %s"%('y', y.shape)

    neff = nJets
    file_out_str = "%s/%s_FC/%s_label%d_run%d.hdf5"%(eosDir, decay, decay, label, i)
    f = h5py.File(file_out_str)
    rowtype = np.dtype('f4,f4,f4,f4') # np.float32
    dt = h5py.special_dtype(vlen=np.dtype(rowtype))
    dset = {}
    dset['X'] = f.create_dataset('events', (neff,), dtype=dt)
    dset['y'] = f.create_dataset('y', (neff,), dtype='f4')
    dset['jetEventId'] = f.create_dataset('jetEventId', (neff,), dtype=int)
    dset['jetRunId'] = f.create_dataset('jetRunId', (neff,), dtype=int)
    count = 0
    for j in range(nJets):
        if j%10000 == 0:
            pass
            print j
            #print(jetRunId[j], r)
        dset['X'][j] = X[j]
        dset['y'][j] = y[j]
        dset['jetEventId'][j] = jetEventId[j]
        dset['jetRunId'][j] = jetRunId[j]
        count += 1
    print('count:',count)

    print " >> Writing to:", file_out_str
