import h5py
import numpy as np
from dask.delayed import delayed
import dask.array as da
import glob

#eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='/eos/uscms/store/user/mba2012/IMG'
eosDir='/eos/cms/store/user/mandrews/IMG'
#decay = 'QCD_Pt_80_170_00000_IMGjet_RH1'
decay = 'QCD_Pt_80_170_00000'
#decay = 'QCD_Pt_80_170_00001'
label = 0
#label = 1
#files = glob.glob('%s/%s_*_label%d_*_list00000_*.hdf5'%(eosDir,decay,label))
#files = glob.glob('%s/%s_n*k_label%d.hdf5'%(eosDir,decay,label))
files = glob.glob('%s/%s_FC/%s_label%d_ijet0_*.hdf5'%(eosDir,decay,decay,label))
#files = glob.glob('QCD_%d.hdf5'%label)
dsets = [h5py.File(f) for f in files]
#tgtNjets = [213626, 294986, 276988]
#tgtNjets = [107778, 148937, 140182]
tgtNjets = [107778, 148990, 140182] #00000
#tgtNjets = [18136, 23770, 27747] #00001

# runId
jetEventId_ = np.concatenate([dset['jetEventId'][:] for dset in dsets])
jetRunId_ = np.concatenate([dset['jetRunId'][:] for dset in dsets])
X_ = np.concatenate([dset['events'][:] for dset in dsets])
y_ = np.concatenate([dset['y'][:] for dset in dsets])
print " >> %s: %s"%('X', X_.shape)
#print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked_.shape)

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
  #X_ECAL_stacked = X_ECAL_stacked_[runMask][:nJets]
  #print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)
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

  #jetMask = (jetIsQuark_==j).compute()
  #n_jets = (len(jetMask[jetMask==True])//1000)*1000
  #n_jets = len(jetMask[jetMask==True])
  #print " >> N jets @ label%d: %d -> %d"%(j, len(jetMask[jetMask==True]), n_jets)
  #print " >> N jets @ label%d: %d"%(j, len(jetMask[jetMask==True]))

  #jetIsQuark = jetIsQuark[jetMask][:n_jets]
  #print " >> %s: %s"%('jetIsQuark', jetIsQuark.shape)
  #jetM = jetM[jetMask][:n_jets]
  #print " >> %s: %s"%('jetM', jetM.shape)
  #jetPt = jetPt[jetMask][:n_jets]
  #print " >> %s: %s"%('jetPt', jetPt.shape)
  #jetSeed_iphi = jetSeed_iphi[jetMask][:n_jets]
  #print " >> %s: %s"%('jetSeed_iphi', jetSeed_iphi.shape)
  #jetSeed_ieta = jetSeed_ieta[jetMask][:n_jets]
  #print " >> %s: %s"%('jetSeed_ieta', jetSeed_ieta.shape)
  #X_jets = X_jets[jetMask][:n_jets]
  #print " >> %s: %s"%('X_jets', X_jets.shape)
  ##X_ECAL_stacked = X_ECAL_stacked[jetMask][:n_jets]
  ##print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)
  #jetEventId = eventId[jetMask][:n_jets]
  #jetRunId = runId[jetMask][:n_jets]

  #jetIsQuark = jetIsQuark_[jetMask]
  #print " >> %s: %s"%('jetIsQuark', jetIsQuark.shape)
  #jetM = jetM_[jetMask]
  #print " >> %s: %s"%('jetM', jetM.shape)
  #jetPt = jetPt_[jetMask]
  #print " >> %s: %s"%('jetPt', jetPt.shape)
  #jetSeed_iphi = jetSeed_iphi_[jetMask]
  #print " >> %s: %s"%('jetSeed_iphi', jetSeed_iphi.shape)
  #jetSeed_ieta = jetSeed_ieta_[jetMask]
  #print " >> %s: %s"%('jetSeed_ieta', jetSeed_ieta.shape)
  #X_jets = X_jets_[jetMask]
  #print " >> %s: %s"%('X_jets', X_jets.shape)
  #jetEventId = eventId[jetMask]
  #jetRunId = runId[jetMask]

  #file_out_str = "%s_label%d_run%d.hdf5"%(decay, label, i)
  #file_out_str = "%s/%s_list00000_n%d_label%d_run%d.hdf5"%(eosDir, decay, nJets, label, i)
  print " >> Writing to:", file_out_str
