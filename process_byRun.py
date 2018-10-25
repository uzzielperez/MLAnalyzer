import h5py
import numpy as np
from dask.delayed import delayed
import dask.array as da
import glob

#eosDir='/eos/uscms/store/user/mba2012/IMGs'
eosDir='/eos/uscms/store/user/mba2012/IMG'
#decay = 'QCD_Pt_80_170_00000_IMGjet_RH1'
list_idx = '00001'
decay = 'QCD_Pt_80_170_%s_IMGjet'%list_idx
label = 0
#label = 1
files = glob.glob('%s/%s/%s_*_label%d_*_list%s_*.hdf5'%(eosDir,decay,decay,label,list_idx))
print " >> N files:",len(files)
#files = glob.glob('%s/%s/%s_RH1_n*_label%d.hdf5'%(eosDir,decay,decay,label))
dsets = [h5py.File(f) for f in files]
#tgtNjets = [213626, 294986, 276988]
#tgtNjets = [107778, 148990, 140182]
tgtNjets = [18136, 23770, 27747]

# runId
jetEventId_ = da.concatenate([da.from_array(dset['jetEventId'], chunks=dset['jetEventId'].chunks) for dset in dsets])
jetRunId_ = da.concatenate([da.from_array(dset['jetRunId'], chunks=dset['jetRunId'].chunks) for dset in dsets])
print " >> %s: %s"%('jetRunId', jetRunId_.shape)
jetM_ = da.concatenate([da.from_array(dset['jetM'], chunks=dset['jetM'].chunks) for dset in dsets])
jetPt_ = da.concatenate([da.from_array(dset['jetPt'], chunks=dset['jetPt'].chunks) for dset in dsets])
X_jets_ = da.concatenate([da.from_array(dset['X_jets'], chunks=dset['X_jets'].chunks) for dset in dsets])
#X_ECAL_stacked_ = da.concatenate([da.from_array(dset['X_ECAL_stacked'], chunks=dset['X_ECAL_stacked'].chunks) for dset in dsets])
y_jets_ = da.concatenate([da.from_array(dset['y_jets'], chunks=dset['y_jets'].chunks) for dset in dsets])
print " >> %s: %s"%('X_jets', X_jets_.shape)
#print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked_.shape)

runs = [194533, 200519, 206859]
for i,r in enumerate(runs):

  print " >> Run[%d]: %d"%(i,r)
  runMask = (jetRunId_==r).compute()
  nJets = len(runMask[runMask==True])
  print " >> nJets: %d -> %d"%(nJets, tgtNjets[i])
  nJets = tgtNjets[i]
  jetEventId = jetEventId_[runMask][:nJets]
  print " >> %s: %s"%('jetEventId', jetEventId.shape)
  jetRunId = jetRunId_[runMask][:nJets]
  print " >> %s: %s"%('jetRunId', jetRunId.shape)
  jetM = jetM_[runMask][:nJets]
  print " >> %s: %s"%('jetM', jetM.shape)
  jetPt = jetPt_[runMask][:nJets]
  print " >> %s: %s"%('jetPt', jetPt.shape)
  X_jets = X_jets_[runMask][:nJets]
  print " >> %s: %s"%('X_jets', X_jets.shape)
  #X_ECAL_stacked = X_ECAL_stacked_[runMask][:nJets]
  #print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)
  y_jets = y_jets_[runMask][:nJets]
  print " >> %s: %s"%('y_jets', y_jets.shape)

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

  #file_out_str = "test_jets.hdf5"
  file_out_str = "%s/%s/%s_n%d_label%d_run%d.hdf5"%(eosDir, decay, decay, nJets, label, i)
  print " >> Writing to:", file_out_str
  da.to_hdf5(file_out_str, {
                            #'runId': runId,
                            #'lumiId': lumiId,
                            #'eventId': eventId,
                            #'X_ECAL_stacked': X_ECAL_stacked,
                            #'y': y,
                            'jetRunId': jetRunId,
                            'jetEventId': jetEventId,
                            #'jetSeed_iphi': jetSeed_iphi,
                            #'jetSeed_ieta': jetSeed_ieta,
                            'jetM': jetM,
                            'jetPt': jetPt,
                            'X_jets': X_jets,
                            'y_jets': y_jets 
                            }, compression='lzf')

  print " >> Done.\n"
