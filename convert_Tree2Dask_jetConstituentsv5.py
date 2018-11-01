import numpy as np
import ROOT
from root_numpy import tree2array, root2array
import glob, os
import h5py

#eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='/eos/uscms/store/user/mba2012/IMG'
eosDir='/eos/cms/store/user/mandrews/IMG'
#eosDir='/uscms/home/mba2012/nobackup'
#eosDir='~/work/MLHEP/CMSSW_8_0_26_patch1/src/ggAnalysis/ggNtuplizer/test'
#eosDir='/eos/cms/store/user/mandrews/OPENDATA/IMGs/MGG90_Eta23'
#decays = ['QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU', 'QCDToQQ_Pt_80_120_13TeV_TuneCUETP8M1_noPU']
#decays = ['dummy']
list_idx = '00000'
#list_idx = '00001'
date_str = '181029_120833'
#date_str = '181025_233834'
decays = ['QCD_Pt_80_170_%s'%list_idx, 'QCD_Pt_80_170_%s'%list_idx]

def load_X(tree, start_, stop_, branches_):
  X = tree2array(tree, start=start_, stop=stop_, branches=branches_)
  #X = root2array(tree, treename='fevt/RHTree', start=start_, stop=stop_, branches=branches_)
  X = [np.concatenate(x).reshape(4,-1).T.astype(np.float32) for x in X]
  return X

def load_single(tree, start_, stop_, branches_):
  X = tree2array(tree, start=start_, stop=stop_, branches=branches_)
  #X = root2array(tree, treename='fevt/RHTree', start=start_, stop=stop_, branches=branches_)
  X = np.array([x[0] for x in X])
  return X

for d,decay in enumerate(decays):

  if d == 0:
    pass
    #continue
  print ' >> decay:',d

  tfile_idxs = glob.glob('%s/%s*_AODSIM/%s/*/output_*.root'%(eosDir,decay,date_str))
  tfile_idxs = [s.replace('.root','').split('_')[-1] for s in tfile_idxs]
  print tfile_idxs

  for n in tfile_idxs:

    #tfile_str = 'output.root'
    #tfile_str = '%s/%s_AODSIM/180919_131137/*/output_1.root'%(eosDir, decay)
    tfile_str = glob.glob('%s/%s*_AODSIM/%s/*/output_%s.root'%(eosDir,decay,date_str,n))[0]
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    print " >> Input file:", tfile_str

    #tfiles = glob.glob('%s/%s*_AODSIM/180919_131137/*/output_*.root'%(eosDir,decay))
    #print " >> %d files found."%len(tfiles)
    #tree = ROOT.TChain("fevt/RHTree")
    #for f in tfiles:
    #  tree.Add(f)
    #nevts = tree.GetEntries()
    #tree = tfiles
    #print " >> Input file[0]:", tfiles[0]

    #neff = (nevts//1000)*1000
    #neff = (nevts//250)*250
    neff = int(nevts)
    print " >> Doing decay:", decay
    print " >> Label:", d
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    for ijet in [0]: # only keep leading jet

      print ' >> ijet:',ijet

      branches = ['subJet%d_E'%ijet, 'subJet%d_Px'%ijet, 'subJet%d_Py'%ijet, 'subJet%d_Pz'%ijet]
      X = load_X(tree, 0, neff, branches)

      branches = ["jetIsQuark%d"%ijet]
      jetIsQuark = load_single(tree, 0, neff, branches)

      branches = ["runId"]
      jetRunId = load_single(tree, 0, neff, branches)

      branches = ["eventId"]
      jetEventId = load_single(tree, 0, neff, branches)

      n_jets = len(jetIsQuark[jetIsQuark==d])
      print '  >> n_jets:', n_jets

      # Initialize output file
      if not os.path.isdir("%s/%s_FC"%(eosDir,decay)):
        os.makedirs("%s/%s_FC"%(eosDir,decay))
      file_out_str = '%s/%s_FC/%s_label%d_ijet%d_%s.hdf5'%(eosDir, decay, decay, d, ijet, n)
      print "  >> Writing to:", file_out_str
      f = h5py.File(file_out_str)
      rowtype = np.dtype('f4,f4,f4,f4') # np.float32
      dt = h5py.special_dtype(vlen=np.dtype(rowtype))
      dset = {}
      dset['X'] = f.create_dataset('events', (n_jets,), dtype=dt)
      dset['y'] = f.create_dataset('y', (n_jets,), dtype='f4')
      dset['jetEventId'] = f.create_dataset('jetEventId', (n_jets,), dtype=int)
      dset['jetRunId'] = f.create_dataset('jetRunId', (n_jets,), dtype=int)

      count = 0
      for i in range(neff):

        if i%10000 == 0:
          pass
          #print '  >> %d/%d'%(i,neff)
        if jetIsQuark[i] != d:
          continue

        # Get jet constituents
        nconst = len(X[i])
        #print "N constituents:", nconst 
        x_p4 = np.zeros((nconst,), dtype=rowtype)
        for k in range(nconst):
          x_p4[k][0] = X[i][k][0]
          x_p4[k][1] = X[i][k][1]
          x_p4[k][2] = X[i][k][2]
          x_p4[k][3] = X[i][k][3]

        # Write to output file
        dset['X'][count] = x_p4 
        dset['y'][count] = jetIsQuark[i] 
        dset['jetEventId'][count] = jetEventId[i]
        dset['jetRunId'][count] = jetRunId[i]
        count += 1

      print '  >> count:', count
      #print f['events'][0]
      #print f['y'][0]
