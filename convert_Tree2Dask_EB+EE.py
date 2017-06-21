import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da

eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs_RAW'
decays = ["DoublePhotonGaussPt55_StdDev20_FEVTDEBUG_HighLumiPileUp"]
#decays = ["H125GGgluonfusion_13TeV_TuneCUETP8M1_HighLumiPileUp"]
#decays = ["SinglePhotonPt50","SingleElectronPt50"]
#decays = ["DoublePhotonFlatPt10To60","DoubleElectronFlatPt10To60"]
chunk_size = 250

@delayed
def load_X(tree, start_, stop_, branches_, readouts):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    # Convert the object array X to a multidim array:
    # 1: for each event x in X, concatenate the object columns (branches) into a flat array of shape (readouts*branches)
    # 2: reshape the flat array into a stacked array: (branches, readouts)
    # 3: embed each stacked array as a single row entry in a list via list comprehension
    # 4: convert this list into an array with shape (events, branches, readouts) 
    X = np.array([np.concatenate(x).reshape(len(branches_),readouts[0]*readouts[1]) for x in X])
    #print "X.shape:",X.shape
    X = X.reshape((-1,len(branches_),readouts[0],readouts[1]))
    X = np.transpose(X, [0,2,3,1])
    X /= 1000. 
    return X

for j,decay in enumerate(decays):

    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_CROPS32.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_pT_CROPS32.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG.root'%(eosDir,decay)
    tfile_str = '%s/%s_n250k_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n350k_IMG.root'%(eosDir,decay)
    #tfile_str = 'output_n10.root'
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    neff = (nevts//chunk_size)*chunk_size
    neff = 250000 
    #neff = 233000
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    # EB
    readouts = [170,360]
    branches = ["EB_energy"]
    X_EB = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EB.shape

    # EE-
    readouts = [100,100]
    branches = ["EEm_energy"]
    X_EEm = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EEm.shape

    # EE+
    readouts = [100,100]
    branches = ["EEp_energy"]
    X_EEp = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EEp.shape

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(X_EB.shape[0], label, dtype=np.float32),\
            chunks=(chunk_size,))

    #file_out_str = "test.hdf5"
    #file_out_str = "%s/%s_IMG_RHraw_n%dk.hdf5"%(eosDir,decay,neff//1000.)
    file_out_str = "%s/%s_IMG_RHv1_n%dk.hdf5"%(eosDir,decay,neff//1000.)
    #file_out_str = "%s/%s_n%dk_IMG_RHraw.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk_pT.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk_DIGI.hdf5"%(eosDir,decay,neff//1000)
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y}, chunks=(chunk_size,s,s,2), compression='lzf')
    da.to_hdf5(file_out_str, {'/X_EB': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, '/y': y}, compression='lzf')

    print " >> Done.\n"
