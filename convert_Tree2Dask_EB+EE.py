import numpy as np
import ROOT
from root_numpy import root2array, tree2array
#from root_pandas.readwrite import convert_to_dataframe
from dask.delayed import delayed
import dask.array as da
import dask.dataframe as df

eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs_RAW'
decays = ["H125"]
#decays = ["SinglePhotonPt50","SingleElectronPt50"]
#decays = ["DoublePhotonFlatPt10To60","DoubleElectronFlatPt10To60"]
chunk_size = 1

@delayed
def load_X(tree, start_, stop_, branch, crop_size):
    X = tree2array(tree, start=start_, stop=stop_, branches=[branch]) 
    X = np.array([np.concatenate(x).reshape(-1,crop_size) for x in X]) # converts the list of list to multidim array
    X = X.reshape((-1,1,crop_size))
    return X

for j,decay in enumerate(decays):

    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_CROPS32.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_pT_CROPS32.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG.root'%(eosDir,decay)
    tfile_str = 'output_EEtest3.root'
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    neff = (nevts//chunk_size)*chunk_size
    #neff = 100 
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    X_EB = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, "EB_energy", 170*360),\
                    shape=(chunk_size,1,170*360),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EB.shape

    X_EEm = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size,"EEm_energy", 100*100),\
                    shape=(chunk_size,1,100*100),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EEm.shape

    print " >> Class label:",j
    y = da.from_array(\
            np.full(X_EB.shape[0], j, dtype=np.float32),\
            chunks=(chunk_size,))

    file_out_str = "test.hdf5"
    #file_out_str = "%s/%s_n%dk_IMG_RHraw.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk_pT.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk_DIGI.hdf5"%(eosDir,decay,neff//1000)
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y}, chunks=(chunk_size,s,s,2), compression='lzf')
    da.to_hdf5(file_out_str, {'/X_EB': X_EB, 'X_EEm': X_EEm, '/y': y}, compression='lzf')

    print " >> Done.\n"
