import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da

#eosDir='/eos/uscms/store/user/mba2012/IMGs/HighLumi_ROOTv2'
eosDir='/eos/uscms/store/user/mba2012/IMGs/h24gamma_eta14'
#decays = ['h22gammaSM_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
#decays = ['SM2gamma_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
#decays = ['SM2gamma_1j_1M_noPU', 'h22gammaSM_1j_1M_noPU']
#decays = ['SM2gamma_1j_1M_noPU', 'h22gammaSM_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
decays = ['dummy']

chunk_size_ = 250
scale = 1.

@delayed
def load_X(tree, start_, stop_, branches_, readouts, scale):
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

    # Rescale
    X /= scale 
    return X

@delayed
def load_single(tree, start_, stop_, branches_):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    if len(branches_) > 1:
      X = np.array([np.concatenate(x).reshape(len(branches_),1) for x in X])
      X = X.reshape((-1,len(branches_)))
    else:
      X = np.array([x[0] for x in X])

    return X

for j,decay in enumerate(decays):

    if j == 0 or j == 1:
        pass
        #continue

    tfile_str = 'output_numEvent5.root'
    #tfile_str = '%s/%s_FEVTDEBUG_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_nXXX_IMG.root'%(eosDir,decay)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    #neff = (nevts//1000)*1000
    neff = int(nevts)
    #neff = 170000
    #chunk_size = chunk_size_
    chunk_size = int(nevts)
    if neff > nevts:
        neff = int(nevts)
        chunk_size = int(nevts)
    #neff = 1000 
    #neff = 233000
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    # EB
    readouts = [170,360]
    branches = ["EB_energy"]
    X = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X.shape

    # SC0
    readouts = [32,32]
    branches = ["SC_energy0"]
    X_crop0 = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_crop0.shape

    ## SC1
    #readouts = [32,32]
    #branches = ["SC_energy1"]
    #X_crop1 = da.concatenate([\
    #            da.from_delayed(\
    #                load_X(tree,i,i+chunk_size, branches, readouts, scale),\
    #                shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
    #                dtype=np.float32)\
    #            for i in range(0,neff,chunk_size)])
    #print " >> Expected shape:", X_crop1.shape

    # pho_pT0 
    branches = ["pho_pT0"]
    pho_pT0 = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", pho_pT0.shape

    # pho_E0 
    branches = ["pho_E0"]
    pho_E0 = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", pho_E0.shape

    # pho_eta0 
    branches = ["pho_eta0"]
    pho_eta0 = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", pho_eta0.shape

    # eventId
    branches = ["eventId"]
    eventId = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.int32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", eventId.shape

    ## Kinematics
    #branches = ["pho_pT", "pho_E", "pho_eta", "pho_phi"]
    #X_p4 = da.concatenate([\
    #            da.from_delayed(\
    #                load_single(tree,i,i+chunk_size, branches),\
    #                shape=(chunk_size,len(branches)),\
    #                dtype=np.float32)\
    #            for i in range(0,neff,chunk_size)])
    #print " >> Expected shape:", X_p4.shape

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(X.shape[0], label, dtype=np.float32),\
            chunks=(chunk_size,))

    #file_out_str = "%s/%s_IMG_RH%d_n%dk_label%d.hdf5"%(eosDir,decay,int(scale),neff//1000.,label)
    file_out_str = "test.hdf5"
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y, 'eventId': eventId, 'X_crop0': X_crop0, 'X_crop1': X_crop1}, compression='lzf')
    da.to_hdf5(file_out_str, {'/X': X, '/y': y,
                              'eventId': eventId,
                              'X_crop0': X_crop0,
                              #'X_crop1': X_crop1
                              #'X_p4': X_p4
                              'pho_pT0': pho_pT0,
                              'pho_E0': pho_E0,
                              'pho_eta0': pho_eta0
                              }, compression='lzf')

    print " >> Done.\n"
