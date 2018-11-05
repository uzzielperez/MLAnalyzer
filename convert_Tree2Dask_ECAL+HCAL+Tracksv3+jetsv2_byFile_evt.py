import numpy as np
import ROOT
from root_numpy import tree2array, root2array
from dask.delayed import delayed
import dask.array as da
from skimage.measure import block_reduce
import glob, os, re

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-l', '--label', required=True, type=int, help='Decay label.')
#parser.add_argument('-n', '--file_start_stop', default=[], nargs='+', type=int, help='List of file indices.')
parser.add_argument('-n', '--file_idx_start', default=1, type=int, help='File index start.')
args = parser.parse_args()

#eosDir='/eos/uscms/store/user/mba2012/IMG'
eosDir='/eos/cms/store/user/mandrews/IMG'
outDir='~/work/IMG'
#xrootd='root://cmsxrootd.fnal.gov' # FNAL
xrootd='root://eoscms.cern.ch' # CERN
#list_idx = '00000'
list_idx = '00001'
#date_str = '181029_120833'
date_str = '181025_233834'
decays = ['QCD_Pt_80_170_%s'%list_idx, 'QCD_Pt_80_170_%s'%list_idx]

#chunk_size_ = 250
chunk_size_ = 200
#chunk_size_ = 250
#scale = [100., 150.]
scale = [1., 1.]
jet_shape = 125

@delayed
def load_X(tree, start_, stop_, branches_, readouts, scale):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    #X = root2array(tree, treename='fevt/RHTree', start=start_, stop=stop_, branches=branches_)
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
    #X = root2array(tree, treename='fevt/RHTree', start=start_, stop=stop_, branches=branches_)
    X = np.array([x[0] for x in X])

    return X

@delayed
def load_single_bool(tree, start_, stop_, branches_):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_)
    #X = root2array(tree, treename='fevt/RHTree', start=start_, stop=stop_, branches=branches_)
    X = np.array([x[0] for x in X]).astype(bool)

    return X

@delayed
def load_X_upsampled(tree, start_, stop_, branches_, readouts, scale, upscale):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    #X = root2array(tree, treename='fevt/RHTree', start=start_, stop=stop_, branches=branches_)
    # Convert the object array X to a multidim array:
    # 1: for each event x in X, concatenate the object columns (branches) into a flat array of shape (readouts*branches)
    # 2: reshape the flat array into a stacked array: (branches, readouts)
    # 3: embed each stacked array as a single row entry in a list via list comprehension
    # 4: convert this list into an array with shape (events, branches, readouts) 
    X = np.array([np.concatenate(x).reshape(len(branches_),readouts[0]*readouts[1]) for x in X])
    #print "X.shape:",X.shape
    X = X.reshape((-1,len(branches_),readouts[0],readouts[1]))

    #print "unsampled.shape",X.shape
    X = np.stack([tile_stacked_array(x, upscale) for x in X])
    #print "upsampled.shape",X.shape
    X = np.transpose(X, [0,2,3,1])

    # Rescale
    X /= scale 
    return X

from numpy.lib.stride_tricks import as_strided

def tile_stacked_array(X, upscale):
    #print "un-tile_stacked.shape",X.shape
    X = np.stack([tile_array(x, upscale, upscale) for x in X])
    #print "tile_stacked.shape",X.shape
    return X
    
def tile_array(x, b0, b1):
    r, c = x.shape                                    # number of rows/columns
    rs, cs = x.strides                                # row/column strides 
    x = as_strided(x, (r, b0, c, b1), (rs, 0, cs, 0)) # view a as larger 4D array
    #return x.reshape(r*b0, c*b1)                      # create new 2D array
    return x.reshape(r*b0, c*b1)/(b0*b1)              # create new 2D array, conserve energy

def block_resample_EE(X):

    return np.array([resample_EE(x) for x in X])

def resample_EE(imgECAL, factor=2):
    
    imgECAL = np.squeeze(imgECAL)
    #print('imgECAL.shape:',imgECAL.shape)
    
    # EE-
    imgEEm = imgECAL[:140-85] # EE- in the first 55 rows
    imgEEm = np.pad(imgEEm, ((1,0),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEm_dn = block_reduce(imgEEm, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEm_dn_up = tile_array(imgEEm_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor**2
    imgECAL[:140-85] = imgEEm_dn_up[1:] ## replace the old EE- rows
    
    # EE+
    imgEEp = imgECAL[140+85:] # EE+ in the last 55 rows
    imgEEp = np.pad(imgEEp, ((0,1),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEp_dn = block_reduce(imgEEp, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEp_dn_up = tile_array(imgEEp_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor*factor
    imgECAL[140+85:] = imgEEp_dn_up[:-1] # replace the old EE+ rows
    
    return np.expand_dims(imgECAL, -1)

def crop_jet(imgECAL, iphi, ieta):
    iphi = int(iphi*5 + 2)
    ieta = int(ieta*5 + 2)
    off = jet_shape//2
    #print('iphi:%d, ieta:%d'%(iphi, ieta))
    if iphi < off:
        #print('1')
        diff = off-iphi
        #diff = diff+150
        #print(diff)
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,-diff:],
                                   imgECAL[ieta-off:ieta+off+1,:iphi+off+1]), axis=1)
    elif 360-iphi < off:
        #print('2')
        diff = off - (360-iphi)
        #diff = diff+150
        #print('diff:',diff)
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,iphi-off:], 
                                   imgECAL[ieta-off:ieta+off+1,:diff+1]), axis=1)
    else:
        #print('0')
        img_crop = imgECAL[ieta-off:ieta+off+1,iphi-off:iphi+off+1]
    #print(img_crop.shape)
    return img_crop

@delayed
def crop_jet_block(Xs, iphis, ietas):
    return np.array([crop_jet(x,iphi,ieta) for x,iphi,ieta in zip(Xs,iphis,ietas)])

def get_chunk_size(i, neff, chunk_size):
    #print neff, i
    if neff-i < chunk_size:
        eff_chunk_size = neff - i
        #eff_chunk_size = chunk_size
    else:
        eff_chunk_size = chunk_size
        #print eff_chunk_size
    return eff_chunk_size

for d,decay in enumerate(decays):

    if d != args.label:
        pass
        continue
    print ' >> decay:',d

    tfile_idxs = glob.glob('%s/%s*_AODSIM/%s/*/output_*.root'%(eosDir,decay,date_str))
    tfile_idxs = [s.replace('.root','').split('_')[-1] for s in tfile_idxs]
    tfile_idxs = [int(i) for i in tfile_idxs]
    tfile_idxs.sort()
    print tfile_idxs

    #for n in tfile_idxs[:1]:
    for n in tfile_idxs:

        if n < args.file_idx_start:
          pass
          continue

        #tfile_str = 'output.root'
        #tfile_str = '%s/%s_AODSIM/180919_131137/*/output_1.root'%(eosDir, decay)
        tfile_str = glob.glob('%s/%s*_AODSIM/%s/*/output_%d.root'%(eosDir,decay,date_str,n))
        assert len(tfile_str) == 1, "More than 1 file of same name found in different dirs: %s"%tfile_str
        tfile_str = '%s/%s'%(xrootd,tfile_str[0])
        print " >> For input file:", tfile_str
        tfile = ROOT.TFile(tfile_str)
        tree = tfile.Get('fevt/RHTree')
        nevts = tree.GetEntries()

        #tfiles = glob.glob('%s/%s*_AODSIM/180919_131137/*/output_*.root'%(eosDir,decay))
        #print " >> %d files found."%len(tfiles)
        #tree = ROOT.TChain("fevt/RHTree")
        #for f in tfiles:
        #  tree.Add(f)
        #nevts = tree.GetEntries()
        #tree = tfiles
        #print " >> Input file[0]:", tfiles[0]

        neff = int(nevts)
        #neff = (nevts//1000)*1000
        #neff = (nevts//chunk_size_)*chunk_size_
        #neff = 84600
        #neff = 175000
        #neff = 135600
        #neff = 500
        #neff = int(nevts)
        chunk_size = chunk_size_
        if neff < chunk_size:
            chunk_size = neff
        if neff > nevts:
            neff = int(nevts)
            chunk_size = int(nevts)
        proc_range = range(0,neff,chunk_size)
        print " >> Doing decay:", decay
        print " >> Label:", d
        print " >> Total events:", nevts
        print " >> Effective events:", neff

        # eventId
        #branches = ["event"]
        branches = ["eventId"]
        eventId = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i, neff, chunk_size), branches),\
                        shape=(get_chunk_size(i, neff, chunk_size),),\
                        dtype=np.int32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],eventId.shape)

        # runId
        branches = ["runId"]
        runId = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i, neff, chunk_size), branches),\
                        shape=(get_chunk_size(i, neff, chunk_size),),\
                        dtype=np.int32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],runId.shape)

        # ECAL
        readouts = [280,360]
        branches = ["ECAL_energy"]
        X_ECAL = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i, neff, chunk_size), branches, readouts, scale[0]),\
                        shape=(get_chunk_size(i, neff, chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        #print " >> %s: %s"%(branches[0],X_ECAL.shape)

        # ECAL with resampled EE
        X_ECAL_EEup = X_ECAL.map_blocks(lambda x: block_resample_EE(x), dtype=np.float32)
        print " >> %s: %s"%('ECAL_EEup_energy',X_ECAL_EEup.shape)

        # Tracks at ECAL
        readouts = [280,360]
        branches = ["ECAL_tracksPt"]
        X_TracksAtECAL = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i, neff, chunk_size), branches, readouts, scale[0]),\
                        shape=(get_chunk_size(i, neff, chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        #print " >> %s: %s"%(branches[0],X_ECAL.shape)

        # HBHE upsample
        readouts = [56,72]
        branches = ["HBHE_energy"]
        upscale = 5
        X_HBHE_up = da.concatenate([\
                    da.from_delayed(\
                        load_X_upsampled(tree,i,i+get_chunk_size(i, neff, chunk_size), branches, readouts, scale[1], upscale),\
                        shape=(get_chunk_size(i, neff, chunk_size), readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        #print " >> %s(upsampled): %s"%(branches[0],X_HBHE_up.shape)

        X_ECAL_stacked_ = da.concatenate([X_TracksAtECAL, X_ECAL_EEup, X_HBHE_up], axis=-1)
        print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked_.shape)

        # Class label
        label = d
        #label = 1
        print " >> Class label:",label
        y = da.from_array(\
                np.full(len(eventId), label, dtype=np.float32),\
                chunks=eventId.chunks)

        #for ijet in [0, 1]: # only keep leading jet
        for ijet in [0]: # only keep leading jet

            print ' >> ijet:',ijet

            # jet is quark
            branches = ["jetIsQuark%d"%ijet]
            jetIsQuark_ = da.concatenate([\
                        da.from_delayed(\
                            load_single(tree,i,i+get_chunk_size(i, neff, chunk_size), branches),\
                            shape=(get_chunk_size(i, neff, chunk_size),),\
                            dtype=np.float32)\
                        for i in proc_range])

            jetMask = (jetIsQuark_==d).compute()
            #n_jets = (len(jetMask[jetMask==True])//1000)*1000
            n_jets = len(jetMask[jetMask==True])
            #print "  >> N jets @ label%d: %d -> %d"%(j, len(jetMask[jetMask==True]), n_jets)
            print "  >> N jets @ label%d: %d"%(d, len(jetMask[jetMask==True]))

            jetIsQuark = jetIsQuark_[jetMask]
            print "  >> %s: %s"%('jetIsQuark', jetIsQuark.shape)
            X_ECAL_stacked = X_ECAL_stacked_[jetMask]
            print "  >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)
            jetEventId = eventId[jetMask]
            jetRunId = runId[jetMask]

            y_jets = da.from_array(\
                    np.full(len(jetEventId), label, dtype=np.float32),\
                    chunks=jetEventId.chunks)

            outPath = '%s/%s_IMG'%(outDir, decay)
            if not os.path.isdir(outPath):
                os.makedirs(outPath)
            file_out_str = "%s/%s_IMG_RH%d_n%d_label%d_list%s_%d.hdf5"%(outPath,decay,int(scale[0]),n_jets,label,list_idx,n)
            print "  >> Writing to:", file_out_str
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

            print "  >> Done.\n"
