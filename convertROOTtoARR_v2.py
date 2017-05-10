import numpy as np
#import ROOT
#from scipy.sparse import csr_matrix
from root_numpy import root2array, tree2array
from root_pandas import read_root
#import tables

pid = "ele"
#pid = "pho"


#df = read_root('histo_%s.root'%(pid), 'fevt/RHTree')
#df.to_sparse(fill_value=0).to_hdf('test.h5', pid, mode='w', complevel=9, complib='blosc')

#for df in read_root('histo_%s.root'%(pid), 'fevt/RHTree', chunksize=100):
#  df.to_sparse(fill_value=0).to_hdf('test.h5', pid, mode='a', complevel=9, complib='blosc')

#TODO: implement as class with class variable iterated at every call to function
# in order to go forward through chunks
#col_names = []
def x_arr(start_,stop_):
  #arr = root2array('histo_%s_.root'%(pid), 'fevt/RHTree', start=start_, stop=stop_) 
  arr = root2array('histo_%s_.root'%(pid), 'fevt/RHTree') 
  print arr.dtype.names
  #np.savez('fulltree_ele',x=arr)
  return np.array([np.concatenate(x).reshape(-1,61200) for x in arr])
  #return arr

x = x_arr(0,10)
#x.dtype.names = col_names
#print col_names
#print x.dtype.names
np.savez_compressed('fulltree_reshape_%s'%pid,x=x)

#np.savez('fulltree_ele',rhtree=arr)
'''
chain_in = ROOT.TChain("fevt/RHTree")
chain_in.Add("histo_%s.root"%(pid))

nEvts = chain_in.GetEntries()
print " >> nEvts:",nEvts


#sw = ROOT.TStopwatch()
#sw.Start()

branch_ = tree2array(chain_in, branches=['EBenergy'])

#f = h5py.File('test.hdf5','w')
#f = f.create_dataset("dset",(nEvts,61200), chunks=True)
#hdf_out = tables.open_file('test.h5',mode='w', title="RHfile")
#hdf_out.create_array(hdf_out.root, "EBenergy", branch_, "title")
#hdf_out.close()
'''

#sw.Stop()
#print " >> Real time:",sw.RealTime()/60.
