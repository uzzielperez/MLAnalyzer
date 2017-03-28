import numpy as np
#import ROOT
#from scipy.sparse import csr_matrix
from root_numpy import root2array, tree2array
from root_pandas import read_root
#import tables

eosDir='/eos/cms/store/user/mandrews/ML'

pid = "Electron"
#pid = "Photon"

#TODO: implement as class with class variable iterated at every call to function
# in order to go forward through chunks
def x_arr(start_,stop_):
  arr = root2array('%s/Single%s_n10000_IMG.root'%(eosDir,pid), 'fevt/RHTree', start=start_, stop=stop_) 
  print arr.dtype.names
  return np.array([np.concatenate(x).reshape(-1,61200) for x in arr])

chunk_size=1000
for i in range(10):
  print "chunk",i
  x = x_arr(i*chunk_size,(i+1)*1000)
  np.savez_compressed('IMGs/Single%s_n1k_IMG_%d'%(pid,i),x=x)
