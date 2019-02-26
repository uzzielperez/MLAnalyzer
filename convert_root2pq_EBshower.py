import pyarrow.parquet as pq
import pyarrow as pa
import ROOT
import numpy as np
import glob, os

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--infile', required=True, type=str, help='Input root file.')
parser.add_argument('-o', '--outdir', required=True, type=str, help='Output pq file dir.')
parser.add_argument('-d', '--decay', required=True, type=str, help='Decay name.')
parser.add_argument('-n', '--idx', required=True, type=int, help='Input root file index.')
parser.add_argument('-w', '--wgt_file', default='', type=str, help='Weight file.')
args = parser.parse_args()

def get_weight_2d(m0, pt, m0_edges, pt_edges, wgts):
    idx_m0 = np.argmax(m0 <= m0_edges)-1
    idx_pt = np.argmax(pt <= pt_edges)-1
    #print(idx_m0, idx_pt)
    return wgts[idx_m0, idx_pt]

w = np.load(args.wgt_file)
hmvpt, m_edges, pt_edges = w['mvpt'], w['m_edges'], w['pt_edges']

rhTreeStr = args.infile 
rhTree = ROOT.TChain("fevt/RHTree")
rhTree.Add(rhTreeStr)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> Input file:",rhTreeStr
print " >> nEvts:",nEvts

data = [
        pa.array([np.zeros((1,32,32)).tolist()]),
        pa.array([np.zeros((2,32,32)).tolist()]),
        pa.array([124.]),
        pa.array([60.]),
        pa.array([360.]),
        pa.array([85.])
        ]
#keys = ['X', 'm', 'pt']
keys = ['X', 'Xtz', 'm', 'pt', 'iphi', 'ieta']
table = pa.Table.from_arrays(data, keys)
outStr = '%s/%s.parquet.%d'%(args.outdir, args.decay, args.idx) 
#outStr = '%s.parquet.%d'%(args.decay, args.idx) 
writer = pq.ParquetWriter(outStr, table.schema, compression='snappy')

##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
iEvtEnd   = 100
iEvtEnd   = nEvts 
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nAcc = 0
sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    # Initialize event
    rhTree.GetEntry(iEvt)

    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    SC_energyT = rhTree.SC_energyT
    SC_energyZ = rhTree.SC_energyZ
    SC_energy = rhTree.SC_energy
    SC_mass = rhTree.SC_mass
    SC_pT = rhTree.SC_pT
    SC_iphi = rhTree.SC_iphi
    SC_ieta = rhTree.SC_ieta
    #EB_energy = rhTree.EB_energy
    rands = np.random.random((len(SC_mass,)))

    for i in range(len(SC_mass)):

        sc_mass = SC_mass[i]
        sc_pT = SC_pT[i]
        sc_iphi = SC_iphi[i]
        sc_ieta = SC_ieta[i]

        if rands[i] < get_weight_2d(sc_mass, sc_pT, m_edges, pt_edges, hmvpt):
            continue

        sc_energy = np.array(SC_energy[i]).reshape(1,32,32).tolist()
        sc_energyT = np.array(SC_energyT[i]).reshape(1,32,32)
        sc_energyZ = np.array(SC_energyT[i]).reshape(1,32,32)
        sc_energyTZ = np.concatenate((sc_energyT, sc_energyZ), axis=0).tolist()
        #eb_energy = np.array(EB_energy).reshape(1,170,360).tolist()

        pqdata = [
                pa.array([sc_energy])
                ,pa.array([sc_energyTZ])
                ,pa.array([sc_mass])
                ,pa.array([sc_pT])
                ,pa.array([sc_iphi])
                ,pa.array([sc_ieta])
                #,pa.array([eb_energy]) 
                ]

        table = pa.Table.from_arrays(pqdata, keys)
        writer.write_table(table)
        nAcc += 1

writer.close()

sw.Stop()
#print " >> nPhos:",nAcc,"/",2*(iEvtEnd-iEvtStart),"(",100.*nAcc/(2*(iEvtEnd-iEvtStart)),"% )"
print " >> nPhos:",nAcc
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"

pqIn = pq.ParquetFile(outStr)
print(pqIn.metadata)
print(pqIn.schema)
#X = pqIn.read_row_group(0, columns=['X.list.item.list.item.list.item', 'm', 'pt'], nthreads=len(keys)).to_pydict()['X']
#X = pqIn.read(['X.list.item.list.item.list.item']).to_pydict()['X']
#X = np.float32(X)
