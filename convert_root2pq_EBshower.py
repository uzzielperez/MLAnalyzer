import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
import glob, os

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('-i', '--infile', default='output.root', type=str, help='Input root file.')
parser.add_argument('-i', '--infile', default=['output.root'], nargs='+', type=str, help='Input root file.')
#parser.add_argument('-i', '--infile', default=['output.root'], type=list, help='Input root file.')
parser.add_argument('-o', '--outdir', default='.', type=str, help='Output pq file dir.')
parser.add_argument('-d', '--decay', default='test', type=str, help='Decay name.')
parser.add_argument('-n', '--idx', default=0, type=int, help='Input root file index.')
parser.add_argument('-w', '--wgt_file', default=None, type=str, help='Weight file.')
args = parser.parse_args()

def crop_EBshower(imgEB, iphi, ieta, window=32):

    # NOTE: image window here should correspond to the one used in RHAnalyzer
    off = window//2
    iphi = int(iphi)+1 # seed positioned at [15,15]
    ieta = int(ieta)+1 # seed positioned at [15,15]

    # Wrap-around on left side
    if iphi < off:
        diff = off-iphi
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,-diff:],
                                   imgEB[:,ieta-off:ieta+off,:iphi+off]), axis=-1)
    # Wrap-around on right side
    elif 360-iphi < off:
        diff = off - (360-iphi)
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,iphi-off:],
                                   imgEB[:,ieta-off:ieta+off,:diff]), axis=-1)
    # Nominal case
    else:
        img_crop = imgEB[:,ieta-off:ieta+off,iphi-off:iphi+off]

    return img_crop

def get_weight_2d(m0, pt, m0_edges, pt_edges, wgts):
    idx_m0 = np.argmax(m0 <= m0_edges)-1
    idx_pt = np.argmax(pt <= pt_edges)-1
    #print(idx_m0, idx_pt)
    return wgts[idx_m0, idx_pt]

if args.wgt_file is not None:
    w = np.load(args.wgt_file)
    hmvpt, m_edges, pt_edges = w['mvpt'], w['m_edges'], w['pt_edges']

rhTreeStr = args.infile 
print " >> Input file:",rhTreeStr
rhTree = ROOT.TChain("fevt/RHTree")
for f in rhTreeStr:
  rhTree.Add(f)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> nEvts:",nEvts
#outStr = '%s/%s.parquet.%d'%(args.outdir, args.decay, args.idx) 
outStr = '%s/%s.reg_2reco.parquet.%d'%(args.outdir, args.decay, args.idx) 
print " >> Output file:",outStr

##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
#iEvtEnd   = 10
iEvtEnd   = nEvts 
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nPhos = 0
data = {} # Arrays to be written to parquet should be saved to data dict
sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    # Initialize event
    rhTree.GetEntry(iEvt)

    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    if rhTree.m0 < 100. or rhTree.m0 > 110.:
      continue

    idx = [rhTree.runId, rhTree.lumiId, rhTree.eventId]

    SC_energyT = rhTree.SC_energyT
    SC_energyZ = rhTree.SC_energyZ
    SC_energy = rhTree.SC_energy
    #SC_mass = rhTree.SC_mass
    #SC_pT = rhTree.SC_pT
    SC_iphi = rhTree.SC_iphi
    SC_ieta = rhTree.SC_ieta

    pho_pTs =             rhTree.pho_pT
    pho_r9s =             rhTree.pho_r9
    pho_sieies =          rhTree.pho_sieie
    pho_phoIsos =         rhTree.pho_phoIso
    pho_chgIsos =         rhTree.pho_chgIso
    pho_chgIsoWrongVtxs = rhTree.pho_chgIsoWrongVtx
    pho_Eraws =           rhTree.pho_Eraw
    pho_phiWidths =       rhTree.pho_phiWidth
    pho_etaWidths =       rhTree.pho_etaWidth
    pho_scEtas =          rhTree.pho_scEta
    pho_sieips =          rhTree.pho_sieip
    pho_s4s =             rhTree.pho_s4

    #EB_time = np.array(rhTree.EB_time).reshape(170,360)
    #TracksAtEB_pt = np.array(rhTree.TracksPt_EB).reshape(170,360)
    #X_EB = np.stack([TracksAtEB_pt, EB_time], axis=0) 
    #X_EB = np.array(rhTree.EB_energy).reshape(1,170,360)

    nphos = len(SC_iphi)
    rands = np.random.random(nphos)
    for i in range(nphos):

        data['idx'] = idx + [i]
        data['m0'] = rhTree.m0

        #data['m'] = SC_mass[i]
        #data['pt'] = SC_pT[i]
        data['iphi'] = SC_iphi[i]
        data['ieta'] = SC_ieta[i]

        #if data['ieta'] >= 170-16:
        #    continue

        if args.wgt_file is not None:
            if rands[i] < get_weight_2d(data['m'], data['pt'], m_edges, pt_edges, hmvpt):
                continue

        data['pt_reco'] = pho_pTs[i]
        data['Xvars'] = [
            pho_r9s[i]
            ,pho_sieies[i]
            ,pho_phoIsos[i]
            ,pho_chgIsos[i]
            ,pho_chgIsoWrongVtxs[i]
            ,pho_Eraws[i]
            ,pho_phiWidths[i]
            ,pho_etaWidths[i]
            ,pho_scEtas[i]
            ,pho_sieips[i]
            ,pho_s4s[i]
        ]

        data['X'] = np.array(SC_energy[i]).reshape(1,32,32)
        sc_energyT = np.array(SC_energyT[i]).reshape(1,32,32)
        sc_energyZ = np.array(SC_energyT[i]).reshape(1,32,32)
        data['Xtz'] = np.concatenate((sc_energyT, sc_energyZ), axis=0)

        #sc_cms = crop_EBshower(X_EB, data['iphi'], data['ieta']) 
        #if sc_cms.shape != data['Xtz'].shape:
        #    print(sc_cms.shape)
        #    print(data['Xtz'].shape)
        #data['Xcms'] = np.concatenate([sc_cms, data['Xtz']], axis=0)

        pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in data.values()]
        table = pa.Table.from_arrays(pqdata, data.keys())

        if nPhos == 0:
            writer = pq.ParquetWriter(outStr, table.schema, compression='snappy')

        writer.write_table(table)

        nPhos += 1

writer.close()

sw.Stop()
print " >> nPhos:",nPhos
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
print " >> ======================================"

pqIn = pq.ParquetFile(outStr)
print(pqIn.metadata)
print(pqIn.schema)
#X = pqIn.read_row_group(0, columns=['m','pt','iphi','ieta','pt_reco']).to_pydict()
X = pqIn.read_row_group(0, columns=['idx.list.item','iphi','ieta','pt_reco']).to_pydict()
X = pqIn.read_row_group(1, columns=['idx.list.item','iphi','ieta','pt_reco']).to_pydict()
print(X)
#X = pqIn.read_row_group(0, columns=['X.list.item.list.item.list.item']).to_pydict()['X']
#X = pqIn.read(['X.list.item.list.item.list.item']).to_pydict()['X']
#X = np.float32(X)
