import ROOT
import numpy as np
np.random.seed(0)
import glob, os

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-d', '--genDR', default=10, type=int, help='gen-level dR.')
parser.add_argument('-p', '--p_drop', default=0.95, type=float, help='p(drop) scale.')
args = parser.parse_args()

genDR = args.genDR
p_drop_scale = args.p_drop

xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'

#pu = 'noPU'
pu = '2016_25ns_Moriond17MC_PoissonOOTPU'
decay = 'DoublePi0Pt15To100_m0To1600_pythia8_%s_mlog_ptexp'%pu
decay = '%s_genDR%d_recoDR16_seedPos_phoVars_IMG'%(decay, genDR)
date_str = '190310_173649' # DR10

# Paths to input files
rhFileList = '%s/%s/%s/*/output_*.root'%(eosDir, decay, date_str)
print(" >> Input file list: %s"%rhFileList)
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print(" >> %d files found"%len(rhFileList))
rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
print(' >> Input File[0]: %s'%rhFileList[0])

rhTree = ROOT.TChain("fevt/RHTree")
for f in rhFileList:
    rhTree.Add(f)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> nEvts:",nEvts

##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
#iEvtEnd   = 10000
iEvtEnd   = nEvts 
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

eb_xtal = 0.0174

nAcc = 0
sc_mass_, sc_pT_, sc_dR_ = [], [], []
hnPho = ROOT.TH2F("hnPho", "hnPho;m_{#pi^{0}};p_{T,#pi^{0}}", 16, 0, 1.6 , 20, 20., 100.)

sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    rhTree.GetEntry(iEvt)

    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    SC_mass = rhTree.SC_mass
    SC_pT = rhTree.SC_pT
    SC_dR = rhTree.SC_DR

    for i in range(len(SC_mass)):

        sc_mass = SC_mass[i] 
        sc_pT = SC_pT[i] 
        sc_dR = SC_dR[i] 

        hnPho.Fill(sc_mass, sc_pT) 
        nAcc += 1

        sc_mass_.append(sc_mass)
        sc_pT_.append(sc_pT)
        sc_dR_.append(sc_dR)

sw.Stop()
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
print " >> nPi0s: ", nAcc

sc_mass_ = np.array(sc_mass_)
sc_pT_ = np.array(sc_pT_)
sc_dR_ = np.array(sc_dR_)

hmvpt, m_edges, pt_edges = np.histogram2d(sc_mass_, sc_pT_, range=((0., 1.6), (20.,100.)), bins=(16, 20))
hmvpt = hmvpt/hmvpt.max()
floor = np.mean(hmvpt.flatten()) - 2.*np.std(hmvpt.flatten())
if floor < 0.:
    print " >> Forcing floor to %f -> 0."%floor
    floor = 0.
print " >> w(m_v_pt) floor:",floor
hmvpt[hmvpt < floor] = floor
hmvpt = hmvpt/hmvpt.max()
hmvpt = hmvpt-hmvpt.min()
print " >> p(drop) scale:",p_drop_scale
hmvpt = p_drop_scale*hmvpt
print " >> w(m_v_pt)*p(drop) min, max:",hmvpt.min(), hmvpt.max()
np.savez('%s_mvpt_weights_pdrop%.2f.npz'%(decay, p_drop_scale),\
        mvpt = hmvpt,
        m_edges = m_edges,
        pt_edges = pt_edges,
        sc_mass = sc_mass_,
        sc_pT = sc_pT_,
        sc_dR = sc_dR_
        )

#hnPho.Draw("COL Z")

def get_weight_2d(m0, pt, m0_edges, pt_edges, wgts):
    idx_m0 = np.argmax(m0 <= m0_edges)-1
    idx_pt = np.argmax(pt <= pt_edges)-1
    #print(idx_m0, idx_pt)
    return wgts[idx_m0, idx_pt]

hnPhow = ROOT.TH2F("hnPhow", "hnPhow;m_{#pi^{0}};p_{T,#pi^{0}}", 16, 0, 1.6 , 20, 20., 100.)
nAcc = 0
for iEvt in range(iEvtStart,iEvtEnd):

    rhTree.GetEntry(iEvt)

    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    SC_mass = rhTree.SC_mass
    SC_pT = rhTree.SC_pT
    #SC_dR = rhTree.SC_DR
    rands = np.random.random((len(SC_mass,)))

    for i in range(len(SC_mass)):

        sc_mass = SC_mass[i] 
        sc_pT = SC_pT[i] 
        #sc_dR = SC_dR[i] 

        if rands[i] < get_weight_2d(sc_mass, sc_pT, m_edges, pt_edges, hmvpt):
            continue

        hnPhow.Fill(sc_mass, sc_pT) 
        nAcc += 1

print " >> nPi0s: ", nAcc
#hnPhow.Draw("COL Z")

hFile = ROOT.TFile("%s_hnPho_pdrop%.2f.root"%(decay, p_drop_scale),"RECREATE")

hnPho.SetMinimum(0.)
hnPho.Write()
hnPhoX = hnPho.ProjectionX()
hnPhoX.GetYaxis().SetRangeUser(0., 1.2*hnPhoX.GetMaximum())
hnPhoX.Write()
hnPhoY = hnPho.ProjectionY()
hnPhoY.GetYaxis().SetRangeUser(0., 1.2*hnPhoY.GetMaximum())
hnPhoY.Write()

hnPhow.SetMinimum(0.)
hnPhow.Write()
hnPhowX = hnPhow.ProjectionX()
hnPhowX.GetYaxis().SetRangeUser(0., 1.2*hnPhowX.GetMaximum())
hnPhowX.Write()
hnPhowY = hnPhow.ProjectionY()
hnPhowY.GetYaxis().SetRangeUser(0., 1.2*hnPhowY.GetMaximum())
hnPhowY.Write()

hFile.Close()
