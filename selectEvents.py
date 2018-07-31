import os
import sys
import numpy as np
import ROOT
#import argparse
import sys

# Register command line options
#parser = argparse.ArgumentParser(description='Run STEALTH selection.')
#parser.add_argument('-s','--inFile', required=True, type=str)
#args = parser.parse_args()

## MAIN ##
def main():

    # Load input TTrees into TChain
    eosDir = "/eos/cms/store/user/mandrews/OPENDATA/IMGs/MGG90_Eta14v2"
    #ggInStr = "%s/DiPhotonBorn_10000_IMG.root"%(eosDir)
    ggInStr = "%s/GluGluHToGG_MGG90_Eta14v2_IMG.root"%(eosDir)
    #ggInStr = args.inFile
    ggInStr = sys.argv[1] 
    ggIn = ROOT.TChain("fevt/RHTree")
    ggIn.Add(ggInStr)
    nEvts = ggIn.GetEntries()
    print " >> nEvts:",nEvts
    nEvts = 57425 
    print " >> Input file(s):",ggInStr
    print " >> nEvts:",nEvts

    runCount = {}
    runCount['194533'] = 0
    runCount['200519'] = 0
    runCount['206859'] = 0

    h_EB_energyom = ROOT.TH1F('h_EB_energyom', 'h_EB_energyom;energy/m_{#gamma#gamma};N_{hits}', 50, 0., 5.)

    sw = ROOT.TStopwatch()
    sw.Start()
    nTotal = 0
    for jEvt in range(nEvts):

        # Initialize event
        if jEvt > nEvts:
            break
        treeStatus = ggIn.LoadTree(jEvt)
        if treeStatus < 0:
            break
        evtStatus = ggIn.GetEntry(jEvt)
        if evtStatus <= 0:
            continue
        if jEvt % 1000 == 0:
            print " .. Processing entry",jEvt
            print ">> 0:", runCount['194533']
            print ">> 1:", runCount['200519']
            print ">> 2:", runCount['206859']

        #run = int(ggIn.runId)
        #runCount['%s'%run] += 1
        nTotal += 1
        #if run == 194533:
        #  print(">>",0)
        #if run == 200519:
        #  print(">>",1)
        #if run == 206859:
        #  print(">>",2)
        m = ggIn.m0
        #for en in ggIn.EB_energy[ggIn.EB_energy > 0.]:
        #continue
        for en in ggIn.EB_energy:
          if en < 0.001:
            continue
          h_EB_energyom.Fill(en/m)


    #for k, val in runCount.iteritems():
    #  print k, val

    print nTotal

    sw.Stop()
    print " >> Real time:",sw.RealTime()/60.,"minutes"
    print " >> CPU time: ",sw.CpuTime() /60.,"minutes"

    #hFile = ROOT.TFile("h_EB_energyom_GluGluH2GG.root","RECREATE")
    hFile = ROOT.TFile("h_EB_energyom_%s.root"%(ggInStr.split('/')[-1].split('_')[0]), "RECREATE")
    h_EB_energyom.Write()
    hFile.Close()


#_____ Call main() ______#
if __name__ == '__main__':
    main()
