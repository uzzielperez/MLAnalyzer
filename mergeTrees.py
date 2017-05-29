import os
import sys
import numpy as np
import argparse
import ROOT

## MAIN ##
def main():

    # Keep time
    sw = ROOT.TStopwatch()
    sw.Start()

    # Load input TTrees into TChain
    ggIn = ROOT.TChain("fevt/RHTree")
    ggIn.Add('/eos/cms/store/user/mandrews/ML/IMGs/SinglePositronPt50_FEVTDEBUG_n125k_IMG.root')
    ggIn.Add('/eos/cms/store/user/mandrews/ML/IMGs/SingleElectronPt50_FEVTDEBUG_n125k_IMG.root')
    nEvts = ggIn.GetEntries()
    print " >> nEvts:",nEvts

    # Initialize output file as empty clone
    outFileStr = "/eos/cms/store/user/mandrews/ML/IMGs/SingleElePosPt50_FEVTDEBUG_n250k_IMG.root"
    outFile = ROOT.TFile(outFileStr, "RECREATE")
    ggOut = ggIn.CloneTree(-1) 
    outFile.Write()
    print " >> Output file:",outFileStr

    sw.Stop()
    print " >> Real time:",sw.RealTime()/60.,"minutes"
    print " >> CPU time: ",sw.CpuTime() /60.,"minutes"


#_____ Call main() ______#
if __name__ == '__main__':
    main()
