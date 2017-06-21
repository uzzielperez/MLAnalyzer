import numpy as np
import ROOT

histname = "h_pT"
#histname = "h_E"
#histname = "h_eta"
hist = "fevt/%s"%histname

#t_Hgg = ROOT.TFile('~/eos/IMGs/TEST_eta14/SingleHiggsPt10to90_eta14_pythia8_FEVTDEBUG_TEST_IMG.root', 'READ')
t_Hgg = ROOT.TFile('~/eos/IMGs/TEST_eta14/SingleHiggsPt10to80_eta14_pythia8_FEVTDEBUG_TEST_IMG.root', 'READ')
#t_Hgg = ROOT.TFile('~/eos/IMGs/TEST_eta14/SingleHiggsPt10to75_eta14_pythia8_FEVTDEBUG_TEST_IMG.root', 'READ')
#t_Hgg = ROOT.TFile('~/eos/IMGs/H125GGgluonfusion_13TeV_TuneCUETP8M1_HighLumiPileUp_FEVTDEBUG_n350k_IMG.root', 'READ')
#t_Hgg = ROOT.TFile('~/eos/IMGs/TEST_eta23/SingleHiggsPt10to100_pythia8_FEVTDEBUG_TEST_IMG.root', 'READ')
#t_Hgg = ROOT.TFile('~/eos/IMGs/H125GGgluonfusion_13TeV_TuneCUETP8M1_FEVTDEBUG_n175k_IMG.root', 'READ')
h_Hgg = ROOT.gDirectory.Get(hist)
print " >> Hgg entries:",h_Hgg.GetEntries()
t_Ggg = ROOT.TFile('~/eos/IMGs/TEST_eta14/DoublePhotonGaussPt55_StdDev20_FEVTDEBUG_TEST_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/IMGs/TEST_eta23/PromptDiPhoton_PtHat10_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_FEVTDEBUG_TEST_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/IMGs/DoublePhotonGaussPt55_StdDev20_FEVTDEBUG_HighLumiPileUp_n250k_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/IMGs/DoublePhotonGaussPt55_StdDev20_FEVTDEBUG_HighLumiPileUp_n175k_IMG_numEvent9216.root', 'READ')
h_Ggg = ROOT.gDirectory.Get(hist)
print " >> gg entries:",h_Ggg.GetEntries()

c = ROOT.TCanvas("c","c",600,600)
c.SetBorderSize(0);
c.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
#ROOT.gPad.SetLogy()
c.cd()

ROOT.gPad.SetLeftMargin(0.15)

h_Ggg.SetLineColor(2)
h_Ggg.Scale(1./h_Ggg.Integral())
h_Ggg.Draw("")
h_Hgg.GetYaxis().SetTitleOffset(2.)
h_Hgg.Scale(1./h_Hgg.Integral())
h_Hgg.Draw("SAME")

leg = ROOT.TLegend(0.6,0.65,0.85,0.8)
leg.AddEntry(h_Hgg,"H#rightarrow#gamma#gamma","LP")
leg.AddEntry(h_Ggg,"#gamma#gamma, GaussPt55","LP")
leg.SetBorderSize(0)
leg.Draw()

c.Print('%s.png'%histname)
