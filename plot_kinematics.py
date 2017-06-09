import numpy as np
import ROOT

#histname = "H_pT"
histname = "H_E"
#histname = "H_eta"
hist = "fevt/%s"%histname

t_Hgg = ROOT.TFile('~/eos/IMGs/H125GGgluonfusion_13TeV_TuneCUETP8M1_FEVTDEBUG_n175k_IMG.root', 'READ')
h_Hgg = ROOT.gDirectory.Get(hist)
print " >> Hgg entries:",h_Hgg.GetEntries()
t_Ggg = ROOT.TFile('~/eos/IMGs/DoublePhotonGaussPt55_StdDev20_FEVTDEBUG_HighLumiPileUp_n175k_IMG_numEvent9216.root', 'READ')
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

h_Hgg.Draw("")
h_Hgg.GetYaxis().SetTitleOffset(2.)
h_Ggg.SetLineColor(2)
h_Ggg.Draw("SAME")

leg = ROOT.TLegend(0.6,0.65,0.85,0.8)
leg.AddEntry(h_Hgg,"H#rightarrow#gamma#gamma","LP")
leg.AddEntry(h_Ggg,"#gamma#gamma, GaussPt55","LP")
leg.SetBorderSize(0)
leg.Draw()

c.Print('%s.png'%histname)
