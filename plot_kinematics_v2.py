import numpy as np
import ROOT

histname = "h_pT"
#histname = "h_m0"
#histname = "h_E"
#histname = "h_eta"
hist = "fevt/%s"%histname

nProc = 2.e5
t_Hgg = ROOT.TFile('~/eos/OPENDATA/IMGs/QCD_00000_IMG_numEvent%d_qq_gg.root'%(nProc), 'READ')
#t_Hgg = ROOT.TFile('~/eos/OPENDATA/IMGs/QCD_00000_IMG_numEvent%d_qg.root'%(nProc), 'READ')
tree_Hgg = t_Hgg.Get('fevt/RHTree')
h_Q = ROOT.TH1F("h_Q", "h_Q", 40, 50., 250.)
h_G = ROOT.TH1F("h_G", "h_G", 40, 50., 250.)
tree_Hgg.Draw("jetPt0 >> h_Q", "jetIsQuark0 == 1. && jetPt0 > 70.")
tree_Hgg.Draw("jetPt0 >> h_G", "jetIsQuark0 == 0. && jetPt0 > 70.")
print " >> Q entries:",h_Q.GetEntries(),float(h_Q.GetEntries())/nProc,29.5e6*float(h_Q.GetEntries())/nProc
print " >> G entries:",h_G.GetEntries(),float(h_G.GetEntries())/nProc,29.5e6*float(h_G.GetEntries())/nProc

c = ROOT.TCanvas("c","c",700,600)
c.SetBorderSize(0);
c.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
#ROOT.gPad.SetLogy()
c.cd()

ROOT.gPad.SetLeftMargin(0.15)

ymax = 0.1
ymax = 1.
#h_G.SetLineColor(2)
h_G.GetYaxis().SetTitleOffset(2.)
h_G.GetYaxis().SetTitle("f_{evts} / 5 GeV")
h_G.GetXaxis().SetTitleOffset(1.2)
h_G.GetXaxis().SetTitle("p_{T} [GeV]")
h_G.Scale(1./h_G.Integral())
h_G.Draw("")
ROOT.gPad.Update()
h_Q.SetLineColor(2)
h_Q.Scale(1./h_Q.Integral())
h_Q.Draw("SAME")

leg = ROOT.TLegend(0.65,0.7,0.85,0.85)
#leg.AddEntry(h_Q,"8 TeV","LP")
#leg.AddEntry(h_G,"13 TeV","LP")
leg.AddEntry(h_G,"G","LP")
leg.AddEntry(h_Q,"Q","LP")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw()

#c.Print('%s_8TeV.png'%histname)
##c.Print('%s_13TeV.png'%histname)
