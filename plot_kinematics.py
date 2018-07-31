import numpy as np
import ROOT

histname = "h_pT"
#histname = "h_m0"
#histname = "h_E"
#histname = "h_eta"
hist = "fevt/%s"%histname

t_Hgg = ROOT.TFile('~/eos/OPENDATA/IMGs/MGG90_Eta14v2/GluGluHToGG_MGG90_Eta14v2_IMG.root', 'READ')
#t_Hgg = ROOT.TFile('~/eos/OPENDATA/IMGs/MGG90_Eta14v2/DiPhotonBorn_MGG90_Eta14v2_IMG.root', 'READ')
#t_Hgg = ROOT.TFile('h_EB_energyom_GluGluHToGG.root', 'READ')
#t_Hgg = ROOT.TFile('h_EB_energyom_DiPhotonBorn.root', 'READ')
#t_Hgg = ROOT.TFile('h_EB_energyom_h22gammaSM.root', 'READ')
#t_Hgg = ROOT.TFile('~/eos/ML/IMGs/h22gammaSM_1j_1M_noPU_FEVTDEBUG_IMG.root', 'READ')
#tree_Hgg = t_Hgg.Get('fevt/RHTree')
#h_Hgg = ROOT.TH1F("h8TeV", "h8TeV", 60, 90., 210.)
#h_Hgg = ROOT.TH1F("h8TeV", "h8TeV", 60, 0., 6.)
#tree_Hgg.Draw("m0 >> h8TeV")
#tree_Hgg.Draw("diPhoPt/m0 >> h8TeV")
#histname = "h_phoEta"
histname = "h_phoPt"
#histname = "h_EB_energyom"
h_Hgg = ROOT.gDirectory.Get("fevt/%s"%histname)
#h_Hgg = ROOT.gDirectory.Get(histname)
print " >> Hgg entries:",h_Hgg.GetEntries()

#t_Ggg = ROOT.TFile('~/eos/ML/IMGs/h22gammaSM_1j_1M_noPU_FEVTDEBUG_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/ML/IMGs/SM2gamma_1j_1M_noPU_FEVTDEBUG_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('h_EB_energyom_h22gammaSM.root', 'READ')
#t_Ggg = ROOT.TFile('h_EB_energyom_SM2gamma.root', 'READ')
#t_Ggg = ROOT.TFile('h_EB_energyom_DiPhotonBorn.root', 'READ')
t_Ggg = ROOT.TFile('~/eos/OPENDATA/IMGs/MGG90_Eta14v2/DiPhotonBorn_MGG90_Eta14v2_IMG.root', 'READ')
#tree_Ggg = t_Ggg.Get('fevt/RHTree')
#h_Ggg = ROOT.TH1F("h13TeV", "h13TeV", 60, 90., 210.)
#h_Ggg = ROOT.TH1F("h13TeV", "h13TeV", 60, 0., 6.)
#tree_Ggg.Draw("diPhoPt/m0 >> h13TeV")
#histname = "h_eta"
#histname = "h_pT"
h_Ggg = ROOT.gDirectory.Get("fevt/%s"%histname)
#h_Ggg = ROOT.gDirectory.Get(histname)
print " >> gg entries:",h_Ggg.GetEntries()

c = ROOT.TCanvas("c","c",600,600)
c.SetBorderSize(0);
c.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetLogy()
c.cd()

ROOT.gPad.SetLeftMargin(0.15)

ymax = 0.1
ymax = 1.
h_Hgg.SetLineColor(2)
h_Hgg.GetYaxis().SetTitleOffset(2.)
h_Hgg.GetYaxis().SetTitle("N_{#gamma} / N_{#gamma,tot}")
h_Hgg.GetXaxis().SetTitleOffset(1.2)
h_Hgg.GetXaxis().SetTitle("p_{T} [GeV]")
h_Hgg.Scale(1./h_Hgg.Integral())
h_Hgg.Draw("")
#h_Hgg.GetXaxis().SetRangeUser(90.,350.)
#h_Hgg.GetXaxis().SetRangeUser(20.,350.)
ROOT.gPad.Update()
h_Ggg.Scale(1./h_Ggg.Integral())

if histname == 'h_eta' or histname == 'h_phoEta':
    ymax = 0.05
    h_Hgg.GetXaxis().SetTitle("#eta")
    h_Hgg.GetXaxis().SetRangeUser(-2.,2.)
    h_Hgg.GetYaxis().SetRangeUser(0.,ymax)
#h_Hgg.GetYaxis().SetRangeUser(0.,ymax)
h_Hgg.GetYaxis().SetRangeUser(5.e-6,3.e-1)
if histname == 'h_EB_energyom':
    ymax = 0.06
    h_Hgg.GetYaxis().SetTitle("N_{hits}")
    h_Hgg.GetXaxis().SetTitle("E/m_{#gamma#gamma}")
    #h_Hgg.GetXaxis().SetRangeUser(-1.5,1.5)
h_Ggg.Draw("SAME")

leg = ROOT.TLegend(0.65,0.7,0.85,0.85)
#leg.AddEntry(h_Hgg,"8 TeV","LP")
#leg.AddEntry(h_Ggg,"13 TeV","LP")
leg.AddEntry(h_Hgg,"H","LP")
leg.AddEntry(h_Ggg,"Prompt","LP")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.Draw()

c.Print('%s_8TeV.png'%histname)
#c.Print('%s_13TeV.png'%histname)
