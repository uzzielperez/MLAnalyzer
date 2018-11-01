import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/h24gamma_1j_1M_100MeV_noPU_FEVTDEBUG/180109_233954/0000/step_full_1.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/VBFHiggs0PHToGG_M-125p6_7TeV-JHUGenV4-pythia6-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/04238D05-CAFF-E311-8B63-848F69FD5027.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/GluGluHToGG_M-125_8TeV-pythia6/AODSIM/PU_RD1_START53_V7N-v2/00000/02590948-234A-E411-87AF-7845C4FC3A61.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/DiPhotonBorn_Pt-25To250_8TeV_ext-pythia6/AODSIM/PU_RD1_START53_V7N-v2/10000/0043B5BE-6ED2-E211-91DA-00266CFFBF50.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/DiPhotonBox_Pt-25To250_8TeV_ext-pythia6/AODSIM/PU_RD1_START53_V7N-v2/10000/04037D44-BFCF-E211-BBF3-AC162DABAF78.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV_ext-pythia6/AODSIM/PU_RD1_START53_V7N-v1/20000/000EDB33-4FD2-E211-8E59-0026189438EB.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/DiPhotonJets_M0_8TeV-madgraph/AODSIM/PU_RD1_START53_V7N-v1/10000/00665E6E-12D6-E211-9202-20CF305616E2.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_RD1_START53_V7N-v1/10000/0037DF27-90D9-E211-9A2E-20CF3056170A.root'
inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/AODSIM/PU_RD1_START53_V7N-v3/00000/001C975B-E099-E311-816C-0025905A610A.root'
maxEvents_=50
skipEvents_=0 #9, 20

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
