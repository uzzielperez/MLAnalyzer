import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/h24gamma_1j_1M_100MeV_noPU_FEVTDEBUG/180109_233954/0000/step_full_1.root'
inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/VBFHiggs0PHToGG_M-125p6_7TeV-JHUGenV4-pythia6-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/04238D05-CAFF-E311-8B63-848F69FD5027.root'
maxEvents_=-1
skipEvents_=0

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
