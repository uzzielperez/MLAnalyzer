import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/FEVTDEBUG/h24gamma_1j_10K_100MeV_FEVTDEBUG_2016_25ns_Moriond17MC_PoissonOOTPU/180109_112606/0000/step_full_1.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/AODSIM/T7WgStealth_800_200_AODSIM_noPU/180227_224514/0000/step_AODSIM_1.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/MINIAODSIM/h2gg.root'
inputFiles_='/store/mc/RunIISummer17MiniAOD/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/30000/FEFDBA0A-74C2-E711-8EB7-00259029E81A.root'

maxEvents_=-1
skipEvents_=0

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
