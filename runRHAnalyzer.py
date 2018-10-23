import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/FEVTDEBUG/h24gamma_1j_10K_100MeV_FEVTDEBUG_2016_25ns_Moriond17MC_PoissonOOTPU/180109_112606/0000/step_full_1.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/AODSIM/T7WgStealth_800_200_AODSIM_noPU/180227_224514/0000/step_AODSIM_1.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/AODSIM/WZToJets_TuneCUETP8M1_13TeV_pythia8_noPU_AODSIM/0/0000/step_full_1.root'
inputFiles_='file:/eos/uscms/store/group/lpcml/QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU_AODSIM/180809_215549/0000/step_full_1.root'
#inputFiles_='file:../step_full.root'

maxEvents_=-1
skipEvents_=0#

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
