import os

cmssw_base = os.getenv("CMSSW_BASE")
cfg= cmssw_base + '/src/MLAnalyzer/RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/FEVTDEBUG/h24gamma_1j_10K_100MeV_FEVTDEBUG_2016_25ns_Moriond17MC_PoissonOOTPU/180109_112606/0000/step_full_1.root'
#inputFiles_='/store/user/mandrews/ML/AODSIM/T7WgStealth_800_200_AODSIM_noPU/180227_224514/0000/step_AODSIM_1.root'
#inputFiles_='/store/user/mandrews/ML/AODSIM/WZToJets_TuneCUETP8M1_13TeV_pythia8_noPU_AODSIM/0/0000/step_full_1.root'
#inputFiles_='/store/user/johnda/AODSIM/Py8PtGun_bb_noPU_AODSIM/180731_210318/0000/step_AODSIM_1.root'
#inputFiles_='/store/user/johnda/AODSIM/Py8PtGun_bb_noPU_AODSIM/180731_210318/0000/step_AODSIM_2.root'

#inputFiles_ ='/store/user/johnda/AODSIM/NonB_Pt_30_70_13TeV_TuneCUETP8M1_noPU_AODSIM/181020_000214/0000/step_AODSIM_9.root'
#inputFiles_ ='/store/user/johnda/AODSIM/QCD_Pt_30_70_13TeV_TuneCUETP8M1_noPU_AODSIM/181019_231226/0000/step_AODSIM_1.root'
#inputFiles_ ='file:/uscms/home/jda102/nobackup/BTaggingML/CMSSW_8_0_30/src/step_AODSIM.root'
#inputFiles_ = 'file:/eos/uscms/store/user/jda102/AODSIM/QCD_Pt_30_70_13TeV_TuneCUETP8M1_noPU_AODSIM/181019_231226/0000/step_AODSIM_83.root'

#inputFiles_='/store/group/lpcml/mandrews/AODSIM/QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU_AODSIM/180809_215549/0000/step_full_1.root'
#inputFiles_='/store/group/lpcml/AODSIM/h24gamma_1j_1M_200MeV_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180907_083212/0000/step_full_1.root'
#inputFiles_='file:../step_AODSIM.root'
#inputFiles_='/store/group/lpcml/mandrews/AODSIM/SM2gammaj_1j_1M_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180906_154222/0000/step_full_98.root'

#inputFiles_='/store/mc/RunIIFall17DRPremix/ADDGravToGG_NegInt-1_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/40000/205BA2DE-1215-E911-A455-001E67DBE79B.root'

#inputFiles_='/store/mc/RunIIFall17DRPremix/ADDGravToGG_NegInt-1_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/40000/105A03BC-DB14-E911-A662-0025905C548A.root'

inputFiles_='/store/mc/RunIIFall17MiniAODv2/ADDGravToGG_NegInt-1_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/58C23153-2A1E-E911-955C-0CC47AFCC372.root'

#inputFiles_='/store/mc/RunIIFall17DRPremix/ADDGravToGG_NegInt-1_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/40000/*.root'

maxEvents_=-1
skipEvents_=0#

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
