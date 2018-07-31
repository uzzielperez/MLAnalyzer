import os

cfg='RecHitAnalyzer/python/SCRegressor_cfg.py'
#inputFiles_='file:../step_full.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/h24gamma_1j_1M_100MeV_noPU_FEVTDEBUG/180109_233954/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m000/180413_215734/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m100/180413_215901/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m000/180416_105657/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/DoublePi0Pt30To90_pythia8_m000_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180604_160025/0004/step_full_4917.root'
inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/DoublePi0Pt50To60_m0To1600_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180620_123051/0000/step_full_filtered_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/DoublePi0Pt50To60_m000_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180622_100748/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m600/180611_092217/0000/step_full_1.root'
#inputFiles_='file:step_full_filtered.root'
maxEvents_=-1
skipEvents_=0
#outputFile_='output.root'
#inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
#inputTag='TEST'

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
#for ievt in range(1):
#if not os.path.isdir(inputTag):
#    os.system('mkdir %s'%(inputTag))
#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,ievt)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
#os.system('mv c*.eps %s/'%(inputTag))
#    os.system('mv cEB*_%d.eps %s/'%(ievt+1,inputTag))

#os.system('scram b -j8')
