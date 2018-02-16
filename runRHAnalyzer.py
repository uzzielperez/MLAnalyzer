import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:/eos/cms/store/user/mandrews/ML/FEVTDEBUG/T7WgStealth_800_200_FEVTDEBUG_noPU/180214_210138/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/h24gamma_1j_1M_100MeV_noPU_FEVTDEBUG/180109_233954/0000/step_full_1.root'
#inputFiles_='root://cmseos.fnal.gov:1094//store/user/mba2012/FEVTDEBUG/step_full_17.root'
maxEvents_=-1
skipEvents_=0
#outputFile_='output.root'
inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
inputTag='TEST'

if not os.path.isdir(inputTag):
    os.system('mkdir %s'%(inputTag))
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
