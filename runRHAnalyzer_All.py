import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:SinglePhotonPt50_FEVTDEBUG.root'
#inputFiles_='file:SingleElectronPt50_FEVTDEBUG.root'
maxEvents_=-1
skipEvents_=0

#inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
inputTag='TEST_AllEvts'

if not os.path.isdir(inputTag):
  os.system('mkdir %s'%(inputTag))
  cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
  print '%s'%cmd
  os.system(cmd)
  os.system('mv cEB*.eps %s/'%(inputTag))

#os.system('scram b -j8')
