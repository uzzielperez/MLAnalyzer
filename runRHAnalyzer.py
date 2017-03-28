import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:SinglePhotonPt50_FEVTDEBUG.root'
#inputFiles_='file:SingleElectronPt50_FEVTDEBUG.root'
maxEvents_=1
skipEvents_=3
#outputFile_='output.root'
inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
for ievt in range(1):
  if not os.path.isdir(inputTag):
    os.system('mkdir %s'%(inputTag))
    cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,ievt)
  print '%s'%cmd
    os.system(cmd)
    os.system('mv cEB*_%d.eps %s/'%(ievt+1,inputTag))

#os.system('scram b -j8')
