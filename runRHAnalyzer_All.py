import os
from glob import glob
import re

eosDir='/eos/cms/store/user/mandrews/ML'
decay='SinglePhotonPt50_FEVTSIM'

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:SinglePhotonPt50_FEVTDEBUG.root'
#inputFiles_='file:SingleElectronPt50_FEVTDEBUG.root'
#inputFiles_='file:../1488894123_SinglePhotonPt50/SinglePhotonPt50_FEVTDEBUG_n5000.root'
#inputFiles_='file:../1488894128_SingleElectronPt50/SingleElectronPt50_FEVTDEBUG_n5000.root'
#inputFiles_='file:~/eos/cms/store/user/mandrews/ML/1489488093_SingleElectronPt50/SingleElectronPt50_FEVTDEBUG_n10.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/1489488767_SinglePhotonPt50/SinglePhotonPt50_FEVTDEBUG_n10000.root'
#inputFiles_='file:/eos/cms/store/user/mandrews/ML/1489488713_SingleElectronPt50/SingleElectronPt50_FEVTDEBUG_n10000.root'
#inputFiles_='file:../1489423343_SinglePhotonPt50/SinglePhotonPt50_FEVTDEBUG_n5000.root'
#inputFiles_ = ['file:%s'%path for path in glob('%s/FEVTSIM/%s/170325_195239/0000/*root'%(eosDir,decay))]
inputFiles_ = ['file:%s'%path for path in glob('%s/FEVTSIM/%s/170325_195323/0000/*root'%(eosDir,decay))]
inputFiles_ = re.sub('[\[ \]]','',str(inputFiles_))

maxEvents_=-1
skipEvents_=0

#inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
#inputTag='TEST_AllEvts_Ele'
#inputTag='TEST_AllEvts_Pho'
inputTag='TEST'

if not os.path.isdir(inputTag):
  os.system('mkdir %s'%(inputTag))
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
#print '%s'%cmd
os.system(cmd)

#os.system('mv cEB*.eps %s/'%(inputTag))
#os.system('mv histo.root /eos/cms/store/user/mandrews/ML/SinglePhoton_n10000_IMG.root')
cmdStr = 'mv histo.root %s/%s_n250k_IMG.root'%(eosDir,decay)
#os.system(cmdStr)
#os.system('mv histo.root /eos/cms/store/user/mandrews/ML/SingleElectron_n10000_IMG.root')

#os.system('scram b -j8')
