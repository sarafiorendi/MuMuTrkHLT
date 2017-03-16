# MuMuTrkHLT  

Package to produce ntuples to study displaced Jpsi + track trigger  
Need some additions to the CMSSW release to have hlt quantities saved 


cmsrel CMSSW_8_0_0  
cd CMSSW_8_0_0/src/  
git cms-addpkg HLTrigger/Configuration  
git cms-addpkg DataFormats/RecoCandidate/  
git cms-addpkg HLTrigger/btau/  

git remote add bph-hlt-cmssw git@github.com:sarafiorendi/bph-hlt-cmssw.git  
git fetch bph-hlt-cmssw  
git checkout bph-hlt-cmssw/master HLTrigger/btau/src/HLTmumutkVtxProducer.cc  
git checkout bph-hlt-cmssw/master DataFormats/RecoCandidate/src/classes_def.xml  
git checkout bph-hlt-cmssw/master DataFormats/RecoCandidate/src/classes.h  
git checkout bph-hlt-cmssw/master DataFormats/RecoCandidate/interface/RecoChargedCandidateDoubleAssociation.h  
git checkout bph-hlt-cmssw/master HLTrigger/Configuration/test/customize_hlt.py  

mkdir MyTools/  
cd MyTools/  
git clone git@github.com:sarafiorendi/MuMuTrkHLT.git  
cd MuMuTrkHLT  
git checkout newntuples  

cd ..  
sc  
cd HLTrigger/Configuration/test  
cp /afs/cern.ch/user/f/fiorendi/public/BPHhlt/hlt_bph.py .  
cmsRun customize_hlt.py  
