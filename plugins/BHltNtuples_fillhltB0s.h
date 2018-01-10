
void BHltNtuples::fillHltMuons(const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands ,
                               const edm::Event                                        & event )
{
  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) 
  {
    HLTMuCand theL3Mu;

    reco::RecoChargedCandidateRef candref(l3cands, il3);
    theL3Mu.pt      = candref -> pt();
    theL3Mu.eta     = candref -> eta();
    theL3Mu.phi     = candref -> phi();
    theL3Mu.charge  = candref -> charge();

    event_.hlt_mu.push_back(theL3Mu);
  }
}

void BHltNtuples::fillHltDiMuons(const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands ,
                                 const edm::Event                                        & event   ,
                                 const edm::EventSetup                                   & eventSetup
                                )
{
  if (l3cands->size() < 2) return;
  HLTDimuonCand theDimuon;
  
  //get the b field
  std::string mfName_ = "";
  edm::ESHandle<MagneticField> bFieldHandle;
  eventSetup.get<IdealMagneticFieldRecord>().get(mfName_, bFieldHandle);  

  TrajectoryStateClosestToPoint mu1TS, mu2TS;
  
  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) 
  {
//     HLTMuCand oneL3Mu;
    reco::RecoChargedCandidateRef oneref(l3cands, il3);
	reco::TrackRef tk1 = oneref->get<reco::TrackRef>();
	reco::TransientTrack mu1TT(*tk1, &(*bFieldHandle));
	mu1TS = mu1TT.impactPointTSCP();
    
    theDimuon.mu1pt  = oneref -> pt();
	theDimuon.mu1eta = oneref -> eta();
	theDimuon.mu1phi = oneref -> phi();
	
	theDimuon.DCA    = -20  ;

	for( unsigned int jl3 = il3+1; jl3 < l3cands->size(); ++jl3) 
	{
      reco::RecoChargedCandidateRef tworef(l3cands, jl3);
//       if ( tworef -> charge() * oneref -> charge() > 0) continue;
	  reco::TrackRef tk2 = tworef->get<reco::TrackRef>();
	  reco::TransientTrack mu2TT(*tk2, &(*bFieldHandle));
	  mu2TS = mu2TT.impactPointTSCP();
       
      theDimuon.mu2pt  = tworef -> pt();
      theDimuon.mu2eta = tworef -> eta();
      theDimuon.mu2phi = tworef -> phi();
      // DCA between the two muons
	  if (mu1TS.isValid() && mu2TS.isValid()) {
	    ClosestApproachInRPhi cApp;
	    cApp.calculate(mu1TS.theState(), mu2TS.theState());
	    if ( cApp.status()) theDimuon.DCA = cApp.distance()  ;
	    else                theDimuon.DCA = -10  ;
	  }
	  
	  theDimuon.charge =  tworef -> charge() * oneref -> charge();

      event_.hlt_dimu.push_back(theDimuon);


    }    
  }
}

void BHltNtuples::fillHltTracks(const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands ,
                                const edm::Event                                        & event    ,
                                edm::ConsumesCollector &&                                 iC       ,
                                edm::InputTag                                           & tkcandTag_,
                                const int ncoll
                                )
{

//   bool d0_valid = false;
//   edm::Handle<reco::RecoChargedCandidateDoubleMap > d0TrkColl; 
//   if   (event.getByToken(d0token_, d0TrkColl)) d0_valid = true;
//   else edm::LogWarning("") << "Online d0 collection not found !!!";
// 
//   bool LS_valid = false;
//   edm::Handle< reco::RecoChargedCandidateDoubleMap > LSigmaColl; 
//   if   (event.getByToken(lsToken_, LSigmaColl)) LS_valid = true;
//   else edm::LogWarning("") << "Online LS collection not found !!!";
// 
//   bool Cos_valid = false;
//   edm::Handle< reco::RecoChargedCandidateDoubleMap > CosineColl; 
//   if   (event.getByToken(cosToken_, CosineColl)) Cos_valid = true;
//   else edm::LogWarning("") << "Online cosine collection not found !!!";
// 
//   bool CL_valid = false;
//   edm::Handle< reco::RecoChargedCandidateDoubleMap > VertexColl; 
//   if   (event.getByToken(vertexToken_, VertexColl)) CL_valid = true;
//   else edm::LogWarning("") << "Online vertex CL collection not found !!!";
// 
// 
//   for( unsigned int itk = 0; itk < trkcands->size(); ++itk) 
//   {
//     HLTTkCand theTk;
// 
//     reco::RecoChargedCandidateRef candref(trkcands, itk);
//     theTk.pt      = candref -> pt();
//     theTk.eta     = candref -> eta();
//     theTk.phi     = candref -> phi();
//     theTk.charge  = candref -> charge();
//     
//     if (d0_valid){
//       for(unsigned int id0 = 0; id0 < (*d0TrkColl).size(); ++id0) {
//         if ((*d0TrkColl).at(id0).first->pt() != candref -> pt()) continue;
//         theTk.d0Sign = d0TrkColl->at(id0).second;
//         break;
//       }
//     }
//     if (LS_valid){
//       for(unsigned int il3 = 0; il3 < (*LSigmaColl).size(); ++il3) {
//         if ((*LSigmaColl).at(il3).first->pt() != candref -> pt()) continue;
//         theTk.LSigma = LSigmaColl->at(il3).second;
//       }
//     }    
//     if (Cos_valid){
//       for(unsigned int ic = 0; ic < (*CosineColl).size(); ++ic) {
//         if ((*CosineColl).at(ic).first->pt() != candref -> pt()) continue;
//         theTk.CosBS = CosineColl->at(ic).second;
//       }
//     }    
//     if (CL_valid){
//       for(unsigned int ic = 0; ic < (*VertexColl).size(); ++ic) {
//         if ((*VertexColl).at(ic).first->pt() != candref -> pt()) continue;
//         theTk.CL = VertexColl->at(ic).second;
//       }
//     }    
//     
//     if      (ncoll ==0) event_.hlt_tk.push_back(theTk);
//     else if (ncoll ==1) event_.hlt_glbtk.push_back(theTk);
//     else if (ncoll ==2) event_.hlt_newtk.push_back(theTk);
//   }
}



void BHltNtuples::fillHltPixTracks(const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands ,
                                   const edm::Event                                        & event    )
{
  for( unsigned int itk = 0; itk < trkcands->size(); ++itk) 
  {
    HLTTkCand theTk;

    reco::RecoChargedCandidateRef candref(trkcands, itk);
    theTk.pt      = candref -> pt();
    theTk.eta     = candref -> eta();
    theTk.phi     = candref -> phi();
    theTk.charge  = candref -> charge();
    
    event_.hlt_pix_tk.push_back(theTk);
  }
}



void BHltNtuples::fillHltTkVtx(const edm::Handle<reco::VertexCollection>               & hltVertexHandle ,
                               const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands,
                               const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands ,
                               const edm::Event                         & event )
{

  // Handle to the online dimuon+trk vtx collection
//   edm::Handle<reco::VertexCollection> hltVertexHandle; 
//   try   { event.getByLabel(MumuVtxProdTag_, hltVertexHandle);}
//   catch (...) { std::cout << "online vtx collection not found" << std::endl; return;}
  const reco::VertexCollection & hltVertices = *hltVertexHandle.product();

  HLTMuMuTkVtxCand theVtx;
  for (unsigned ivtx = 0 ;  ivtx < hltVertices.size(); ivtx++){
    theVtx.x  = hltVertices.at(ivtx).position().x();
    theVtx.y  = hltVertices.at(ivtx).position().y();
    theVtx.z  = hltVertices.at(ivtx).position().z();

    theVtx.ex = hltVertices.at(ivtx).xError();
    theVtx.ey = hltVertices.at(ivtx).yError();
    theVtx.ez = hltVertices.at(ivtx).zError();
    
    theVtx.CL   = TMath::Prob(hltVertices.at(ivtx).chi2(), hltVertices.at(ivtx).ndof() );
    theVtx.normchi2 = hltVertices.at(ivtx).normalizedChi2();
    
    theVtx.mu1pt = -1;
    theVtx.mu2pt = -1;
    theVtx.tkpt  = -1;

    // not working, missing RCCRef somehow
//     bool foundMu0 = false;
//     bool foundMu1 = false;
//     bool foundTk  = false;
//     reco::Vertex::trackRef_iterator trki;
//     for (trki  = hltVertices.at(ivtx).tracks_begin(); trki != hltVertices.at(ivtx).tracks_end(); ++trki) 
//     {
// 	   reco::RecoChargedCandidateRef tmp1Ref(l3cands, (*trki).key());
// 	   if (! (tmp1Ref -> track().isNull()) && !foundMu0 )  {
// 	     theVtx.mu1pt = tmp1Ref -> track() -> pt();
//          foundMu0 = true;
//        }
//        else if (! (tmp1Ref -> track().isNull()) && !foundMu1 ){
// 	     theVtx.mu2pt = tmp1Ref -> track() -> pt();
//          foundMu1 = true;
//        }
// 
// 	   reco::RecoChargedCandidateRef tmp2Ref(trkcands, (*trki).key());
// 	   if (! (tmp2Ref -> track().isNull()) && !foundTk )  {
// 	     theVtx.tkpt = tmp2Ref -> track() -> pt();
//          foundTk = true;
//        }
//     }
    event_.hlt_tkvtx.push_back(theVtx);
  }
    
}



void BHltNtuples::fillHltMuVtx(const edm::Handle<reco::VertexCollection>               & hltVertexHandle ,
                               const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands,
                               const edm::Event                                        & event,
                               const edm::EventSetup                                   & eventSetup)
{

  const reco::VertexCollection & hltVertices = *hltVertexHandle.product();

  HLTMuMuVtxCand theVtx;
  for (unsigned ivtx = 0 ;  ivtx < hltVertices.size(); ivtx++){
    theVtx.x  = hltVertices.at(ivtx).position().x();
    theVtx.y  = hltVertices.at(ivtx).position().y();
    theVtx.z  = hltVertices.at(ivtx).position().z();

    theVtx.ex = hltVertices.at(ivtx).xError();
    theVtx.ey = hltVertices.at(ivtx).yError();
    theVtx.ez = hltVertices.at(ivtx).zError();
    
    theVtx.CL = TMath::Prob(hltVertices.at(ivtx).chi2(), hltVertices.at(ivtx).ndof() );
    
//     bool foundMu0 = false;
//     bool foundMu1 = false;
// 
//     reco::Vertex::trackRef_iterator trki;
//     for (trki  = hltVertices.at(ivtx).tracks_begin(); trki != hltVertices.at(ivtx).tracks_end(); ++trki) 
//     {
// 	   reco::RecoChargedCandidateRef tmp1Ref(l3cands, (*trki).key());
// 	   if (! (tmp1Ref -> track().isNull()) && !foundMu0 )  {
// 	     theVtx.mu1pt = tmp1Ref -> track() -> pt();
//          foundMu0 = true;
//        }
//        else if (! (tmp1Ref -> track().isNull()) && !foundMu1 ){
// 	     theVtx.mu2pt = tmp1Ref -> track() -> pt();
//          foundMu1 = true;
//        }
//     }
    event_.hlt_muvtx.push_back(theVtx);
  }
    
}


void BHltNtuples::fillL1Muons(const edm::Handle<l1t::MuonBxCollection> & l1cands ,
                              const edm::Event                         & event    
                              )
{

  for (int ibx = l1cands->getFirstBX(); ibx <= l1cands->getLastBX(); ++ibx) {
    if (ibx != 0) continue;
    for (auto it = l1cands->begin(ibx); it != l1cands->end(ibx); it++){

      l1t::MuonRef muon(l1cands, distance(l1cands->begin(l1cands->getFirstBX()),it) );

      L1MuonCand theL1Mu;

      theL1Mu.pt       = muon -> pt();
      theL1Mu.eta      = muon -> eta();
      theL1Mu.phi      = muon -> phi();
      theL1Mu.etaAtVtx = -999;
      theL1Mu.phiAtVtx = -999;
//       theL1Mu.etaAtVtx = muon -> etaAtVtx();
//       theL1Mu.phiAtVtx = muon -> phiAtVtx();
      theL1Mu.charge   = muon -> charge();
      theL1Mu.quality  = muon -> hwQual();

      event_.L1muons.push_back(theL1Mu);
    }
  }
}


