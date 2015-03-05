/** \class JpsiTrkTreeMaker
 *  Class to measure trigger efficiencies (very rough)
 *
 KStar     = K+ pi-
 antiKStar = K- pi+

 */

#include "JpsiTrkTreeMaker.h"

/// default constructor
JpsiTrkTreeMaker::JpsiTrkTreeMaker(const edm::ParameterSet& cfg): 
  vertexes_           (cfg.getParameter<edm::InputTag>("vertexes")), 
  muons_              (cfg.getParameter<edm::InputTag>("muons")), 
  tracks_             (cfg.getParameter<edm::InputTag>("tracks")), 
  triggerProcess_     (cfg.getParameter<std::string>("triggerProcess")), 
  tagTriggerName_     (cfg.getParameter<std::string>("tagTriggerName")), 
  triggerName_        (cfg.getParameter<std::string>("triggerName")), 
  probeFilterDen_     (cfg.getParameter<std::string>("probeFilterDen")),
  probeFilterNum_     (cfg.getParameter<std::string>("probeFilterNum")),
  nMaxMuons_          (cfg.getUntrackedParameter<unsigned int>("maxNumberMuons", 999999)),
  l3candTag_          (cfg.getParameter<edm::InputTag>("L3CandidatesLabel")),
  VtxTag_             (cfg.getParameter<edm::InputTag>("VtxLabel")), 
  MumuVtxTag_         (cfg.getParameter<edm::InputTag>("MuMuVtxLabel")), 
  MumuVtxProdTag_     (cfg.getParameter<edm::InputTag>("MuMuVtxProdLabel")), 
  beamSpot_           (cfg.getUntrackedParameter<edm::InputTag>("beamspot")), 
  maxEta_             (cfg.getUntrackedParameter<double>("maxEta")), 
  minPtTrk_           (cfg.getUntrackedParameter<double>("minPtTrk")), 
  thirdTrackMass_     (cfg.getUntrackedParameter<double>("thirdTrkMass")),  //kaon mass
  fourthTrackMass_    (cfg.getUntrackedParameter<double>("fourthTrkMass")), //pion mass
  minInvMass_         (cfg.getUntrackedParameter<double>("minInvMass")), 
  maxInvMass_         (cfg.getUntrackedParameter<double>("maxInvMass")), 
  minTrkTrkMass_      (cfg.getUntrackedParameter<double>("minTrkTrkMass")), 
  maxTrkTrkMass_      (cfg.getUntrackedParameter<double>("maxTrkTrkMass")), 
  minCandPt_          (cfg.getUntrackedParameter<double>("minCandPt")), 
  minCandLxy_         (cfg.getUntrackedParameter<double>("minCandLxy")), 
  minCandCos_         (cfg.getUntrackedParameter<double>("minCandCos")), 
  maxNormChi2_        (cfg.getUntrackedParameter<double>("maxNormChi2")), 
  minVtxProb_         (cfg.getUntrackedParameter<double>("minVtxProb")), 
  mind0Sign_          (cfg.getUntrackedParameter<double>("mind0Sign")), 
  skimJpsiDisplaced_  (cfg.getUntrackedParameter<bool>("DisplacedJpsiRequirements")),
  debug_              (cfg.getUntrackedParameter<bool>("Debug")) 
{}


void JpsiTrkTreeMaker::endJob() {
}


void JpsiTrkTreeMaker::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {

  bool changed = true;
  if( hltConfig_.init(run, eventSetup, triggerProcess_, changed) ) {
  }
  else {
    std::cout << "Warning, didn't find process " << triggerProcess_.c_str() << std::endl;
    assert(false);
  }

  triggerIndex_    = -1; 
  tagTriggerIndex_ = -1; 

  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    std::string tempName = hltConfig_.triggerName(iHltPath);
    if(tempName.compare(triggerName_   ) == 0)      triggerIndex_    = int(iHltPath);
    if(tempName.compare(tagTriggerName_) == 0)      tagTriggerIndex_ = int(iHltPath);
    if( triggerIndex_>-1 && tagTriggerIndex_>-1 )     break; 
  } // end for each path

  if( triggerIndex_ == -1 )     std::cout << "Warning, didn't find trigger "     <<  triggerName_.c_str()    << std::endl;
  if( tagTriggerIndex_ == -1 )  std::cout << "Warning, didn't find tag trigger " <<  tagTriggerName_.c_str() << std::endl;
}

void JpsiTrkTreeMaker::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void JpsiTrkTreeMaker::analyze(const edm::Event &event, const edm::EventSetup &eventSetup) {

  using reco::Muon;
  clearVar();

  const double thirdTrackMass2  = thirdTrackMass_ *thirdTrackMass_ ; // kaon
  const double fourthTrackMass2 = fourthTrackMass_*fourthTrackMass_; // pion
  
  if (skimJpsiDisplaced_){
    minJpsiPt_  = 6.9;
    minJpsiCos_ = 0.9; 
    minJpsiLS_  =   3;
    minJpsiCL_  = 0.1;
  }

  // Handle to the offline PV collection
  edm::Handle<reco::VertexCollection> pvHandle; 
  event.getByLabel(vertexes_, pvHandle);
  const reco::VertexCollection & vertices = *pvHandle.product();

  // Get trigger results
  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(edm::InputTag("TriggerResults", "", triggerProcess_), triggerResults);
  if(!triggerResults.isValid()) {
    std::cout << "Trigger results not valid" << std::endl;
    return;
  } 

  //Print trigger in the event
  if (debug_){
    edm::TriggerNames triggerNames_ = event.triggerNames(*triggerResults);
    int ntrigs=triggerResults->size();
    std::vector<std::string> triggernames = triggerNames_.triggerNames();
    for(int itrig = 0; itrig != ntrigs; ++itrig)
    {
      std::cout << triggerNames_.triggerName(itrig) << std::endl;
    }
  }

  if( !triggerResults->accept(tagTriggerIndex_) ) return; // there are no tags
  hists_["countEvents"] -> Fill(counter++);

  // Handle to the offline muon collection
  edm::Handle<std::vector<Muon> > muons;
  event.getByLabel(muons_, muons);
  if( nMaxMuons_>0 && muons->size()>nMaxMuons_ ) return; 

  // Handle to the offline track collection
  edm::Handle<reco::TrackCollection > tracks;
  event.getByLabel(tracks_, tracks);

  // Handle to the online track collection
//   edm::Handle<reco::TrackCollection > onlineTracks;
//   event.getByLabel(onlineTracks_, onlineTracks);

  //get the transient track builder:
  edm::ESHandle<TransientTrackBuilder> theB;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  //get offline beamspot position
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  event.getByLabel(beamSpot_,recoBeamSpotHandle);
  const reco::BeamSpot& vertexBeamSpot = *recoBeamSpotHandle;

  //get the b field
  std::string mfName_ = "";
  edm::ESHandle<MagneticField> bFieldHandle;
  eventSetup.get<IdealMagneticFieldRecord>().get(mfName_, bFieldHandle);  
  const MagneticField* magField = bFieldHandle.product();
  TSCBLBuilderNoMaterial blsBuilder;
  
  //Handle to the online track collection
  edm::Handle<reco::RecoChargedCandidateCollection> trkcands;
  event.getByLabel("hltJpsiTkAllConeTracksIter", trkcands);
  
  // Handle to the online dimuon vtx collection
  edm::Handle<reco::VertexCollection> mumuVtx; 
  event.getByLabel(MumuVtxTag_, mumuVtx);

  // Handle to the online dimuon+trk vtx collection
  edm::Handle<reco::VertexCollection> hltVertexHandle; 
  try   { event.getByLabel(MumuVtxProdTag_, hltVertexHandle);}
  catch (...) { std::cout << "online vtx collection not found" << std::endl; return;}
  const reco::VertexCollection & hltVertices = *hltVertexHandle.product();
   
  // Handle to the online d0 vtx collection
  edm::Handle< reco::RecoChargedCandidateDoubleMap > d0Trk1Coll; 
  event.getByLabel("hltLowMassNonResonantTkVertexProducer", "d0firstTrk",  d0Trk1Coll);
  edm::Handle< reco::RecoChargedCandidateDoubleMap > LSigmaColl; 
  event.getByLabel("hltLowMassNonResonantTkVertexProducer", "LSigma"    , LSigmaColl);
  edm::Handle< reco::RecoChargedCandidateDoubleMap > CosineColl; 
  event.getByLabel("hltLowMassNonResonantTkVertexProducer", "Cosine"    , CosineColl);
  edm::Handle< reco::RecoChargedCandidateDoubleMap > VertexColl; 
  event.getByLabel("hltLowMassNonResonantTkVertexProducer", "VertexCL"  , VertexColl);
  
  T_Run    = event.id().run();
  T_Lumi   = event.id().luminosityBlock();
  T_Event  = event.id().event();
  T_NTrks  = tracks->size(); 
  T_Nprim  = vertices.size();

  MonteCarloStudies(event);
  if (!boolGEN) return;
  
  // Get trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", triggerProcess_), triggerEvent);
  if(!triggerEvent.isValid()) { 
    std::cout << "TriggerEvent not valid" << std::endl;
    return;
  }

  // Sanity check
  assert(triggerResults->size()==hltConfig_.size());

  // Get trigger objects from trigger summary
  const trigger::TriggerObjectCollection & toc = triggerEvent->getObjects();

  // Modules in tag trigger path
  const std::vector<std::string>& tagModuleLabels(hltConfig_.moduleLabels(tagTriggerIndex_));
  assert( tagModuleLabels.size()==hltConfig_.size(tagTriggerIndex_) );
  const unsigned int tagModuleIndex( hltConfig_.size(tagTriggerIndex_)-2 ); // index of last filter (excluding HLTEndBool)
  const unsigned int tagFilterIndex( triggerEvent->filterIndex( edm::InputTag( tagModuleLabels[tagModuleIndex], "", triggerProcess_) ) );
  assert( tagFilterIndex < triggerEvent->sizeFilters() );
  const trigger::Vids & tagVids( triggerEvent->filterIds(tagFilterIndex) );
  const trigger::Keys & tagKeys( triggerEvent->filterKeys(tagFilterIndex) );
  assert( tagVids.size()==tagKeys.size() );
  const unsigned int nTagTrig(tagVids.size());

  // Modules in probe trigger path
  const std::vector<std::string>& probeModuleLabels(hltConfig_.moduleLabels(triggerIndex_));
  assert( probeModuleLabels.size()==hltConfig_.size(triggerIndex_) );
  const unsigned int lastModuleIndex(triggerResults->index(triggerIndex_));
  //Find filter different from the last one
  unsigned int probeFilterIndex = triggerEvent->sizeFilters();
  for(unsigned int j=0; j<=lastModuleIndex; ++j) { 
    if( probeFilterNum_.compare(probeModuleLabels[j])!=0 ) continue; 
    probeFilterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterNum_, "", triggerProcess_));
    break; 
   }

  assert( probeFilterIndex < triggerEvent->sizeFilters() );
  const trigger::Vids & probeVids( triggerEvent->filterIds(probeFilterIndex) );
  const trigger::Keys & probeKeys( triggerEvent->filterKeys(probeFilterIndex) );
  assert( probeVids.size()==probeKeys.size() );
  const unsigned int nProbeTrig(probeVids.size());

  // Loop muon collection
  for(std::vector<Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) 
  { 
    if( muon::isLooseMuon( (*mu1) ) && (*mu1).pt()>4 && fabs( (*mu1).eta())< 2.2) 
    {
       T_Mu1Pt  = (*mu1).pt( );
       T_Mu1Eta = (*mu1).eta();
       T_Mu1Phi = (*mu1).phi();
       
      // Is this matched to an hlt object? 
      double finTagDeltaR = 10.0; 
      for(unsigned int i=0; i!=nTagTrig; ++i) 
      {
        const trigger::TriggerObject & tagTo = toc[tagKeys[i]];
        double tmpTagDeltaR = deltaR( (*mu1), tagTo ); 
        if( tmpTagDeltaR<finTagDeltaR ) 
        {
          finTagDeltaR    = tmpTagDeltaR;
          T_hlt_Mu1Pt     = tagTo.pt(  );
          T_hlt_Mu1Eta    = tagTo.eta( );
          T_hlt_Mu1Phi    = tagTo.phi( );
          T_Mu1_hltRecoDR = finTagDeltaR; 
        }
      }

      // Go on and look for the second tag muon
      for(std::vector<Muon>::const_iterator mu2=mu1; mu2!=muons->end(); ++mu2) {
        if( mu2==mu1 ) continue; 
        if( muon::isLooseMuon( (*mu2)) && (*mu2).pt()>4 && fabs( (*mu2).eta())< 2.2) 
        {
          if(!( mu1->charge()*mu2->charge()<0 ))         continue; 
          
          T_MuMuMass = (mu1->p4()+mu2->p4()).mass();
          hists_["mumuMass_all"]->Fill( T_MuMuMass );
          T_Mu2Pt  = (*mu2).pt( );
          T_Mu2Eta = (*mu2).eta();
          T_Mu2Phi = (*mu2).phi();

//           if(!( mumuMass > 2.85 && mumuMass < 3.35) )    continue;
          
          //Second muon found
          double finTag2DeltaR = 10.0; 
          for(unsigned int i=0; i!=nTagTrig; ++i) 
          {
            const trigger::TriggerObject & tagTo = toc[tagKeys[i]];
            double tmpTag2DeltaR = deltaR( (*mu2), tagTo ); 
            if( tmpTag2DeltaR < finTag2DeltaR ) 
            {
              finTag2DeltaR   = tmpTag2DeltaR;
              T_hlt_Mu2Pt     = tagTo.pt(  );
              T_hlt_Mu2Eta    = tagTo.eta( );
              T_hlt_Mu2Phi    = tagTo.phi( );
              T_Mu2_hltRecoDR = finTag2DeltaR; 
            }
          }

          // do jpsi vertex fit
          vector<reco::TransientTrack> j_tks;
          j_tks.push_back((*theB).build(mu1->track().get()));
          j_tks.push_back((*theB).build(mu2->track().get()));
          if (j_tks.size()!=2) continue;
       
          KalmanVertexFitter jkvf;
          TransientVertex jtv = jkvf.vertex(j_tks);
          if (!jtv.isValid()) continue;
        
          reco::Vertex jpsivertex = jtv;
          if( (jpsivertex.chi2()>=0.0) && (jpsivertex.ndof()>0) ) T_MuMuCL = TMath::Prob(jpsivertex.chi2(), jpsivertex.ndof() );

          // calculate three-track transverse momentum
          math::XYZVector jpperp(mu1->px() + mu2->px() ,
                                 mu1->py() + mu2->py() ,
                                 0.);
         
          GlobalPoint jVertex = jtv.position();
          GlobalError jerr    = jtv.positionError();
          //calculate decay length  significance w.r.t. the beamspot
          GlobalPoint displacementFromBeamspotJpsi( -1*((vertexBeamSpot.x0() -jVertex.x()) +  (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), -1*((vertexBeamSpot.y0() - jVertex.y())+ (jVertex.z() -vertexBeamSpot.z0()) * vertexBeamSpot.dydz()), 0);
          reco::Vertex::Point vperpj(displacementFromBeamspotJpsi.x(),displacementFromBeamspotJpsi.y(),0.);

          T_JpsiL          = displacementFromBeamspotJpsi.perp();
          T_JpsiSigma      = sqrt(jerr.rerr(displacementFromBeamspotJpsi));
          T_JpsiPt         = jpperp.R();
          T_JpsiCosBS      = vperpj.Dot(jpperp)/(vperpj.R()*jpperp.R());
          T_JpsiPosition_x = jVertex.x();
          T_JpsiPosition_y = jVertex.y();
          T_JpsiPosition_z = jVertex.z();
         
          T_JpsiVtx_hltReco_dx = 10;
          T_JpsiVtx_hltReco_dy = 10;
          T_JpsiVtx_hltReco_dz = 10;
          for (unsigned vtxIt =0 ;  vtxIt < mumuVtx->size(); vtxIt++){
            float idx = mumuVtx->at(vtxIt).position().x() -  jtv.position().x();
            float idy = mumuVtx->at(vtxIt).position().y() -  jtv.position().y();
            float idz = mumuVtx->at(vtxIt).position().z() -  jtv.position().z();
            if (fabs(idx) < fabs(T_JpsiVtx_hltReco_dx)  &&  fabs(idy) < fabs(T_JpsiVtx_hltReco_dy)){
              T_JpsiVtx_hltReco_dx      = idx;
              T_JpsiVtx_hltReco_dy      = idy;
              T_JpsiVtx_hltReco_dz      = idz;
              T_hlt_JpsiPosition_x      =  mumuVtx->at(vtxIt).position().x();
              T_hlt_JpsiPosition_y      =  mumuVtx->at(vtxIt).position().y();
              T_hlt_JpsiPosition_z      =  mumuVtx->at(vtxIt).position().z();
              T_hlt_JpsiPositionError_x =  mumuVtx->at(vtxIt).xError();
              T_hlt_JpsiPositionError_y =  mumuVtx->at(vtxIt).yError();
              T_hlt_JpsiPositionError_z =  mumuVtx->at(vtxIt).zError();
            }
          }
          
          hists_["pullX_mumu"  ]->Fill( T_JpsiVtx_hltReco_dx/T_hlt_JpsiPositionError_x );
          hists_["pullY_mumu"  ]->Fill( T_JpsiVtx_hltReco_dy/T_hlt_JpsiPositionError_y );
          hists_["pullZ_mumu"  ]->Fill( T_JpsiVtx_hltReco_dz/T_hlt_JpsiPositionError_z );

          hists_["jpsiPt"  ]->Fill( T_JpsiPt            );
          hists_["JpsiLS"  ]->Fill( T_JpsiL/T_JpsiSigma );
          hists_["JpsiCos" ]->Fill( T_JpsiCosBS         );
          hists_["JpsiCL"  ]->Fill( T_MuMuCL            );
         
          if (T_JpsiPt            < minJpsiPt_    ) continue;
          if (T_JpsiL/T_JpsiSigma < minJpsiLS_    ) continue;
          if (T_JpsiCosBS         < minJpsiCos_   ) continue;
          if (T_MuMuCL            < minJpsiCL_    ) continue;

		  if ( mu1->charge() > 0 ){
			T_MatchMu1 = deltaR(mu1->eta(),mu1->phi(),p_muP.Eta(),p_muP.Phi());
			T_MatchMu2 = deltaR(mu2->eta(),mu2->phi(),p_muM.Eta(),p_muM.Phi());
		  }
		  else {
			T_MatchMu1 = deltaR(mu2->eta(),mu2->phi(),p_muP.Eta(),p_muP.Phi());
			T_MatchMu2 = deltaR(mu1->eta(),mu1->phi(),p_muM.Eta(),p_muM.Phi());
		  }

          // Loop on track collection - trk 1
          for (unsigned tracksIt =0 ;  tracksIt < tracks->size(); tracksIt++)
          {
            reco::Track itrk1       = tracks->at(tracksIt) ;                                                
            
            //charge = +1 
            if (itrk1.charge() != 1)    continue;
            if (overlap(*mu1,itrk1))    continue;
            if (overlap(*mu2,itrk1))    continue;
            
			if      (itrk1.charge() == T_PiCharge)   T_MatchTrkP = deltaR(itrk1.eta(), itrk1.phi(), p_pi.Eta(), p_pi.Phi()); 
			else if (itrk1.charge() == T_KCharge)    T_MatchTrkP = deltaR(itrk1.eta(), itrk1.phi(), p_k.Eta(),  p_k.Phi() ); 
			if (T_MatchTrkP > 0.2 ) continue;

            T_TrkPPt  = itrk1.pt( ) ;
            T_TrkPEta = itrk1.eta() ;
            T_TrkPPhi = itrk1.phi() ;
            T_TrkPHQ  = 0;
            if (itrk1.quality(reco::TrackBase::highPurity)) T_TrkPHQ = 1;
            
            FreeTrajectoryState InitialFTS = initialFreeState(itrk1, magField);
            TrajectoryStateClosestToBeamLine tscb( blsBuilder(InitialFTS, *recoBeamSpotHandle) );
            T_TrkPd0Sign = tscb.transverseImpactParameter().significance();

            hists_["trkPt"] -> Fill( T_TrkPPt   );
            hists_["D0sig"] -> Fill(T_TrkPd0Sign);
 
            // eta and pt cut
            if (fabs(itrk1.eta()) > maxEta_   )                  continue;
            if (itrk1.pt()        < minPtTrk_ )                  continue;
            if (T_TrkPd0Sign      < mind0Sign_)                  continue;
            if (! itrk1.quality(reco::TrackBase::highPurity))    continue;

            for (unsigned tracksIt2 = tracksIt + 1 ;  tracksIt2 < tracks->size(); tracksIt2++)
            {
              reco::Track itrk2       = tracks->at(tracksIt2) ;                                                
              //charge = -1 
              if (itrk2.charge() != -1)                continue;
              if (itrk2.charge()*itrk1.charge() != -1) continue;
              if (overlap(*mu1,itrk2))                 continue;
              if (overlap(*mu2,itrk2))                 continue;

              if      (itrk2.charge() == T_PiCharge)   T_MatchTrkM = deltaR(itrk2.eta(), itrk2.phi(), p_pi.Eta(), p_pi.Phi()); 
              else if (itrk2.charge() == T_KCharge)    T_MatchTrkM = deltaR(itrk2.eta(), itrk2.phi(), p_k.Eta(),  p_k.Phi() ); 
              if (T_MatchTrkM > 0.2 ) continue;

        	  T_TrkMPt  = itrk2.pt( ) ;
			  T_TrkMEta = itrk2.eta() ;
			  T_TrkMPhi = itrk2.phi() ;
			  T_TrkMHQ  = 0;
			  if (itrk2.quality(reco::TrackBase::highPurity)) T_TrkMHQ = 1;

			  if ( (*d0Trk1Coll).size() != (*LSigmaColl).size())  std::cout << "different sizes" << std::endl;

              FreeTrajectoryState InitialFTS2 = initialFreeState(itrk2, magField);
              TrajectoryStateClosestToBeamLine tscb2( blsBuilder(InitialFTS2, *recoBeamSpotHandle) );
              T_TrkMd0Sign = tscb2.transverseImpactParameter().significance();
              hists_["trkPt"] -> Fill(T_TrkMPt    );
              hists_["D0sig"] -> Fill(T_TrkMd0Sign);
              // eta and pt cut
              if (fabs(itrk2.eta()) > maxEta_   )                 continue;
              if (itrk2.pt()        < minPtTrk_ )                 continue;
              if (T_TrkMd0Sign      < mind0Sign_)                 continue;
              if (! itrk2.quality(reco::TrackBase::highPurity))   continue;

              double e1,e2;
              double e3_k, e4_p; //kstar hp 
              double e3_p, e4_k; //kstar b hp
              reco::Particle::LorentzVector pB, pbarB, p1, p2, p3_k, p4_p, p3_p, p4_k, pKstar, pKstarBar, pJ;

              // Combined system
              e1   = sqrt(mu1->momentum().Mag2()  + MuMass2          );
              e2   = sqrt(mu2->momentum().Mag2()  + MuMass2          );
              e3_k = sqrt(itrk1.momentum().Mag2() + thirdTrackMass2  );
              e3_p = sqrt(itrk1.momentum().Mag2() + fourthTrackMass2 );
              e4_k = sqrt(itrk2.momentum().Mag2() + thirdTrackMass2  );
              e4_p = sqrt(itrk2.momentum().Mag2() + fourthTrackMass2 );
            
              p1   = reco::Particle::LorentzVector(mu1->px() , mu1->py() , mu1->pz() , e1  );
              p2   = reco::Particle::LorentzVector(mu2->px() , mu2->py() , mu2->pz() , e2  );
              p3_k = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3_k);
              p3_p = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3_p);
              p4_k = reco::Particle::LorentzVector(itrk2.px(), itrk2.py(), itrk2.pz(), e4_k);
              p4_p = reco::Particle::LorentzVector(itrk2.px(), itrk2.py(), itrk2.pz(), e4_p);
            
              pB        = p1 + p2 + p3_k + p4_p;
              pbarB     = p1 + p2 + p3_p + p4_k;
              pJ        = p1 + p2;
              pKstar    = p3_k + p4_p;
              pKstarBar = p3_p + p4_k;
            
              T_KStarMass    = pKstar.mass()    ;
              T_barKStarMass = pKstarBar.mass() ;
              T_KStarPt      = pKstar.pt()      ;
              T_KStarEta     = pKstar.eta()     ;
              T_KStarPhi     = pKstar.phi()     ;
              T_B0Mass       = pB.mass()        ;
              T_barB0Mass    = pbarB.mass()     ;
              T_B0StarPt     = pB.pt()          ;
              T_B0StarEta    = pB.eta()         ;
              T_B0StarPhi    = pB.phi()         ;
              
              
              // do the vertex fit
              std::vector<reco::TransientTrack> t_tks;
              t_tks.push_back((*theB).build(mu1->track().get()));
              t_tks.push_back((*theB).build(mu2->track().get()));
              t_tks.push_back((*theB).build(&itrk1));
              t_tks.push_back((*theB).build(&itrk2));
              if (t_tks.size()!=4) continue;
            
              reco::TransientTrack trkTrn = (*theB).build(&itrk1); 

              KalmanVertexFitter kvf;
              TransientVertex tv  = kvf.vertex(t_tks);
              reco::Vertex vertex = tv;
              if (!tv.isValid()) continue;
              hists_["B0InvMass"]->Fill(T_B0Mass);

              if ((vertex.chi2()>=0.0) && (vertex.ndof()>0) )   T_JpsiTkTkCL = TMath::Prob(vertex.chi2(), vertex.ndof() );

			  T_JpsiTkTkVtx_hltReco_dx = 10;
			  T_JpsiTkTkVtx_hltReco_dy = 10;
			  T_JpsiTkTkVtx_hltReco_dz = 10;
			  GlobalPoint hlt_secondaryVertex;
			  GlobalError hlt_secondaryVertexErr;
              for (unsigned vtxIt = 0 ;  vtxIt < hltVertices.size(); vtxIt++){
                float idx = hltVertices.at(vtxIt).position().x() - tv.position().x();
                float idy = hltVertices.at(vtxIt).position().y() - tv.position().y();
                float idz = hltVertices.at(vtxIt).position().z() - tv.position().z();
                if (fabs(idx) < fabs(T_JpsiTkTkVtx_hltReco_dx)  &&  fabs(idy) < fabs(T_JpsiTkTkVtx_hltReco_dy)){
				  T_JpsiTkTkVtx_hltReco_dx      = idx;
				  T_JpsiTkTkVtx_hltReco_dy      = idy;
				  T_JpsiTkTkVtx_hltReco_dz      = idz;
				  T_hlt_JpsiTkPosition_x      = hltVertices.at(vtxIt).position().x();
				  T_hlt_JpsiTkPosition_y      = hltVertices.at(vtxIt).position().y();
				  T_hlt_JpsiTkPosition_z      = hltVertices.at(vtxIt).position().z();
				  T_hlt_JpsiTkPositionError_x = hltVertices.at(vtxIt).xError();
				  T_hlt_JpsiTkPositionError_y = hltVertices.at(vtxIt).yError();
				  T_hlt_JpsiTkPositionError_z = hltVertices.at(vtxIt).zError();
                  T_hlt_JpsiTkCL = TMath::Prob(hltVertices.at(vtxIt).chi2(), hltVertices.at(vtxIt).ndof() );
//                   hlt_secondaryVertex    = hltVertices.at(vtxIt).position();
//                   hlt_secondaryVertexErr = hltVertices.at(vtxIt).positionError();
                }
              }

              hists_["pullX_4trks" ]->Fill( T_JpsiTkTkVtx_hltReco_dx/T_hlt_JpsiTkPositionError_x );
              hists_["pullY_4trks" ]->Fill( T_JpsiTkTkVtx_hltReco_dy/T_hlt_JpsiTkPositionError_y );
              hists_["pullZ_4trks" ]->Fill( T_JpsiTkTkVtx_hltReco_dz/T_hlt_JpsiTkPositionError_z );
              hists_["hltVtxProb"  ]->Fill( T_hlt_JpsiTkCL      );
              

              // calculate three-track transverse momentum
              math::XYZVector pperp(mu1->px() + mu2->px() + itrk1.px()+ itrk2.px(),
                                    mu1->py() + mu2->py() + itrk1.py()+ itrk2.py(),
                                    0.);
              // get vertex position and error to calculate the decay length significance
              GlobalPoint secondaryVertex = tv.position();
              GlobalError err = tv.positionError();
              GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0() - secondaryVertex.x()) + (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), 
                                                    -1*((vertexBeamSpot.y0() - secondaryVertex.y()) + (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()), 0);
              reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
              
			  T_JpsiTkTkL          = displacementFromBeamspot.perp();
			  T_JpsiTkTkSigma      = sqrt(err.rerr(displacementFromBeamspot));
              T_JpsiTkTkCosBS      = vperp.Dot(pperp)/(vperp.R()*pperp.R());
			  T_JpsiTkTkPosition_x = secondaryVertex.x();
			  T_JpsiTkTkPosition_y = secondaryVertex.y();
			  T_JpsiTkTkPosition_z = secondaryVertex.z();

              hists_["MatchTrkP" ] -> Fill(T_MatchTrkP           );            
              hists_["MatchTrkM" ] -> Fill(T_MatchTrkM           );            
              
              double finProbe1DeltaR = 1.0; 
              for(unsigned int i=0; i!=nProbeTrig; ++i) 
              {
                const trigger::TriggerObject & probeTo = toc[probeKeys[i]];
                double tmpProbeDeltaR = deltaR( (itrk1), probeTo ); 
                if( tmpProbeDeltaR < finProbe1DeltaR ) 
                {
                  finProbe1DeltaR  = tmpProbeDeltaR ;
                  T_TrkP_hltRecoDR = finProbe1DeltaR;
                  T_hlt_TrkPPt     = probeTo.pt(   );
                  T_hlt_TrkPEta    = probeTo.eta(  );
                  T_hlt_TrkPPhi    = probeTo.phi(  );
                }
              }
   			  if (T_hlt_TrkPPt!=0){
   			    for(unsigned int il3=0; il3< (*d0Trk1Coll).size(); ++il3) {
   			      if ((*d0Trk1Coll).at(il3).first->pt() != T_hlt_TrkPPt) continue;
   			      T_hlt_TrkPd0Sign = d0Trk1Coll->at(il3).second;
   			    }
   			    for(unsigned int il3=0; il3< (*LSigmaColl).size(); ++il3) {
   			      if ((*LSigmaColl).at(il3).first->pt() != T_hlt_TrkPPt) continue;
   			      T_hlt_JpsiTkPLSigma = LSigmaColl->at(il3).second;
   			    }
   			    for(unsigned int il3=0; il3< (*CosineColl).size(); ++il3) {
   			      if ((*CosineColl).at(il3).first->pt() != T_hlt_TrkPPt) continue;
   			      T_hlt_JpsiTkPCosBS = CosineColl->at(il3).second;
   			    }
   			    for(unsigned int il3=0; il3< (*VertexColl).size(); ++il3) {
   			      if ((*VertexColl).at(il3).first->pt() != T_hlt_TrkPPt) continue;
   			      T_hlt_JpsiTkPCLTest = VertexColl->at(il3).second;
   			    }
   			  }  

              double finProbe2DeltaR = 1.0; 
              for(unsigned int i=0; i!=nProbeTrig; ++i) 
              {
                const trigger::TriggerObject & probeTo = toc[probeKeys[i]];
                double tmpProbeDeltaR = deltaR( (itrk2), probeTo ); 
                if( tmpProbeDeltaR < finProbe2DeltaR ) 
                {
                  finProbe2DeltaR  = tmpProbeDeltaR ;
                  T_TrkM_hltRecoDR = finProbe2DeltaR;
                  T_hlt_TrkMPt     = probeTo.pt(   );
                  T_hlt_TrkMEta    = probeTo.eta(  );
                  T_hlt_TrkMPhi    = probeTo.phi(  );
                }
              }
   			  if (T_hlt_TrkMPt!=0){
   			    for(unsigned int il3=0; il3< (*d0Trk1Coll).size(); ++il3) {
   			      if ((*d0Trk1Coll).at(il3).first->pt() != T_hlt_TrkMPt) continue;
   			      T_hlt_TrkMd0Sign = d0Trk1Coll->at(il3).second;
   			    }
   			    for(unsigned int il3=0; il3< (*LSigmaColl).size(); ++il3) {
   			      if ((*LSigmaColl).at(il3).first->pt() != T_hlt_TrkMPt) continue;
   			      T_hlt_JpsiTkMLSigma = LSigmaColl->at(il3).second;
   			    }
   			    for(unsigned int il3=0; il3< (*CosineColl).size(); ++il3) {
   			      if ((*CosineColl).at(il3).first->pt() != T_hlt_TrkMPt) continue;
   			      T_hlt_JpsiTkMCosBS = CosineColl->at(il3).second;
   			    }
   			    for(unsigned int il3=0; il3< (*VertexColl).size(); ++il3) {
   			      if ((*VertexColl).at(il3).first->pt() != T_hlt_TrkMPt) continue;
   			      T_hlt_JpsiTkMCLTest = VertexColl->at(il3).second;
   			    }
   			  }  
              
              outTree_ -> Fill();            
            }
          }
        }
      }
    }
  }
}
void JpsiTrkTreeMaker::MonteCarloStudies(const edm::Event& iEvent)
{
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  int BcplusId =    511;
  int JpsiId   =    443;
  int kstar    =    313;
  int kplus    =    321;
  int piplus   =    211;

  for ( size_t i=0; i< genParticles->size(); ++i) 
  { 
    const reco::GenParticle &p = (*genParticles)[i];
    int id = p.pdgId();
    hists_["PDG_id"] -> Fill(p.pdgId());
    if(fabs(id) != BcplusId )     continue; 

    // return all daughters of B 
    if (p.numberOfDaughters()!=2 ) continue;  
    bool boolJpsi = false;
    bool boolpi   = false;
    bool boolk    = false;
    for ( size_t ides=0; ides < p.numberOfDaughters(); ++ides ) 
    {
      const reco::Candidate *des = p.daughter(ides);
      int dauId = des->pdgId();
      hists_["dau_id"] -> Fill(des->pdgId());
      if( dauId == JpsiId ) 
      {
        boolJpsi = true;
        T_gen_JpsiPt = des->pt();
        // return all daughters of J/psi    
        if (des->numberOfDaughters()!=2)          continue;  
        for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
        {
          const reco::Candidate *mumu = des->daughter(imumu);
          if      ( mumu->pdgId() == -13 )   p_muP.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
          else if ( mumu->pdgId() ==  13 )   p_muM.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
          else                               continue;         
        }
      } 
      if( fabs(dauId) == kstar)  
      {
        for ( size_t ikstar=0; ikstar < des->numberOfDaughters(); ++ikstar ) 
        {
          const reco::Candidate *pi = des->daughter(ikstar);
          if( fabs(pi->pdgId()) == kplus)  
          {
            p_k.SetPxPyPzE(pi->px(),pi->py(),pi->pz(),pi->energy());
            boolk = true;
            T_KCharge  = pi->charge();
            T_gen_KPt  = pi->pt(    );
            T_gen_KEta = pi->eta(   );
            hists_["kaonPt"    ] -> Fill(T_gen_KPt     );
            hists_["kaonCharge"] -> Fill(T_KCharge     );
          }
          else if( fabs(pi->pdgId()) == piplus)  
          {
            p_pi.SetPxPyPzE(pi->px(),pi->py(),pi->pz(),pi->energy());
            boolpi = true;
            T_PiCharge = pi->charge();
            T_gen_KPt  = pi->pt(    );
            T_gen_KEta = pi->eta(   );
            hists_["pionPt"    ] -> Fill(T_gen_PiPt  );
            hists_["pionCharge"] -> Fill(T_PiCharge  );
          }
        }
      }//end if dau id = kstar                
    } // end for des
    if (!(boolJpsi && boolpi && boolk)) continue;
    boolGEN = true;
    p_B.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
    p_jpsi = p_muP + p_muM;
    T_gen_BPt = p_B.Pt();    
  }
}





void JpsiTrkTreeMaker::beginJob() {

  TH1::SetDefaultSumw2() ;
  
  outTree_ = outfile_->make<TTree>("outTree_", "outTree_");
  outTree_->Branch("Run"                           , &T_Run                         );
  outTree_->Branch("Lumi"                          , &T_Lumi                        );
  outTree_->Branch("Event"                         , &T_Event                       );
  outTree_->Branch("NTrks"                         , &T_NTrks                       );
  outTree_->Branch("Nprim"                         , &T_Nprim                       );
  outTree_->Branch("TrueNI"                        , &T_TrueNI                      );
            
  outTree_->Branch("T_Mu1Pt"                       , &T_Mu1Pt                       );
  outTree_->Branch("T_Mu1Eta"                      , &T_Mu1Eta                      );
  outTree_->Branch("T_Mu1Phi"                      , &T_Mu1Phi                      );
  outTree_->Branch("T_Mu2Pt"                       , &T_Mu2Pt                       );
  outTree_->Branch("T_Mu2Eta"                      , &T_Mu2Eta                      );
  outTree_->Branch("T_Mu2Phi"                      , &T_Mu2Phi                      );
            
  outTree_->Branch("T_MuMuMass"                    , &T_MuMuMass                    );
  outTree_->Branch("T_MuMuCL"                      , &T_MuMuCL                      );
  outTree_->Branch("T_JpsiPosition_x"              , &T_JpsiPosition_x              );
  outTree_->Branch("T_JpsiPosition_y"              , &T_JpsiPosition_y              );
  outTree_->Branch("T_JpsiPosition_z"              , &T_JpsiPosition_z              );
  outTree_->Branch("T_JpsiPt"                      , &T_JpsiPt                      );
  outTree_->Branch("T_JpsiL"                       , &T_JpsiL                       );
  outTree_->Branch("T_JpsiSigma"                   , &T_JpsiSigma                   );
  outTree_->Branch("T_JpsiCosBS"                   , &T_JpsiCosBS                   );
                 
  outTree_->Branch("T_hlt_Mu1Pt"                   , &T_hlt_Mu1Pt                   );
  outTree_->Branch("T_hlt_Mu1Eta"                  , &T_hlt_Mu1Eta                  );
  outTree_->Branch("T_hlt_Mu1Phi"                  , &T_hlt_Mu1Phi                  );
  outTree_->Branch("T_hlt_Mu2Pt"                   , &T_hlt_Mu2Pt                   );
  outTree_->Branch("T_hlt_Mu2Eta"                  , &T_hlt_Mu2Eta                  );
  outTree_->Branch("T_hlt_Mu2Phi"                  , &T_hlt_Mu2Phi                  );
  outTree_->Branch("T_Mu1_hltRecoDR"               , &T_Mu1_hltRecoDR               );
  outTree_->Branch("T_Mu2_hltRecoDR"               , &T_Mu2_hltRecoDR               );
            
  outTree_->Branch("T_hlt_MuMuCL"                  , &T_hlt_MuMuCL                  );
  outTree_->Branch("T_hlt_JpsiPosition_x"          , &T_hlt_JpsiPosition_x          );
  outTree_->Branch("T_hlt_JpsiPosition_y"          , &T_hlt_JpsiPosition_y          );
  outTree_->Branch("T_hlt_JpsiPosition_z"          , &T_hlt_JpsiPosition_z          );
  outTree_->Branch("T_hlt_JpsiPositionError_x"     , &T_hlt_JpsiPositionError_x     );
  outTree_->Branch("T_hlt_JpsiPositionError_y"     , &T_hlt_JpsiPositionError_y     );
  outTree_->Branch("T_hlt_JpsiPositionError_z"     , &T_hlt_JpsiPositionError_z     );
  outTree_->Branch("T_JpsiVtx_hltReco_dx"          , &T_JpsiVtx_hltReco_dx          );
  outTree_->Branch("T_JpsiVtx_hltReco_dy"          , &T_JpsiVtx_hltReco_dy          );
  outTree_->Branch("T_JpsiVtx_hltReco_dz"          , &T_JpsiVtx_hltReco_dz          );
    
  outTree_->Branch("T_TrkPPt"                      , &T_TrkPPt                      );
  outTree_->Branch("T_TrkPEta"                     , &T_TrkPEta                     );
  outTree_->Branch("T_TrkPPhi"                     , &T_TrkPPhi                     );
  outTree_->Branch("T_TrkPHQ"                      , &T_TrkPHQ                      );
  outTree_->Branch("T_TrkPd0Sign"                  , &T_TrkPd0Sign                  ); 
  outTree_->Branch("T_TrkMPt"                      , &T_TrkMPt                      );
  outTree_->Branch("T_TrkMEta"                     , &T_TrkMEta                     );
  outTree_->Branch("T_TrkMPhi"                     , &T_TrkMPhi                     );
  outTree_->Branch("T_TrkMHQ"                      , &T_TrkMHQ                      );
  outTree_->Branch("T_TrkMd0Sign"                  , &T_TrkMd0Sign                  ); 
            
  outTree_->Branch("T_KStarMass"                   , &T_KStarMass                   ); 
  outTree_->Branch("T_barKStarMass"                , &T_barKStarMass                ); 
  outTree_->Branch("T_KStarPt"                     , &T_KStarPt                     ); 
  outTree_->Branch("T_KStarEta"                    , &T_KStarEta                    ); 
  outTree_->Branch("T_KStarPhi"                    , &T_KStarPhi                    ); 
  outTree_->Branch("T_B0Mass"                      , &T_B0Mass                      ); 
  outTree_->Branch("T_barB0Mass"                   , &T_barB0Mass                   ); 
  outTree_->Branch("T_B0StarPt"                    , &T_B0StarPt                    ); 
  outTree_->Branch("T_B0StarEta"                   , &T_B0StarEta                   ); 
  outTree_->Branch("T_B0StarPhi"                   , &T_B0StarPhi                   ); 
            
  outTree_->Branch("T_hlt_TrkCharge"               , &T_hlt_TrkCharge               );
  outTree_->Branch("T_hlt_TrkPPt"                  , &T_hlt_TrkPPt                  );
  outTree_->Branch("T_hlt_TrkMPt"                  , &T_hlt_TrkMPt                  );
  outTree_->Branch("T_hlt_TrkPEta"                 , &T_hlt_TrkPEta                 );
  outTree_->Branch("T_hlt_TrkMEta"                 , &T_hlt_TrkMEta                 );
  outTree_->Branch("T_hlt_TrkPd0Sign"              , &T_hlt_TrkPd0Sign              ); 
  outTree_->Branch("T_hlt_TrkMd0Sign"              , &T_hlt_TrkMd0Sign              ); 


  outTree_->Branch("T_JpsiTkTkPosition_x"          , &T_JpsiTkTkPosition_x          );
  outTree_->Branch("T_JpsiTkTkPosition_y"          , &T_JpsiTkTkPosition_y          );
  outTree_->Branch("T_JpsiTkTkPosition_z"          , &T_JpsiTkTkPosition_z          );
  outTree_->Branch("T_JpsiTkTkL"                   , &T_JpsiTkTkL                   );
  outTree_->Branch("T_JpsiTkTkSigma"               , &T_JpsiTkTkSigma               );
  outTree_->Branch("T_JpsiTkTkCosBS"               , &T_JpsiTkTkCosBS               );
  outTree_->Branch("T_JpsiTkTkCL"                  , &T_JpsiTkTkCL                  );

  outTree_->Branch("T_hlt_JpsiTkPosition_x"        , &T_hlt_JpsiTkPosition_x        );
  outTree_->Branch("T_hlt_JpsiTkPosition_y"        , &T_hlt_JpsiTkPosition_y        );
  outTree_->Branch("T_hlt_JpsiTkPosition_z"        , &T_hlt_JpsiTkPosition_z        );
  outTree_->Branch("T_hlt_JpsiTkPositionError_x"   , &T_hlt_JpsiTkPositionError_x   ) ;
  outTree_->Branch("T_hlt_JpsiTkPositionError_y"   , &T_hlt_JpsiTkPositionError_y   );
  outTree_->Branch("T_hlt_JpsiTkPositionError_z"   , &T_hlt_JpsiTkPositionError_z   );
  outTree_->Branch("T_JpsiTkTkVtx_hltReco_dx"      , &T_JpsiTkTkVtx_hltReco_dx      );
  outTree_->Branch("T_JpsiTkTkVtx_hltReco_dy"      , &T_JpsiTkTkVtx_hltReco_dy      );
  outTree_->Branch("T_JpsiTkTkVtx_hltReco_dz"      , &T_JpsiTkTkVtx_hltReco_dz      );

  outTree_->Branch("T_hlt_JpsiTkCL"                , &T_hlt_JpsiTkCL             	);
  outTree_->Branch("T_hlt_JpsiTkPCLTest"           , &T_hlt_JpsiTkPCLTest          	);
  outTree_->Branch("T_hlt_JpsiTkMCLTest"           , &T_hlt_JpsiTkMCLTest          	);
  outTree_->Branch("T_hlt_JpsiTkPLSigma"           , &T_hlt_JpsiTkPLSigma        	);
  outTree_->Branch("T_hlt_JpsiTkPCosBS"            , &T_hlt_JpsiTkPCosBS          	);
  outTree_->Branch("T_hlt_JpsiTkMLSigma"           , &T_hlt_JpsiTkMLSigma        	);
  outTree_->Branch("T_hlt_JpsiTkMCosBS"            , &T_hlt_JpsiTkMCosBS          	);

  outTree_->Branch("T_TrkP_hltRecoDR"              , &T_TrkP_hltRecoDR          	);
  outTree_->Branch("T_TrkM_hltRecoDR"              , &T_TrkM_hltRecoDR          	);

  outTree_->Branch("T_MatchMu1"                    , &T_MatchMu1           		    );
  outTree_->Branch("T_MatchMu2"                    , &T_MatchMu2                    );
  outTree_->Branch("T_MatchTrkP"                   , &T_MatchTrkP                   );
  outTree_->Branch("T_MatchTrkM"                   , &T_MatchTrkM                   );

  outTree_->Branch("T_KCharge"                     , &T_KCharge                     );
  outTree_->Branch("T_PiCharge"                    , &T_PiCharge                    );
  outTree_->Branch("T_gen_PiPt"                    , &T_gen_PiPt                    );
  outTree_->Branch("T_gen_KPt"                     , &T_gen_KPt                     );
  outTree_->Branch("T_gen_PiEta"                   , &T_gen_PiEta                   );
  outTree_->Branch("T_gen_KEta"                    , &T_gen_KEta                    );
  outTree_->Branch("T_gen_JpsiPt"                  , &T_gen_JpsiPt                  );
  outTree_->Branch("T_gen_BPt"                     , &T_gen_BPt                     );


  hists_["countEvents"  ] = outfile_->make<TH1F>("countEvents" , "countEvents"              ,    4,     0.,    4 );
  hists_["mumuMass_all" ] = outfile_->make<TH1F>("mumuMass_all", "mass"                     , 2000,     0.,   20 ); 

  hists_["jpsiPt"  ] = outfile_->make<TH1F>("jpsiPt", "mass"                     ,  400,     0.,   40 ); 
  hists_["JpsiLS"  ] = outfile_->make<TH1F>("JpsiLS", "mass"                     ,  400,     0.,   40 ); 
  hists_["JpsiCos" ] = outfile_->make<TH1F>("JpsiCos", "mass"                    , 2000,    -1.,   1 ); 
  hists_["JpsiCL"  ] = outfile_->make<TH1F>("JpsiCL", "mass"                     , 1000,     0.,   1 ); 


  hists_["trkPt"        ] = outfile_->make<TH1F>("trkPt"       , "pt trk"                   ,  150,     0.,   15 );
  hists_["onlineTrkPt"  ] = outfile_->make<TH1F>("onlineTrkPt" , "pt onl trk"               ,  150,     0.,   15 );
  hists_["D0sig"        ] = outfile_->make<TH1F>("D0sig"       , ""                         ,  600,    -1.,    5 ); 
  hists_["Trk1D0sigCut" ] = outfile_->make<TH1F>("Trk1D0sigCut", ""                         ,  600,    -1.,    5 ); 
  hists_["Trk2D0sigCut" ] = outfile_->make<TH1F>("Trk2D0sigCut", ""                         ,  600,    -1.,    5 ); 
  hists_["B0InvMass"    ] = outfile_->make<TH1F>("B0InvMass"   , "B0InvMass"                , 2000,     0.,   20.);

  hists_["pullX_mumu"   ] = outfile_->make<TH1F>("pullX_mumu"  , "Dimuon delta positionX"   ,  500,     -5,    5 ); 
  hists_["pullY_mumu"   ] = outfile_->make<TH1F>("pullY_mumu"  , "Dimuon delta positionY"   ,  500,     -5,    5 ); 
  hists_["pullZ_mumu"   ] = outfile_->make<TH1F>("pullZ_mumu"  , "Dimuon delta positionZ"   ,  500,   -10.,   10 ); 
  hists_["pullX_4trks"  ] = outfile_->make<TH1F>("pullX_4trks" , "mumutktk delta positionX" ,  500,     -5,    5 ); 
  hists_["pullY_4trks"  ] = outfile_->make<TH1F>("pullY_4trks" , "mumutktk delta positionX" ,  500,     -5,    5 ); 
  hists_["pullZ_4trks"  ] = outfile_->make<TH1F>("pullZ_4trks" , "mumutktk delta positionX" ,  500,     -5,    5 ); 
  hists_["hltVtxProb"   ] = outfile_->make<TH1F>("hltVtxProb"  , "hltVtxProb"               ,  100,      0,    1 ); 

  hists_["PDG_id"       ] = outfile_->make<TH1F>("PDG_id"      , "PDG_id"                   , 2000, -1000., 1000 ); 
  hists_["dau_id"       ] = outfile_->make<TH1F>("dau_id"      , "dau_id"                   , 2000, -1000., 1000 ); 
  hists_["kaonPt"       ] = outfile_->make<TH1F>("kaonPt"      , "kaonPt"                   ,  300,     0.,   30 ); 
  hists_["pionPt"       ] = outfile_->make<TH1F>("pionPt"      , "pionPt"                   ,  300,     0.,   30 ); 
  hists_["pionCharge"   ] = outfile_->make<TH1F>("pionCharge"  , "pionCharge"               ,    3,    -1.,    2 ); 
  hists_["kaonCharge"   ] = outfile_->make<TH1F>("kaonCharge"  , "kaonCharge"               ,    3,    -1.,    2 ); 
  hists_["MatchTrkP"    ] = outfile_->make<TH1F>("MatchTrkP"   , "MatchTrkP"                , 1000,     0.,   50 ); 
  hists_["MatchTrkM"    ] = outfile_->make<TH1F>("MatchTrkM"   , "MatchTrkM"                , 1000,     0.,   50 ); 

}    


void JpsiTrkTreeMaker::clearVar()
{

 T_Run                         = 9999;  
 T_Lumi                        = 9999;      
 T_Event                       = 9999;      
 T_NTrks                       = 9999;      
 T_Nprim                       = 9999;      
 T_TrueNI                      = 9999;      
 
 T_Mu1Pt                       = 9999;      
 T_Mu2Pt                       = 9999;      
 T_Mu1Eta                      = 9999;      
 T_Mu1Phi                      = 9999;      
 T_Mu2Eta                      = 9999;      
 T_Mu2Phi                      = 9999;      
 
 T_MuMuMass                    = 9999;      
 T_MuMuCL                      = 9999;      
 T_JpsiPosition_x              = 9999;      
 T_JpsiPosition_y              = 9999;      
 T_JpsiPosition_z              = 9999;      
 T_JpsiPt                      = 9999;      
 T_JpsiL                       = 9999;      
 T_JpsiSigma                   = 9999;      
 T_JpsiCosBS                   = 9999;      
 
 T_hlt_Mu1Pt                   = 9999;      
 T_hlt_Mu1Eta                  = 9999;      
 T_hlt_Mu1Phi                  = 9999;      
 T_hlt_Mu2Pt                   = 9999;      
 T_hlt_Mu2Eta                  = 9999;      
 T_hlt_Mu2Phi                  = 9999;      
 T_Mu1_hltRecoDR               = 9999;      
 T_Mu2_hltRecoDR               = 9999;      
 T_hlt_TrkPd0Sign              = 9999;
 T_hlt_TrkMd0Sign              = 9999;
 
 T_hlt_MuMuCL                  = 9999;      
 T_hlt_JpsiPosition_x          = 9999;      
 T_hlt_JpsiPosition_y          = 9999;      
 T_hlt_JpsiPosition_z          = 9999;      
 T_hlt_JpsiPositionError_x     = 9999;      
 T_hlt_JpsiPositionError_y     = 9999;      
 T_hlt_JpsiPositionError_z     = 9999;      
 T_JpsiVtx_hltReco_dx          = 9999;      
 T_JpsiVtx_hltReco_dy          = 9999;      
 T_JpsiVtx_hltReco_dz          = 9999;      
 
 T_TrkPPt                      = 9999;      
 T_TrkPEta                     = 9999;      
 T_TrkPPhi                     = 9999;      
 T_TrkPHQ                      = 9999;      
 T_TrkPd0Sign                  = 9999;      
 T_TrkMPt                      = 9999;      
 T_TrkMEta                     = 9999;      
 T_TrkMPhi                     = 9999;      
 T_TrkMHQ                      = 9999;      
 T_TrkMd0Sign                  = 9999;      
 
 T_KStarMass                   = 9999;      
 T_barKStarMass                = 9999;      
 T_KStarPt                     = 9999;      
 T_KStarEta                    = 9999;      
 T_KStarPhi                    = 9999;      
 T_B0Mass                      = 9999;      
 T_barB0Mass                   = 9999;      
 T_B0StarPt                    = 9999;      
 T_B0StarEta                   = 9999;      
 T_B0StarPhi                   = 9999;      
 
 T_hlt_TrkPPt                  = 9999;      
 T_hlt_TrkMPt                  = 9999;      
 T_hlt_TrkPEta                 = 9999;      
 T_hlt_TrkMEta                 = 9999;      
 
 T_JpsiTkTkPosition_x          = 9999;      
 T_JpsiTkTkPosition_y          = 9999;      
 T_JpsiTkTkPosition_z          = 9999;      
 T_JpsiTkTkL                   = 9999;      
 T_JpsiTkTkSigma               = 9999;      
 T_JpsiTkTkCosBS               = 9999;      
 
 T_hlt_TrkCharge               = 9999;
 
 T_hlt_JpsiTkPosition_x        = 9999;      
 T_hlt_JpsiTkPosition_y        = 9999;      
 T_hlt_JpsiTkPosition_z        = 9999;      
 T_hlt_JpsiTkPositionError_x   = 9999;      
 T_hlt_JpsiTkPositionError_y   = 9999;      
 T_hlt_JpsiTkPositionError_z   = 9999;      
 T_JpsiTkTkVtx_hltReco_dx      = 9999;      
 T_JpsiTkTkVtx_hltReco_dy      = 9999;      
 T_JpsiTkTkVtx_hltReco_dz      = 9999;      
 T_hlt_JpsiTkCL                = 9999;      
 T_hlt_JpsiTkPLSigma           = 9999;      
 T_hlt_JpsiTkPCosBS            = 9999; 
 T_hlt_JpsiTkPCLTest           = 9999;     
 T_hlt_JpsiTkMCLTest           = 9999;     
 T_hlt_JpsiTkMLSigma           = 9999;      
 T_hlt_JpsiTkMCosBS            = 9999;      
 T_TrkP_hltRecoDR              = 9999;      
 T_TrkM_hltRecoDR              = 9999;      
 T_MatchMu1                    = 9999;      
 T_MatchMu2                    = 9999;      
 T_MatchTrkP                   = 9999;      
 T_MatchTrkM                   = 9999;      
 T_KCharge                     = 9999;      
 T_PiCharge                    = 9999;      
 T_gen_PiPt                    = 9999;      
 T_gen_KPt                     = 9999;      
 T_gen_PiEta                   = 9999;      
 T_gen_KEta                    = 9999;      
 T_gen_JpsiPt                  = 9999;      
 T_gen_BPt                     = 9999;      
}




// define this as a plug-in
DEFINE_FWK_MODULE(JpsiTrkTreeMaker);
