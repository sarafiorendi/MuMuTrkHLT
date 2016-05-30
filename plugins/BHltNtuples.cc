/** \class BHltNtuples
 */
      
#include "BHltNtuples.h"
#include "BHltNtuples_utils.h"
#include "BHltNtuples_fillHlt.h"
#include "BHltNtuples_fillhltB0s.h"
#include "BHltNtuples_fillGen.h"
#include "BHltNtuples_fillB0s.h"


/// default constructor
BHltNtuples::BHltNtuples(const edm::ParameterSet& cfg): 

  triggerResultTag_       (cfg.getUntrackedParameter<edm::InputTag>("triggerResult")), 
    triggerResultToken_     (consumes<edm::TriggerResults>(triggerResultTag_)),
  triggerSummTag_         (cfg.getUntrackedParameter<edm::InputTag>("triggerSummary")), 
    triggerSummToken_       (consumes<trigger::TriggerEvent>(triggerSummTag_)),

  tagTriggerResultTag_    (cfg.getUntrackedParameter<edm::InputTag>("tagTriggerResult")), 
    tagTriggerResultToken_  (consumes<edm::TriggerResults>(tagTriggerResultTag_)),
  tagTriggerSummTag_      (cfg.getUntrackedParameter<edm::InputTag>("tagTriggerSummary")), 
    tagTriggerSummToken_    (consumes<trigger::TriggerEvent>(tagTriggerSummTag_)),

  puTag_                  (cfg.getUntrackedParameter<edm::InputTag>("puInfoTag")),
    puToken_                (consumes<std::vector< PileupSummaryInfo>>(puTag_)), 
  offlinePVTag_           (cfg.getParameter<edm::InputTag>("offlineVtx")), 
    offlinePVToken_         (consumes<reco::VertexCollection>(offlinePVTag_)), 
  beamspotTag_           (cfg.getParameter<edm::InputTag>("beamspot")), 
    beamspotToken_         (consumes<reco::BeamSpot>(beamspotTag_)), 

  genTag_                 (cfg.getUntrackedParameter<edm::InputTag>("genParticlesTag")),
    genToken_               (consumes<reco::GenParticleCollection>(genTag_)), 

  l3candTag_              (cfg.getParameter<edm::InputTag>("L3CandidatesTag")),
    l3candToken_            (consumes<reco::RecoChargedCandidateCollection>(l3candTag_)), 
  tkcandTag_              (cfg.getParameter<edm::InputTag>("TkCandidatesTag")),
    tkcandToken_            (consumes<reco::RecoChargedCandidateCollection>(tkcandTag_)), 

  mumuVtxTag_             (cfg.getUntrackedParameter<edm::InputTag>("MuMuVtxTag")),
    mumuVtxToken_           (consumes<reco::VertexCollection>(mumuVtxTag_)), 
  tkVtxTag_               (cfg.getUntrackedParameter<edm::InputTag>("TkVtxTag")),
    tkVtxToken_             (consumes<reco::VertexCollection>(tkVtxTag_)), 

  offlineMuonsTag_        (cfg.getUntrackedParameter<edm::InputTag>("OfflineMuonsTag")),
    offlineMuonsToken_      (consumes<reco::MuonCollection>(offlineMuonsTag_)), 
  offlineTksTag_          (cfg.getUntrackedParameter<edm::InputTag>("OfflineTkTag")),
    offlineTksToken_        (consumes<reco::TrackCollection>(offlineTksTag_)), 

  thirdTrackMass_         (cfg.getUntrackedParameter<double>("thirdTrkMass")),  //kaon mass
  fourthTrackMass_        (cfg.getUntrackedParameter<double>("fourthTrkMass")), //pion mass
  skimJpsiDisplaced_      (cfg.getUntrackedParameter<bool>("displacedJpsi")),
  doOffline_              (cfg.getUntrackedParameter<bool>("doOffline")),

  maxEta_                 (cfg.getUntrackedParameter<double>("maxEta")), 
  minPtTrk_               (cfg.getUntrackedParameter<double>("minPtTrk")),
  mind0Sign_              (cfg.getUntrackedParameter<double>("mind0Sign")) 
{
}


void BHltNtuples::beginJob() {

  TH1::SetDefaultSumw2() ;
  outTree_["ntupleTree"] = outfile_-> make<TTree>("ntupleTree","ntupleTree");
  outTree_["ntupleTree"] -> Branch("event" ,&event_, 64000,2);

  hists_["countEvents"  ] = outfile_->make<TH1F>("countEvents" , "countEvents"              ,    4,     0.,    4 );
  hists_["mumuMass_all" ] = outfile_->make<TH1F>("mumuMass_all", "mass"                     , 2000,     0.,   20 ); 

  hists_["JpsiPt"       ] = outfile_->make<TH1F>("JpsiPt"      , "mass"                     ,  400,     0.,   40 ); 
  hists_["JpsiLS"       ] = outfile_->make<TH1F>("JpsiLS"      , "mass"                     ,  400,     0.,   40 ); 
  hists_["JpsiCos"      ] = outfile_->make<TH1F>("JpsiCos"     , "mass"                     , 2000,    -1.,    1 ); 
  hists_["JpsiCL"       ] = outfile_->make<TH1F>("JpsiCL"      , "mass"                     , 1000,     0.,    1 ); 

  hists_["trkPt"        ] = outfile_->make<TH1F>("trkPt"       , "pt trk"                   ,  150,     0.,   15 );
  hists_["onlineTrkPt"  ] = outfile_->make<TH1F>("onlineTrkPt" , "pt onl trk"               ,  150,     0.,   15 );
  hists_["D0sig"        ] = outfile_->make<TH1F>("D0sig"       , ""                         ,  600,    -1.,    5 ); 
  hists_["B0InvMass"    ] = outfile_->make<TH1F>("B0InvMass"   , "B0InvMass"                , 2000,     0.,   20.);

}    



 
void BHltNtuples::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) {

  beginEvent();
  
  if (skimJpsiDisplaced_){
    minJpsiPt_  = 6.9;
    minJpsiCos_ = 0.9; 
    minJpsiLS_  =   3;
    minJpsiCL_  = 0.1;
  }

  // Fill general info
  event_.runNumber             = event.id().run();
  event_.luminosityBlockNumber = event.id().luminosityBlock();
  event_.eventNumber           = event.id().event();

  // Fill PU info
  if (!event.isRealData()) {
    edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
    if ( event.getByToken(puToken_,puInfo)){
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) 
      {
        if(PVI->getBunchCrossing()==0){
          event_.trueNI   = PVI->getTrueNumInteractions();
          continue;
        }
      }
    } 
    else  
      edm::LogError("") << "PU collection not found !!!";
  }

  // Fill trigger information for probe muon
  edm::Handle<edm::TriggerResults>   triggerResults;
  edm::Handle<trigger::TriggerEvent> triggerEvent;

  if (event.getByToken(triggerResultToken_, triggerResults) &&
      event.getByToken(triggerSummToken_  , triggerEvent)) {
      
    edm::TriggerNames triggerNames_ = event.triggerNames(*triggerResults);
    fillHlt(triggerResults, triggerEvent, triggerNames_, event, false);
  }
  else 
    edm::LogError("") << "Trigger collection for probe muon not found !!!";

  // Fill trigger information for tag muon
//   edm::Handle<edm::TriggerResults>   tagTriggerResults;
//   edm::Handle<trigger::TriggerEvent> tagTriggerEvent;
      
//   if (event.getByToken(tagTriggerResultToken_, tagTriggerResults) &&
//       event.getByToken(tagTriggerSummToken_  , tagTriggerEvent)) {
//       
//     edm::TriggerNames tagTriggerNames_ = event.triggerNames(*tagTriggerResults);
//     fillHlt(tagTriggerResults, tagTriggerEvent, tagTriggerNames_, event, true);
//   }
//   else 
//     edm::LogError("") << "Trigger collection for tag muon not found !!!";


  if (doOffline_){
    
    // Fill vertex info
    edm::Handle<reco::VertexCollection> vtxColl; 
    event.getByToken(offlinePVToken_, vtxColl);
    for(reco::VertexCollection::const_iterator it = vtxColl->begin(); it != vtxColl->end(); ++it) {
      if( !it->isValid())  continue;
      nGoodVtx++;
    }
    event_.nVtx = nGoodVtx;

    // Handle the offline collections and fill offline b0s
    edm::Handle<reco::MuonCollection> muons;
    event.getByToken(offlineMuonsToken_,  muons);
    edm::Handle<reco::TrackCollection >   tracks;
    event.getByToken(offlineTksToken_,    tracks);
    fillB0s(muons, tracks, event, eventSetup );

  }
  
  // Fill MC GEN info
  if (!event.isRealData()) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    event.getByLabel("genParticles", genParticles);
      fillGen(genParticles, event);
  }  
  

  // Handle the online muon collection and fill online muons
  edm::Handle<reco::RecoChargedCandidateCollection> l3cands;
  if (event.getByToken(l3candToken_, l3cands))
    fillHltMuons(l3cands, event);
  else
    edm::LogWarning("") << "Online muon collection not found !!!";


  // Handle the online track collection and fill online tracks
  edm::Handle<reco::RecoChargedCandidateCollection> tkcands;
  if (event.getByToken(tkcandToken_, tkcands))
    fillHltTracks(tkcands, event);
  else
    edm::LogWarning("") << "Online track collection not found !!!";


  // Handle the online mumu vtx collection and fill container
  edm::Handle<reco::VertexCollection> hlt_muvtx;
  if (event.getByToken(mumuVtxToken_, hlt_muvtx) && 
      event.getByToken( l3candToken_, l3cands  ) )
    fillHltMuVtx(hlt_muvtx, l3cands, event);
  else
    edm::LogWarning("") << "Online dimuon vertex collection not found !!!";


  // Handle the online mumutk vtx collection and fill container
  edm::Handle<reco::VertexCollection> hlt_tkvtx;
  if (event.getByToken( tkVtxToken_ , hlt_tkvtx) && 
      event.getByToken( l3candToken_, l3cands  )  && 
      event.getByToken( tkcandToken_, tkcands  ) )
    fillHltTkVtx(hlt_tkvtx, l3cands, tkcands, event);
  else
    edm::LogWarning("") << "Online muon-trk vertex collection not found !!!";
  
  outTree_["ntupleTree"] -> Fill();
}



// define this as a plug-in
DEFINE_FWK_MODULE(BHltNtuples);
