/** \class BHltNtuples
 */      
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateDoubleAssociation.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
// #include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"



#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"
#include "TLorentzVector.h"

#include "MyTools/MuMuTrkHLT/src/ntupleTree.h"


class BHltNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

 public:
  explicit BHltNtuples(const edm::ParameterSet& cfg);
  virtual ~BHltNtuples() {};

  virtual void beginJob();
  /// everything that needs to be done after the event loop
  virtual void endJob();
  /// everything that needs to be done before each run
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done after each run
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  virtual void beginEvent();

  typedef std::map<uint, reco::Track  >  selTracksDef ;  

 private:

  void fillHlt      (const edm::Handle<edm::TriggerResults> &, 
                     const edm::Handle<trigger::TriggerEvent> &,
                     const edm::TriggerNames &,
                     const edm::Event &,
                     bool 
                    );

  void fillB0s      (const edm::Handle<reco::MuonCollection> &,
                     const edm::Handle<reco::TrackCollection> &,
                     const edm::Handle<reco::VertexCollection >& ,
                     const edm::Event      &, 
                     const edm::EventSetup & 
                    );

  void fillBp       (const edm::Handle<reco::MuonCollection> &,
                     const edm::Handle<reco::TrackCollection> &,
                     const edm::Handle<reco::VertexCollection >& ,
                     const edm::Event      &, 
                     const edm::EventSetup & 
                    );

  void fillL1Muons  (const edm::Handle<l1t::MuonBxCollection> &,
                     const edm::Event   &
                    );

  void fillHltMuons (const edm::Handle<reco::RecoChargedCandidateCollection> &,
                     const edm::Event   & 
                    );
  void fillHltTracks(const edm::Handle<reco::RecoChargedCandidateCollection> &,
                     const edm::Event   &                                     ,
                     edm::ConsumesCollector &&,
                     edm::InputTag &
                    );

  void fillHltDiMuons(const edm::Handle<reco::RecoChargedCandidateCollection> & ,
                      const edm::Event                                        & ,
                      const edm::EventSetup                                   & 
                     );

  void fillHltTkVtx (const edm::Handle<reco::VertexCollection> &,
                     const edm::Handle<reco::RecoChargedCandidateCollection> &,
                     const edm::Handle<reco::RecoChargedCandidateCollection> &,
                     const edm::Event   & 
                    );

  void fillHltMuVtx (const edm::Handle<reco::VertexCollection> &,
                     const edm::Handle<reco::RecoChargedCandidateCollection> &,
                     const edm::Event   & ,
                     const edm::EventSetup & 
                    );

  void fillGen      (const edm::Handle<reco::GenParticleCollection> &  ,
                     const edm::Event                               & 
                    );


  int                       overlap            ( const reco::Candidate&, const reco::Track    ) ;
  FreeTrajectoryState       initialFreeState   ( const reco::Track&,     const MagneticField* ) ;
  std::pair<double,double>  pionIPBeamSpot     ( reco::TransientTrack,   GlobalPoint          ) ;
  std::pair<double,double>  pionImpactParameter(reco::TransientTrack,    TransientVertex      ) ;
  
  // Trigger process
  edm::InputTag triggerResultTag_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_;
  edm::InputTag triggerSummTag_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummToken_;
  edm::InputTag tagTriggerResultTag_;
  edm::EDGetTokenT<edm::TriggerResults>   tagTriggerResultToken_;
  edm::InputTag tagTriggerSummTag_;
  edm::EDGetTokenT<trigger::TriggerEvent> tagTriggerSummToken_;

  edm::InputTag puTag_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;
  edm::InputTag offlinePVTag_;
  edm::EDGetTokenT<reco::VertexCollection> offlinePVToken_;
  edm::InputTag beamspotTag_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  edm::InputTag lumiScalerTag_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalerToken_; 
 
  edm::InputTag genTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

  edm::InputTag l1candTag_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1candToken_; 
  edm::InputTag l3candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l3candToken_;
  edm::InputTag tkcandTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> tkcandToken_;
  
  edm::InputTag mumuVtxTag_;
  edm::EDGetTokenT<reco::VertexCollection> mumuVtxToken_;
  edm::InputTag  tkVtxTag_;
  edm::EDGetTokenT<reco::VertexCollection> tkVtxToken_;

  edm::InputTag offlineMuonsTag_;
  edm::EDGetTokenT<reco::MuonCollection> offlineMuonsToken_;
  edm::InputTag offlineTksTag_;
  edm::EDGetTokenT<reco::TrackCollection> offlineTksToken_;

  edm::EDGetTokenT<reco::RecoChargedCandidateDoubleMap> d0token_;
  edm::EDGetTokenT<reco::RecoChargedCandidateDoubleMap> lsToken_;
  edm::EDGetTokenT<reco::RecoChargedCandidateDoubleMap> cosToken_;
  edm::EDGetTokenT<reco::RecoChargedCandidateDoubleMap> vertexToken_;

  /// file service
  edm::Service<TFileService> outfile_;
  ntupleEvent event_;
  TTree* outTree_;

  /// histograms
//   std::map<std::string, TH1*> hists_;
  
  unsigned int nGoodVtx; 

  const double MuMass = 0.106;
  const double MuMass2 = MuMass*MuMass;

  //displaced Jpsi requirements 
  double minJpsiPt_  = -20;
  double minJpsiCos_ = -20;
  double minJpsiLS_  = -20;
  double minJpsiCL_  = -20;

  // Services
  edm::ESHandle<MagneticField>         magneticField_;
  edm::ESHandle<Propagator>            propagator_;
  
  bool debug_;
  int counter = 0;

  double thirdTrackMass_   ;
  double fourthTrackMass_  ;
  bool   skimJpsiDisplaced_;
  bool   doOffline_;
  bool   doBplus_;
  bool   doBzero_;

  double maxEta_;
  double minPtTrk_;
  double mind0Sign_;
  
};


void BHltNtuples::endJob() {}

void BHltNtuples::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}

void BHltNtuples::endRun  (const edm::Run & run, const edm::EventSetup & eventSetup) {}

//---------------------------------------------
void BHltNtuples::beginEvent()
{

  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();

  event_.hltTag.triggers.clear();
  event_.hltTag.objects.clear();

  event_.genParticles.clear();
  event_.b0cands.clear();
  event_.bpcands.clear();

  event_.L1muons.clear();
  event_.hlt_mu.clear();
  event_.hlt_tk.clear();
  event_.hlt_dimu.clear();
  event_.hlt_muvtx.clear();
  event_.hlt_tkvtx.clear();
  
//   for (unsigned int ix=0; ix<3; ++ix) {
//     event_.primaryVertex[ix] = 0.;
//     for (unsigned int iy=0; iy<3; ++iy) {
//       event_.cov_primaryVertex[ix][iy] = 0.;
//     }
//   }
  event_.nVtx       = -1;
  event_.trueNI     = -1;
  
  nGoodVtx = 0; 
}

