/** \class JpsiTrkTreeMaker
 *  Class to measure trigger efficiencies (very rough)
 *
 KStar     = K+ pi-
 antiKStar = K- pi+

 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"

#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateDoubleAssociation.h"

#include <DataFormats/Math/interface/deltaPhi.h>
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <map>
#include <string>
#include <memory>
#include <iomanip>

// using namespace reco;
// using namespace edm;
// using namespace trigger;
using namespace std;


class JpsiTrkTreeMaker : public edm::EDAnalyzer {

 public:
  /// default constructor
  JpsiTrkTreeMaker(const edm::ParameterSet& cfg);
  /// default destructor
  virtual ~JpsiTrkTreeMaker() {};
  /// everything that needs to be done before the event loop


 private:
  virtual void beginJob();
  /// everything that needs to be done after the event loop
  virtual void endJob();
  /// everything that needs to be done before each run
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done after each run
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  void clearVar();

  TLorentzVector p_muP,p_muM, p_pi, p_k, p_jpsi,p_B;
  int PICH = 0;
  int KCH  = 0;
  bool boolGEN = false;

  void MonteCarloStudies(const edm::Event&);



  
  /// input tag for vertices     
  edm::InputTag vertexes_;
  /// input tag for muons
  edm::InputTag muons_;
  /// input tag for tracks
  edm::InputTag tracks_;
  /// file service
  edm::Service<TFileService> outfile_;
  /// histograms
  std::map<std::string, TH1*> hists_;
  std::map<std::string, TH2*> hist2D_;
  
  TTree *outTree_; 

  // Trigger process
  std::string triggerProcess_;
  // Trigger names
  std::string tagTriggerName_;
  std::string triggerName_;
  std::string probeFilterDen_;
  std::string probeFilterNum_;

  // Trigger indexes
  int tagTriggerIndex_;
  int triggerIndex_;
  // HLTConfig
  HLTConfigProvider hltConfig_;

  unsigned int nMaxMuons_;

  //Deposits tag & l3 muons
  edm::InputTag l3candTag_;
  edm::InputTag VtxTag_;
  edm::InputTag MumuVtxTag_;
  edm::InputTag MumuVtxProdTag_;
  edm::InputTag AllVtxTag_;
  edm::InputTag beamSpot_;

  /// Cuts
  double maxEta_;
  double minPtTrk_;
  
  double thirdTrackMass_;
  double fourthTrackMass_;
  double minInvMass_;
  double maxInvMass_;
  double minTrkTrkMass_;
  double maxTrkTrkMass_;

  double minCandPt_;
  double minCandLxy_;
  double minCandCos_;
  double maxNormChi2_;
  double minVtxProb_;
  double mind0Sign_;
  
  bool   skimJpsiDisplaced_;


  const double MuMass = 0.106;
  const double MuMass2 = MuMass*MuMass;

  //Tree variables

  int T_Run,  T_Lumi, T_Event, T_NTrks, T_Nprim, T_TrueNI; 

  //muons reco
  double T_Mu1Pt, T_Mu1Eta, T_Mu1Phi;
  double T_Mu2Pt, T_Mu2Eta, T_Mu2Phi;
  double T_MuMuMass;
  double T_MuMuCL, T_JpsiPosition_x, T_JpsiPosition_y, T_JpsiPosition_z;  
  double T_JpsiPt;
  float  T_JpsiSigma, T_JpsiL, T_JpsiCosBS;
  
  //muon hlt
  double T_hlt_Mu1Pt, T_hlt_Mu1Eta, T_hlt_Mu1Phi;
  double T_hlt_Mu2Pt, T_hlt_Mu2Eta, T_hlt_Mu2Phi;
  double T_Mu1_hltRecoDR, T_Mu2_hltRecoDR;
  double T_hlt_MuMuCL;
  double T_hlt_JpsiPosition_x,      T_hlt_JpsiPosition_y,      T_hlt_JpsiPosition_z;
  double T_hlt_JpsiPositionError_x, T_hlt_JpsiPositionError_y, T_hlt_JpsiPositionError_z;  
  double T_JpsiVtx_hltReco_dx, T_JpsiVtx_hltReco_dy, T_JpsiVtx_hltReco_dz;  

  //trks reco
  double T_TrkPPt, T_TrkPEta, T_TrkPPhi, T_TrkPd0Sign;
  int    T_TrkPHQ;
  double T_TrkMPt, T_TrkMEta, T_TrkMPhi, T_TrkMd0Sign;
  int    T_TrkMHQ;

  //trks hlt
  double T_hlt_TrkPPt, T_hlt_TrkPEta, T_hlt_TrkPPhi, T_hlt_TrkPd0Sign;
  double T_hlt_TrkMPt, T_hlt_TrkMEta, T_hlt_TrkMPhi, T_hlt_TrkMd0Sign;
  int    T_hlt_TrkCharge;
  
  //mumjutrktrk cand offline
  double T_KStarMass,  T_barKStarMass, T_KStarPt,  T_KStarEta,  T_KStarPhi;
  double T_B0Mass,     T_barB0Mass,    T_B0StarPt, T_B0StarEta, T_B0StarPhi;   
  double T_JpsiTkTkCL, T_JpsiTkTkPosition_x, T_JpsiTkTkPosition_y, T_JpsiTkTkPosition_z; 
  float  T_JpsiTkTkSigma, T_JpsiTkTkL, T_JpsiTkTkCosBS;

  //mumutrktrk cand hlt  
  double T_hlt_JpsiTkCL, T_hlt_JpsiTkPCLTest, T_hlt_JpsiTkMCLTest;
  double T_hlt_JpsiTkPosition_x,      T_hlt_JpsiTkPosition_y,      T_hlt_JpsiTkPosition_z;
  double T_hlt_JpsiTkPositionError_x, T_hlt_JpsiTkPositionError_y, T_hlt_JpsiTkPositionError_z;  
  float  T_hlt_JpsiTkPLSigma, T_hlt_JpsiTkPCosBS;
  float  T_hlt_JpsiTkMLSigma, T_hlt_JpsiTkMCosBS;
  double T_JpsiTkTkVtx_hltReco_dx, T_JpsiTkTkVtx_hltReco_dy, T_JpsiTkTkVtx_hltReco_dz;  
  double T_TrkP_hltRecoDR, T_TrkM_hltRecoDR;
  
  //gen 
  int    T_KCharge, T_PiCharge;
  float  T_MatchMu1, T_MatchMu2, T_MatchTrkP, T_MatchTrkM;  
  double T_gen_PiPt, T_gen_KPt, T_gen_PiEta, T_gen_KEta, T_gen_JpsiPt, T_gen_BPt;
  
  
  //displaced Jpsi requirements 
  double minJpsiPt_  = -20;
  double minJpsiCos_ = -20;
  double minJpsiLS_  = -20;
  double minJpsiCL_  = -20;


  // Services
  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<MuonDetLayerGeometry> detLayerGeometry_;
  
  bool debug_;
  int counter = 0;
  
  

  //useful function definition
  std::pair<double,double> pionImpactParameter(reco::TransientTrack piTT, TransientVertex jpsiVtx)
  {
    std::pair<double,double> measure;
    std::pair<bool,Measurement1D>  piIP_pair = IPTools::absoluteImpactParameter3D(piTT, jpsiVtx);
    if (piIP_pair.first)
    {
      measure.first  = piIP_pair.second.value();
      measure.second = piIP_pair.second.significance();
    }
    else 
    {
      measure.first  = 0;
      measure.second = 0;
    }  
    return measure;
  }


  std::pair<double,double> pionIPBeamSpot(reco::TransientTrack piTT, GlobalPoint BsGp)
  {
    std::pair<double,double> measureBS;
    TrajectoryStateClosestToPoint pion_BeamSpot = piTT.trajectoryStateClosestToPoint(BsGp);
    if(pion_BeamSpot.isValid())
    {
      measureBS.first = pion_BeamSpot.perigeeParameters().transverseImpactParameter();
      if(pion_BeamSpot.hasError() && !(pion_BeamSpot.hasError()==0)) 
      {
        measureBS.second = measureBS.first/pion_BeamSpot.perigeeError().transverseImpactParameterError();
      }
    }
    return measureBS;       
  }
  
  int overlap(const reco::Candidate &a, const reco::Track b) 
  {
    double eps(1.e-3);
    double dphi = deltaPhi(a.phi(), b.phi());
    dphi *= dphi;
    double deta = a.eta() - b.eta();
    deta *= deta;
    if ((dphi + deta) < eps) {
      return 1;
    }
    return 0;
  }

  FreeTrajectoryState initialFreeState( const reco::Track& tk, const MagneticField* field)
  {
    Basic3DVector<float> pos( tk.vertex());
    GlobalPoint gpos( pos);
    Basic3DVector<float> mom( tk.momentum());
    GlobalVector gmom( mom);
    GlobalTrajectoryParameters par( gpos, gmom, tk.charge(), field);
    CurvilinearTrajectoryError err( tk.covariance());
    return FreeTrajectoryState( par, err);
  }



  
  
  
  
};

