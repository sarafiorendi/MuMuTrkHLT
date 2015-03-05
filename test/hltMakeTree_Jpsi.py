import FWCore.ParameterSet.Config as cms

process = cms.Process("RATE")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#                                 inFile
                            "file:/afs/cern.ch/work/f/fiorendi/private/BPH-HLT/CMSSW_7_2_0_pre6/src/HLTrigger/Configuration/test/redoTheLMNRPathOnFile6_noCuts_alsoReco.root"
#                            "/store/group/phys_bphys/fiorendi/720pre6/FromTemplateMenu_onMC_etaPhi117/HLT_DisplacedPaths__phys_bphys_0.root",
                            ),
                            secondaryFileNames = cms.untracked.vstring()  
                            )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "PRE_R_71_V3::All"

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryIdeal_cff")#sara
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

NAMEPLOT = [ "Full", "L1", "L2overL1", "L3overL2", "TrkOverL3"]


TRIGGERS = [ "HLT_DoubleMu3p5_LowMass_Displaced_v6"
           ]


j = 4
for i in range(len(TRIGGERS)):
  PROCMODU = TRIGGERS[i]+NAMEPLOT[j]
  PROCPATH = 'Path'+PROCMODU
  module =cms.EDAnalyzer("JpsiTreeMaker",
                         vertexes           = cms.InputTag("offlinePrimaryVertices"),
                         muons              = cms.InputTag("muons"),
                         tracks             = cms.InputTag("generalTracks"),
                         triggerProcess     = cms.string("HLT"),
                         tagTriggerProcess  = cms.string("TRKHLT"),
                         tagTriggerName     = cms.string("HLT_DoubleMu3p5_LowMassNonResonant_Displaced_v1"),    
                         triggerName        = cms.string("HLT_DoubleMu3p5_LowMassNonResonant_Displaced_v6"),
                         probeFilterDen     = cms.string("hltDisplacedmumuFilterDoubleMu3p5LowMassNonResonant"),#not used
                         probeFilterNum     = cms.string("hltDisplacedmumuFilterDoubleMu3p5LowMassNonResonant"),
                         maxNumberMuons     = cms.untracked.uint32(999999),
                         L3CandidatesLabel  = cms.InputTag("hltL3MuonCandidates"),
                         TrkCandidatesLabel = cms.InputTag("hltIter2L3MuonMerged"),
                         VtxLabel           = cms.InputTag("hltOnlinePrimaryVertices"),
                         MuMuVtxLabel       = cms.InputTag("hltDisplacedmumuVtxProducerDoubleMu3p5LowMassNonResonant"),
                         MuMuVtxProdLabel   = cms.InputTag("hltLowMassNonResonantTkVertexProducer"),
                         maxEta             = cms.untracked.double(2.2),
                         minInvMass         = cms.untracked.double(4.),
                         maxInvMass         = cms.untracked.double(7.),
                         minCandPt          = cms.untracked.double(0),
                         minCandLxy         =  cms.untracked.double(0),
                         minCandCos         = cms.untracked.double(-1),
                         maxNormChi2        = cms.untracked.double(10.),
                         DisplacedJpsiRequirements = cms.untracked.bool(False),# to be set true!
                         Debug              = cms.untracked.bool(False),
                         beamspot           = cms.untracked.InputTag("offlineBeamSpot")
                         )
  setattr(process, PROCMODU, module)
  setattr(process, PROCPATH, cms.Path(module))
  process.mypath  = cms.Path(getattr(process,PROCMODU))

process.TFileService = cms.Service("TFileService",
#                                    fileName = cms.string(outFile),
                                   fileName = cms.string("myTestTreeJpsi_eventsFailingTrigger.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.options   = cms.untracked.PSet(
                       SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
