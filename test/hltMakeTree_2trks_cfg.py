import FWCore.ParameterSet.Config as cms

process = cms.Process("RATE")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                inFile
#                          	"/store/group/phys_bphys/fiorendi/720pre6/FromTemplateMenu_onMC_etaPhi117_addBMass/HLT_DisplacedPaths__phys_bphys_10.root",
#                               'file:/afs/cern.ch/work/f/fiorendi/private/BPH-HLT/CMSSW_7_2_0_pre6/src/HLTrigger/Configuration/test/cfgs_2014_10_16_12_7/HLT_DisplacedPaths__phys_bphys_14.root',
                            ),
                            secondaryFileNames = cms.untracked.vstring()
                            )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "POSTLS172_V3::All" 

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryIdeal_cff")#sara
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

NAMEPLOT = [ "Full", "L1", "L2overL1", "L3overL2", "TrkOverL3"]

FILTERS  = [ "hltJpsiTkTkVertexFilterOS", 
           ]

TRIGGERS = [ "HLT_DoubleMu4_JpsiTrkTrk_Displaced_OS"]
# TRIGGERS = [ "HLT_DoubleMu4_JpsiTrkTrk_Displaced_v1"           ]


j = 4
for i in range(len(TRIGGERS)):
  PROCMODU = TRIGGERS[i]+NAMEPLOT[j]
  PROCPATH = 'Path'+PROCMODU
  module =cms.EDAnalyzer("JpsiTrkTrkTreeMaker",
						 vertexes           = cms.InputTag("offlinePrimaryVertices"),
						 muons              = cms.InputTag("muons"),
						 tracks             = cms.InputTag("generalTracks"),
						 triggerProcess     = cms.string("TRKHLT"), 
						 tagTriggerName     = cms.string("HLT_DoubleMu4_Jpsi_Displaced_v1"),    
                         triggerName        = cms.string(TRIGGERS[i]),
                         probeFilterDen     = cms.string("hltJpsiTkTkVertexFilter"),#not used
                         probeFilterNum     = cms.string(FILTERS[i]),
						 maxNumberMuons     = cms.untracked.uint32(999999),
						 L3CandidatesLabel  = cms.InputTag("hltL3MuonCandidates"), 
						 TrkCandidatesLabel = cms.InputTag("hltIter2L3MuonMerged"), 
						 VtxLabel           = cms.InputTag("hltOnlinePrimaryVertices"), 
						 MuMuVtxLabel       = cms.InputTag("hltDisplacedmumuVtxProducerDoubleMu4Jpsi"), 
                         MuMuVtxProdLabel   = cms.InputTag("hltJpsiTkTkVertexProducer"),
						 maxEta             = cms.untracked.double(2.2),
						 minPtTrk           = cms.untracked.double(0.5),
						 thirdTrkMass       = cms.untracked.double(0.493677),# kaon
						 fourthTrkMass      = cms.untracked.double(0.139570),# pion
						 minInvMass         = cms.untracked.double(4.),
						 maxInvMass         = cms.untracked.double(7),
						 minTrkTrkMass      = cms.untracked.double(0.692),
						 maxTrkTrkMass      = cms.untracked.double(1.092),
						 minCandPt          = cms.untracked.double(0),
						 minCandLxy         = cms.untracked.double(0),
						 minCandCos         = cms.untracked.double(-1),
						 maxNormChi2        = cms.untracked.double(10.),
						 minVtxProb         = cms.untracked.double(0.005),
						 mind0Sign          = cms.untracked.double(0.),
						 DisplacedJpsiRequirements = cms.untracked.bool(True), 
						 Debug              = cms.untracked.bool(False), 
						 beamspot           = cms.untracked.InputTag("offlineBeamSpot") 
						 )   
  setattr(process, PROCMODU, module)
  setattr(process, PROCPATH, cms.Path(module))
  process.mypath  = cms.Path(getattr(process,PROCMODU))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outFile),
#                                    fileName = cms.string("Efficiency_2trks_test.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))   

process.options   = cms.untracked.PSet( 
                       SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
