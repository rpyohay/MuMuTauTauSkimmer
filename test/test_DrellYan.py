import FWCore.ParameterSet.Config as cms
from subprocess import *
import FWCore.Utilities.FileUtils as FileUtils
mylist=FileUtils.loadListFromFile('/afs/cern.ch/user/m/mshi/CMSSW_7_6_3/src/GGHAA2Mu2TauAnalysis/testDrellYan.txt')
process = cms.Process("testSKIM")

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
W_PDGID = 24
TAU_PDGID = 15
MU_PDGID = 13
NUMU_PDGID = 14
D_PDGID = 1
U_PDGID = 2
S_PDGID = 3
C_PDGID = 4
B_PDGID = 5
T_PDGID = 6
G_PDGID = 21
ANY_PDGID = 0

#tau decay types
TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3
ANY_PT_RANK = -1
#tau hadronic decay types
TAU_ALL_HAD = -1
TAU_1PRONG_0NEUTRAL = 0
TAU_1PRONG_1NEUTRAL = 1
TAU_1PRONG_2NEUTRAL = 2
TAU_1PRONG_3NEUTRAL = 3
TAU_1PRONG_NNEUTRAL = 4
TAU_2PRONG_0NEUTRAL = 5
TAU_2PRONG_1NEUTRAL = 6
TAU_2PRONG_2NEUTRAL = 7
TAU_2PRONG_3NEUTRAL = 8
TAU_2PRONG_NNEUTRAL = 9
TAU_3PRONG_0NEUTRAL = 10
TAU_3PRONG_1NEUTRAL = 11
TAU_3PRONG_2NEUTRAL = 12
TAU_3PRONG_3NEUTRAL = 13
TAU_3PRONG_NNEUTRAL = 14
TAU_RARE = 15


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
                SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*mylist))

process.source.inputCommands = cms.untracked.vstring("keep *")

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V7F::All')

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
#for mu-less jets
process.load('Configuration.StandardSequences.MagneticField_cff') #I changed it from: process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') # Kyle Added this
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi') # Kyle Added this
process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15') # Kyle added this
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoTauTag.RecoTau.RecoTauPiZeroProducer_cfi")
process.load('Tools/CleanJets/cleanjets_cfi')
ZMuMuPSet = cms.PSet(momPDGID = cms.vint32(Z_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))
process.genZMuSelector=cms.EDFilter(
 'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(MU_PDGID),     #choose a gen muon...
    sisterAbsMatchPDGID = cms.uint32(MU_PDGID), #...whose sister is another gen muon...
    genTauDecayIDPSet = ZMuMuPSet,              #...and whose mother is a pseudoscalar a
    primaryTauDecayType = cms.uint32(TAU_ALL),  #choose TAU_ALL when the gen particle is not a tau
    sisterTauDecayType = cms.uint32(TAU_ALL),   #choose TAU_ALL when the gen particle sister is not a tau
    primaryTauPTRank = cms.int32(ANY_PT_RANK),  #should always be ANY_PT_RANK
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD), #choose TAU_ALL_HAD when the gen particle is not a tau
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),     #choose TAU_ALL_HAD when the gen particle sister is not a tau
    primaryTauAbsEtaMax = cms.double(-1.0), #no cut on gen particle |eta|
    primaryTauPTMin = cms.double(-1.0),     #no cut on gen particle pT
    countSister = cms.bool(True),           #True if you want to put both muons in the output collection, False if just one
    applyPTCuts = cms.bool(False),          #should always be False
    countKShort = cms.bool(False),          #should always be False
    minNumGenObjectsToPassFilter = cms.uint32(2), #EDFilter only returns true if >=2 muons are found
    makeAllCollections = cms.bool(False) #should always be False
    )
process.NearestMuonFromTau=cms.EDFilter(
'NearestRecoObject',
 tauTag=cms.InputTag('muHadTauSelector','','SKIM'),
 muonTag=cms.InputTag('tauMuonSelector','','SKIM'),
 minNumObjsToPassFilter=cms.uint32(1)
)
process.NearestMuonRef=cms.EDFilter('MuonRefSelector',
                                 src = cms.InputTag('NearestMuonFromTau'),
                                 cut = cms.string('pt > 0.0'),
                                 filter = cms.bool(True)
)
process.genMatchZMuSelector = cms.EDFilter(
    'GenMatchedMuonProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genZMuSelector'), #must be a reco::GenParticleRefVector
    recoObjTag = cms.InputTag('NearestMuonRef'),              #must be a reco::MuonRefVector
    baseRecoObjTag = cms.InputTag('muons'),
    genTauDecayIDPSet = ZMuMuPSet,      #need to know the pseudoscalar a mother
    applyPTCuts = cms.bool(False),        #should always be false
    countKShort = cms.bool(False),        #should always be false
    dR = cms.double(0.1),                 #dR criteria for matching
   Bins=cms.vdouble(0.0, 5.0,10.0, 15.0, 20.0, 40.0, 60.0, 80.0, 100.0),
    minNumGenObjectsToPassFilter = cms.uint32(0) #EDFilter returns true if >=1 gen-matched reco muon is found
    )
#require event to fire IsoMu24_eta2p1
process.RECOAnalyze=cms.EDAnalyzer(
'MuMuTauTauRecoAnalyzer',
  tauTag=cms.InputTag('muHadTauSelector','','SKIM'),
 jetMuonMapTag=cms.InputTag('CleanJets','muonValMap','SKIM'),
  genParticleTag=cms.InputTag('genParticles'),
  muHadMassBins=cms.vdouble(0.0, 2.0, 4.0, 6.0, 8.0, 10.0,12.0,14.0),
  outFileName=cms.string('testDrellYan_RECOAnalyzer.root')
)

process.noSelectionSequence = cms.Sequence(
process.genZMuSelector*
process.NearestMuonFromTau*
process.NearestMuonRef*
process.genMatchZMuSelector*
process.RECOAnalyze
)

process.TFileService = cms.Service("TFileService",
    fileName =  cms.string('DrellYan_Tfile.root')
)
#no selection path
process.p = cms.Path(process.noSelectionSequence)
