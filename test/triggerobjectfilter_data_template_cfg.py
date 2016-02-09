import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
mylist=FileUtils.loadListFromFile('/afs/cern.ch/user/m/mshi/CMSSW_7_4_12_patch4/src/GGHAA2Mu2TauAnalysis/heavy750light9Reco.txt')
process = cms.Process("MUHADANALYSIS")

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
TAU_PDGID = 15
MU_PDGID = 13
ANY_PDGID = 0
#tau decay types
TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3

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

#no consideration of pT rank
ANY_PT_RANK = -1

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
#process.load('trigger_match.TriggerObjectFilter.MuonTriggerObjectFilter_cfi')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
                                     SkipEvent = cms.untracked.vstring('ProductNotFound'))


process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(*mylist)
)

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9::All')

#for HLT selection
#process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#process.Mu16Mu0OniaSelector = process.hltHighLevel.clone()
#process.Mu16Mu0OniaSelector.HLTPaths = cms.vstring('HLT_Mu16_TkMu0_dEta18_Onia_v1')

process.MuonIWant = cms.EDFilter('MuonRefSelector',
                                 src = cms.InputTag('muons'),
                                 cut = cms.string('pt > 0.0'),
                                 filter = cms.bool(False)
)
process.SingleMuon = cms.EDFilter('MuonRefSelector',
                                 src = cms.InputTag('muons'),
                                 cut = cms.string('pt > 45.0 & abs( eta ) < 2.1'),
                                 filter = cms.bool(True)
)
process.muonTriggerObjectFilter1 = cms.EDFilter(
    'MuonTriggerObjectFilter',
    recoObjTag = cms.InputTag('SingleMuon'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu45_eta2p1_v1", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_Mu45_eta2p1_v1"),
    theRightHLTSubFilter1 = cms.InputTag("hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter1= cms.uint32(1),
)
process.SingleMuonLooseID=cms.EDFilter('LooseMuon',
                                       muonTag=cms.InputTag('muonTriggerObjectFilter1'),
                                       minNumObjsToPassFilter=cms.uint32(1)
)
process.afterVetoSingleMuon = cms.EDFilter('VetoMuon',
                              muonTag=cms.InputTag('MuonIWant'),
                              vetoMuonTag=cms.InputTag('muonTriggerObjectFilter1'),
                              minNumObjsToPassFilter=cms.uint32(1)
)
process.SingleMuonsPartnerSelector=cms.EDFilter('MuonPartner',
                                                muonTag=cms.InputTag('afterVetoSingleMuon'),
                                                minNumObjsToPassFilter=cms.uint32(1)
)
process.SingleMuonPartnerLooseID=cms.EDFilter('LooseMuon',
                                              muonTag=cms.InputTag('SingleMuonsPartnerSelector'),
                                              minNumObjsToPassFilter=cms.uint32(1)
)
process.OppositeSign=cms.EDFilter('OppositeSign',
                                   muonTag=cms.InputTag('SingleMuonPartnerLooseID'),
                                   SingleMuonTag=cms.InputTag('SingleMuonLooseID'),
                                   minNumObjsToPassFilter=cms.uint32(1)
)

process.SeventeenGeVMuSelector = cms.EDFilter('MuonRefSelector',
                                        src =cms.InputTag('muons'),
                                        cut =cms.string('pt>17.0'),
                                        filter=cms.bool(True)
)
process.EightGeVMuSelector = cms.EDFilter('MuonRefSelector',
                                       src=cms.InputTag('muons'),
                                       cut=cms.string('pt>8.0'),
                                       filter=cms.bool(True)
)
process.highestSeventeen=cms.EDFilter('HighestPtSelector',
                                muonTag=cms.InputTag('SeventeenGeVMuSelector')
)
process.afterVetoHighestSeventeen=cms.EDFilter('VetoMuon',
                                  muonTag=cms.InputTag('EightGeVMuSelector'),
                                  vetoMuonTag=cms.InputTag('highestSeventeen'),
                                  minNumObjsToPassFilter= cms.uint32(1)
) 
process.highestEight=cms.EDFilter('HighestPtSelector',
                               muonTag=cms.InputTag('afterVetoHighestSeventeen')
)

commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0)
                                   )

process.genTauMuSelector0 = cms.EDFilter('GenObjectProducer',
   genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    primaryTauDecayType = cms.uint32(TAU_MU),
    sisterTauDecayType = cms.uint32(TAU_HAD),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )
process.genAMuSelector1 = cms.EDFilter('GenObjectProducer1',
                                       genParticleTag = cms.InputTag('genParticles'),
                                       minNumGenObjectsToPassFilter = cms.uint32(1)
)
process.genMatchSeventeen0=cms.EDFilter('GenMatchedMuonProducer',
                                genParticleTag = cms.InputTag('genParticles'),
                                selectedGenParticleTag=cms.InputTag('genTauMuSelector0'),
                                recoObjTag=cms.InputTag('highestSeventeen'),
                                baseRecoObjTag=cms.InputTag('muons'),
                                genTauDecayIDPSet=commonGenTauDecayIDPSet,
                                applyPTCuts=cms.bool(False),
                                countKShort=cms.bool(True),
                                pTRank=cms.int32(0),
                                makeAllCollections=cms.bool(False),
                                useGenObjPTRank=cms.bool(False),
                                nOutputColls=cms.uint32(1),
                                dR=cms.double(0.1),
                                minNumGenObjectsToPassFilter=cms.uint32(1)
)
process.genMatchEight0=cms.EDFilter('GenMatchedMuonProducer',
                                genParticleTag = cms.InputTag('genParticles'),
                                selectedGenParticleTag=cms.InputTag('genTauMuSelector0'),
                                recoObjTag=cms.InputTag('highestEight'),
                                baseRecoObjTag=cms.InputTag('muons'),
                                genTauDecayIDPSet=commonGenTauDecayIDPSet,
                                applyPTCuts=cms.bool(False),
                                countKShort=cms.bool(True),
                                pTRank=cms.int32(0),
                                makeAllCollections=cms.bool(False),
                                useGenObjPTRank=cms.bool(False),
                                nOutputColls=cms.uint32(1),
                                dR=cms.double(0.1),
                                minNumGenObjectsToPassFilter=cms.uint32(1)
)

process.genMatchSeventeen1=cms.EDFilter('GenMatchedMuonProducer1',
                               selectedGenParticleTag=cms.InputTag('genAMuSelector1'),
                               recoObjTag=cms.InputTag('highestSeventeen'),
                               baseRecoObjTag=cms.InputTag('muons'),
                               dR=cms.double(0.1),
                               minNumGenObjectsToPassFilter=cms.uint32(1)
)
process.genMatchEight1=cms.EDFilter('GenMatchedMuonProducer1',
                               selectedGenParticleTag=cms.InputTag('genAMuSelector1'),
                               recoObjTag=cms.InputTag('highestEight'),
                               baseRecoObjTag=cms.InputTag('muons'),
                               dR=cms.double(0.1),
                               minNumGenObjectsToPassFilter=cms.uint32(1)
)

                           
#muon trigger object filter
process.muonTriggerObjectFilter = cms.EDFilter(
    'MuonTriggerObjectFilter',
    recoObjTag = cms.InputTag('MuonIWant'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_TripleMu_12_10_5_v1", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_TripleMu_12_10_5_v1"),
    theRightHLTSubFilter1 = cms.InputTag("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter1= cms.uint32(1),
    )
process.muonTriggerObjectFilter1 = cms.EDFilter(
    'MuonTriggerObjectFilter',
    recoObjTag = cms.InputTag('SingleMuon'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu45_eta2p1_v1", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_Mu45_eta2p1_v1"),
    theRightHLTSubFilter1 = cms.InputTag("hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter1= cms.uint32(1),
    )
process.muonTriggerObjectFilter2 = cms.EDFilter(
    'MuonTriggerObjectFilter',
    recoObjTag = cms.InputTag('MuonIWant'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu17_Mu8_SameSign_v1", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_Mu17_Mu8_SameSign_v1"),
    theRightHLTSubFilter1 = cms.InputTag("hltL3fL1sMu5L1f0L2f5L3Filtered8"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter1= cms.uint32(1),
    )
process.muonTriggerObjectFilter3 = cms.EDFilter(
    'MuonTriggerObjectFilter',
    recoObjTag = cms.InputTag('MuonIWant'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu17_Mu8_SameSign_v1", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_Mu17_Mu8_SameSign_v1"),
    theRightHLTSubFilter1 = cms.InputTag("hltL3fL1sMu12L1f0L2f12L3Filtered17"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter1= cms.uint32(1),
    )

process.muonTriggerObjectFilter8 = cms.EDFilter(
    'MuonTriggerObjectFilter8',
    recoObjTag = cms.InputTag('MuonIWant'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu17_Mu8_SameSign_v1", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_Mu17_Mu8_SameSign_v1"),
    theRightHLTSubFilter0 = cms.InputTag("hltL3fL1sMu5L1f0L2f5L3Filtered8"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter0= cms.uint32(1),
    )

process.muonTriggerObjectFilter17 = cms.EDFilter(
    'MuonTriggerObjectFilter17',
    recoObjTag = cms.InputTag('MuonIWant'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu17_Mu8_SameSign_v1", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_Mu17_Mu8_SameSign_v1"),
    theRightHLTSubFilter1 = cms.InputTag("hltL3fL1sMu12L1f0L2f12L3Filtered17"),
    HLTSubFilters = cms.untracked.VInputTag(""),
    minNumObjsToPassFilter1= cms.uint32(1)
    )
process.muonPtCorrelation=cms.EDFilter(
    'MuonPtCorrelation',
    recoObjTag17=cms.InputTag('muonTriggerObjectFilter17'),
    recoObjTag8=cms.InputTag('muonTriggerObjectFilter8')
)
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
process.filter_all_explicit = hlt.hltHighLevel.clone(
    HLTPaths = ['HLT_Mu45_eta2p1_v1'],
    andOr = False,
    throw = False
)
#output
process.TFileService = cms.Service("TFileService",
    fileName =  cms.string('testtrigger.root')
)
process.SingleMuonSelection=cms.Sequence(
process.SingleMuon*
process.MuonIWant*
process.muonTriggerObjectFilter1*
process.SingleMuonLooseID*
process.afterVetoSingleMuon*
process.SingleMuonsPartnerSelector*
process.SingleMuonPartnerLooseID*
process.OppositeSign
)

process.SeventeenTauMu = cms.Sequence(
    process.SeventeenGeVMuSelector*
    process.highestSeventeen*
    process.genTauMuSelector0*
    process.genMatchSeventeen1 
)
process.SeventeenAMu = cms.Sequence(
    process.SeventeenGeVMuSelector*
    process.highestSeventeen*
    process.genAMuSelector1*
    process.genMatchSeventeen1*
    process.muonTriggerObjectFilter17 
)
process.EightTauMu = cms.Sequence(
    process.SeventeenGeVMuSelector*
    process.highestSeventeen*
    process.EightGeVMuSelector*
    process.afterVetoHighestSeventeen*
    process.highestEight*
    process.genTauMuSelector0*
    process.genMatchEight0
)

process.p = cms.Path(process.SingleMuonSelection)
                           
