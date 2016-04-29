import FWCore.ParameterSet.Config as cms
from subprocess import *

process = cms.Process("MERGE")

#verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

#input
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('NMSSM_ggH_a9_H1125_H2500_H3500_2mu2tau.txt')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(*mylist)
    )

#output
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('heavyHiggs_125_light_9_2mu2tau_STEP2.root')
    )

#path
process.e = cms.EndPath(process.output)
