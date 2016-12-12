import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
import sys

# ******* reconstructed or/and generated *********
reco = 1
gen  = 0
mc   = 1
if gen == 1 and mc == 0: 
  sys.exit("Error: gen = 1 requires mc = 1")
# ************************************************

inFile1MuEG = 'file:root://eospublic.cern.ch//eos/opendata/cms/Run2011A/MuEG/AOD/12Oct2013-v1/20001/00440023-843E-E311-A760-02163E008D8E.root'
inFile1DoubleMu = 'file:root://eospublic.cern.ch//eos/opendata/cms/Run2011A/DoubleMu/AOD/12Oct2013-v1/10000/000D143E-9535-E311-B88B-002618943934.root'
inFile1DoubleEl = 'file:root://eospublic.cern.ch//eos/opendata/cms/Run2011A/DoubleElectron/AOD/12Oct2013-v1/20000/0014CE62-9C3E-E311-8FCC-00261894389F.root'

#inFileTest = 'file:/home/cms-opendata/CMSSW_5_3_32/src/ttbar/Analyzer/small.root'
#inFileTest = 'file:root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00003/B47F2161-0DB6-E311-B0FB-00304867902E.root'
#inFileTest = 'file:/home/cms-opendata/cmsOpenDataFiles/00440023-843E-E311-A760-02163E008D8E.root'
#inFileTest = 'file:/home/cms-opendata/cmsOpenDataFiles/825B6645-B8CE-E311-A85A-0025B3E05CBC.root'
inFileTest = 'file:/home/cms-opendata/cmsOpenDataFiles/00345E78-7CC5-E311-A9C7-001E6739730A.root'

outFileTest = 'ttbarSelTmp.root'

# arguments
if len(sys.argv) < 4:
  print("Usage: cmsRun analyzer_cfg.py <input list> <output file>")
  inputList = inFileTest;outFile = outFileTest
  #sys.exit("Wrong usage!")
else:                 
  inputList = FileUtils.loadListFromFile(sys.argv[2])
  outFile = sys.argv[3]

process = cms.Process("Demo")
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# global tag from http://opendata.cern.ch/getting-started/CMS?ln=en
if mc == 0:
# DATA
# Before should be done:
#    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1
#    ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
  process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
  process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
# MC
# Before should be done:
#ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db START53_LV6A1.db
#ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
else:
  process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db')
  process.GlobalTag.globaltag = 'START53_LV6A1::All'
   
# intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

if len(sys.argv) < 4:
  process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inputList))
else:
  process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(*inputList))

# JSON
goodJSON = '/home/cms-opendata/CMSSW_5_3_32/src/ttbar/Analyzer/data/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',') 
if mc == 0:
  process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
  process.source.lumisToProcess.extend(myLumis) 

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")

process.demo = cms.EDAnalyzer('Analyzer', outFile = cms.string(outFile), mc = cms.int32(mc), reco = cms.int32(reco), gen = cms.int32(gen))
process.p = cms.Path(process.demo)

# print MC generated particles
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                   src = cms.InputTag("genParticles"),                                                                 
#                                   printP4 = cms.untracked.bool(False),
#                                   printPtEtaPhi = cms.untracked.bool(True),
#                                   printVertex = cms.untracked.bool(False),
#                                   printStatus = cms.untracked.bool(True),
#                                   printIndex = cms.untracked.bool(False),
#                                   #status = cms.untracked.vint32( 3 )
#                                   )
#process.p = cms.Path( process.printTree)
