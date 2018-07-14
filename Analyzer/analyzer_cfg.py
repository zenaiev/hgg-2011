import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
import sys

########################################################################
#################### Passed arguments ##################################
########################################################################
#
# can be invoked with no parameters passed, in this case use default values
#
# input file name
#inFileTest = 'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/Photon/AOD/12Oct2013-v1/20000/04F7597C-9D38-E311-A47C-0025904B2C4E.root'
#inFileTest = 'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/VBFHiggs0PToGG_M-125p6_7TeV-JHUGenV4-pythia6-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/385469CA-1C01-E411-92A7-002590AC4C24.root'
# 2012
inFileTest = 'root://eospublic.cern.ch//eos/opendata/cms/Run2012B/DoublePhoton/AOD/22Jan2013-v1/20000/02ACA618-0F6D-E211-B9A1-003048FFD744.root'
# for fast tests, you can copy input ROOT files to the local machine
#inFileTest = 'file:/home/cms-opendata/cmsOpenDataFiles/00A4E1AF-B3C7-E311-BA6D-002590200808.root'
#
# output file name
outFileTest = 'hggTmp.root'
#
# flags which determine what will be done
flag_reco = 1   # process reconstruction level
flag_gen  = 0   # process generated level
flag_mc   = 0   # 1 for mc, 0 for data
flag_year = 1   # 0 for 2011, 1 for 2012
#
# number of events to be processed (-1 for all)
maxEvents = -1
#
# process passed arguments, if any
#
if len(sys.argv) < 5:
  print("Usage: cmsRun analyzer_cfg.py <input list> <output file> <reco flag> <gen flag> <mc flag> <year>")
  inputList = inFileTest
  outFile = outFileTest
  maxEvents = 1000
  # do not stop execution at this point, run with default arguments
  #sys.exit("Wrong usage!")
else:                 
  inputList = FileUtils.loadListFromFile(sys.argv[2])
  outFile   = sys.argv[3]
  flag_reco = int(sys.argv[4])
  flag_gen  = int(sys.argv[5])
  flag_mc   = int(sys.argv[6])
  flag_year = int(sys.argv[7])
# consistency check
if flag_gen == 1 and flag_mc == 0: 
  sys.exit("Error: gen = 1 requires mc = 1")
#
########################################################################
#
########################################################################
#################### Prepare and run Analyzer ##########################
########################################################################
#
process = cms.Process("Demo")
# load needed tools
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#
# global tag as described at http://opendata.cern.ch/getting-started/CMS?ln=en
if flag_mc == 0:
# DATA
  if flag_year == 0:
    # Before should be done:
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
    process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
  elif flag_year == 1:
    # Before should be done:
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL FT53_V21A_AN6
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL.db FT53_V21A_AN6_FULL.db
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL FT53_V21A_AN6_FULL
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL.db')
    process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'
# MC
else:
  if flag_year == 0:
    # Before should be done:
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db START53_LV6A1.db
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db')
    process.GlobalTag.globaltag = 'START53_LV6A1::All'
  elif flag_year == 1:
    # Before should be done:
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27.db START53_V27.db
    # ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27 START53_V27
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_V27.db')
    process.GlobalTag.globaltag = 'START53_V27::All'
#
# intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('Demo')
# change the value below if you want more or less status output
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#
# Change this to a positive value to limit the number of processed events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
#
# supply processor with input files
if len(sys.argv) < 4:
  # only one file in the list
  process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(inputList))
else:
  # many files in the list
  process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(*inputList))
#
# JSON (good luminosity sections), only if processing data
if flag_mc == 0:
  if flag_year == 0:
    goodJSON = 'data/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
  elif flag_year == 1:
    goodJSON = 'data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'
  myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',') 
  process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
  process.source.lumisToProcess.extend(myLumis) 
#
# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
#
# all is ready: pass all arguments to Analyzer (C++ code in src/Analyzer.cc)
process.demo = cms.EDAnalyzer('Analyzer', outFile = cms.string(outFile), mc = CfgTypes.int32(flag_mc), reco = CfgTypes.int32(flag_reco), gen = CfgTypes.int32(flag_gen))
process.p = cms.Path(process.demo)
#
########################################################################
#
########################################################################
#################### ParticleTreeDrawer ################################
########################################################################
#
# this is another processor which prints MC generated particles
# (if you wanr this, disable Analyzer above)
#
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
#
########################################################################
