// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc ttbar/Analyzer/src/Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Sun Jul  3 15:22:14 CEST 2016
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//------ EXTRA HEADER FILES--------------------//
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

// electrons
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// jets
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/BTauReco/interface/JetTag.h>
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// MET
#include "DataFormats/METReco/interface/PFMET.h"

//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// MC
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//
// class declaration
//

#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>

#include <TTree.h>
#include <TDirectory.h>


class Analyzer : public edm::EDAnalyzer {
   public:
      explicit Analyzer(const edm::ParameterSet&);
      ~Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
      // OZ
      int SelectEvent(const edm::Event& iEvent);
      int SelectMu(const edm::Handle<reco::MuonCollection>& muons, const reco::VertexCollection::const_iterator& pv);
      int SelectEl(const edm::Handle<reco::GsfElectronCollection>& electrons, const reco::VertexCollection::const_iterator& pv);
      int SelectJet(const edm::Handle<reco::PFJetCollection>& jets, const reco::JetTagCollection& bTags, const edm::Event& iEvent, const edm::EventSetup& iSetup);
      int SelectMET(const edm::Handle<edm::View<reco::PFMET> >& pfmets);
      reco::PFJetCollection::const_iterator MatchJet(const edm::Handle<reco::PFJetCollection>& jets, const edm::RefToBase<reco::Jet>& jettag, double& diff, double& nextdiff);
      void FindTriggerBits(const HLTConfigProvider& trigConf);
      void SelectTriggerBits(const edm::Handle<edm::TriggerResults>& HLTR);
      void PrintTriggerBits();
      int SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex);
      //std::vector<const reco::Candidate*> FindGenDecay(const reco::Candidate* particle, const std::vector<int>& daugterIDs);
      //std::vector<const reco::Candidate*> FindGenDecayW(const reco::Candidate* particle, const std::vector<int>& leptonID);
      const reco::Candidate* GetFinalState(const reco::Candidate* particle, const int id);
      void FillFourMomentum(const reco::Candidate* particle, float* p);
      void SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles);
      void DisableUnnededBranches();
      void InitBranchVars();

      // ----------member data ---------------------------
      edm::InputTag _inputTagMuons;
      edm::InputTag _inputTagElectrons;
      edm::InputTag _inputTagBtags;
      edm::InputTag _inputTagJets;
      edm::InputTag _inputTagMet;
      edm::InputTag _inputTagTriggerResults;
      edm::InputTag _inputTagPrimaryVertex;
      edm::InputTag _inputTagMCgen;

      // Jet correction label
      std::string mJetCorr;

      // general
      int _flagMC;
      int _flagRECO;
      int _flagGEN;
      int _nevents;
      int _neventsSelected;
      int _signLeptonP;
      int _signLeptonM;
      // storage
      TFile* _file;
      TTree* _tree;
      // event
      int _evRunNumber;
      int _evLumiBlock;
      int _evEventNumber;
      // muons
      static const int _maxNmu = 10;
      int _Nmu;
      float _muPt[_maxNmu];
      float _muEta[_maxNmu];
      float _muPhi[_maxNmu];
      float _muIso03[_maxNmu];
      float _muIso04[_maxNmu];
      int _muHitsValid[_maxNmu];
      int _muHitsPixel[_maxNmu];
      float _muDistPV0[_maxNmu];
      float _muDistPVz[_maxNmu];
      float _muTrackChi2NDOF[_maxNmu];
      // electrons
      static const int _maxNel = 10;
      int _Nel;
      float _elPt[_maxNel];
      float _elEta[_maxNel];
      float _elPhi[_maxNel];
      float _elIso03[_maxNel];
      float _elIso04[_maxNel];
      int _elConvFlag[_maxNel];
      float _elConvDist[_maxNel];
      float _elConvDcot[_maxNel];
      float _elMissHits[_maxNel];
      float _elDistPV0[_maxNel];
      float _elDistPVz[_maxNel];
      // jets
      static const int _maxNjet = 25;
      int _Njet;
      float _jetPt[_maxNjet];
      float _jetEta[_maxNjet];
      float _jetPhi[_maxNjet];
      float _jetMass[_maxNjet];
      float _jetMuEn[_maxNjet];
      float _jetElEn[_maxNjet];
      float _jetBTagDiscr[_maxNjet];
      float _jetBTagMatchDiff1[_maxNjet];
      float _jetBTagMatchDiff2[_maxNjet];
      // MET
      float _metPx;
      float _metPy;
      // triggers
      int _triggers;
      std::vector<std::vector<int> > _vecTriggerBits;
      std::vector<std::string> _vecTriggerNames;
      // primary vertex
      int _Npv;
      int _pvNDOF;
      float _pvZ;
      float _pvRho;
      
      // MC generated info
      int _mcNTTbar;
      int _mcNTWb;
      int _mcNTtbarDilepton;
      int _mcNTtbarDileptonEE;
      int _mcNTtbarDileptonMuMu;
      int _mcNTtbarDileptonTauTau;
      int _mcNTtbarDileptonEMu;
      // MC generated info to store
      int _mcEventType;
      float _mcT[4];
      float _mcTbar[4];
      float _mcWp[4];
      float _mcWm[4];
      float _mcB[4];
      float _mcBbar[4];
      float _mcLp[4];
      float _mcNu[4];
      float _mcLm[4];
      float _mcNubar[4];
};

//
// constants, enums and typedefs
//
double _massMu = 0.105658;
double _massEl = 0.000511;

//
// static data member definitions
//

//
// constructors and destructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig)
{
  setbuf(stdout, NULL);
  //DisableUnnededBranches();
  
  // input tags
  _inputTagMuons = edm::InputTag("muons");
  _inputTagElectrons = edm::InputTag("gsfElectrons");
  _inputTagBtags = edm::InputTag("combinedSecondaryVertexBJetTags");
  _inputTagJets = edm::InputTag("ak5PFJets");
  _inputTagMet = edm::InputTag("pfMet");
  _inputTagTriggerResults = edm::InputTag("TriggerResults", "", "HLT");
  //_inputTagTriggerResults = edm::InputTag("TriggerResults");
  _inputTagPrimaryVertex = edm::InputTag("offlinePrimaryVerticesWithBS");
  _inputTagMCgen = edm::InputTag("genParticles");
  
  // jet correction label
  mJetCorr = "ak5PFL1FastL2L3Residual";

  //now do what ever initialization is needed
  _flagMC = iConfig.getParameter<int>("mc");
  _flagRECO = iConfig.getParameter<int>("reco");
  _flagGEN = iConfig.getParameter<int>("gen");
  _nevents = 0;
  _neventsSelected = 0;
  std::string fileout = iConfig.getParameter<std::string>("outFile");
  _file = new TFile(fileout.c_str(), "recreate");
  _tree = new TTree("tree", "ttbar");
  // event
  _tree->Branch("evRunNumber", &_evRunNumber, "evRunNumber/I");
  _tree->Branch("evLumiBlock", &_evLumiBlock, "evLumiBlock/I");
  _tree->Branch("evEventNumber", &_evEventNumber, "evEventNumber/I");

  if(_flagRECO)
  {
    // muons
    _tree->Branch("Nmu", &_Nmu, "Nmu/I");
    _tree->Branch("muPt", _muPt, "muPt[Nmu]/F");
    _tree->Branch("muEta", _muEta, "muEta[Nmu]/F");
    _tree->Branch("muPhi", _muPhi, "muPhi[Nmu]/F");
    _tree->Branch("muIso03", _muIso03, "muIso03[Nmu]/F");
    _tree->Branch("muIso04", _muIso04, "muIso04[Nmu]/F");
    _tree->Branch("muHitsValid", _muHitsValid, "muHitsValid[Nmu]/I");
    _tree->Branch("muHitsPixel", _muHitsPixel, "muHitsPixel[Nmu]/I");
    _tree->Branch("muDistPV0", _muDistPV0, "muDistPV0[Nmu]/F");
    _tree->Branch("muDistPVz", _muDistPVz, "muDistPVz[Nmu]/F");
    _tree->Branch("muTrackChi2NDOF", _muTrackChi2NDOF, "muTrackChi2NDOF[Nmu]/F");
    // electrons
    _tree->Branch("Nel", &_Nel, "Nel/I");
    _tree->Branch("elPt", _elPt, "elPt[Nel]/F");
    _tree->Branch("elEta", _elEta, "elEta[Nel]/F");
    _tree->Branch("elPhi", _elPhi, "elPhi[Nel]/F");
    _tree->Branch("elIso03", _elIso03, "elIso03[Nel]/F");
    _tree->Branch("elIso04", _elIso04, "elIso04[Nel]/F");
    _tree->Branch("elConvFlag", _elConvFlag, "elConvFlag[Nel]/I");
    _tree->Branch("elConvDist", _elConvDist, "elConvDist[Nel]/F");
    _tree->Branch("elConvDcot", _elConvDcot, "elConvDcot[Nel]/F");
    _tree->Branch("elMissHits", _elMissHits, "elMissHits[Nel]/F");
    _tree->Branch("elDistPV0", _elDistPV0, "elDistPV0[Nel]/F");
    _tree->Branch("elDistPVz", _elDistPVz, "elDistPVz[Nel]/F");
    // jets
    _tree->Branch("Njet", &_Njet, "Njet/I");
    _tree->Branch("jetPt", _jetPt, "jetPt[Njet]/F");
    _tree->Branch("jetEta", _jetEta, "jetEta[Njet]/F");
    _tree->Branch("jetPhi", _jetPhi, "jetPhi[Njet]/F");
    _tree->Branch("jetMass", _jetMass, "jetMass[Njet]/F");
    _tree->Branch("jetMuEn", _jetMuEn, "jetMuEn[Njet]/F");
    _tree->Branch("jetElEn", _jetElEn, "jetElEn[Njet]/F");
    _tree->Branch("jetBTagDiscr", _jetBTagDiscr, "jetBTagDiscr[Njet]/F");
    _tree->Branch("jetBTagMatchDiff1", _jetBTagMatchDiff1, "jetBTagMatchDiff1[Njet]/F");
    _tree->Branch("jetBTagMatchDiff2", _jetBTagMatchDiff2, "jetBTagMatchDiff2[Njet]/F");
    // MET
    _tree->Branch("metPx", &_metPx, "metPx/F");
    _tree->Branch("metPy", &_metPy, "metPy/F");
    // triggers
    _tree->Branch("Triggers", &_triggers, "Triggers/I");
    // mumu
    _vecTriggerNames.push_back("HLT_DoubleMu7");
    _vecTriggerNames.push_back("HLT_Mu13_Mu8");
    _vecTriggerNames.push_back("HLT_Mu17_Mu8");
    _vecTriggerNames.push_back("HLT_DoubleMu6");
    _vecTriggerNames.push_back("HLT_DoubleMu45");
    _vecTriggerNames.push_back("HLT_Mu10_Ele10_CaloIdl");
    // ee
    _vecTriggerNames.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL");
    _vecTriggerNames.push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL");
    _vecTriggerNames.push_back("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL");
    _vecTriggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_");
    _vecTriggerNames.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL");
    _vecTriggerNames.push_back("HLT_DoubleEle45_CaloIdL");
    // emu
    _vecTriggerNames.push_back("HLT_Mu10_Ele10_CaloIdL");
    _vecTriggerNames.push_back("HLT_Mu8_Ele17_CaloIdL");
    _vecTriggerNames.push_back("HLT_Mu17_Ele8_CaloIdL");
    _vecTriggerNames.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL");
    _vecTriggerNames.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL");
    _vecTriggerNames.push_back("HLT_Mu10_Ele10_CaloIdL");
    // primary vertex
    _tree->Branch("Npv", &_Npv, "Npv/I");
    _tree->Branch("pvNDOF", &_pvNDOF, "pvNDOF/I");
    _tree->Branch("pvZ", &_pvZ, "pvZ/F");
    _tree->Branch("pvRho", &_pvRho, "pvRho/F");
  }

  // MC generated info
  if(_flagGEN)
  {
    _mcNTTbar = 0;
    _mcNTWb = 0;
    _mcNTtbarDilepton = 0;
    _mcNTtbarDileptonEE = 0;
    _mcNTtbarDileptonMuMu = 0;
    _mcNTtbarDileptonTauTau = 0;
    _mcNTtbarDileptonEMu = 0;
    // MC generated info to store
    _tree->Branch("mcEventType", &_mcEventType, "mcEventType/I");
    _tree->Branch("mcT", _mcT, "mcT[4]/F");
    _tree->Branch("mcTbar", _mcTbar, "mcTbar[4]/F");
    _tree->Branch("mcWp", _mcWp, "mcWp[4]/F");
    _tree->Branch("mcWm", _mcWm, "mcWm[4]/F");
    _tree->Branch("mcB", _mcB, "mcB[4]/F");
    _tree->Branch("mcBbar", _mcBbar, "mcBbar[4]/F");
    _tree->Branch("mcLp", _mcLp, "mcLp[4]/F");
    _tree->Branch("mcNu", _mcNu, "mcNu[4]/F");
    _tree->Branch("mcLm", _mcLm, "mcLm[4]/F");
    _tree->Branch("mcNubar", _mcNubar, "mcNubar[4]/F");
  }
}


Analyzer::~Analyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  _file->cd();
  _tree->Write();
  _file->Close();

  if(_flagGEN)
  {
    // MC gen info
    printf("%25s = %10d\n", "mcTTbar", _mcNTTbar);
    printf("%25s = %10d(%5.2f%%)\n", "_mcNTWb", _mcNTWb, 100. * _mcNTWb / _mcNTTbar / 2);
    printf("%25s = %10d(%5.2f%%)\n", "_mcNTtbarDilepton", _mcNTtbarDilepton, 100. * _mcNTtbarDilepton / _mcNTTbar);
    printf("%25s = %10d(%5.2f%%)\n", "_mcNTtbarDileptonEE", _mcNTtbarDileptonEE, 100. * _mcNTtbarDileptonEE / _mcNTTbar);
    printf("%25s = %10d(%5.2f%%)\n", "_mcNTtbarDileptonMuMu", _mcNTtbarDileptonMuMu, 100. * _mcNTtbarDileptonMuMu / _mcNTTbar);
    printf("%25s = %10d(%5.2f%%)\n", "_mcNTtbarDileptonTauTau", _mcNTtbarDileptonTauTau, 100. * _mcNTtbarDileptonTauTau / _mcNTTbar);
    printf("%25s = %10d(%5.2f%%)\n", "_mcNTtbarDileptonEMu", _mcNTtbarDileptonEMu, 100. * _mcNTtbarDileptonEMu / _mcNTTbar);
  }

  printf("Processed %d events, selected %d\n", _nevents, _neventsSelected);
}


//
// member functions
//

void Analyzer::InitBranchVars()
{
  _evRunNumber = 0;
  _evLumiBlock = 0;
  _evEventNumber = 0;
  _Nmu = 0;
  _Nel = 0;
  _Njet = 0;
  _metPx = 0;
  _metPy = 0;
  _triggers = 0;
  _Npv= 0;
  _pvNDOF = 0;
  _pvZ = 0;
  _pvRho = 0;
  _mcEventType = 0;
}

// faster tree processing
void Analyzer::DisableUnnededBranches()
{
   std::vector<TString> names;
   names.push_back("Events");
   names.push_back("MetaData");
   names.push_back("Parentage");
   names.push_back("ParameterSets");
   names.push_back("LuminosityBlocks");
   names.push_back("Runs");
   
   //gDirectory->ls();
   for(unsigned int name = 0; name < names.size(); name++)
   {
     TTree *oldtree = (TTree*)gDirectory->Get(names[name]);
     if(oldtree)
     {
       if(name == 0)
       {
         printf("disabling tree %s\n", names[name].Data());
         std::vector<TString> enable;
         enable.push_back("BranchListIndexes");
         enable.push_back("Event*");
         if(_flagRECO)
         {
           enable.push_back("recoMuons_muons__RECO.*");
           enable.push_back("recoGsfElectrons_gsfElectrons__RECO.*");
           enable.push_back("recoPFJets_ak5PFJets__RECO.*");
           enable.push_back("recoJetedmRefToBaseProdTofloatsAssociationVector_combinedSecondaryVertexBJetTags__RECO.*");
           enable.push_back("recoPFMETs_pfMet__RECO.*");
           enable.push_back("recoVertexs_offlinePrimaryVerticesWithBS__RECO.*");
           enable.push_back("edmTriggerResults_TriggerResults__HLT.*");
           enable.push_back("recoTracks_globalMuons__RECO.*");
           enable.push_back("recoGsfElectronCores_gsfElectronCores__RECO.*");
           enable.push_back("recoGsfTracks_electronGsfTracks__RECO.*");
           enable.push_back("recoCaloJets_ak5CaloJets__RECO.*");
         }
         if(_flagGEN)
         {
           enable.push_back("recoGenParticles_genParticles__SIM.*");
         }

         std::vector<TString> disable;
         
         oldtree->SetBranchStatus("*",0);
         for(std::vector<TString>::iterator it = enable.begin(); it != enable.end(); it++)
         {
           oldtree->SetBranchStatus(*it, 1);
           oldtree->AddBranchToCache(*it, 1);
         }
         for(std::vector<TString>::iterator it = disable.begin(); it != disable.end(); it++)
         {
           oldtree->SetBranchStatus(*it, 0);
         }
       }
       else
       {
         oldtree->AddBranchToCache("*", 1);
       }
       oldtree->StopCacheLearningPhase();
     }
     else
       printf("no tree %s\n", names[name].Data());
   }
}

// Event info
int Analyzer::SelectEvent(const edm::Event& iEvent)
{
  _evRunNumber = iEvent.id().run();
  _evLumiBlock = iEvent.id().luminosityBlock();
  _evEventNumber = iEvent.id().event();
  return 0;
}

// MET selection
int Analyzer::SelectMET(const edm::Handle<edm::View<reco::PFMET> >& pfmets)
{
  _metPx = (pfmets->front()).px();
  _metPy = (pfmets->front()).py();
  return 0;
}

// muon selection
int Analyzer::SelectMu(const edm::Handle<reco::MuonCollection>& muons, const reco::VertexCollection::const_iterator& pv)
{
  _Nmu = 0;
  // loop over muons
  for (reco::MuonCollection::const_iterator it = muons->begin(); it != muons->end(); it++)
  {
    if(_Nmu == _maxNmu)
    {
      printf("Maximum number of muons %d reached, skipping the rest\n", _maxNmu);
      return 0;
    }
    // selection
    if(it->pt() < 20)
      continue;
    if(TMath::Abs(it->eta()) > 2.4)
      continue;
    if((it->globalTrack()).isNull())
      continue;
    // fill isolation
    const reco::MuonPFIsolation& iso03 = it->pfIsolationR03();
    _muIso03[_Nmu] = ( iso03.sumChargedParticlePt + iso03.sumNeutralHadronEt ) / it->pt();
    const reco::MuonPFIsolation& iso04 = it->pfIsolationR04();
    _muIso04[_Nmu] = ( iso04.sumChargedParticlePt + iso04.sumNeutralHadronEt ) / it->pt();
    if(_muIso03[_Nmu] > 0.20 && _muIso04[_Nmu] > 0.125)
      continue;
    // fill hits
    _muHitsValid[_Nmu] = 0;
    _muHitsPixel[_Nmu] = 0;
    const reco::HitPattern& p = (it->globalTrack())->hitPattern();
    //const reco::HitPattern& p = it->hitPattern();
    for (int i = 0; i < p.numberOfHits(); i++) 
    {
      uint32_t hit = p.getHitPattern(i);
      if (p.validHitFilter(hit) && p.pixelHitFilter(hit))
        _muHitsPixel[_Nmu]++;
      if (p.validHitFilter(hit))
        _muHitsValid[_Nmu]++;
    }
    // fill momentum
    //TLorentzVector vMu;
    //vMu.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(), _massMu);
    //_muMom[0][_Nmu] = vMu.Px();
    //_muMom[1][_Nmu] = vMu.Py();
    //_muMom[2][_Nmu] = vMu.Pz();
    _muPt[_Nmu] = it->pt() * it->charge();
    _muEta[_Nmu] = it->eta();
    _muPhi[_Nmu] = it->phi();
    // chi2 / ndof
    _muTrackChi2NDOF[_Nmu] = it->globalTrack()->chi2() / it->globalTrack()->ndof();
    // distance to primary vertex
    _muDistPV0[_Nmu] = TMath::Sqrt(TMath::Power(pv->x() - it->globalTrack()->vx(), 2.0) + TMath::Power(pv->y() - it->globalTrack()->vy(), 2.0));
    _muDistPVz[_Nmu] = TMath::Abs(pv->z() - it->globalTrack()->vz());
    // store muon
    _Nmu++;
    if(it->charge() == +1)
        _signLeptonP = 1;
    if(it->charge() == -1)
        _signLeptonM = 1;
  }
  return 0;
}

// electron selection
int Analyzer::SelectEl(const edm::Handle<reco::GsfElectronCollection>& electrons, const reco::VertexCollection::const_iterator& pv)
{
  _Nel = 0;
  for (reco::GsfElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); it++)
  {
    if(_Nel == _maxNel)
    {
      printf("Maximum number of electrons %d reached, skipping the rest\n", _maxNel);
      return 0;
    }
    // selection
    if(it->pt() < 20)
      continue;
    if(TMath::Abs(it->eta()) > 2.4)
      continue;
    // fill isolation
    /*_elIso03TkSumPt[_Nel] = it->dr03TkSumPt();
    _elIso03EcalRecHitSumEt[_Nel] = it->dr03EcalRecHitSumEt();
    _elIso03HcalTowerSumEt[_Nel] = it->dr03HcalTowerSumEt();
    _elIso04TkSumPt[_Nel] = it->dr04TkSumPt();
    _elIso04EcalRecHitSumEt[_Nel] = it->dr04EcalRecHitSumEt();
    _elIso04HcalTowerSumEt[_Nel] = it->dr04HcalTowerSumEt();*/
    _elIso03[_Nel] = (it->dr03TkSumPt() + it->dr03EcalRecHitSumEt() + it->dr03HcalTowerSumEt()) / it->pt();
    _elIso04[_Nel] = (it->dr04TkSumPt() + it->dr04EcalRecHitSumEt() + it->dr04HcalTowerSumEt()) / it->pt();
    if(_elIso03[_Nel] > 0.17 && _elIso04[_Nel] > 0.125)
      continue;
    // fill momentum
    //TLorentzVector vEl;
    //vEl.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(), _massMu);
    //_elMom[0][_Nel] = vEl.Px();
    //_elMom[1][_Nel] = vEl.Py();
    //_elMom[2][_Nel] = vEl.Pz();
    _elPt[_Nel] = it->pt() * it->charge();
    _elEta[_Nel] = it->eta();
    _elPhi[_Nel] = it->phi();
    // fill conversion
    _elConvFlag[_Nel] = it->convFlags();
    _elConvDist[_Nel] = it->convDist();
    _elConvDcot[_Nel] = it->convDcot();
    _elMissHits[_Nel] = ((it->gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
    // distance to primary vertex
    _elDistPV0[_Nel] = TMath::Sqrt(TMath::Power(pv->x() - it->vx(), 2.0) + TMath::Power(pv->y() - it->vy(), 2.0));
    _elDistPVz[_Nel] = TMath::Abs(pv->z() - it->vz());
    // store electron
    _Nel++;
    if(it->charge() == +1)
        _signLeptonP = 1;
    if(it->charge() == -1)
        _signLeptonM = 1;
  }
  return 0;
}

// match jet (for b-tagging)
reco::PFJetCollection::const_iterator Analyzer::MatchJet(const edm::Handle<reco::PFJetCollection>& jets, const edm::RefToBase<reco::Jet>& jettag, double& diff, double& nextdiff)
{
  diff = nextdiff = 1000.0;
  reco::PFJetCollection::const_iterator bestit = jets->begin();
  for(reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); it++)
  {
    double ldiff = TMath::Power(it->eta() - jettag->eta(), 2.0) + TMath::Power(it->phi() - jettag->phi(), 2.0);
    if(ldiff < nextdiff && ldiff > diff)
      nextdiff = ldiff;
    if(ldiff > diff)
      continue;
    nextdiff = diff;
    diff = ldiff;
    bestit = it;
  }
  return bestit;
}

// jet selection
int Analyzer::SelectJet(const edm::Handle<reco::PFJetCollection>& jets, const reco::JetTagCollection& bTags, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  _Njet = 0;
    
  // Loop over jets and study b tag info
  std::vector<std::pair<reco::PFJetCollection::const_iterator, std::vector<double> > > vItBTagMatched;
  for (unsigned int i = 0; i != bTags.size(); ++i) 
  {
    double discriminator = bTags[i].second;
    if(discriminator < 0.2)
      continue;
    //printf("************** New jet: ************** \n");
    std::vector<double> btag(3);
    btag[1] = btag[2] = -1.0;
    reco::PFJetCollection::const_iterator itBTagMatched = MatchJet(jets, bTags[i].first, btag[1], btag[2]);
    //printf("%.3f  %.3f\n", btag[1], btag[2]);
    btag[0] = discriminator;
    vItBTagMatched.push_back(std::pair<reco::PFJetCollection::const_iterator, std::vector<double> >(itBTagMatched, btag));
  }
  
  // Load jet energy correction service
  const JetCorrector* corrector = JetCorrector::getJetCorrector(mJetCorr, iSetup);

  // Loop over needed jets
  int status = 1;
  for (reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); it++)
  {
    if(_Njet == _maxNjet)
    {
      printf("Maximum number of jets %d reached, skipping the rest\n", _maxNjet);
      return 0;
    }
    //printf("%f\n", it->pt());
    // Apply jet energy correction (JEC)
    double jec = corrector->correction(*it, iEvent, iSetup);
    // copy original (uncorrected) jet;
    reco::PFJet corjet = *it;
    // apply JEC
    corjet.scaleEnergy(jec);
    // store jet
    if(corjet.pt() < 30)
      continue;
    if(TMath::Abs(corjet.eta()) > 2.4)
      continue;
    //TLorentzVector vJet;
    //vJet.SetPtEtaPhiM(corjet.pt(), corjet.eta(), corjet.phi(), corjet.mass());
    //_jetMom[0][_Njet] = vJet.Px();
    //_jetMom[1][_Njet] = vJet.Py();
    //_jetMom[2][_Njet] = vJet.Pz();
    //_jetMom[3][_Njet] = vJet.M();
    _jetPt[_Njet] = corjet.pt();
    _jetEta[_Njet] = corjet.eta();
    _jetPhi[_Njet] = corjet.phi();
    _jetMass[_Njet] = corjet.mass();
    //printf("%f %f %f %f\n", _jetMom[0][_Njet], _jetMom[1][_Njet], _jetMom[2][_Njet], _jetMom[3][_Njet]);
    _jetMuEn[_Njet] = corjet.muonEnergy();
    _jetElEn[_Njet] = corjet.electronEnergy();
    _jetBTagDiscr[_Njet] = _jetBTagMatchDiff1[_Njet] = _jetBTagMatchDiff2[_Njet] = -1.0;
    for(std::vector<std::pair<reco::PFJetCollection::const_iterator, std::vector<double> > >::const_iterator ittag = vItBTagMatched.begin(); ittag != vItBTagMatched.end(); ittag++)
    {
      if(it == ittag->first)
      {
        status = 0;
        _jetBTagDiscr[_Njet] = (ittag->second)[0];
        _jetBTagMatchDiff1[_Njet] = (ittag->second)[1];
        _jetBTagMatchDiff2[_Njet] = (ittag->second)[2];
        break;
      }
    }
    _Njet++;
  }
  if(_Njet < 2)
    status = 1;
  return status;
}

// returns vector of integers which are needed trigger bits
// (called in the beginning of each run)
void Analyzer::FindTriggerBits(const HLTConfigProvider& trigConf)
{
  _vecTriggerBits.clear();
  _vecTriggerBits.resize(_vecTriggerNames.size());
  std::vector<std::string> trigNames;
  trigNames = trigConf.triggerNames();
  for(unsigned int i = 0; i < trigNames.size(); i++)
  {
    std::string currentName = trigConf.triggerNames()[i];
    //printf("%5d  %s\n", i, currentName.c_str());
    for(unsigned int n = 0; n < _vecTriggerNames.size(); n++)
    {
      if(currentName.find(_vecTriggerNames[n]) != std::string::npos)
      {
        _vecTriggerBits[n].push_back(i);
      }
    }
  }
}

// print found trigger indices
void Analyzer::PrintTriggerBits()
{
  printf("********* Trigger Bits: **********\n");
  for(unsigned int n = 0; n < _vecTriggerNames.size(); n++)
  {
    printf("%s: ", _vecTriggerNames[n].c_str());
    for(unsigned int i = 0; i < _vecTriggerBits[n].size(); i++)
      printf(" %d ", _vecTriggerBits[n][i]);
    printf("\n");
  }
}

// fill trigger bits
void Analyzer::SelectTriggerBits(const edm::Handle<edm::TriggerResults>& HLTR)
{
  for(unsigned int i = 0; i < _vecTriggerBits.size(); i++)
  {
    int status = 0;
    for(unsigned int j = 0; j < _vecTriggerBits[i].size(); j++)
    {
      status = status || HLTR->accept(_vecTriggerBits[i][j]);
    }
    //printf("TRIGGER:  %d %s\n", status, _vecTriggerNames[i].c_str());
    _triggers ^= (-status ^ _triggers) & (1 << i);
  }
  //printf("*************\n");
}

// select primary vertex
int Analyzer::SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex)
{
  bool sel = false;
  if(primVertex->size() == 0)
    return false;
  reco::VertexCollection::const_iterator pv = primVertex->begin();
  _pvZ = pv->z();
  _pvRho = TMath::Sqrt(TMath::Power(pv->x(), 2.0) + TMath::Power(pv->y(), 2.0));
  _Npv = primVertex->size();
  _pvNDOF = pv->ndof();
  sel = true;
  return sel;
}

// get final-state stable particle
const reco::Candidate* Analyzer::GetFinalState(const reco::Candidate* particle, const int id)
{
  for(unsigned int i = 0; i < particle->numberOfDaughters(); i++)
  {
    const reco::Candidate* daughter = particle->daughter(i);
    if(daughter->pdgId() == id && daughter->status() == 1)
      return daughter;
    const reco::Candidate* result = GetFinalState(daughter, id);
    if(result)
      return result;
  }
  return NULL;
}

/*
// find W leptonic decay
std::vector<const reco::Candidate*> Analyzer::FindGenDecayW(const reco::Candidate* particle, const std::vector<int>& leptonID)
{
  // vector to store the results
  std::vector<const reco::Candidate*>  matchedDaughters;
  
  // check provided id
  //if(leptonID != 11 && leptonID != 13 && leptonID != 15)
  //{
  //  printf("Error in FindGenDecayW: wrong provided id %d\n", leptonID);
  //  return mathcedDaughters;
  //}

  // status of W must be 3: hard-scattering particle
  if(particle->status() != 3)
    return matchedDaughters;
  
  // parent sign
  int sign = TMath::Sign(1, particle->pdgId());
  
  // get number of W daughters  
  int nDaughters = particle->numberOfDaughters();
  if(nDaughters < 2)
    return matchedDaughters;
  
  // indices of found daughters
  const reco::Candidate* daughterIndices[2] = {NULL, NULL};
  for(int i = 0; i < nDaughters; i++)
  {
    const reco::Candidate* daughter = particle->daughter(i);
    // lepton and nu must be hard-scattering particles
    if(daughter->status() != 3)
      continue;
    for(unsigned int j = 0; j < leptonID.size(); j++)
    {
      
      if( (daughter->pdgId() * sign) == (leptonID[j] * -1) )
        daughterIndices[0] = daughter; // lepton found
      if( (daughter->pdgId() * sign) == (leptonID[j] + 1) )
        daughterIndices[1] = daughter; // nu found
    }
  }
  
  // check whether all daughters were found
  if(!daughterIndices[0] || !daughterIndices[1])
    return matchedDaughters;
  
  // now get final state lepton and mu
  matchedDaughters.resize(2);
  matchedDaughters[0] = GetFinalState(daughterIndices[0], daughterIndices[0]->pdgId());
  if(!matchedDaughters[0])
  {
    return std::vector<const reco::Candidate*>();
    printf("Error in FindGenDecayW(): no stable daughter for lepton %d\n", daughterIndices[0]->pdgId());
  }
  matchedDaughters[1] = GetFinalState(daughterIndices[1], daughterIndices[1]->pdgId());
  if(!matchedDaughters[1])
  {
    return std::vector<const reco::Candidate*>();
    printf("Error in FindGenDecayW(): no stable daughter for neutrino %d\n", daughterIndices[1]->pdgId());
  }
  return matchedDaughters;
}

// find exact decay
std::vector<const reco::Candidate*> Analyzer::FindGenDecay(const reco::Candidate* particle, const std::vector<int>& daugterIDs)
{
  // vector to store the results
  std::vector<const reco::Candidate*> mathcedDaughters;
  
  // number of daughters of the provided particle
  unsigned int nDaughters = particle->numberOfDaughters();
  
  // check if the number of daughters matches the number of provided IDs
  if(nDaughters == daugterIDs.size())
  {
    // now check if all daughters match provided IDs
    bool matched = true;
    // vector to store daighter indices
    std::vector<int> daughterIDFound(nDaughters, -1);
    // loop over daugters of the provided particle
    for(unsigned int i = 0; i < nDaughters; i++)
    {
      bool thisDaughterFound = false;
      // loop over provided daugter IDs
      for(unsigned int j = 0; j < nDaughters; j++)
      {
        // skip those which are already found
        if(daughterIDFound[j] != -1)
          continue;
        // check if the two match
        if(particle->daughter(i)->pdgId() * TMath::Sign(1, particle->pdgId()) == daugterIDs[i])
        {
          // store index, mark this daughter as "found", break the loop over IDs
          daughterIDFound[j] = j;
          thisDaughterFound = true;
          break;
        }
      }
      // if the current daughter does not match any ID, break the loop over daughters
      if(!thisDaughterFound)
      {
        matched = false;
        break;
      }
    }
    
    // check whether all daughters have matched provided IDs
    if(matched)
    {
      // all found, fill pointers and return
      mathcedDaughters.resize(nDaughters);
      for(unsigned int i = 0; i < daughterIDFound.size(); i++)
      {
        mathcedDaughters[i] = particle->daughter(daughterIDFound[i]);
      }
      return mathcedDaughters;
    }
  }
    
  // daughters of the provided particle do not match all IDs: 
  // now call this function iteratively for each daughter
  for(unsigned int i = 0; i < nDaughters; i++)
  {
    const reco::Candidate* daughter = particle->daughter(i);
    std::vector<const reco::Candidate*> mathcedDaughters = FindGenDecay(daughter, daugterIDs);
    // if the returned vector is not empty, this is the needed decay: return it
    if(mathcedDaughters.size() > 0)
      return mathcedDaughters;
  }
    
  // nothing found: return empty vector
  return mathcedDaughters;
}
*/

// fill 4-momentum with provided particle pointer
void Analyzer::FillFourMomentum(const reco::Candidate* particle, float* p)
{
  if(particle == NULL)
  {
    p[0] = p[1] = p[2] = p[3] = 0.0;
    return;
  }
  
  p[0] = particle->px();
  p[1] = particle->py();
  p[2] = particle->pz();
  p[3] = particle->mass();
}

void Analyzer::SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles)
{
  const reco::Candidate* genT = NULL;
  const reco::Candidate* genTbar = NULL;
  const reco::Candidate* genB = NULL;
  const reco::Candidate* genBbar = NULL;
  const reco::Candidate* genWp = NULL;
  const reco::Candidate* genWm = NULL;
  const reco::Candidate* genLp = NULL;
  const reco::Candidate* genNu = NULL;
  const reco::Candidate* genLm = NULL;
  const reco::Candidate* genNubar = NULL;
  for(unsigned int p = 0; p < genParticles->size(); p++)
  {
    const reco::Candidate* particle = &genParticles->at(p);
    const bool sign = (particle->pdgId() > 0); // true for top, false for antitop
    
    // find hard-scattering top
    if(particle->pdgId() == (sign ? 6 : -6) && particle->status() == 3) 
    {
      const reco::Candidate*& t = (sign ? genT : genTbar);
      const reco::Candidate*& b = (sign ? genB : genBbar);
      const reco::Candidate*& W = (sign ? genWp : genWm);
      const reco::Candidate*& l = (sign ? genLp : genLm);
      const reco::Candidate*& nu = (sign ? genNu : genNubar);
      if(t) printf("Error: multiple hard-scattering t\n");
      t = particle;
      // find t -> bW decay
      for(unsigned int d = 0; d < t->numberOfDaughters(); d++)
      {
        const reco::Candidate* daughter = t->daughter(d);
        if(daughter->status() != 3)
          continue;
        if(daughter->pdgId() == (sign ? 5 : -5))
        {
          if(b) printf("Error: multiple hard-scattering b\n");
          b = daughter;
        }
        if(daughter->pdgId() == (sign ? 24 : -24))
        {
          if(W) printf("Error: multiple hard-scattering W\n");
          W = daughter;
        }
      }
      if(!b || !W)
        continue;
      // find W -> lnu decay
      for(unsigned int d = 0; d < W->numberOfDaughters(); d++)
      {
        const reco::Candidate* daughter = W->daughter(d);
        if(daughter->status() != 3)
          continue;
        // electron
        if(daughter->pdgId() == (sign ? -11 : 11))
        {
          if(l) printf("Error: multiple hard-scattering l\n");
          l = daughter;
        }
        if(daughter->pdgId() == (sign ? 12 : -12))
        {
          if(nu) printf("Error: multiple hard-scattering nu\n");
          nu = daughter;
        }
        // muon
        if(daughter->pdgId() == (sign ? -13 : 13))
        {
          if(l) printf("Error: multiple hard-scattering l\n");
          l = daughter;
        }
        if(daughter->pdgId() == (sign ? 14 : -14))
        {
          if(nu) printf("Error: multiple hard-scattering nu\n");
          nu = daughter;
        }
        // tau
        if(daughter->pdgId() == (sign ? -15 : 15))
        {
          if(l) printf("Error: multiple hard-scattering l\n");
          l = daughter;
        }
        if(daughter->pdgId() == (sign ? 16 : -16))
        {
          if(nu) printf("Error: multiple hard-scattering nu\n");
          nu = daughter;
        }
      }
      if(!l || !nu)
        continue;
      // get final states
      l = GetFinalState(l, l->pdgId());
      nu = GetFinalState(nu, nu->pdgId());
    }
  }

  // fill branch variables
  FillFourMomentum(genT, _mcT);
  FillFourMomentum(genTbar, _mcTbar);
  FillFourMomentum(genB, _mcB);
  FillFourMomentum(genBbar, _mcBbar);
  FillFourMomentum(genWp, _mcWp);
  FillFourMomentum(genWm, _mcWm);
  FillFourMomentum(genLp, _mcLp);
  FillFourMomentum(genNu, _mcNu);
  FillFourMomentum(genLm, _mcLm);
  FillFourMomentum(genNubar, _mcNubar);
  
  // now classify this generated event
  _mcEventType = 0;
  if(genT && genTbar)
    _mcNTTbar++;
  if(genWp && genB)
    _mcNTWb++;
  if(genWm && genBbar)
    _mcNTWb++;
  if(genLp && genLm && genNu && genNubar)
  {
    _mcNTtbarDilepton++;
    if(genLp->pdgId() == -11 && genLm->pdgId() == 11 && genNu->pdgId() == 12 && genNubar->pdgId() == -12)
    {
      _mcEventType = 1;
      _mcNTtbarDileptonEE++;
    }
    if(genLp->pdgId() == -13 && genLm->pdgId() == 13 && genNu->pdgId() == 14 && genNubar->pdgId() == -14)
    {
      _mcEventType = 2;
      _mcNTtbarDileptonMuMu++;
    }
    if(genLp->pdgId() == -11 && genLm->pdgId() == 13 && genNu->pdgId() == 12 && genNubar->pdgId() == -14)
    {
      _mcEventType = 3;
      _mcNTtbarDileptonEMu++;
    }
    if(genLp->pdgId() == -13 && genLm->pdgId() == 11 && genNu->pdgId() == 14 && genNubar->pdgId() == -12)
    {
      _mcEventType = 3;
      _mcNTtbarDileptonEMu++;
    }
    if(genLp->pdgId() == -15 && genLm->pdgId() == 15 && genNu->pdgId() == 16 && genNubar->pdgId() == -16)
      _mcNTtbarDileptonTauTau++;
  }
}

// ------------ method called for each event  ------------
void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
	using namespace reco;
	using namespace std;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //Handle<ExampleData> pIn;
   //iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   //ESHandle<SetupData> pSetup;
   //iSetup.get<SetupRecord>().get(pSetup);
#endif

  // event counting
  _nevents++;
  //printf("*** EVENT %6d ***\n", _nevents);
  if( (_nevents % 1000) == 0)
  {
    //printf("*****************************************************************\n");
    printf("************* NEVENTS = %d K, selected = %d *************\n", _nevents / 1000, _neventsSelected);
    //printf("*****************************************************************\n");
  }
  //return;
  
  // declare event contents
  edm::Handle<reco::GenParticleCollection> genParticles;
  Handle<reco::VertexCollection> primVertex;
  edm::Handle<reco::MuonCollection> muons;
  edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::JetTagCollection> bTagHandle;
  edm::Handle<edm::View<reco::PFMET> > pfmets;
  Handle<TriggerResults> HLTR;

  // read event contents immediately (may be needed because of strange perfomance issues)
  /*if(_flagGEN)
    iEvent.getByLabel(_inputTagMCgen, genParticles);
  if(_flagRECO)
  {
    iEvent.getByLabel(_inputTagPrimaryVertex, primVertex);
    iEvent.getByLabel(_inputTagElectrons, electrons);
    iEvent.getByLabel(_inputTagMuons, muons);
    iEvent.getByLabel(_inputTagJets, jets);
    iEvent.getByLabel(_inputTagBtags, bTagHandle);
    iEvent.getByLabel(_inputTagMet,pfmets);
    iEvent.getByLabel(_inputTagTriggerResults, HLTR);
  }*/

  // event selection
  InitBranchVars();
  bool selGEN = false;
  if(_flagGEN)
  {
    iEvent.getByLabel(_inputTagMCgen, genParticles);
    SelectMCGen(genParticles);
    if(_mcEventType != 0)
      selGEN = true;
    if(!selGEN && !_flagRECO)
      return;
  }
  bool selRECO = false;
  if(_flagRECO)
  {
    // primary vertex
    iEvent.getByLabel(_inputTagPrimaryVertex, primVertex);
    reco::VertexCollection::const_iterator pv = primVertex->begin();
    // electrons
    _signLeptonP = _signLeptonM = 0;
    iEvent.getByLabel(_inputTagElectrons, electrons);
    SelectEl(electrons, pv);
    // muons
    iEvent.getByLabel(_inputTagMuons, muons);
    SelectMu(muons, pv);
    // require pair of opposite signed leptons
    if( _signLeptonP && _signLeptonM )
      selRECO = true;
    if(!selRECO && !selGEN)
      return;
    // jets and b-tagging
    iEvent.getByLabel(_inputTagJets, jets);
    iEvent.getByLabel(_inputTagBtags, bTagHandle);
    const reco::JetTagCollection& bTags = *(bTagHandle.product());
    SelectJet(jets, bTags, iEvent, iSetup);
    // require two jets
    if( _Njet >= 2)
      selRECO = true;
    if(!selRECO && !selGEN)
      return;
    // MET
    iEvent.getByLabel(_inputTagMet, pfmets);
    SelectMET(pfmets);
    // primary vertex
    SelectPrimaryVertex(primVertex);
    // triggers
    iEvent.getByLabel(_inputTagTriggerResults, HLTR);
    SelectTriggerBits(HLTR);
  }
  
  // event info
  SelectEvent(iEvent);
  // store event
  _tree->Fill();
  _neventsSelected++;
}


// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob() {;}


// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() {;}

// ------------ method called when starting to processes a run  ------------
void Analyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  // trigger stuff
  HLTConfigProvider triggerConfig;
  bool changed = true;
  triggerConfig.init(iRun, iSetup, _inputTagTriggerResults.process(), changed);
  FindTriggerBits(triggerConfig);
}

// ------------ method called when ending the processing of a run  ------------
void Analyzer::endRun(edm::Run const& run, edm::EventSetup const& setup) {;}

// ------------ method called when starting to processes a luminosity block  ------------
void Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {;}

// ------------ method called when ending the processing of a luminosity block  ------------
void Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {;}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);
