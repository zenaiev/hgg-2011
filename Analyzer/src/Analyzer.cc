// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc ttbar/Analyzer/src/Analyzer.cc

 Description: ttbar ntuple production

 Implementation:
     for Open Data 2011
*/


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
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for photons
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"

// for electrons
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// for jets
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/BTauReco/interface/JetTag.h>
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// for MET
#include "DataFormats/METReco/interface/PFMET.h"

//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// for MC generator level
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ROOT
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>

#include <TTree.h>
#include <TDirectory.h>

//
// class declaration
//
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

    // user routines (detailed description given with the method implementations)
    int SelectEvent(const edm::Event& iEvent);
    int SelectPhotons(const edm::Handle<reco::PhotonCollection>& photons,const edm::Handle<reco::GsfElectronCollection>& electrons, const reco::VertexCollection::const_iterator& pv, const Handle<reco::VertexCollection>& Vertices, const edm::Handle<reco::PFCandidateCollection>& PF);
    int SelectJet(const edm::Handle<reco::PFJetCollection>& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup);
    int SelectMET(const edm::Handle<edm::View<reco::PFMET> >& pfmets);
    void FindTriggerBits(const HLTConfigProvider& trigConf);
    void SelectTriggerBits(const edm::Handle<edm::TriggerResults>& HLTR);
    void PrintTriggerBits();
    int SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex);
    //const reco::Candidate* GetFinalState(const reco::Candidate* particle, const int id);
    void FillFourMomentum(const reco::Candidate* particle, float* p);
    double MatchPhoton(const double etaGen, const double phiGen, const double etaRec, const double phiRec, const float max = 0.1);
    void SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles);
    void InitBranchVars();

    // input tags
    edm::InputTag _inputTagPhotons;
    edm::InputTag _inputTagElectrons;
    edm::InputTag _inputTagJets;
    edm::InputTag _inputTagMet;
    edm::InputTag _inputTagPF;
    edm::InputTag _inputTagTriggerResults;
    edm::InputTag _inputTagPrimaryVertex;
    edm::InputTag _inputTagRho;
    edm::InputTag _inputTagMCgen;

    // jet correction label
    std::string mJetCorr;
    
    // photon isolation
    PFIsolationEstimator mIsolator;

    // general flags and variables
    int _flagMC;
    int _flagRECO;
    int _flagGEN;
    int _flagYEAR;
    int _nevents;
    int _neventsSelected;

    // storage
    TFile* _file;
    TTree* _tree;

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>> event variables >>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (their description given when tree branches are created)
    // event
    int _evRunNumber;
    int _evEventNumber;
    double _rho;
    // photons
    static const int _maxNph = 50;
    int _Nph;
    float _phE[_maxNph];
    float _phPt[_maxNph];
    float _phEta[_maxNph];
    float _phPhi[_maxNph];
    float _phSCx[_maxNph];
    float _phSCy[_maxNph];
    float _phSCz[_maxNph];
    float _phR9[_maxNph];
    float _phTrkSumPtHollowConeDR04[_maxNph];
    float _phEcalRecHitSumEtConeDR04[_maxNph];
    float _phHcalTowerSumEtConeDR04[_maxNph];
    //float _phHcalDepth1TowerSumEtConeDR04[_maxNph];
    //float _phHcalDepth2TowerSumEtConeDR04[_maxNph];
    float _phTrkSumPtHollowConeDR03[_maxNph];
    float _phEcalRecHitSumEtConeDR03[_maxNph];
    float _phHcalTowerSumEtConeDR03[_maxNph];
    //float _phHcalDepth1TowerSumEtConeDR03[_maxNph];
    //float _phHcalDepth2TowerSumEtConeDR03[_maxNph];
    
    //PFlow isolation
    float _phChargedHadronIso[_maxNph];
    float _phNeutralHadronIso[_maxNph];
    float _phPhotonIso[_maxNph];
    float _phIsolationSum[_maxNph];
    float _phIsolationSumWrongVtx[_maxNph];
    
    //electrons in same SuperCluster
    int _phNumElectronsSuperCluster[_maxNph];
    int _elMissingHits[_maxNph];
    float _phElectronDR[_maxNph];
    int _phHasConversionTracks[_maxNph];
    
    // H/E
    float _phHadronicOverEm[_maxNph];
    // covietaieta
    float _phSigmaIetaIeta[_maxNph];
    // matching
    float _phMatch[_maxNph];
    float _phMatchDeltaE[_maxNph];
    // jets
    static const int _maxNjet = 50;
    int _Njet;
    float _jetPt[_maxNjet];
    float _jetEta[_maxNjet];
    float _jetPhi[_maxNjet];
    float _jetMass[_maxNjet];
    float _jetMuEn[_maxNjet];
    float _jetElEn[_maxNjet];
    // MET
    float _metPx;
    float _metPy;
    // triggers
    int _triggers;
    std::vector<std::vector<int> > _vecTriggerBits;
    std::vector<std::string> _vecTriggerNames;
    std::vector<std::vector<std::string> > _vecTriggerNamesFull;
    // primary vertex
    int _Npv;
    int _pvNDOF;
    float _pvZ;
    float _pvRho;
    // MC generated info
    int _mcEventType;
    int _mcNHgg;
    // generator level four vectors
    float _mcH[4];
    float _mcPh[2][4];
    TLorentzVector _mcPhLV[2];
};

//
// constructor
//
Analyzer::Analyzer(const edm::ParameterSet& iConfig)
{
  // for proper log files writing (immediate output)
  setbuf(stdout, NULL);
  
  // input tags
  _inputTagPhotons = edm::InputTag("photons");
  _inputTagElectrons = edm::InputTag("gsfElectrons");
  _inputTagJets = edm::InputTag("ak5PFJets");
  _inputTagMet = edm::InputTag("pfMet");
  _inputTagPF = edm::InputTag("particleFlow");
  _inputTagTriggerResults = edm::InputTag("TriggerResults", "", "HLT");
  _inputTagPrimaryVertex = edm::InputTag("offlinePrimaryVerticesWithBS");
  _inputTagRho = edm::InputTag("fixedGridRhoAll");
  _inputTagMCgen = edm::InputTag("genParticles");
  
  // jet correction label
  mJetCorr = "ak5PFL1FastL2L3Residual";
  
  //photon isolator 
  mIsolator.initializePhotonIsolation(kTRUE);
  mIsolator.setConeSize(0.3);

  // read configuration parameters
  _flagMC = iConfig.getParameter<int>("mc"); // true for MC, false for data
  _flagRECO = iConfig.getParameter<int>("reco"); // if true, RECO level processed
  _flagGEN = iConfig.getParameter<int>("gen"); // if true, generator level processed (works only for MC)
  _flagYEAR = iConfig.getParameter<int>("year"); // 0 for 2011, 1 for 2012
  _nevents = 0; // number of processed events
  _neventsSelected = 0; // number of selected events
  std::string fileout = iConfig.getParameter<std::string>("outFile"); // output file name
  _file = new TFile(fileout.c_str(), "recreate"); // output file
  _tree = new TTree("tree", "hgg"); // output tree

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>> tree branches >>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // event
  _tree->Branch("evRunNumber", &_evRunNumber, "evRunNumber/I"); // run number
  _tree->Branch("evEventNumber", &_evEventNumber, "evEventNumber/I"); // event number
  _tree->Branch("rho", &_rho, "rho/F"); // rho

  if(_flagRECO)
  {
    // muons
    _tree->Branch("Nph", &_Nph, "Nph/I"); // number of photons
    _tree->Branch("phE", _phE, "phE[Nph]/F"); // photon E
    _tree->Branch("phPt", _phPt, "phPt[Nph]/F"); // photon pT
    _tree->Branch("phEta", _phEta, "phEta[Nph]/F"); // photon pseudorapidity
    _tree->Branch("phPhi", _phPhi, "phPhi[Nph]/F"); // photon phi
    _tree->Branch("phPt", _phSCx, "phSCx[Nph]/F"); // photon supercluster x
    _tree->Branch("phPt", _phSCy, "phSCy[Nph]/F"); // photon supercluster y
    _tree->Branch("phEta", _phSCz, "phSCz[Nph]/F"); // photon supercluster z
    _tree->Branch("phR9", _phR9, "phR9[Nph]/F"); // photon R9
    _tree->Branch("phTrkSumPtHollowConeDR04", _phTrkSumPtHollowConeDR04, "phTrkSumPtHollowConeDR04[Nph]/F"); // photon R9
    _tree->Branch("phEcalRecHitSumEtConeDR04", _phR9, "phEcalRecHitSumEtConeDR04[Nph]/F"); // photon R9
    _tree->Branch("phHcalTowerSumEtConeDR04", _phHcalTowerSumEtConeDR04, "phHcalTowerSumEtConeDR04[Nph]/F"); // photon R9
    //_tree->Branch("phHcalDepth1TowerSumEtConeDR04", _phHcalDepth1TowerSumEtConeDR04, "phHcalDepth1TowerSumEtConeDR04[Nph]/F"); // photon R9
    //_tree->Branch("phHcalDepth2TowerSumEtConeDR04", _phHcalDepth2TowerSumEtConeDR04, "phHcalDepth2TowerSumEtConeDR04[Nph]/F"); // photon R9
    _tree->Branch("phTrkSumPtHollowConeDR03", _phTrkSumPtHollowConeDR03, "phTrkSumPtHollowConeDR03[Nph]/F"); // photon R9
    _tree->Branch("phEcalRecHitSumEtConeDR03", _phEcalRecHitSumEtConeDR03, "phEcalRecHitSumEtConeDR03[Nph]/F"); // photon R9
    _tree->Branch("phHcalTowerSumEtConeDR03", _phHcalTowerSumEtConeDR03, "phHcalTowerSumEtConeDR03[Nph]/F"); // photon R9
    //_tree->Branch("phHcalDepth1TowerSumEtConeDR03", _phHcalDepth1TowerSumEtConeDR03, "phHcalDepth1TowerSumEtConeDR03[Nph]/F"); // photon R9
    //_tree->Branch("phHcalDepth2TowerSumEtConeDR03", _phHcalDepth2TowerSumEtConeDR03, "phHcalDepth2TowerSumEtConeDR03[Nph]/F"); // photon R9
    _tree->Branch("phHadronicOverEm", _phHadronicOverEm, "phHadronicOverEm[Nph]/F"); // photon R9
    _tree->Branch("phSigmaIetaIeta", _phSigmaIetaIeta, "phSigmaIetaIeta[Nph]/F"); // photon R9
	  //PFlow isolation
    _tree->Branch("phChargedHadronIso",_phChargedHadronIso , "phChargedHadronIso[Nph]/F");
    _tree->Branch("phNeutralHadronIso", _phNeutralHadronIso , "phNeutralHadronIso[Nph]/F");
	  _tree->Branch("phPhotonIso",_phPhotonIso,"phPhotonIso[Nph]/F");
    _tree->Branch("phIsolationSum", _phIsolationSum, "phIsolationSum[Nph]/F");
    _tree->Branch("phIsolationSumWrongVtx",_phIsolationSumWrongVtx, "phIsolationSumWrongVtx[Nph]/F");
    
    //electrons in same superCluster (sc)
    _tree->Branch("phNumElectronsSuperCluster",_phNumElectronsSuperCluster,"phNumElectronsSuperCluster[Nph]/I");
    _tree->Branch("elMissingHits", _elMissingHits, "elMissingHits[Nph]/I");
    _tree->Branch("phElectronDR", _phElectronDR, "phElectronDR[Nph]/F");
    _tree->Branch("phHasConversionTracks" , _phHasConversionTracks, "phHasConversionTracks[Nph]/I");
    
    
    // jets
    _tree->Branch("Njet", &_Njet, "Njet/I"); // number of jets
    _tree->Branch("jetPt", _jetPt, "jetPt[Njet]/F"); // jet pT
    _tree->Branch("jetEta", _jetEta, "jetEta[Njet]/F"); // jet pseudorapidity
    _tree->Branch("jetPhi", _jetPhi, "jetPhi[Njet]/F"); // jet phi
    _tree->Branch("jetMass", _jetMass, "jetMass[Njet]/F"); // jet mass
    _tree->Branch("jetMuEn", _jetMuEn, "jetMuEn[Njet]/F"); // jet muon energy
    _tree->Branch("jetElEn", _jetElEn, "jetElEn[Njet]/F"); // jet electron energy
    // MET
    //_tree->Branch("metPx", &_metPx, "metPx/F"); // missing transverse energy x component
    //_tree->Branch("metPy", &_metPy, "metPy/F"); // missing transverse energy y component
    // triggers
    _tree->Branch("Triggers", &_triggers, "Triggers/I"); // trigger bits (see trigger names below)
    //_vecTriggerNames.push_back("HLT_Photon26_CaloID_Iso_Photon18_CaloID_Iso");
    //_vecTriggerNames.push_back("HLT_Photon26_R9ID_Photon18_R9ID");
    //_vecTriggerNames.push_back("HLT_Photon26_R9ID_Photon18_CaloID_Iso");
    //_vecTriggerNames.push_back("HLT_Photon26_CaloID_Iso_Photon18_R9ID");
    _vecTriggerNames.push_back("HLT_Photon20");
    _vecTriggerNames.push_back("HLT_Photon26");
    _vecTriggerNames.push_back("HLT_Photon36");
    // primary vertex
    _tree->Branch("Npv", &_Npv, "Npv/I"); // total number of primary vertices
    _tree->Branch("pvNDOF", &_pvNDOF, "pvNDOF/I"); // number of degrees of freedom of the primary vertex
    _tree->Branch("pvZ", &_pvZ, "pvZ/F"); // z component of the primary vertex
    _tree->Branch("pvRho", &_pvRho, "pvRho/F"); // rho of the primary vertex (projection on transverse plane)
  }

  // MC generated info
  if(_flagGEN)
  {
    // matching
    _tree->Branch("phMatch", _phMatch, "phMatch[Nph]/F"); // photon R9
    _tree->Branch("phMatchDeltaE", _phMatchDeltaE, "phMatchDeltaE[Nph]/F"); // photon R9
    // true level info
    _mcNHgg = 0; // number of H->gg events
    // MC generated info to store
    _tree->Branch("mcEventType", &_mcEventType, "mcEventType/I"); // MC generator level event type: 1 signal, 0 anything else
    _tree->Branch("mcH", _mcH, "mcH[4]/F"); // generator level Higgs four vector
    _tree->Branch("mcPh", _mcPh, "mcPh[2][4]/F"); // generator level photons four vectors
  }
}


// destructor
Analyzer::~Analyzer()
{
  // close files, deallocate resources etc.
  _file->cd();
  _tree->Write();
  _file->Close();

  if(_flagGEN)
  {
    // MC generator level info printout
    printf("%25s = %10d\n", "mcNHgg", _mcNHgg);
  }

  // print total number of processed and selected events
  printf("Processed %d events, selected %d\n", _nevents, _neventsSelected);
}


//
// member functions
//

// initialise event variables with needed default (zero) values; called in the beginning of each event
void Analyzer::InitBranchVars()
{
  _evRunNumber = 0;
  _evEventNumber = 0;
  _Nph = 0;
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

// Store event info (fill corresponding tree variables)
int Analyzer::SelectEvent(const edm::Event& iEvent)
{
  _evRunNumber = iEvent.id().run();
  _evEventNumber = iEvent.id().event();

  // pileup
  edm::Handle<double> rho;
  iEvent.getByLabel(_inputTagRho, rho);
  _rho = rho.failedToGet() ? -999. : *rho;

  return 0;
}

// Missing transverse energy (MET) selection
int Analyzer::SelectMET(const edm::Handle<edm::View<reco::PFMET> >& pfmets)
{
  _metPx = (pfmets->front()).px();
  _metPy = (pfmets->front()).py();
  return 0;
}

// photon selection
int Analyzer::SelectPhotons(const edm::Handle<reco::PhotonCollection>& photons,const edm::Handle<reco::GsfElectronCollection>& electrons, const reco::VertexCollection::const_iterator& pv, const Handle<reco::VertexCollection>& Vertices, const edm::Handle<reco::PFCandidateCollection>& PF)
{
  _Nph = 0;
  // loop over photons
  for (reco::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); it++)
  {
    if(_Nph == _maxNph)
    {
      printf("Maximum number of photons %d reached, skipping the rest\n", _maxNph);
      return 0;
    }
    
    //get the supercluster of the photon
    float scEta = it->superCluster()->position().eta();
    float scPhi = it->superCluster()->position().phi();
    float scEtaWidth = it->superCluster()->etaWidth();
    float scPhiWidth = it->superCluster()->phiWidth();
   
    _phNumElectronsSuperCluster[_Nph] = 0;
    _elMissingHits[_Nph] = -1;
    //check if there are any electrons within the supercluster 
    //maybe there are other methods for PFlow in 2012 to find supercluster
    for (reco::GsfElectronCollection::const_iterator itEl = electrons->begin(); itEl != electrons->end(); itEl++)
    {
      float scEtaEl = itEl->superCluster()->position().eta();
      float scPhiEl = itEl->superCluster()->position().phi();
      if(scEtaEl > scEta + scEtaWidth / 2 || scEtaEl < scEta - scEtaWidth / 2)
        continue;
      if(scPhiEl > scPhi + scPhiWidth / 2 || scPhiEl < scPhi - scPhiWidth / 2)
        continue;
      //energy cut on electron
      if(itEl->pt() < 2.5)
        continue;
      //increase electron counter
      _phNumElectronsSuperCluster[_Nph] += 1; //should never be higher than 1
      //store number of (missing) hits 
      _elMissingHits[_Nph] = itEl->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      //deltaR for 2011 cut
      _phElectronDR[_Nph] = TMath::Sqrt(TMath::Power(it->eta()-itEl->eta(),2.0)+TMath::Power(it->phi()-itEl->phi(),2.0));
      
    }
    
    
    // selection for 2011 data
    if(_flagYEAR == 0)
    {
      if(it->pt() < 25.0)
        continue;
      if(TMath::Abs(it->eta()) > 2.5 || (TMath::Abs(it->eta()) > 1.44 && TMath::Abs(it->eta()) < 1.57))
        continue;
      if(it->r9() < 0.32)
        continue;
      if(it->trkSumPtHollowConeDR03() > 3.5)
        continue;
      if(it->hadronicOverEm() > 0.082)
        continue;
      if(_phNumElectronsSuperCluster[_Nph] > 0 && _elMissingHits[_Nph] == 0)
      {
        //not rejected if it is in class two
        if(TMath::Abs(it->eta()) > 1.44 || it->r9() > 0.94)
          continue;
      }
    } 
    // selection for 2012 data
    if(_flagYEAR == 1)
    {
      if(it->pt() < 25.0)
        continue;
      if(TMath::Abs(it->eta()) > 2.5 || (TMath::Abs(it->eta()) > 1.44 && TMath::Abs(it->eta()) < 1.57))
        continue;
      if(it->r9() < 0.24)
        continue;
      if(it->hadronicOverEm() > 0.142)
        continue;
      //preselection on PFlow
      reco::VertexRef myVtx(Vertices, 0); //chosen vertex is first
      mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
      if(mIsolator.getIsolationCharged() > 3.8)
        continue;
      if(_phNumElectronsSuperCluster[_Nph] > 0 && _elMissingHits[_Nph] == 0 && it->hasConversionTracks() == true)
        continue;
    } 
    //conversion track
    _phHasConversionTracks[_Nph] = it->hasConversionTracks();
    // fill four momentum (pT, eta, phi, E)
    _phE[_Nph] = it->energy();
    // enum P4type { undefined=-1, ecal_standard=0, ecal_photons=1, regression1=2, regression2= 3 } ;
    //_phE[_Nph] = it->getCorrectedEnergy(reco::Photon::ecal_standard);
    //_phE[_Nph] = it->getCorrectedEnergy(reco::Photon::ecal_photons);
    //double corr = _phE[_Nph] / it->energy();
    //printf("corr = %f %f %f %f %f\n", it->energy(), it->getCorrectedEnergy(reco::Photon::ecal_standard), it->getCorrectedEnergy(reco::Photon::ecal_photons),
    //       it->getCorrectedEnergy(reco::Photon::regression1), it->getCorrectedEnergy(reco::Photon::regression2));
    //printf("%f\n", it->pt() * TMath::CosH(it->eta()));
    _phPt[_Nph] = it->pt();
    _phEta[_Nph] = it->eta();
    _phPhi[_Nph] = it->phi();
    // fill distance to primary vertex
    //_muDistPV0[_Nph] = TMath::Sqrt(TMath::Power(pv->x() - it->globalTrack()->vx(), 2.0) + TMath::Power(pv->y() - it->globalTrack()->vy(), 2.0));
    //_muDistPVz[_Nph] = TMath::Abs(pv->z() - it->globalTrack()->vz());
    // store photon

    // supercluster
    _phSCx[_Nph] = it->superCluster()->position().x();
    _phSCy[_Nph] = it->superCluster()->position().y();
    _phSCz[_Nph] = it->superCluster()->position().z();

    // r9 and other vars
    _phR9[_Nph] = it->r9();
    _phTrkSumPtHollowConeDR04[_Nph] = it->trkSumPtHollowConeDR04();
    _phEcalRecHitSumEtConeDR04[_Nph] = it->ecalRecHitSumEtConeDR04();
    _phHcalTowerSumEtConeDR04[_Nph] = it->hcalTowerSumEtConeDR04();
    //_phHcalDepth1TowerSumEtConeDR04[_Nph] = it->hcalDepth1TowerSumEtConeDR04();
    //_phHcalDepth2TowerSumEtConeDR04[_Nph] = it->hcalDepth2TowerSumEtConeDR04();
    _phTrkSumPtHollowConeDR03[_Nph] = it->trkSumPtHollowConeDR03();
    _phEcalRecHitSumEtConeDR03[_Nph] = it->ecalRecHitSumEtConeDR03();
    _phHcalTowerSumEtConeDR03[_Nph] = it->hcalTowerSumEtConeDR03();
    //_phHcalDepth1TowerSumEtConeDR03[_Nph] = it->hcalDepth1TowerSumEtConeDR03();
    //_phHcalDepth2TowerSumEtConeDR03[_Nph] = it->hcalDepth2TowerSumEtConeDR03();
    
	  //PFlow isolation (not in PhotonCollection)
    //use selected vertex, which was first 
    reco::VertexRef myVtx(Vertices, 0);
    mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
	  _phChargedHadronIso[_Nph] = mIsolator.getIsolationCharged();
    _phNeutralHadronIso[_Nph] = mIsolator.getIsolationNeutral();
    _phPhotonIso[_Nph] = mIsolator.getIsolationPhoton();
    _phIsolationSum[_Nph] = _phChargedHadronIso[_Nph] + _phNeutralHadronIso[_Nph] + _phPhotonIso[_Nph];
    
    //worst ChargedHadronIsolation
    _phIsolationSumWrongVtx[_Nph] = _phIsolationSum[_Nph];
    //loop over all possible vertices
    for(int unsigned long i = 0; i < Vertices->size(); i++)
    {
      reco::VertexRef myVtx(Vertices, i);
      mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
      float curVal = mIsolator.getIsolationCharged() + mIsolator.getIsolationPhoton() + mIsolator.getIsolationNeutral();
      if( _phIsolationSumWrongVtx[_Nph] < curVal) 
        _phIsolationSumWrongVtx[_Nph] = curVal;
    }
    
    // H/E
    _phHadronicOverEm[_Nph] = it->hadronicOverEm();
    _phSigmaIetaIeta[_Nph] = it->sigmaIetaIeta();

    // matching
    if(_flagMC && _mcEventType == 1)
    {
      double delta1 = MatchPhoton(_mcPhLV[0].Eta(), _mcPhLV[0].Phi(), _phEta[_Nph], _phPhi[_Nph]);
      double delta2 = MatchPhoton(_mcPhLV[1].Eta(), _mcPhLV[1].Phi(), _phEta[_Nph], _phPhi[_Nph]);
      double delta = std::min(delta1, delta2);
      _phMatch[_Nph] = delta;
      if(delta < 0.1)
        _phMatchDeltaE[_Nph] = _phE[_Nph] / (_mcPhLV[(delta1 < delta2) ? 0 : 1].E());
      else
        _phMatchDeltaE[_Nph] = 0.0;
        
    }
    
    _Nph++;
  }
  return _Nph;
}

// jet selection
int Analyzer::SelectJet(const edm::Handle<reco::PFJetCollection>& jets, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Load jet energy correction service
  const JetCorrector* corrector = JetCorrector::getJetCorrector(mJetCorr, iSetup);

  _Njet = 0;
  // Loop over jets
  for (reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); it++)
  {
    if(_Njet == _maxNjet)
    {
      printf("Maximum number of jets %d reached, skipping the rest\n", _maxNjet);
      return 0;
    }
    // Apply jet energy correction (JEC)
    double jec = corrector->correction(*it, iEvent, iSetup);
    // copy original (uncorrected) jet;
    reco::PFJet corjet = *it;
    // apply JEC
    corjet.scaleEnergy(jec);
    // select jet: pT > 20 GeV, |eta| < 4.7
    if(corjet.pt() < 20)
      continue;
    if(TMath::Abs(corjet.eta()) > 4.7)
      continue;
    // fill jet four momentum (pT, eta, phi, mass)
    _jetPt[_Njet] = corjet.pt();
    _jetEta[_Njet] = corjet.eta();
    _jetPhi[_Njet] = corjet.phi();
    _jetMass[_Njet] = corjet.mass();
    // fill jet muon and electron energy fractions
    _jetMuEn[_Njet] = corjet.muonEnergy();
    _jetElEn[_Njet] = corjet.electronEnergy();
    _Njet++;
  }
  return _Njet;
}

// returns vector of integers which are trigger bits needed in the analysis
// (called in the beginning of each run)
void Analyzer::FindTriggerBits(const HLTConfigProvider& trigConf)
{
  // reset containers
  _vecTriggerBits.clear();
  _vecTriggerBits.resize(_vecTriggerNames.size());
  _vecTriggerNamesFull.clear();
  _vecTriggerNamesFull.resize(_vecTriggerNames.size());
  std::vector<std::string> trigNames;
  trigNames = trigConf.triggerNames();
  // for interesting trigger names find corresponding trigger bits
  for(unsigned int i = 0; i < trigNames.size(); i++)
  {
    std::string currentName = trigConf.triggerNames()[i];
    //printf("%5d  %s\n", i, currentName.c_str());
    for(unsigned int n = 0; n < _vecTriggerNames.size(); n++)
    {
      if(currentName.find(_vecTriggerNames[n]) != std::string::npos)
      {
        _vecTriggerBits[n].push_back(i);
        _vecTriggerNamesFull[n].push_back(currentName);
      }
    }
  }
  PrintTriggerBits();
}

// print found trigger indices
// (for debugging purpose)
void Analyzer::PrintTriggerBits()
{
  printf("********* Trigger Bits: **********\n");
  for(unsigned int n = 0; n < _vecTriggerNames.size(); n++)
  {
    printf("%s: ", _vecTriggerNames[n].c_str());
    for(unsigned int i = 0; i < _vecTriggerBits[n].size(); i++)
    {
      printf(" %s ", _vecTriggerNamesFull[n][i].c_str());
      printf(" %d ", _vecTriggerBits[n][i]);
    }
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
    // set ith bit of _triggers integer to the current trigger bit status
    _triggers ^= (-status ^ _triggers) & (1 << i);
  }
  //printf("*************\n");
}

// select primary vertex
int Analyzer::SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex)
{
  // if no primary vertices in the event, return false status
  if(primVertex->size() == 0)
    return false;
  // take the first primary vertex
  reco::VertexCollection::const_iterator pv = primVertex->begin();
  // fill z and rho (projection on transverse plane)
  _pvZ = pv->z();
  _pvRho = TMath::Sqrt(TMath::Power(pv->x(), 2.0) + TMath::Power(pv->y(), 2.0));
  // fill number of primary veritces
  _Npv = primVertex->size();
  // fill number of degrees of freedom
  _pvNDOF = pv->ndof();
  // return true status
  return true;
}

// get final-state stable generator level particle with required id
// (if not found, return NULL pointer)
/*const reco::Candidate* Analyzer::GetFinalState(const reco::Candidate* particle, const int id)
{
  // loop over daughters
  for(unsigned int i = 0; i < particle->numberOfDaughters(); i++)
  {
    const reco::Candidate* daughter = particle->daughter(i);
    // if this daughter has required id, return its pointer
    if(daughter->pdgId() == id && daughter->status() == 1)
      return daughter;
    // otherwise call itself recursively
    const reco::Candidate* result = GetFinalState(daughter, id);
    // return the result of the recursive call
    if(result)
      return result;
  }
  // if gets here, there are no daughter with required id: return NULL pointer
  return NULL;
}*/

// fill 4-momentum (p) with provided particle pointer
void Analyzer::FillFourMomentum(const reco::Candidate* particle, float* p)
{
  // if NULL pointer provided, initialise with default (zero)
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

double Analyzer::MatchPhoton(const double etaGen, const double phiGen, const double etaRec, const double phiRec, const float max/* = 0.1*/)
{
  double deltaEta = etaGen - etaRec;
  double deltaPhi = phiGen - phiRec;
  if(deltaPhi > TMath::Pi())  deltaPhi -= 2 * TMath::Pi();
  if(deltaPhi < -TMath::Pi()) deltaPhi += 2 * TMath::Pi();
  double delta = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
  if(delta > max)
    return 10.0;
  else
    return delta;
}

// select MC generator level information
// (analysis specific ttbar dileptonic decay)
void Analyzer::SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles)
{
  // reset all particle momenta
  FillFourMomentum(NULL, _mcH);
  FillFourMomentum(NULL, _mcPh[0]);
  FillFourMomentum(NULL, _mcPh[1]);

  // initialise all particle pointers with NULL
  const reco::Candidate* genH = NULL;
  const reco::Candidate* genPh[2] = { NULL, NULL };
  // loop over generated particeles
  for(unsigned int p = 0; p < genParticles->size(); p++)
  {
    const reco::Candidate* particle = &genParticles->at(p);
    
    // find Higgs
    if(particle->pdgId() == 25 && particle->status() == 3)
    {
      // determine sign, use it to set references (t, b, W, l, nu) to either particles or antiparticles
      // there should be no more than 1 Higgs in the event but who knows
      if(genH)
        printf("Error: multiple hard-scattering Higgs\n");
      // determine top decay
      genH = particle;
      // find H -> gg decay
      for(unsigned int d = 0; d < genH->numberOfDaughters(); d++)
      {
        const reco::Candidate* daughter = genH->daughter(d);
        if(daughter->status() != 3)
          continue;
        if(daughter->pdgId() == 22)
        {
          // again there should be no more than 2 photons from Higgs
          if(genPh[1])
            printf("Error: more than 2 photons from Higgs\n");
          else if(genPh[0])
            genPh[1] = daughter;
          else
            genPh[0] = daughter;
        }
      }
    }
  }

  // fill branch variables
  FillFourMomentum(genH, _mcH);
  FillFourMomentum(genPh[0], _mcPh[0]);
  FillFourMomentum(genPh[1], _mcPh[1]);
  for(int g = 0; g < 2; g++)
    _mcPhLV[g].SetXYZM(_mcPh[g][0], _mcPh[g][1], _mcPh[g][2], _mcPh[g][3]);

  // now classify this generated event
  _mcEventType = 0;
  if(genH && genPh[0] && genPh[1])
  {
    _mcEventType = 1;
    _mcNHgg++;
  }
}

// ------------ method called for each event  ------------
void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  // event counting, printout after each 1K processed events
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
  edm::Handle<reco::VertexCollection> primVertex;
  edm::Handle<reco::PFCandidateCollection> PF;
  edm::Handle<reco::PhotonCollection> photons;
  edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<reco::PFJetCollection> jets;
  //edm::Handle<edm::View<reco::PFMET> > pfmets;
  Handle<TriggerResults> HLTR;

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>> event selection >>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // initialise event variables with default values
  InitBranchVars();
  // process generator level, if needed
  bool selGEN = false;
  if(_flagGEN)
  {
    iEvent.getByLabel(_inputTagMCgen, genParticles);
    SelectMCGen(genParticles);
    if(_mcEventType != 0)
      selGEN = true;
    // if nothing interesting at generator level and not required to process reco level, return here
    if(!selGEN && !_flagRECO)
      return;
  }
  // process reco level, if needed
  bool selRECO = false;
  if(_flagRECO)
  {
    // primary vertex
    iEvent.getByLabel(_inputTagPrimaryVertex, primVertex);
    reco::VertexCollection::const_iterator pv = primVertex->begin();
    // photons
    iEvent.getByLabel(_inputTagPhotons, photons);
    //electrons
    iEvent.getByLabel(_inputTagElectrons, electrons);
    //particle Flow
    iEvent.getByLabel(_inputTagPF, PF);
    SelectPhotons(photons,electrons, pv, primVertex, PF);
    // require pair of opposite signed leptons
    if( _Nph >= 2 )
      selRECO = true;
    if(!selRECO && !selGEN)
      return;
    // jets
    iEvent.getByLabel(_inputTagJets, jets);
    SelectJet(jets, iEvent, iSetup);
    // fill MET
    //iEvent.getByLabel(_inputTagMet, pfmets);
    //SelectMET(pfmets);
    // fill primary vertex
    SelectPrimaryVertex(primVertex);
    // fill triggers
    iEvent.getByLabel(_inputTagTriggerResults, HLTR);
    SelectTriggerBits(HLTR);
  }
  // fill event info
  SelectEvent(iEvent);
  // all done: store event
  _tree->Fill();
  _neventsSelected++;
}


// ------------ method called when starting to processes a run  ------------
void Analyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  // trigger stuff
  HLTConfigProvider triggerConfig;
  bool changed = true;
  triggerConfig.init(iRun, iSetup, _inputTagTriggerResults.process(), changed);
  FindTriggerBits(triggerConfig);
}

// below is some default stuff, was not modified

// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob() {;}

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() {;}

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
