// -*- C++ -*-
//
// Package:    Analyzer
// Class:      Analyzer
// 
/**\class Analyzer Analyzer.cc hgg-2011/Analyzer/src/Analyzer.cc

 Description: hgg ntuple production

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
    void FindTriggerBits(const HLTConfigProvider& trigConf);
    void SelectTriggerBits(const edm::Handle<edm::TriggerResults>& HLTR);
    void PrintTriggerBits();
    int SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex);
    void FillFourMomentum(const reco::Candidate* particle, float* p);
    double MatchPhoton(const double etaGen, const double phiGen, const double etaRec, const double phiRec, const float max = 0.1);
    void SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles);
    void InitBranchVars();

    // input tags
    edm::InputTag _inputTagPhotons;
    edm::InputTag _inputTagElectrons;
    edm::InputTag _inputTagMet;
    edm::InputTag _inputTagPF;
    edm::InputTag _inputTagTriggerResults;
    edm::InputTag _inputTagPrimaryVertex;
    edm::InputTag _inputTagRho;
    edm::InputTag _inputTagMCgen;
    
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
    float _phTrkSumPtHollowConeDR03[_maxNph];
    float _phEcalRecHitSumEtConeDR03[_maxNph];
    float _phHcalTowerSumEtConeDR03[_maxNph];

    //PFlow isolation
    float _phChargedHadronIsoDR02[_maxNph];
    float _phChargedHadronIsoDR03[_maxNph];
    float _phChargedHadronIsoDR04[_maxNph];
    float _phIsolationSumDR04[_maxNph];
    float _phIsolationSumDR03[_maxNph];
    float _phIsolationSumWrongVtxDR04[_maxNph];
    float _phIsolationSumWrongVtxDR03[_maxNph];
    
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
  _inputTagMet = edm::InputTag("pfMet");
  _inputTagPF = edm::InputTag("particleFlow");
  _inputTagTriggerResults = edm::InputTag("TriggerResults", "", "HLT");
  _inputTagPrimaryVertex = edm::InputTag("offlinePrimaryVerticesWithBS");
  _inputTagRho = edm::InputTag("fixedGridRhoAll");
  _inputTagMCgen = edm::InputTag("genParticles");

  //photon isolator 
  mIsolator.initializePhotonIsolation(kTRUE);
  mIsolator.setRingSize(0.3);

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
    // photons
    _tree->Branch("Nph", &_Nph, "Nph/I"); // number of photons
    _tree->Branch("phE", _phE, "phE[Nph]/F"); // photon E
    _tree->Branch("phPt", _phPt, "phPt[Nph]/F"); // photon pT
    _tree->Branch("phEta", _phEta, "phEta[Nph]/F"); // photon pseudorapidity
    _tree->Branch("phPhi", _phPhi, "phPhi[Nph]/F"); // photon phi
    _tree->Branch("phPt", _phSCx, "phSCx[Nph]/F"); // photon supercluster x
    _tree->Branch("phPt", _phSCy, "phSCy[Nph]/F"); // photon supercluster y
    _tree->Branch("phEta", _phSCz, "phSCz[Nph]/F"); // photon supercluster z
    _tree->Branch("phR9", _phR9, "phR9[Nph]/F"); // photon R9
    _tree->Branch("phTrkSumPtHollowConeDR04", _phTrkSumPtHollowConeDR04, "phTrkSumPtHollowConeDR04[Nph]/F");
    _tree->Branch("phEcalRecHitSumEtConeDR04", _phR9, "phEcalRecHitSumEtConeDR04[Nph]/F");
    _tree->Branch("phHcalTowerSumEtConeDR04", _phHcalTowerSumEtConeDR04, "phHcalTowerSumEtConeDR04[Nph]/F");
    _tree->Branch("phTrkSumPtHollowConeDR03", _phTrkSumPtHollowConeDR03, "phTrkSumPtHollowConeDR03[Nph]/F");
    _tree->Branch("phEcalRecHitSumEtConeDR03", _phEcalRecHitSumEtConeDR03, "phEcalRecHitSumEtConeDR03[Nph]/F");
    _tree->Branch("phHcalTowerSumEtConeDR03", _phHcalTowerSumEtConeDR03, "phHcalTowerSumEtConeDR03[Nph]/F");
    _tree->Branch("phHadronicOverEm", _phHadronicOverEm, "phHadronicOverEm[Nph]/F");
    _tree->Branch("phSigmaIetaIeta", _phSigmaIetaIeta, "phSigmaIetaIeta[Nph]/F");

	//PFlow isolation
    _tree->Branch("phChargedHadronIsoDR04",_phChargedHadronIsoDR04 , "phChargedHadronIsoDR04[Nph]/F");
    _tree->Branch("phChargedHadronIsoDR02", _phChargedHadronIsoDR02, "phChargedHadronIsoDR02[Nph]/F");
    _tree->Branch("phChargedHadronIsoDR03", _phChargedHadronIsoDR03, "phChargedHadronIsoDR03[Nph]/F");
    _tree->Branch("phIsolationSumDR04", _phIsolationSumDR04, "phIsolationSumDR04[Nph]/F");
    _tree->Branch("phIsolationSumDR03", _phIsolationSumDR03, "phIsolationSumDR03[Nph]/F");
    _tree->Branch("phIsolationSumWrongVtxDR04",_phIsolationSumWrongVtxDR04, "phIsolationSumWrongVtxDR04[Nph]/F");
    _tree->Branch("phIsolationSumWrongVtxDR03", _phIsolationSumWrongVtxDR03, "phIsolationSumWrongVtxDR03[Nph]/F");
    
    //electrons in same superCluster (sc)
    _tree->Branch("phNumElectronsSuperCluster",_phNumElectronsSuperCluster,"phNumElectronsSuperCluster[Nph]/I");
    _tree->Branch("elMissingHits", _elMissingHits, "elMissingHits[Nph]/I");
    _tree->Branch("phElectronDR", _phElectronDR, "phElectronDR[Nph]/F");
    _tree->Branch("phHasConversionTracks" , _phHasConversionTracks, "phHasConversionTracks[Nph]/I");

    // Triggers
    _tree->Branch("Triggers", &_triggers, "Triggers/I"); // trigger bits (see trigger names below)
    //define interesting trigger bits (All triggers with same beginning string will be chosen)
    _vecTriggerNames.push_back("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v");
    _vecTriggerNames.push_back("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v");
    _vecTriggerNames.push_back("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v");
    _vecTriggerNames.push_back("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v");
    _vecTriggerNames.push_back("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v");
    _vecTriggerNames.push_back("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v");
    _vecTriggerNames.push_back("HLT_Photon26_Photon18_v");
    _vecTriggerNames.push_back("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v");
    _vecTriggerNames.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
    _vecTriggerNames.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
    _vecTriggerNames.push_back("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
    _vecTriggerNames.push_back("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
    _vecTriggerNames.push_back("HLT_Photon36_R9Id85_Photon22_R9Id85_v");
    _vecTriggerNames.push_back("HLT_Photon36_Photon22_v");
    _vecTriggerNames.push_back("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v");

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
      _phNumElectronsSuperCluster[_Nph] += 1;
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
      //matching electron
      if(_phNumElectronsSuperCluster[_Nph] > 0 && _elMissingHits[_Nph] == 0)
        continue;

      //charged hadron Isolation
      reco::VertexRef myVtx(Vertices, 0);
      mIsolator.setRingSize(0.2);
      mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
      _phChargedHadronIsoDR02[_Nph] = mIsolator.getIsolationCharged();
      //remove pile-up
      if((_phChargedHadronIsoDR02[_Nph] - 0.017 * _rho) > 4)
        continue;  
      
    } 

    //conversion track
    _phHasConversionTracks[_Nph] = it->hasConversionTracks();
    //four momentum
    _phE[_Nph] = it->energy();
    _phPt[_Nph] = it->pt();
    _phEta[_Nph] = it->eta();
    _phPhi[_Nph] = it->phi();

    // supercluster
    _phSCx[_Nph] = it->superCluster()->position().x();
    _phSCy[_Nph] = it->superCluster()->position().y();
    _phSCz[_Nph] = it->superCluster()->position().z();

    // r9 and other vars
    _phR9[_Nph] = it->r9();
    _phTrkSumPtHollowConeDR04[_Nph] = it->trkSumPtHollowConeDR04();
    _phEcalRecHitSumEtConeDR04[_Nph] = it->ecalRecHitSumEtConeDR04();
    _phHcalTowerSumEtConeDR04[_Nph] = it->hcalTowerSumEtConeDR04();
    _phTrkSumPtHollowConeDR03[_Nph] = it->trkSumPtHollowConeDR03();
    _phEcalRecHitSumEtConeDR03[_Nph] = it->ecalRecHitSumEtConeDR03();
    _phHcalTowerSumEtConeDR03[_Nph] = it->hcalTowerSumEtConeDR03();
    _phHadronicOverEm[_Nph] = it->hadronicOverEm();
    _phSigmaIetaIeta[_Nph] = it->sigmaIetaIeta();
    
	//PFlow isolation (not in PhotonCollection)
    //use selected vertex
    reco::VertexRef myVtx(Vertices, 0);
    //Set cone size
    mIsolator.setRingSize(0.4);
    mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
    
	_phChargedHadronIsoDR04[_Nph] = mIsolator.getIsolationCharged();
    double neutralIsoDR04 = mIsolator.getIsolationNeutral();
    double photonIsoDR04 = mIsolator.getIsolationPhoton();
    _phIsolationSumDR04[_Nph] = _phChargedHadronIsoDR04[_Nph] + neutralIsoDR04 + photonIsoDR04;
    
    mIsolator.setRingSize(0.3);
    mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
    double neutralIsoDR03 = mIsolator.getIsolationNeutral();
    double photonIsoDR03 = mIsolator.getIsolationPhoton();
    _phChargedHadronIsoDR03[_Nph] = mIsolator.getIsolationCharged();
    _phIsolationSumDR03[_Nph] = _phChargedHadronIsoDR03[_Nph] + neutralIsoDR03 + photonIsoDR03;

    //worst ChargedHadronIsolation
    _phIsolationSumWrongVtxDR04[_Nph] = _phIsolationSumDR04[_Nph];
    _phIsolationSumWrongVtxDR03[_Nph] = _phIsolationSumDR03[_Nph];
    //loop over all possible vertices and save highest value
    for(int unsigned long i = 0; i < Vertices->size(); i++)
    {
      reco::VertexRef myVtx(Vertices, i);
      //Cone Size of 0.3
      mIsolator.setRingSize(0.3);
      mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
      float curVal = mIsolator.getIsolationCharged() + neutralIsoDR03 + photonIsoDR03;
      if( _phIsolationSumWrongVtxDR03[_Nph] < curVal) 
        _phIsolationSumWrongVtxDR03[_Nph] = curVal;
      //Cone Size of 0.4
      mIsolator.setRingSize(0.4);
      mIsolator.fGetIsolation(&(*it), &(*PF), myVtx, Vertices);
      curVal = mIsolator.getIsolationCharged() + neutralIsoDR04 + photonIsoDR04;
      if (_phIsolationSumWrongVtxDR04[_Nph] < curVal)
        _phIsolationSumWrongVtxDR04[_Nph] = curVal;
    }
    
    
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

// returns vector of integers which are trigger bits needed in the analysis
// Implementation:
// Interesting triggers are added in the beginning of this file
//, all triggers with the same starting string are then stored and printed
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
    for(unsigned int n = 0; n < _vecTriggerNames.size(); n++)
    {
      //if trigger is found in _vecTriggerNames
      if(currentName.find(_vecTriggerNames[n]) != std::string::npos)
      {
        //store position of trigger and name
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
      printf(" :  %d \n", _vecTriggerBits[n][i]);
    }
    printf("\n");
  }
  printf("****************\n");
}

// fill trigger bits and store them in _triggers variable
void Analyzer::SelectTriggerBits(const edm::Handle<edm::TriggerResults>& HLTR)
{
  for(unsigned int i = 0; i < _vecTriggerBits.size(); i++)
  {
    int status = 0;
    for(unsigned int j = 0; j < _vecTriggerBits[i].size(); j++)
    {
      //check if one of the triggers is accepted
      status = status || HLTR->accept(_vecTriggerBits[i][j]);
    }
    // set ith bit of _triggers integer to the current trigger bit status
    _triggers ^= (-status ^ _triggers) & (1 << i);
  }
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
void Analyzer::SelectMCGen(const edm::Handle<reco::GenParticleCollection>& genParticles)
{
  // reset all particle momenta
  FillFourMomentum(NULL, _mcH);
  FillFourMomentum(NULL, _mcPh[0]);
  FillFourMomentum(NULL, _mcPh[1]);

  // Initialize all particle pointers with NULL
  const reco::Candidate* genH = NULL;
  const reco::Candidate* genPh[2] = { NULL, NULL };
  // loop over generated particles
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
  if( (_nevents % 1000) == 0)
  {
    printf("************* NEVENTS = %d K, selected = %d *************\n", _nevents / 1000, _neventsSelected);
  }
  
  // declare event contents
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<reco::VertexCollection> primVertex;
  edm::Handle<reco::PFCandidateCollection> PF;
  edm::Handle<reco::PhotonCollection> photons;
  edm::Handle<reco::GsfElectronCollection> electrons;
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
    //Select Photons by applying soft Cuts
    SelectPhotons(photons,electrons, pv, primVertex, PF);
    // require at least two photons
    if( _Nph >= 2 )
      selRECO = true;
    if(!selRECO && !selGEN)
      return;

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
  // triggers
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
