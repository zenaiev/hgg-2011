// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This class contains all needed variables to read ROOT ntuples for ttbar analysis
// (automaticlly produced by ROOT, then slightly tuned manually)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifndef ZTree_h
#define ZTree_h

#include <TROOT.h>
#include <TChain.h>

// Class which gives access to all information in each event stored in ntuples
class ZTree {
public :
   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   // public members (for direct access outside the class)
   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
   // pointer to the analyzed TTree or TChain
   TTree *fChain;  

   // MC flag (true for MC, false for data)
   bool _flagMC;

   // variable array max sizes
   static const int maxNph = 50; // electrons
   static const int maxNjet = 59; // jets

   // Ntuple variables (their description can be found also in Analyzer/src/Analyzer.cc)
   //[N] means that this is fixed size array with N elements
   Int_t           evRunNumber; // run number
   Int_t           evEventNumber; // event number
   Float_t         rho;
   
   Int_t           Nph;
   Float_t         phPt[maxNph];   //[Nph]
   Float_t         phEta[maxNph];   //[Nph]
   Float_t         phPhi[maxNph];   //[Nph]
   Float_t         phR9[maxNph];   //[Nph]
   Float_t         phTrkSumPtHollowConeDR04[maxNph];   //[Nph]
   Float_t         phEcalRecHitSumEtConeDR04[maxNph];   //[Nph]
   Float_t         phHcalTowerSumEtConeDR04[maxNph];   //[Nph]
   //Float_t         phHcalDepth1TowerSumEtConeDR04[maxNph];   //[Nph]
   //Float_t         phHcalDepth2TowerSumEtConeDR04[maxNph];   //[Nph]
   Float_t         phTrkSumPtHollowConeDR03[maxNph];   //[Nph]
   Float_t         phEcalRecHitSumEtConeDR03[maxNph];   //[Nph]
   Float_t         phHcalTowerSumEtConeDR03[maxNph];   //[Nph]
   //Float_t         phHcalDepth1TowerSumEtConeDR03[maxNph];   //[Nph]
   //Float_t         phHcalDepth2TowerSumEtConeDR03[maxNph];   //[Nph]
   Float_t         phHadronicOverEm[maxNph];   //[Nph]
   Float_t         phChargedHadronIso[maxNph];
   Float_t         phPhotonIsoWrongVtx[maxNph];
   Float_t         phNeutralHadronIso[maxNph];
   Float_t         phPhotonIso[maxNph];
   Int_t           phNumElectronsSuperCluster[maxNph];
   Int_t           elMissingHits[maxNph];
   Float_t         phElectronDR[maxNph];
   Float_t         phSigmaIetaIeta[maxNph];   //[Nph]
   Float_t         phMatch[maxNph];   //[Nph]
   Int_t           Njet;
   Float_t         jetPt[maxNjet];   //[Njet]
   Float_t         jetEta[maxNjet];   //[Njet]
   Float_t         jetPhi[maxNjet];   //[Njet]
   Float_t         jetMass[maxNjet];   //[Njet]
   Float_t         jetMuEn[maxNjet];   //[Njet]
   Float_t         jetElEn[maxNjet];   //[Njet]
   Int_t           Triggers;
   Int_t           Npv;
   Int_t           pvNDOF;
   Float_t         pvZ;
   Float_t         pvRho;
   
   // variables for MC only
   Int_t           mcEventType; // type of event: 1 ttbar decay into ee, 2 ttbar decay into mumu, 3 ttbar decay into emu, 0 anything else
   float mcH[4];    // top quark four momentum
   float mcPh[2][4]; // antitop quark four momentum

   // List of branches (their names follow variable names with prefix b_)
   TBranch        *b_evRunNumber;   //!
   TBranch        *b_evEventNumber;   //!
   TBranch        *b_rho; //!
   TBranch        *b_Nph;   //!
   TBranch        *b_phPt;   //!
   TBranch        *b_phEta;   //!
   TBranch        *b_phPhi;   //!
   TBranch        *b_phR9;   //!
   TBranch        *b_phTrkSumPtHollowConeDR04;   //!
   TBranch        *b_phEcalRecHitSumEtConeDR04;   //!
   TBranch        *b_phHcalTowerSumEtConeDR04;   //!
   //TBranch        *b_phHcalDepth1TowerSumEtConeDR04;   //!
   //TBranch        *b_phHcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_phTrkSumPtHollowConeDR03;   //!
   TBranch        *b_phEcalRecHitSumEtConeDR03;   //!
   TBranch        *b_phHcalTowerSumEtConeDR03;   //!
   //TBranch        *b_phHcalDepth1TowerSumEtConeDR03;   //!
   //TBranch        *b_phHcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_phHadronicOverEm;
   TBranch        *b_phChargedHadronIso;
   TBranch        *b_phPhotonIsoWrongVtx;
   TBranch        *b_phNeutralHadronIso;
   TBranch        *b_phPhotonIso;
   TBranch        *b_phNumElectronsSuperCluster;
   TBranch        *b_elMissingHits;
   TBranch        *b_phElectronDR;
   TBranch        *b_phSigmaIetaIeta;
   TBranch        *b_phMatch;
   TBranch        *b_Njet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetMuEn;   //!
   TBranch        *b_jetElEn;   //!
   TBranch        *b_Triggers;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_pvNDOF;   //!
   TBranch        *b_pvZ;   //!
   TBranch        *b_pvRho;   //!
   // for MC only
   TBranch        *b_mcEventType;   //!
   TBranch        *b_mcH;   //!
   TBranch        *b_mcPh;   //!

   // constructor
   // argument: true for MC, false (default) for data
   ZTree(bool flagMC = false) : fChain(0), _flagMC(flagMC) { }
   
   // destructor
   virtual ~ZTree() { }
   
   // initialise with provided tree pointer
   virtual void    Init(TTree *tree);
};

// initialise with provided tree pointer
void ZTree::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evRunNumber", &evRunNumber, &b_evRunNumber);
   fChain->SetBranchAddress("evEventNumber", &evEventNumber, &b_evEventNumber);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("Nph", &Nph, &b_Nph);
   fChain->SetBranchAddress("phPt", phPt, &b_phPt);
   fChain->SetBranchAddress("phEta", phEta, &b_phEta);
   fChain->SetBranchAddress("phPhi", phPhi, &b_phPhi);
   fChain->SetBranchAddress("phR9", phR9, &b_phR9);
   fChain->SetBranchAddress("phTrkSumPtHollowConeDR04", phTrkSumPtHollowConeDR04, &b_phTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("phEcalRecHitSumEtConeDR04", phEcalRecHitSumEtConeDR04, &b_phEcalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("phHcalTowerSumEtConeDR04", phHcalTowerSumEtConeDR04, &b_phHcalTowerSumEtConeDR04);
   //fChain->SetBranchAddress("phHcalDepth1TowerSumEtConeDR04", phHcalDepth1TowerSumEtConeDR04, &b_phHcalDepth1TowerSumEtConeDR04);
   //fChain->SetBranchAddress("phHcalDepth2TowerSumEtConeDR04", phHcalDepth2TowerSumEtConeDR04, &b_phHcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("phTrkSumPtHollowConeDR03", phTrkSumPtHollowConeDR03, &b_phTrkSumPtHollowConeDR03);
   fChain->SetBranchAddress("phEcalRecHitSumEtConeDR03", phEcalRecHitSumEtConeDR03, &b_phEcalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("phHcalTowerSumEtConeDR03", phHcalTowerSumEtConeDR03, &b_phHcalTowerSumEtConeDR03);
   //fChain->SetBranchAddress("phHcalDepth1TowerSumEtConeDR03", phHcalDepth1TowerSumEtConeDR03, &b_phHcalDepth1TowerSumEtConeDR03);
   //fChain->SetBranchAddress("phHcalDepth2TowerSumEtConeDR03", phHcalDepth2TowerSumEtConeDR03, &b_phHcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("phHadronicOverEm", phHadronicOverEm, &b_phHadronicOverEm);
   fChain->SetBranchAddress("phChargedHadronIso", phChargedHadronIso, &b_phChargedHadronIso);
   fChain->SetBranchAddress("phPhotonIsoWrongVtx", phPhotonIsoWrongVtx, &b_phPhotonIsoWrongVtx);
   fChain->SetBranchAddress("phNeutralHadronIso", phNeutralHadronIso, &b_phNeutralHadronIso);
   fChain->SetBranchAddress("phPhotonIso", phPhotonIso , &b_phPhotonIso);
   fChain->SetBranchAddress("phNumElectronsSuperCluster", phNumElectronsSuperCluster, &b_phNumElectronsSuperCluster);
   fChain->SetBranchAddress("elMissingHits", elMissingHits, &b_elMissingHits);
   fChain->SetBranchAddress("phElectronDR", phElectronDR , &b_phElectronDR);
   fChain->SetBranchAddress("phSigmaIetaIeta", phSigmaIetaIeta, &b_phSigmaIetaIeta);
   fChain->SetBranchAddress("Njet", &Njet, &b_Njet);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMass", jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetMuEn", jetMuEn, &b_jetMuEn);
   fChain->SetBranchAddress("jetElEn", jetElEn, &b_jetElEn);
   fChain->SetBranchAddress("Triggers", &Triggers, &b_Triggers);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("pvNDOF", &pvNDOF, &b_pvNDOF);
   fChain->SetBranchAddress("pvZ", &pvZ, &b_pvZ);
   fChain->SetBranchAddress("pvRho", &pvRho, &b_pvRho);
   // MC
   if(_flagMC)
   {
     fChain->SetBranchAddress("mcEventType", &mcEventType, &b_mcEventType);
     fChain->SetBranchAddress("mcH", mcH, &b_mcH);
     fChain->SetBranchAddress("mcPh", mcPh, &b_mcPh);
     fChain->SetBranchAddress("phMatch", phMatch, &b_phMatch);
   }
}

#endif // #ifdef ZTree_h
