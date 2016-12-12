//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug  4 23:45:53 2016 by ROOT version 6.04/06
// from TTree tree/ttbar
// found on file: ttbarSel_emu10p.root
//////////////////////////////////////////////////////////

#ifndef ZTree_h
#define ZTree_h

#include <TROOT.h>
#include <TChain.h>
//#include <TFile.h>
//#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

class ZTree {//: public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // MC flag
   bool _flagMC;

   // variable array max sizes
   static const int maxNel = 10;
   static const int maxNmu = 10;
   static const int maxNjet = 25;

   // Declaration of leaf types
   Int_t           evRunNumber;
   Int_t           evLumiBlock;
   Int_t           evEventNumber;
   Int_t           Nmu;
   Float_t         muPt[maxNmu];   //[Nmu]
   Float_t         muEta[maxNmu];   //[Nmu]
   Float_t         muPhi[maxNmu];   //[Nmu]
   Float_t         muIso03[maxNmu];   //[Nmu]
   Float_t         muIso04[maxNmu];   //[Nmu]
   Int_t           muHitsValid[maxNmu];   //[Nmu]
   Int_t           muHitsPixel[maxNmu];   //[Nmu]
   Float_t         muDistPV0[maxNmu];   //[Nmu]
   Float_t         muDistPVz[maxNmu];   //[Nmu]
   Float_t         muTrackChi2NDOF[maxNmu];   //[Nmu]
   Int_t           Nel;
   Float_t         elPt[maxNel];   //[Nel]
   Float_t         elEta[maxNel];   //[Nel]
   Float_t         elPhi[maxNel];   //[Nel]
   /*Float_t         elIso03TkSumPt[4];   //[Nel]
   Float_t         elIso03EcalRecHitSumEt[4];   //[Nel]
   Float_t         elIso03HcalTowerSumEt[4];   //[Nel]
   Float_t         elIso04TkSumPt[4];   //[Nel]
   Float_t         elIso04EcalRecHitSumEt[4];   //[Nel]
   Float_t         elIso04HcalTowerSumEt[4];   //[Nel]*/
   Float_t         elIso03[maxNel];   //[Nel]
   Float_t         elIso04[maxNel];   //[Nel]
   Int_t           elConvFlag[maxNel];   //[Nel]
   Float_t         elConvDist[maxNel];   //[Nel]
   Float_t         elConvDcot[maxNel];   //[Nel]
   Float_t         elMissHits[maxNel];   //[Nel]
   Int_t           Njet;
   Float_t         jetPt[maxNjet];   //[Njet]
   Float_t         jetEta[maxNjet];   //[Njet]
   Float_t         jetPhi[maxNjet];   //[Njet]
   Float_t         jetMass[maxNjet];   //[Njet]
   Float_t         jetMuEn[maxNjet];   //[Njet]
   Float_t         jetElEn[maxNjet];   //[Njet]
   Float_t         jetBTagDiscr[maxNjet];   //[Njet]
   Float_t         jetBTagMatchDiff1[maxNjet];   //[Njet]
   Float_t         jetBTagMatchDiff2[maxNjet];   //[Njet]
   Float_t         metPx;
   Float_t         metPy;
   Int_t           Npv;
   Int_t           pvNDOF;
   Float_t         pvZ;
   Float_t         pvRho;
   Int_t           Triggers;
   // MC
   Int_t           mcEventType;
   float mcT[4];
   float mcTbar[4];

   // List of branches
   TBranch        *b_evRunNumber;   //!
   TBranch        *b_evLumiBlock;   //!
   TBranch        *b_evEventNumber;   //!
   TBranch        *b_Nmu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muIso03;   //!
   TBranch        *b_muIso04;   //!
   TBranch        *b_muHitsValid;   //!
   TBranch        *b_muHitsPixel;   //!
   TBranch        *b_muDistPV0;   //!
   TBranch        *b_muDistPVz;   //!
   TBranch        *b_muTrackChi2NDOF;   //!
   TBranch        *b_Nel;   //!
   TBranch        *b_elPt;   //!
   TBranch        *b_elEta;   //!
   TBranch        *b_elPhi;   //!
   /*TBranch        *b_elIso03TkSumPt;   //!
   TBranch        *b_elIso03EcalRecHitSumEt;   //!
   TBranch        *b_elIso03HcalTowerSumEt;   //!
   TBranch        *b_elIso04TkSumPt;   //!
   TBranch        *b_elIso04EcalRecHitSumEt;   //!
   TBranch        *b_elIso04HcalTowerSumEt;   //!*/
   TBranch        *b_elIso03;   //!
   TBranch        *b_elIso04;   //!
   TBranch        *b_elConvFlag;   //!
   TBranch        *b_elConvDist;   //!
   TBranch        *b_elConvDcot;   //!
   TBranch        *b_elMissHits;   //!
   TBranch        *b_Njet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetMuEn;   //!
   TBranch        *b_jetElEn;   //!
   TBranch        *b_jetBTagDiscr;   //!
   TBranch        *b_jetBTagMatchDiff1;   //!
   TBranch        *b_jetBTagMatchDiff2;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_pvNDOF;   //!
   TBranch        *b_pvZ;   //!
   TBranch        *b_pvRho;   //!
   TBranch        *b_Triggers;   //!
   // MC
   TBranch        *b_mcEventType; //!
   TBranch        *b_mcT; //!
   TBranch        *b_mcTbar; //!

   ZTree(bool flagMC = false) : fChain(0), _flagMC(flagMC) { }
   virtual ~ZTree() { }
   //virtual Int_t   Version() const { return 2; }
   //virtual void    Begin(TTree *tree);
   //virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   /*virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();*/

   //ClassDef(ZTree,0);
};

//#endif

//#ifdef ZTree_cxx
void ZTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evRunNumber", &evRunNumber, &b_evRunNumber);
   //fChain->SetBranchAddress("evLumiBlock", &evLumiBlock, &b_evLumiBlock);
   fChain->SetBranchAddress("evEventNumber", &evEventNumber, &b_evEventNumber);
   fChain->SetBranchAddress("Nmu", &Nmu, &b_Nmu);
   fChain->SetBranchAddress("muPt", muPt, &b_muPt);
   fChain->SetBranchAddress("muEta", muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
   fChain->SetBranchAddress("muIso03", muIso03, &b_muIso03);
   fChain->SetBranchAddress("muIso04", muIso04, &b_muIso04);
   fChain->SetBranchAddress("muHitsValid", muHitsValid, &b_muHitsValid);
   fChain->SetBranchAddress("muHitsPixel", muHitsPixel, &b_muHitsPixel);
   fChain->SetBranchAddress("muDistPV0", muDistPV0, &b_muDistPV0);
   fChain->SetBranchAddress("muDistPVz", muDistPVz, &b_muDistPVz);
   fChain->SetBranchAddress("muTrackChi2NDOF", muTrackChi2NDOF, &b_muTrackChi2NDOF);
   fChain->SetBranchAddress("Nel", &Nel, &b_Nel);
   fChain->SetBranchAddress("elPt", elPt, &b_elPt);
   fChain->SetBranchAddress("elEta", elEta, &b_elEta);
   fChain->SetBranchAddress("elPhi", elPhi, &b_elPhi);
   /*fChain->SetBranchAddress("elIso03TkSumPt", elIso03TkSumPt, &b_elIso03TkSumPt);
   fChain->SetBranchAddress("elIso03EcalRecHitSumEt", elIso03EcalRecHitSumEt, &b_elIso03EcalRecHitSumEt);
   fChain->SetBranchAddress("elIso03HcalTowerSumEt", elIso03HcalTowerSumEt, &b_elIso03HcalTowerSumEt);
   fChain->SetBranchAddress("elIso04TkSumPt", elIso04TkSumPt, &b_elIso04TkSumPt);
   fChain->SetBranchAddress("elIso04EcalRecHitSumEt", elIso04EcalRecHitSumEt, &b_elIso04EcalRecHitSumEt);
   fChain->SetBranchAddress("elIso04HcalTowerSumEt", elIso04HcalTowerSumEt, &b_elIso04HcalTowerSumEt);*/
   fChain->SetBranchAddress("elIso03", elIso03, &b_elIso03);
   fChain->SetBranchAddress("elIso04", elIso04, &b_elIso04);
   fChain->SetBranchAddress("elConvFlag", elConvFlag, &b_elConvFlag);
   fChain->SetBranchAddress("elConvDist", elConvDist, &b_elConvDist);
   fChain->SetBranchAddress("elConvDcot", elConvDcot, &b_elConvDcot);
   fChain->SetBranchAddress("elMissHits", elMissHits, &b_elMissHits);
   fChain->SetBranchAddress("Njet", &Njet, &b_Njet);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMass", jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetMuEn", jetMuEn, &b_jetMuEn);
   fChain->SetBranchAddress("jetElEn", jetElEn, &b_jetElEn);
   fChain->SetBranchAddress("jetBTagDiscr", jetBTagDiscr, &b_jetBTagDiscr);
   fChain->SetBranchAddress("jetBTagMatchDiff1", jetBTagMatchDiff1, &b_jetBTagMatchDiff1);
   fChain->SetBranchAddress("jetBTagMatchDiff2", jetBTagMatchDiff2, &b_jetBTagMatchDiff2);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("pvNDOF", &pvNDOF, &b_pvNDOF);
   fChain->SetBranchAddress("pvZ", &pvZ, &b_pvZ);
   fChain->SetBranchAddress("pvRho", &pvRho, &b_pvRho);
   fChain->SetBranchAddress("Triggers", &Triggers, &b_Triggers);
   // MC
   if(_flagMC) fChain->SetBranchAddress("mcEventType", &mcEventType, &b_mcEventType);
   if(_flagMC) fChain->SetBranchAddress("mcT", mcT, &b_mcT);
   if(_flagMC) fChain->SetBranchAddress("mcTbar", mcTbar, &b_mcTbar);
}

/*Bool_t ZTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}*/

#endif // #ifdef ZTree_h
