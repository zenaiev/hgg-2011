// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>> Helper for ttbar event reconstruction >>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#ifndef TTBAR_EVENTRECO_H
#define TTBAR_EVENTRECO_H

// additional files from this analysis 
#include "tree.h"
#include "kinReco.h"
#include "selection.h"
#include "settings.h"
// C++ library or ROOT header files
#include <map>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>> ZVarHisto class >>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Class for control plot and cross section histograms: 
// it stores variable name and histogram. 
// A bunch of histograms can be filled using proper ttbar kinematics 
// input with just one line (see void FillHistos() below).
// Also see void StoreHistos() for histogram storage.
//
class ZVarHisto
{
  private:
    TH1D* zHisto; // histogram
    TString zVar; // variable name (see void FillHistos() fot its usage)
    int zEventClass; // 1 all, 2 dijet, ...

  
  public:
    // constructor
    ZVarHisto(const TString& str, TH1D* h, int eventClass = 1)
    {
      zHisto = h;
      zVar = str;
      zEventClass = eventClass;
    }

    // copy constructor
    ZVarHisto(const ZVarHisto& old)
    {
      zHisto = new TH1D(*(old.zHisto));
      zVar = old.zVar;
      zEventClass = old.zEventClass;
    }

    // access histogram
    TH1D*& H() {return zHisto;}
    
    // access variable name
    TString& V() {return zVar;}

    // access event class
    int& EventClass() {return zEventClass;}
};
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>> Fill histogram from ZVarHisto class >>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Arguments:
//   std::vector<ZVarHisto>& VecVarHisto: vector of objects to be filled
//   double w: weight
//   TLorentzVector* t: top quark momentum
//   TLorentzVector* tbar: antitop quark momentum
//   TLorentzVector* vecLepM: leptoni momentum (if needed, can be omitted)
//   TLorentzVector* vecLepP: lepton+ momentum (if needed, can be omitted)
//
void FillHistos(std::vector<ZVarHisto>& VecVarHisto, double w, TLorentzVector* ph1, TLorentzVector* ph2, int eventClass = 1)
{
  // momentum of higgs
  TLorentzVector higgs = *ph1 + *ph2;
  // loop over provided histograms to be filled
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    // check event class
    if(VecVarHisto[h].EventClass() != 1 && VecVarHisto[h].EventClass() != eventClass)
      continue;
    // retrieve variable name
    TString var = VecVarHisto[h].V();
    // retrieve histogram
    TH1* histo = VecVarHisto[h].H();
    // now fill histograms depending on the variable:
    // top pT
    //ph1->Print("all");
    //ph2->Print("all");
    //printf("mgg = %f\n", higgs.M());
    if(var == "mgg")
      histo->Fill(higgs.M(), w);
    // uknown (not implemented) variable
    // (you can implement more variables here if needed)
    else
    {
      //printf("Error: unknown variable %s\n", var.Data());
      //exit(1);
      continue;
    } // end variable for histo
  } // end loop over histos
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>> Store hostogram from ZVarHisto class >>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Store bunch of histogram (argument std::vector<ZVarHisto>& VecVarHisto)
// (stores a copy of histogram)
//
void StoreHistos(std::vector<ZVarHisto>& VecVarHisto)
{
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    TH1* histo = VecVarHisto[h].H();
    TString name = histo->GetName() + std::to_string(VecVarHisto[h].EventClass());
    TString title = histo->GetTitle() + TString("class = ") + std::to_string(VecVarHisto[h].EventClass());
    histo->SetNameTitle(name, title);
    histo->Write();
  }
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>> Input parameters for eventreco routine (see below)  >>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class ZEventRecoInput
{
  public:
    TString Name; // name pattern (to be used in output histograms)
    std::vector<ZVarHisto> VecVarHisto; // container with needed histograms
    int Type; // 1 data, 2 MC signal, 3 MC ttbar other, 4 MC background
    bool Gen; // if true, the histogram is filled at true level
    std::vector<TString> VecInFile; // container with input files
    double Weight; // weight for histogram filling
    long MaxNEvents; // maximum number of processed events
    
    // contstructor
    ZEventRecoInput()
    {
      // set default values
      Weight = 1.0;
      MaxNEvents = 100e10;
      Gen = false;
    }
    
    // add one more input file (str) to the chain
    void AddToChain(const TString& str)
    {
      VecInFile.push_back(str);
    }

    // erase all input files form the chain
    void ClearChain()
    {
      VecInFile.clear();
    }
};

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>> Basic routine for ttbar event reconstruction >>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void eventreco(ZEventRecoInput in)
{ 
  printf("****** EVENTRECO ******\n");
  printf("input sample: %s type: %d\n", in.Name.Data(), in.Type);
  
  // steering
  const bool flagReqPV = 0;
  const bool flagReqTrig = 0;
  // directory for output ROOT files with histograms
  TString outDir = gHistDir; 

  // this flag determines whether generator level information is available
  // (should be available for signal MC)
  bool flagMC = (in.Type == 2 || in.Type == 3);
  
  // output file
  TFile* fout = TFile::Open(TString::Format("%s/%s.root", outDir.Data(), in.Name.Data()), "recreate");
  
  // input tree
  TChain* chain = new TChain("tree");
  for(int f = 0; f < in.VecInFile.size(); f++)
    chain->Add(in.VecInFile[f]);
  ZTree* preselTree = new ZTree(flagMC);
  preselTree->Init(chain);

  // process generator level, if needed
  if(in.Gen)
  {
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("mcEventType", 1);
    chain->SetBranchStatus("mcH", 1);
    chain->SetBranchStatus("mcPh", 1);
  }
    
  // event counters
  long nSel = 0;
  int nGen = 0;
  
  // determine number of events
  long nEvents = chain->GetEntries();
  //limit it if exceeds the specified maximum number
  if(nEvents > in.MaxNEvents)
    nEvents = in.MaxNEvents;
  printf("nEvents: %ld\n", nEvents);
  // event loop
  for(int e = 0; e < nEvents; e++)
  {
    chain->GetEntry(e);
    if(flagMC)
    {
      // skip background events for MC signal
      if(in.Type == 2 && preselTree->mcEventType != 1) continue;
      // skip signal events for MC 'ttbar other' (background)
      if(in.Type == 3 && preselTree->mcEventType == 1) continue;
    }
    // process generator level if needed
    if(in.Gen)
    {
      // prepare four vectors for top and antitop
      //TLorentzVector h, ph[2];
      //h.SetXYZM(preselTree->mcH[0], preselTree->mcH[1], preselTree->mcH[2], preselTree->mcH[3]);
      TLorentzVector ph[2];
      for(int p = 0; p < 2; p++)
        ph[p].SetXYZM(preselTree->mcPh[p][0], preselTree->mcPh[p][1], preselTree->mcPh[p][2], preselTree->mcPh[p][3]);
      //ph[0].Print("all");
      //ph[1].Print("all");
      // fill histos
      double w = in.Weight;
      FillHistos(in.VecVarHisto, w, &ph[0], &ph[1]);
      nGen++;
      continue;
    }
    if(in.Type > 1)
      nGen++;
    
    // process reco level if needed
    //
    // primary vertex selection
    if(flagReqPV)
      if(preselTree->Npv < 1 || preselTree->pvNDOF < 4 || preselTree->pvRho > 2.0 || TMath::Abs(preselTree->pvZ) > 24.0)
        continue;
    // select dilepton pair
    TLorentzVector ph[2];
    bool trig = false;

    // *****************************************
    // ***************** all *******************
    // *****************************************
    // trigger bits (accept the event if at least one needed trigger bit is fired)
    if(flagReqTrig)
    {
      for(int bit = 0; bit < 3; bit++)
        if((preselTree->Triggers >> bit) & 1)
        {
          trig = true;
          break;
        }
    }
    else
      trig = true;

    // call selection routine (see selection.h for description)
    if(!trig)
      continue;

    int eventClass = SelectHgg(preselTree, 1, ph[0], ph[1]);
    if(eventClass == 0)
      continue;

    // event selection done: increment the counter of selected events
    nSel++;
      
    // fill histograms
    double w = in.Weight;
    double mgg = (ph[0] + ph[1]).M();
    //ph[0].Print("all");
    //ph[1].Print("all");
    FillHistos(in.VecVarHisto, w, &ph[0], &ph[1], eventClass);

    /*for(int p = 0; p < 2; p++)
      ph[p].SetXYZM(preselTree->mcPh[p][0], preselTree->mcPh[p][1], preselTree->mcPh[p][2], preselTree->mcPh[p][3]);
    printf("rec gen --- %.1f %.1f\n", mgg, (ph[0] + ph[1]).M());
    ph[0].Print("all");
    ph[1].Print("all");*/
  } // end event loop
  
  // print the numbers of selected events and events with successfull kinematic reconstruction
  printf("nSel  : %ld\n", nSel);
  // for signal MC, print the number of signal events at generator level and detector efficiency
  // (with and without kinematic reconstruction)
  if(in.Type == 2) 
  {
    printf("nGen  : %ld\n", nGen);
    printf("C = %.2f%%\n", 100. * nSel / nGen);
  }

  // store histograms, close output file
  fout->cd();
  StoreHistos(in.VecVarHisto);
  fout->Close();
}

#endif
