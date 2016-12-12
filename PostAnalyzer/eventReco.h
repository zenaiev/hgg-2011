#ifndef TTBAR_EVENTRECO_H
#define TTBAR_EVENTRECO_H

#include "tree.h"
#include "kinReco.h"
#include "selection.h"
#include "settings.h"
#include <map>
#include <TChain.h>
//#include <mystaff.cxx>
#include <TCanvas.h>
#include <TFile.h>


//const double massEl = 0.000511;
//const double massMu = 0.105658;

class ZVarHisto
{
  private:
    TH1D* zHisto;
    TString zVar;
  
  public:
    ZVarHisto(const TString& str, TH1D* h)
    {
      zHisto = h;
      zVar = str;
    }

    ZVarHisto(const ZVarHisto& old)
    {
      zHisto = new TH1D(*(old.zHisto));
      zVar = old.zVar;
    }

    TH1* H() {return zHisto;}
    TString V() {return zVar;}
};

void FillHistos(std::vector<ZVarHisto>& VecVarHisto, double w, TLorentzVector* t, TLorentzVector* tbar, TLorentzVector* vecLepM = NULL, TLorentzVector* vecLepP = NULL)
{
  TLorentzVector vecttbar = *t + *tbar;
  TLorentzVector* ttbar = &vecttbar;
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    TString var = VecVarHisto[h].V();
    TH1* histo = VecVarHisto[h].H();
    if(var == "ptt") 
      histo->Fill(t->Pt(), w);
    else if(var == "ptat") 
      histo->Fill(tbar->Pt(), w);
    else if(var == "pttat") 
    {
      histo->Fill(t->Pt(), w);
      histo->Fill(tbar->Pt(), w);
    }
    else if(var == "pttt") 
      histo->Fill(ttbar->Pt(), w);
    else if(var == "yt") 
      histo->Fill(t->Rapidity(), w);
    else if(var == "yat") 
      histo->Fill(tbar->Rapidity(), w);
    else if(var == "ytat") 
    {
      histo->Fill(t->Rapidity(), w);
      histo->Fill(tbar->Rapidity(), w);
    }
    else if(var == "ytt") 
      histo->Fill(ttbar->Rapidity(), w);
    else if(var == "mtt") 
      histo->Fill(ttbar->M(), w);
    else if(var == "ptl") 
    {
      histo->Fill(vecLepM->Pt(), w);
      histo->Fill(vecLepP->Pt(), w);
    }
    else
    {
      //printf("Error: unknown variable %s\n", var.Data());
      //exit(1);
      continue;
    } // end variable for histo
  } // end loop over histos
}

void StoreHistos(std::vector<ZVarHisto>& VecVarHisto)
{
  for(int h = 0; h < VecVarHisto.size(); h++)
  {
    //printf("h = %d\n", h);
    TH1* histo = VecVarHisto[h].H();
    //TString name = TString::Format("%s_t%d_ch%d+%s\n", histo->GetName(), in.Type, in.Channel, in.Name.Data());
    //TString title = TString::Format("%s-t%d-ch%d-%s\n", histo->GetTitle(), in.Type, in.Channel, in.Name.Data());
    TString name = histo->GetName();
    TString title = histo->GetTitle();
    histo->SetNameTitle(name, title);
    histo->Write();
  }
}

class ZEventRecoInput
{
  public:
    TString Name;
    std::vector<ZVarHisto> VecVarHisto;
    int Channel; // 1 ee, 2 mumu, 3 emu
    int Type; // 1 data, 2 MC signal, 3 MC ttbar other, 4 MC background
    bool Gen;
    //TChain* Chain;
    std::vector<TString> VecInFile;
    //TString NameFout;
    double Weight;
    long MaxNEvents;
    
    ZEventRecoInput()
    {
      Weight = 1.0;
      MaxNEvents = 100e10;
      Gen = false;
    }
    
    void AddToChain(const TString& str)
    {
      VecInFile.push_back(str);
    }

    void ClearChain()
    {
      VecInFile.clear();
    }
};

void eventreco(ZEventRecoInput in)
{ 
  printf("****** EVENTRECO ******\n");
  printf("input sample: %s\n", in.Name.Data());
  printf("type: %d   channel: %d\n", in.Type, in.Channel);
  
  // some steerings
  const double bTagDiscrL = 0.244;
  TString outDir = gHistDir;

  // some flags
  bool flagMC = (in.Type == 2 || in.Type == 3);
  
  // output file
  TFile* fout = TFile::Open(TString::Format("%s/%s-c%d.root", outDir.Data(), in.Name.Data(), in.Channel), "recreate");
  
  // input tree
  TChain* chain = new TChain("tree");
  for(int f = 0; f < in.VecInFile.size(); f++)
    chain->Add(in.VecInFile[f]);
  ZTree* preselTree = new ZTree(flagMC);
  preselTree->Init(chain);
  // gen level
  if(in.Gen)
  {
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("mcEventType", 1);
    chain->SetBranchStatus("mcT", 1);
    chain->SetBranchStatus("mcTbar", 1);
  }
    
  // event count
  long nSel = 0;
  long nReco = 0;
  int nGen = 0;
  
  // special histos
  TH1D* hInacc = new TH1D("hInacc", "KinReco inaccuracy", 1000, 0.0, 100.0);
  TH1D* hAmbig = new TH1D("hAmbig", "KinReco ambiguity", 100, 0.0, 100.0);

  // event loop
  long nEvents = chain->GetEntries();
  if(nEvents > in.MaxNEvents)
    nEvents = in.MaxNEvents;
  printf("nEvents: %ld\n", nEvents);
  for(int e = 0; e < nEvents; e++)
  {
    chain->GetEntry(e);
    if(flagMC)
    {
      //preselTree->b_mcEventType->GetEntry(e);
      //printf("%d", preselTree->mcEventType);
      // MC signal?
      if(in.Type == 2 && preselTree->mcEventType != in.Channel) continue;
      // MC ttbar other?
      if(in.Type == 3 && preselTree->mcEventType == in.Channel) continue;
      //if(in.Type == 3 && preselTree->mcEventType != 0) continue;
    }
    // gen level
    if(in.Gen)
    {
      //preselTree->b_mcT->GetEntry(e);
      //preselTree->b_mcTbar->GetEntry(e);
      TLorentzVector t, tbar;
      t.SetXYZM(preselTree->mcT[0], preselTree->mcT[1], preselTree->mcT[2], preselTree->mcT[3]);
      tbar.SetXYZM(preselTree->mcTbar[0], preselTree->mcTbar[1], preselTree->mcTbar[2], preselTree->mcTbar[3]);
      // fill histos
      double w = in.Weight;
      FillHistos(in.VecVarHisto, w, &t, &tbar);
      nGen++;
      continue;
    }
    
    // continue reco level
    // load event
    //chain->GetEntry(e);
    if(in.Type > 1)
      nGen++;
    
    // primary vertex selection
    if(preselTree->Npv < 1 || preselTree->pvNDOF < 4 || preselTree->pvRho > 2.0 || TMath::Abs(preselTree->pvZ) > 24.0)
      continue;
    // primary dataset
    TString inFile = chain->GetCurrentFile()->GetName();
    // select dilepton pair
    TLorentzVector vecLepM, vecLepP;
    double maxPtDiLep = -1.0;
    bool trig = false;
    // *****************************************
    // ***************** emu *******************
    // *****************************************
    if(in.Channel == 3)
    {
      // trigger: 12th to 17th bits
      for(int bit = 12; bit < 17; bit++)
        if((preselTree->Triggers >> bit) & 1)
        {
          trig = true;
          break;
        }
      //trig = true;
      if(trig)
        SelectDilepEMu(preselTree, vecLepM, vecLepP, maxPtDiLep);
    }
    // *****************************************
    // ***************** ee ********************
    // *****************************************
    if(in.Channel == 1)
    {
      // trigger: 6th to 11th bits
      for(int bit = 6; bit < 11; bit++)
        if((preselTree->Triggers >> bit) & 1)
        {
          trig = true;
          break;
        }
      //trig = true;
      double met = TMath::Sqrt(TMath::Power(preselTree->metPx, 2.0) + TMath::Power(preselTree->metPy, 2.0));
      if(trig && met > 30.0)
        SelectDilepEE(preselTree, vecLepM, vecLepP, maxPtDiLep);
    }
    // *****************************************
    // **************** mumu *******************
    // *****************************************
    if(in.Channel == 2)
    {
      // trigger: 0th to 5th bits
      for(int bit = 0; bit < 5; bit++)
        if((preselTree->Triggers >> bit) & 1)
        {
          trig = true;
          break;
        }
      //trig = true;
      double met = TMath::Sqrt(TMath::Power(preselTree->metPx, 2.0) + TMath::Power(preselTree->metPy, 2.0));
      if(trig && met > 30.0)
        SelectDilepMuMu(preselTree, vecLepM, vecLepP, maxPtDiLep);
    }
    // check if there is dilepton pair
    if(maxPtDiLep < 0.0)
      continue;
    // dilepton pair found, select jets: fill array
    std::vector<TLorentzVector> vecJets;
    bool oneBTagJet = false;
    for(int j = 0; j < preselTree->Njet; j++)
    {
      if(TMath::Abs(preselTree->jetEta[j]) > 2.4)
        continue;
      TLorentzVector vecJet;
      vecJet.SetPtEtaPhiM(preselTree->jetPt[j], preselTree->jetEta[j], preselTree->jetPhi[j], preselTree->jetMass[j]);
      double corrE = vecJet.E() - preselTree->jetMuEn[j] - preselTree->jetElEn[j];
      //double corrE = vecJet.E();
      double corrPt = preselTree->jetPt[j] * corrE / vecJet.E();
      if(corrPt < 30.0)
        continue;
      TLorentzVector corrVec;
      corrVec.SetPtEtaPhiE(corrPt, preselTree->jetEta[j], preselTree->jetPhi[j], corrE);
      // b-tagging
      if(preselTree->jetBTagDiscr[j] > bTagDiscrL)
      {
        corrVec.SetPtEtaPhiM(corrVec.Pt(), corrVec.Eta(), corrVec.Phi(), -1 * corrVec.M());
        oneBTagJet = true;
      }
      vecJets.push_back(corrVec);
    }
    if(vecJets.size() < 2)
      continue;
    // require b-tagged jet
    if(!oneBTagJet)
      continue;
    // jet array filled
    
    // kinematic reconstruction
    nSel++;
    TLorentzVector t, tbar;
    int status = KinRecoDilepton(vecLepM, vecLepP, vecJets, preselTree->metPx, preselTree->metPy, t, tbar, hInacc, hAmbig);
    //printf("STATUS: %d\n", status);
    if(status > 0) // successfull kinreco
    {
      //printf("top:      (%8.3f  %8.3f  %8.3f  %8.3f)\n", t.X(), t.Y(), t.Z(), t.M());
      //printf("antitop:  (%8.3f  %8.3f  %8.3f  %8.3f)\n", tbar.X(), tbar.Y(), tbar.Z(), tbar.M());
      nReco++;
      
      // fill histos
      double w = in.Weight;
      FillHistos(in.VecVarHisto, w, &t, &tbar, &vecLepM, &vecLepP);
    } // end kinreco
  } // end event loop
  
  // print out
  printf("nReco : %ld\n", nReco);
  printf("nSel  : %ld\n", nSel);
  if(in.Type == 2) 
  {
    printf("nGen  : %ld\n", nGen);
    printf("C = %.2f%% (no KINRECO %.2f%%)\n", 100. * nReco / nGen, 100. * nSel / nGen);
  }

  // store histos
  fout->cd();
  StoreHistos(in.VecVarHisto);
  fout->Close();
}

#endif
