#include "eventReco.h"
#include "settings.h"


TString dataDir = gDataDir;
TString mcDir = gMcDir;

// flags what to run
bool flagData = 1;
bool flagMCsig = 1;
bool flagMCother = 1;
//bool flagMCback = 1;
bool flagMCstop = 1;
bool flagMCwjets = 1;
bool flagMCdy = 1;

int main(int argc, char** argv)
{
  // common vars
  std::vector<TString> nameInFile;
  
  // histograms
  TH1::SetDefaultSumw2();
  std::vector<ZVarHisto> vecVH, vecVHGen;
  // control plots
  vecVHGen.push_back(ZVarHisto("ptt", new TH1D("h_ptt", "pT top", 21, 0.0, 420.0)));
  vecVHGen.push_back(ZVarHisto("ptat", new TH1D("h_ptat", "pT atop", 21, 0.0, 420.0)));
  vecVHGen.push_back(ZVarHisto("pttat", new TH1D("h_pttat", "pT tatop", 21, 0.0, 420.0)));
  vecVHGen.push_back(ZVarHisto("pttt", new TH1D("h_pttt", "pT ttbar", 30, 0.0, 300.0)));
  vecVHGen.push_back(ZVarHisto("yt", new TH1D("h_yt", "y top", 26, -2.6, 2.6)));
  vecVHGen.push_back(ZVarHisto("yat", new TH1D("h_yat", "y atop", 26, -2.6, 2.6)));
  vecVHGen.push_back(ZVarHisto("ytat", new TH1D("h_ytat", "y tatop", 26, -2.6, 2.6)));
  vecVHGen.push_back(ZVarHisto("ytt", new TH1D("h_ytt", "y ttbar", 26, -2.6, 2.6)));
  // cross sections
  {
    double bins[] = {0.,80.,130.,200.,300.,400.};
    vecVHGen.push_back(ZVarHisto("ptt", new TH1D("h_ptt_cs", "pT top", 5, bins)));
    vecVHGen.push_back(ZVarHisto("ptat", new TH1D("h_ptat_cs", "pT atop", 5, bins)));
    vecVHGen.push_back(ZVarHisto("pttat", new TH1D("h_pttat_cs", "pT tatop", 5, bins)));
  }
  {
    double bins[] = {-2.5,-1.3,-0.8,-0.4,0.0,0.4,0.8,1.3,2.5};
    vecVHGen.push_back(ZVarHisto("yt", new TH1D("h_yt_cs", "y top", 8, bins)));
    vecVHGen.push_back(ZVarHisto("yat", new TH1D("h_yat_cs", "y atop", 8, bins)));
    vecVHGen.push_back(ZVarHisto("ytat", new TH1D("h_ytat_cs", "y tatop", 8, bins)));
  }
  {
    double bins[] = {0.,20.,60.,120.,300.};
    vecVHGen.push_back(ZVarHisto("pttt", new TH1D("h_pttt_cs", "pT ttbar", 4, bins)));
  }
  {
    double bins[] = {-2.5,-1.5,-0.7,0.0,0.7,1.5,2.5};
    vecVHGen.push_back(ZVarHisto("ytt", new TH1D("h_ytt_cs", "y ttbar", 6, bins)));
  }
  {
    double bins[] = {0.,345.,400.,470.,550.,650.,800.,1100.,1600.};
    vecVHGen.push_back(ZVarHisto("mtt", new TH1D("h_mtt_cs", "M ttbar", 8, bins)));
  }

  vecVH = vecVHGen;
  vecVH.push_back(ZVarHisto("ptl", new TH1D("h_ptl", "pT leptons", 23, 30.0, 260.0)));
  //vecVH.push_back(ZVarHisto("ptj", new TH1D("h-ptj", "pT jet", 37, 30.0, 400.0)));
  
  for(int ch = 1; ch <= 3; ch++)
  {
    //if(ch != 3) continue;
    // *****************************************
    // **************** DATA *******************
    // *****************************************
    if(flagData)
    {
      ZEventRecoInput in;
      //in.MaxNEvents = 100;
      in.Name = "data";
      in.Type = 1;
      in.Channel = ch;
      in.VecVarHisto = vecVH;
      //in.NameFout = TString::Format("%s/data-c%d.root", outDir.Data(), ch);
      if(ch == 1) // ee
      {
        in.AddToChain(dataDir + "/DoubleElectron/*.root");
        //in.AddToChain(dataDir + "/DoubleElectron/part0/*.root");
        //in.AddToChain(dataDir + "/DoubleElectron/part0_2/*.root");
        //in.AddToChain(dataDir + "/DoubleElectron/part1/*.root");
      }
      else if(ch == 2) // mumu
      {
        in.AddToChain(dataDir + "/DoubleMu/*.root");
        //in.AddToChain(dataDir + "/DoubleMu/2/*.root");
      }
      else if(ch == 3) // emu
      {
        in.AddToChain(dataDir + "/MuEG/*.root");
        //in.AddToChain(dataDir + "/MuEG/2/*.root");
        //in.AddToChain(dataDir + "/MuEG/3/*.root");
        //in.AddToChain(dataDir + "/MuEG/4/*.root");
      }
      eventreco(in);
    }
    
    // *****************************************
    // ************** MC signal ****************
    // *****************************************
    // events: 54990752
    // sigma: 95.43 -> 165.6
    // weight: 0.5 * 5000.0 / (54990752. / 95.43) * (165.6 / 95.43) = 0.007529
    if(flagMCsig)
    {
      ZEventRecoInput in;
      //in.MaxNEvents = 1000;
      in.Weight = 0.007529;
      in.Name = "mcSigReco";
      in.Type = 2;
      in.Channel = ch;
      //in.NameFout = TString::Format("%s/mcSigReco-c%d.root", outDir.Data(), ch);
      in.VecVarHisto = vecVH;
      //in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_00001/ttbarSel_1.root");
      //in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_00001/ttbarSel_2.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/00001/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010000/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010003/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010002/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010001/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/00000/*.root");
      // signal reco level
      eventreco(in);
      // MC ttbar other
      in.Name = "mcSigOtherReco";
      //in.NameFout = TString::Format("%s/mcSigReco-c%d.root", outDir.Data(), ch);
      in.Type = 3;
      eventreco(in);
      // signal gen level
      in.Name = "mcSigGen";
      in.Type = 2;
      //in.NameFout = TString::Format("%s/mcSigGen-c%d.root", outDir.Data(), ch);
      in.VecVarHisto = vecVHGen;
      in.Gen = true;
      eventreco(in);
    }
    // *****************************************
    // ************ MC single top **************
    // *****************************************
    // Events: 744859 + 801626
    // Sigma: 7.475 -> 7.87
    // weight: 0.5 * 5000. / ((744859. + 801626.) / (7.475 * 2.)) * (7.87 / 7.475) = 0.02544
    if(flagMCstop)
    {
      ZEventRecoInput in;
      //in.MaxNEvents = 1000;
      in.Name = "mcSingleTopReco";
      in.Type = 4;
      in.Weight = 0.02544;
      in.Channel = ch;
      //in.NameFout = TString::Format("%s/mcSingleTopReco-c%d.root", outDir.Data(), ch);
      in.VecVarHisto = vecVH;
      in.AddToChain(mcDir + "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*.root");
      in.AddToChain(mcDir + "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*.root");
      eventreco(in);
    }
    // *****************************************
    // ************** MC W+jets ****************
    // *****************************************
    // Events: 78347691
    // Sigma: 25430 -> 31314
    // weight: 0.5 * 5000.0 / (78347691. / (25430. * 0.32)) * (31314. / 25430.) = 0.3197
    if(flagMCstop)
    {
      ZEventRecoInput in;
      //in.MaxNEvents = 1000;
      in.Name = "mcWjetsReco";
      in.Weight = 0.3197;
      in.Type = 4;
      in.Channel = ch;
      //in.NameFout = TString::Format("%s/mcWjetsReco-c%d.root", outDir.Data(), ch);
      in.VecVarHisto = vecVH;
      in.AddToChain(mcDir + "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*.root");
      eventreco(in);
    }
    // *****************************************
    // **************** MC DY ******************
    // *****************************************
    if(flagMCstop)
    {
      ZEventRecoInput in;
      //in.MaxNEvents = 1000;
      in.Type = 4;
      in.Channel = ch;
      in.VecVarHisto = vecVH;
      // low mass
      // Events: 39909640
      // Sigma: 9487 -> 11908
      // weight: 0.5 * 5000. / (39909640. / (9487. * 0.1)) * (11908. / 9487.) = 0.07459
      in.Name = "mcDYlmReco";
      in.Weight = 0.07459;
      //in.NameFout = TString::Format("%s/mcDYlmReco-c%d.root", outDir.Data(), ch);
      in.AddToChain(mcDir + "/DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6/*.root");
      eventreco(in);
      // high mass
      // Events: 36408225
      // Sigma: 2513 -> 3048
      // weight: 0.5 * 5000. / (36408225. / (2513. * 0.1)) * (3048. / 2513.) = 0.02093
      in.Name = "mcDYhmReco";
      in.Weight = 0.2093;
      //in.NameFout = TString::Format("%s/mcDYhmReco-c%d.root", outDir.Data(), ch);
      in.ClearChain();
      in.AddToChain(mcDir + "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/*.root");
      eventreco(in);
    }
  }

  return 0;
}
