// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This code processes ROOT ntuples for ttbar analysis (see 
// Analyzer/src/Analyzer.cc) and produces histograms, which are 
// further used to make final plots (see ttbarMakePlots.cxx).
// Run: ./ttbarMakeHist
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// additional files from this analysis (look there for description) 
#include "eventReco.h"
#include "settings.h"
//
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>> Main function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char** argv)
{
  //
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>> Settings >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // set directories to data and MC ntuples
  TString dataDir = gDataDir;
  TString mcDir = gMcDir;
  //
  // flags what to run
  bool flagData    = 1; // if 1, data will be processed
  bool flagMCsig   = 1; // if 1, signal MC will be processed
  //bool flagMCstop  = 1; // if 1, MC single top (background) will be processed
  //
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // common purpose variables
  std::vector<TString> nameInFile; // container to store input file names
  
  // histograms
  TH1::SetDefaultSumw2(); // keep histogram weights by default
  // ZVarHisto is a simple class which incorporates a histogram and a variable name. 
  // This class is used to store needed input settings (variable names, binning) 
  // for control plots and cross sections (as in TOP-11-013)
  std::vector<ZVarHisto> vecVH, vecVHGen; // vecVH for reconstruction level, vecVHGen for generator level
  // histograms and variables for control plots
  vecVHGen.push_back(ZVarHisto("mgg", new TH1D("h_mgg", "m_{#gamma#gamma}", 80, 100.0, 180.0))); // m(gammagamma)
  // histograms and variables for cross sections
  /*{
    double bins[] = {0.,80.,130.,200.,300.,400.};
    vecVHGen.push_back(ZVarHisto("ptt", new TH1D("h_ptt_cs", "pT top", 5, bins)));
    vecVHGen.push_back(ZVarHisto("ptat", new TH1D("h_ptat_cs", "pT atop", 5, bins)));
    vecVHGen.push_back(ZVarHisto("pttat", new TH1D("h_pttat_cs", "pT tatop", 5, bins)));
  }*/

  // for reconstruction level the same binning is needed, but more event classes
  vecVH = vecVHGen;
  vecVH.push_back(vecVHGen[0]); // m(gammagamma)
  vecVH.back().EventClass() = 2;
  vecVH.push_back(vecVHGen[0]); // m(gammagamma)
  vecVH.back().EventClass() = 3;
  vecVH.push_back(vecVHGen[0]); // m(gammagamma)
  vecVH.back().EventClass() = 4;
  vecVH.push_back(vecVHGen[0]); // m(gammagamma)
  vecVH.back().EventClass() = 5;
  vecVH.push_back(vecVHGen[0]); // m(gammagamma)
  vecVH.back().EventClass() = 6;
  // (here you can add more reconstruction level histograms)
  
  // *****************************************
  // **************** DATA *******************
  // *****************************************
  if(flagData)
  {
    // ZEventRecoInput is a class for event reconstruction, see its description in eventReco.h
    ZEventRecoInput in;
    //in.MaxNEvents = 100; // if you need to limit the number of processed events
    in.Name = "data"; // name pattern for output histograms
    in.Type = 1; // type = 1 for data
    in.VecVarHisto = vecVH; // need to copy it, because further will be changed
    // input ROOT ntuples
    in.AddToChain(dataDir + "/Photon/*.root");
    // main part: event reconstruction call
    eventreco(in);
  }

  // *****************************************
  // ************** MC signal ****************
  // *****************************************
  //
  // MC event weights need to be changed to most precise theoretical predictions,
  // the formula is:
  // weight = lumi / (nevents / sigma_MC) * (sigma_theory / sigma_MC) = lumi * nevents / sigma_theory
  //
  // Number of events can be obtained from webpage (see http://opendata.cern.ch/collection/CMS-Simulated-Datasets),
  // but it should be checked that all events have been processed at the Analyzer step (see end of log files)
  //
  // number of events: 54990752
  // MC cross section -> theory: 95.43 -> 165.6
  // weight: 2500.0 / (54990752. / 95.43) * (165.6 / 95.43) = 0.007529
  if(flagMCsig)
  {
    // MC signal reco level
    ZEventRecoInput in;
    //in.MaxNEvents = 1000;
    in.Weight = 1.0; // weight (see above)
    in.Name = "mcSigReco";
    in.Type = 2;
    in.VecVarHisto = vecVH;
    in.AddToChain(mcDir + "/VBFHiggs0PToGG_M-125p6_7TeV-JHUGenV4-pythia6-tauola/*.root");
    eventreco(in);
    /*// MC other (background): re-use existing ZEventRecoInput, just change type
    in.Name = "mcSigOtherReco";
    in.Type = 3;
    eventreco(in);*/
    // MC signal, generator level: again re-use existing ZEventRecoInput, change type and set proper flag (see below)
    in.Name = "mcSigGen";
    in.Type = 2;
    in.VecVarHisto = vecVHGen;
    in.Gen = true; // flag to notify that generator level should be processed
    eventreco(in);
  }
  // *****************************************
  // ************ MC single top **************
  // *****************************************
  // number of events: 744859 + 801626
  // MC cross section -> theory: 7.475 -> 7.87
  // weight: 2500.0 / ((744859. + 801626.) / (7.475 * 2.)) * (7.87 / 7.475) = 0.02544
  /*if(flagMCstop)
  {
    ZEventRecoInput in;
    //in.MaxNEvents = 1000;
    in.Name = "mcSingleTopReco";
    in.Type = 4;
    in.Weight = 0.02544;
    in.Class = ch;
    in.VecVarHisto = vecVH;
    in.AddToChain(mcDir + "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*.root");
    in.AddToChain(mcDir + "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*.root");
    eventreco(in);
  }*/

  return 0;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
