// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This code processes ROOT ntuples for hgg analysis (see
// Analyzer/src/Analyzer.cc) and produces histograms, which are 
// further used to make final plots (see hggMakePlots.cxx).
// Run: ./hggMakeHist
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// additional files from this analysis (look there for description) 
#include "eventReco.h"
#include "settings.h"
#include "TROOT.h"

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
  bool flagData2011_10bins = 1; // if 1, 2011 data will be processed into 1.0 GeV bins
  bool flagData2011_15bins = 1; // if 1, 2011 data will be processed into 1.5 GeV bins
  bool flagData2012_10bins = 1; // if 1, 2012 data will be processed into 1.0 GeV bins
  bool flagData2012_15bins = 1; // if 1, 2012 data will be processed into 1.5 GeV bins
  bool flagMCsig    = 1; // if 1, signal MC will be processed (only 2012)
  //
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // common purpose variables
  std::vector<TString> nameInFile; // container to store input file names
  
  // histograms
  TH1::SetDefaultSumw2(); // keep histogram weights by default
  // ZVarHisto is a simple class which incorporates a histogram and a variable name. 
  // This class is used to store needed input settings (variable names, binning)
  // **********************************************
  // **************** DATA 2011 *******************
  // **********************************************
  if(flagData2011_10bins)
  {
    //setup histograms
	//ZVarHisto is a class that stores histogram and variable name, see its definition in eventReco.h
    std::vector<ZVarHisto> vecVH, vecVHGen;
    vecVHGen.push_back(ZVarHisto("mgg", new TH1D("h_mgg2011_1_", "m_{#gamma#gamma}", 80, 100.0, 180.0)));
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
    // ZEventRecoInput is a class for event reconstruction, see its description in eventReco.h
    ZEventRecoInput in;
    //in.MaxNEvents = 100; // if you need to limit the number of processed events
    in.Name = "data2011_10GeV"; // name pattern for output histograms
    in.Type = 1; // type = 1 for data
    in.VecVarHisto = vecVH; // need to copy it, because further will be changed
    // input ROOT ntuples
    in.AddToChain(dataDir + "/Photon/*.root"); // 2011
    // main part: event reconstruction call
    eventreco(in);
  }
  
  if(flagData2011_15bins)
  {
    //setup histograms
    std::vector<ZVarHisto> vecVH, vecVHGen;
    vecVHGen.push_back(ZVarHisto("mgg", new TH1D("h_mgg2011_2_", "m_{#gamma#gamma}", 60, 100.0, 190.0))); // m(gammagamma)
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
    // ZEventRecoInput is a class for event reconstruction, see its description in eventReco.h
    ZEventRecoInput in;
    //in.MaxNEvents = 100; // if you need to limit the number of processed events
    in.Name = "data2011_15GeV"; // name pattern for output histograms
    in.Type = 1; // type = 1 for data
    in.VecVarHisto = vecVH; // need to copy it, because further will be changed
    // input ROOT ntuples
    in.AddToChain(dataDir + "/Photon/*.root"); // 2011
    // main part: event reconstruction call
    eventreco(in);
  }

  // **********************************************
  // **************** DATA 2012 *******************
  // **********************************************
  if(flagData2012_10bins)
  {
    //setup histograms
    std::vector<ZVarHisto> vecVH, vecVHGen;
    vecVHGen.push_back(ZVarHisto("mgg", new TH1D("h_mgg2012_1_", "m_{#gamma#gamma}", 80, 100.0, 180.0))); // m(gammagamma)
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
    // ZEventRecoInput is a class for event reconstruction, see its description in eventReco.h
    ZEventRecoInput in;
    //in.MaxNEvents = 100; // if you need to limit the number of processed events
    in.Name = "data2012_10GeV"; // name pattern for output histograms
    in.Type = 1; // type = 1 for data
    in.VecVarHisto = vecVH; // need to copy it, because further will be changed
    // input ROOT ntuples
    in.AddToChain(dataDir + "/DoublePhoton/*.root"); // 2012
    // main part: event reconstruction call
    eventreco(in);
  }

  if(flagData2012_15bins)
  {
    //setup histograms
    std::vector<ZVarHisto> vecVH, vecVHGen;
    vecVHGen.push_back(ZVarHisto("mgg", new TH1D("h_mgg2012_2_", "m_{#gamma#gamma}", 60, 100.0, 190.0))); // m(gammagamma)
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
    // ZEventRecoInput is a class for event reconstruction, see its description in eventReco.h
    ZEventRecoInput in;
    //in.MaxNEvents = 100; // if you need to limit the number of processed events
    in.Name = "data2012_15GeV"; // name pattern for output histograms
    in.Type = 1; // type = 1 for data
    in.VecVarHisto = vecVH; // need to copy it, because further will be changed
    // input ROOT ntuples
    in.AddToChain(dataDir + "/DoublePhoton/*.root"); // 2012
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
  // MC cross section can be obtained from any ROOT file in the mC sample: open the ROOT file, create TBrowser, navigate to
  // Runs -> GenRunInfoProduct_generator__SIM. -> GenRunInfoProduct_generator__SIM.obj -> InternalXSec -> value_
  // (nevertheless sigma_MC cancels)
  //
  // sigma_theory accounts for H->gaga branching ratio (0.000229)
  //
  // number of events: 292178
  // MC cross section -> theory: 12.93 -> 19.5
  // weight: 9850.0 / (292178. / 12.93) * (19.5 * 0.00229 / 12.93) = 0.0015054
  if(flagMCsig)
  {
    //setup histograms
    std::vector<ZVarHisto> vecVH, vecVHGen;
    vecVHGen.push_back(ZVarHisto("mgg", new TH1D("h_mgg", "m_{#gamma#gamma}", 80, 100.0, 180.0))); // m(gammagamma)
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
    // MC signal reco level
    ZEventRecoInput in;
    in.Weight = 0.0015054; // weight (see above)
    in.Name = "mcSigReco";
    in.Type = 2;
    in.VecVarHisto = vecVH;
    in.AddToChain(mcDir + "/GluGluToHToGG_M-125_8TeV-powheg15-pythia6/*.root");
    eventreco(in);
    // MC signal, generator level: again re-use existing ZEventRecoInput, change type and set proper flag (see below)
    in.Name = "mcSigGen";
    in.Type = 2;
    in.VecVarHisto = vecVHGen;
    in.Gen = true; // flag to notify that generator level should be processed
    eventreco(in);
  }

  return 0;
}
