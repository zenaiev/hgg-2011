// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>> General settings (directories) >>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifndef TTBAR_SETTINGS_H
#define TTBAR_SETTINGS_H

#include <TString.h>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// main directory (by default current directory)
TString gBaseDir  = "./"; 
//
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// these are "normal" settings 
//TString gDataDir  = gBaseDir + "/data"; // directory with data ntuples
//TString gMcDir    = gBaseDir + "/mc"; // directory with MC ntuples
//TString gHistDir  = gBaseDir + "/hist"; // directory with histograms
//TString gPlotsDir = gBaseDir + "/plots"; // directory with final plots
//
// For exercises, you could use existing "reference" histograms 
// (they are provided at git) to produce final plots, or even existing 
// ROOT ntuples to produce histograms, for this uncomment proper lines 
// below (and comment out above, of course)
TString gDataDir  = gBaseDir + "/data-REF";
TString gMcDir    = gBaseDir + "/mc-REF";
TString gHistDir  = gBaseDir + "/hist-REF";
TString gPlotsDir = gBaseDir + "/plots";
//
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#endif
