This part processes ROOT ntuples (output of Analyzer, see 
Analyzer/README.txt) to produce ROOT histograms (ttbarMakeHist.cxx) 
and final plots and numbers from them (ttbarMakePlots.cxx).
This is pure C++ and ROOT code: does not require CMSSW, works outside VM 
(although it can work on VM also, of course).

General description of contents (find further description inside the files):
   ttbarMakeHist.cxx: master file to produce histograms
   eventReco.h: ttbar event reconstruction
   selection.h: ttbar event selection
   kinReco.h: kinematic reconstruction
   tree.h: tree structure of input ROOT ntuples
   settings.h: global settings (directory names)
   ttbarMakePlots.cxx: master file to produce final plots and numbers
   plots.h: helper file for plotting

To run the analysis, compile the code:
./compile.sh
and run two commands:
./ttbarMakeHist
./ttbarMakePlots
