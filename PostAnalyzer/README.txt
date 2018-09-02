This part processes ROOT ntuples (output of Analyzer, see 
Analyzer/README.txt) to produce ROOT histograms (ttbarMakeHist.cxx) 
and final plots and numbers from them (ttbarMakePlots.cxx).
This is pure C++ and ROOT code: does not require CMSSW, works outside VM 
(although it can work on VM also, of course).

General description of contents (find further description inside the files):
   hggMakeHist.cxx: master file to produce histograms
   eventReco.h: ttbar event reconstruction
   selection.h: ttbar event selection
   kinReco.h: kinematic reconstruction
   tree.h: tree structure of input ROOT ntuples
   settings.h: global settings (directory names)
   hggMakePlots.cxx: master file to produce final plots and numbers
   plots.h: helper file for plotting

To run the analysis, make sure input ntuples are in place, for default 
directory structure you need to run from the root analysis directory:
mv Analyzer/ntuples-data PostAnalyzer/ntuples-data
mv Analyzer/ntuples-mc PostAnalyzer/ntuples-mc
then compile the code:
./compile.sh
and run three commands:
./hggMakeHist
./hggMakePlots
./ pvalPlot

Also you could do only the last step (plotting) by using "reference" 
histograms produced with the full samples and available with the code 
(PostAnalyzerhist-REF directory), for this modify settings.h. 
Another application of the "reference" histograms could be for 
validation (produce new histograms and compare to the reference ones).

To calculate luminosity, run:
python lumicalc.py
