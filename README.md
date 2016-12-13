# CMS measurements of top quark pair production in dilepton channel at 7 TeV

Relevant CMS publication:
 * normalised cross sections: EPJ C73 (2013) 2339 [arXiv:1211.2220, TOP-11-013]
 * total cross section: JHEP 1608 (2016) 029 [arXiv:1603.02303,TOP-13-004]

Further relevant information can be found also in DESY-THESIS-2012-037

There are two parts in this analysis:
 * Analyzer: ntuple production, requires CMSSW (the instructions assume that you will work on a VM properly contextualized for CMS, available from http://opendata.cern.ch/VM/CMS) and network connection; will take ~ 2 weeks to process the full data + MC samples and ~ 3GB free space for the produced ntuples
 * PostAnalyzer: ntuple processing, produces final numbers and plots, standalone code (requires only gcc and ROOT); will take about 15 minutes

## Creating the working area

This step is only needed the first time you run this program:
```
mkdir WorkingArea
cd ./WorkingArea
cmsrel CMSSW_5_3_32
cd ./CMSSW_5_3_32/src
cmsenv
git clone https://github.com/zenaiev/2011-ttbar.git
scram b
cd 2011-ttbar/Analyzer
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1     
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
```
(no need to download data/MC input file lists and JSON: provided with the git)

## Analyzer
Look at Analyzer/README.txt

## PostAnalyzer
Look at Postanalyzer/README.txt
