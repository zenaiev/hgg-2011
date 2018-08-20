# CMS analysis of Higgs -> gamma gamma at 7 and 8 TeV

Relevant CMS publications:
 * 7 TeV: PLB710 (2012) 403 [arXiv:1202.1487, HIG-11-033]
 * 8 TeV: EPJ C74 (2014) 3076 [arXiv:1407.0558,HIG-13-001]

For the general description of the analysis see also attached description-hgg.pdf (to be added)

There are two parts in this analysis:
 * Analyzer: ntuple production, requires CMSSW (the instructions assume that you will work on a VM properly contextualized for CMS, available from http://opendata.cern.ch/VM/CMS) and network connection; takes ~ 2 weeks to process the full data + MC samples and ~ 3GB free space for the produced ntuples
 * PostAnalyzer: ntuple processing, produces final numbers and plots, standalone code (requires only gcc and ROOT); takes about 2 minutes

## Instructions how to run the analysis

First you need to create the working area (this step is only needed the first time you setup this program). You can create the working area for this analysis on the VM which has other instances of CMSSW, just keep them in different directories.
```
mkdir WorkingArea
cd ./WorkingArea
cmsrel CMSSW_5_3_32
cd ./CMSSW_5_3_32/src
cmsenv
git clone https://github.com/christian512/hgg-2011.git
scram b
cd 2011-doubleelectron-doublemu-mueg-ttbar/Analyzer
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL FT53_V21A_AN6
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL.db FT53_V21A_AN6_FULL.db
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT53_V21A_AN6_FULL FT53_V21A_AN6_FULL
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27.db START53_V27.db
ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_V27 START53_V27
```
(no need to download data/MC input file lists and JSON: provided with the code)

Also the code in PostAnalyzer should be compiled:
```
cd PostAnalyzer
./compile.sh
```

## Running the analysis
Generally, the analysis steps are:
 * run Analyzer/run.sh (look inside first), this processes AOD files (CMS data stored at CERN server, several TB) and produces plain ROOT ntuple files (~...GB), takes ~ ... weeks, extensive network access
 * move produced ntuples to PostAnalyzer directories (this step is manual on purpose, in order not to overwrite accidentally ntuples produced taking long time etc.)
 * run PostAnalyzer/ttbarMakeHist to process ROOT ntuples to create histograms (~2 mins)
 * run PostAnalyzer/ttbarMakePlots to produce final plots from created histograms (few seconds)

Further description of these steps you can find Analyzer/README.txt and Postanalyzer/README.txt
