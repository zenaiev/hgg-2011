# CMS analysis of Higgs -> gamma gamma at 7 and 8 TeV

Relevant CMS publications:
 * 7 TeV: PLB710 (2012) 403 [arXiv:1202.1487, HIG-11-033]
 * 8 TeV: EPJ C74 (2014) 3076 [arXiv:1407.0558,HIG-13-001]

For the general description of the analysis see also attached description_hgg.pdf

There are two parts in this analysis:
 * Analyzer: ntuple production, requires CMSSW (the instructions assume that you will work on a VM properly contextualized for CMS, available from http://opendata.cern.ch/VM/CMS) and network connection; takes ~ 4.5 month to process the full data + MC samples and ~ 3 GB free space for the produced ntuples
 * PostAnalyzer: ntuple processing, produces final numbers and plots, standalone code (requires only gcc and ROOT); takes about 2 minutes

## Instructions how to run the analysis

First you need to create the working area (this step is only needed the first time you setup this program). You can create the working area for this analysis on the VM which has other instances of CMSSW, just keep them in different directories.
```
mkdir WorkingArea
cd ./WorkingArea
cmsrel CMSSW_5_3_32
cd ./CMSSW_5_3_32/src
cmsenv
git clone git://github.com/christian512/2011-photon-2012-doublephoton-higgs-hgaga.git
scram b
cd 2011-photon-2012-doublephoton-higgs-hgaga/Analyzer
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
 * run Analyzer/run.sh (look inside first), this processes AOD files (CMS data stored at CERN server, ~ 25 TB) and produces plain ROOT ntuple files (~ 3 GB), takes ~ 4.5 month, extensive network access
 * move produced ntuples to PostAnalyzer directories (this step is manual on purpose, in order not to overwrite accidentally ntuples produced taking long time etc.)
 * run PostAnalyzer/hggMakeHist to process ROOT ntuples to create histograms (~5 mins)
 * run PostAnalyzer/hggMakePlots to produce final plots from created histograms (few seconds)
 * run PostAnalyzer/pvalPlot to create a simplified significance plot (few seconds) 
 * run PostAnalyzer/lumicalc.py to calculate the integrated luminosity (~5 mins)
 
 ```
 cd hgg-2011/Analyzer
 ./run.sh 1 #run 2011 data
 ./run.sh 2 #run 2011 MC not available
 ./run.sh 3 #run 2012 data
 ./run.sh 4 #run 2012 MC
 cd ..
 mv Analyzer/ntuples-data PostAnalyzer
 mv Analyzer/ntuples-mc PostAnalyzer
 cd PostAnalyzer
 ./compile.sh
 ./hggMakeHist
 ./hggMakePlots
 ./pvalPlot
 python lumicalc.py
 ```
Further description of these steps can be found in the desription-hgg.pdf.

## Git Pushing
If you connect your local git repository to a remote repository you need to push your changes. Due to several problems, the used version of the CMS-Software (CMSSW_5_3_32) and the included version of git are not able to push repositories. Therefore you can change to a newer version of CMSSW for pushing only. Steps are described below:

```
cd WorkingArea/ 
cmsrel CMSSW_10_1_9
cd CMSSW_10_1_9/src/
cmsenv
cd ../../CMSSW_5_3_32/src/2011-photon-2012-doublephoton-higgs-hgaga
#perform your git actions here, e.g.:
git push
```
After you executed your git commands you can switch to the old version of CMSSW_5_3_32 (to execute code) by
```
#somewhere in CMSSW_5_3_32/src/
cmsenv 
```
If you have to push any new changes, it is sufficient to switch to the new `CMSSW_10_1_9/src` folder and execute `cmsenv`

