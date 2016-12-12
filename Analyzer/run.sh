#/bin/bash

#INPUTLIST='data/CMS_Run2011A_MuEG_AOD_12Oct2013-v1-first10.txt'
#INPUTLIST='data/CMS_Run2011A_MuEG_AOD_12Oct2013-v1-10percent_index.txt'

#INPUTLIST='data/CMS_Run2011A_MuEG_AOD_12Oct2013-v1-all_file_index.txt'
#INPUTLIST='data/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1-all_file_index.txt'
INPUTLIST='data/CMS_Run2011A_DoubleElectron_AOD_12Oct2013-v1-all_file_index.txt'

#INPUTLIST='data/CMS_Run2011A_DoubleElectron_AOD_12Oct2013-v1_20000_file_index.txt'

# MC TTJets_TuneZ2_7TeV-madgraph-tauola
#PROCESSING 
#DONE INPUTLIST='mc/TTJets_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt'
#DONE INPUTLIST='mc/TTJets_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_00001_file_index.txt'
#DONE INPUTLIST='mc/TTJets_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_010000_file_index.txt'
#DONE INPUTLIST='mc/TTJets_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_010001_file_index.txt'
#DONE INPUTLIST='mc/TTJets_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_010002_file_index.txt'
#DONE INPUTLIST='mc/TTJets_TuneZ2_7TeV-madgraph-tauola/CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_010003_file_index.txt'
# Single t
#INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_AODSIM_PU_S13_START53_LV6_file_index.txt'
INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt'
# Diboson
#INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_WW_TuneZ2_7TeV_pythia6_tauola_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt'
#INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_WZ_TuneZ2_7TeV_pythia6_tauola_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt'
#INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_ZZ_TuneZ2_7TeV_pythia6_tauola_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt'
# WZjets
#DONE INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6_file_index.txt'
#DONE INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt'
#DONE INPUTLIST='mc/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6_file_index.txt'

#INPUTLIST='failed.txt'

NP=5
outrootsuffix=''

#OUTPUTDIR='ttbarSelected_TEST'

#OUTPUTDIR='data22112016_MuEG'
#OUTPUTDIR='data_DoubleMu'
#OUTPUTDIR='data_DoubleElectron'

#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_TTJets_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6-v1_010002'
# Single t
#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_AODSIM_PU_S13_START53_LV6/111'
OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_AODSIM_PU_S13_START53_LV6'
# Diboson
#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_WW_TuneZ2_7TeV_pythia6_tauola_AODSIM_PU_S13_START53_LV6'
#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_WZ_TuneZ2_7TeV_pythia6_tauola_AODSIM_PU_S13_START53_LV6'
#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_ZZ_TuneZ2_7TeV_pythia6_tauola_AODSIM_PU_S13_START53_LV6'
# WZjets
#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6'
#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6_AODSIM_PU_S13_START53_LV6'
#OUTPUTDIR='ttbarSelected_CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_AODSIM_PU_S13_START53_LV6'

if [ ! -f $INPUTLIST ]
then
  echo "Error: no input file $INPUTLIST"
  exit 1
fi
if [ -d $OUTPUTDIR ]
then
  echo "Error: output directory $OUTPUTDIR exists"
  exit 1
else
  mkdir -p $OUTPUTDIR
  p=0
  for inputFile in `cat $INPUTLIST`
  do
    p=$[$p+1]
    echo $inputFile >> "$OUTPUTDIR/inputList"${outrootsuffix}"_$p.txt"
    if [ $p == $NP ]; then p=0; fi
  done
fi


for p in `seq 1 $NP`
#for p in `echo 3 6`
do
  command='time cmsRun analyzer_cfg.py '${OUTPUTDIR}'/inputList'${outrootsuffix}'_'${p}'.txt '${OUTPUTDIR}'/ttbarSel'${outrootsuffix}'_'${p}'.root'
  nohup $command >& ${OUTPUTDIR}/log${outrootsuffix}_${p}.txt&
done

exit 0

