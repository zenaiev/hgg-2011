# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>> CALCULATOR FOR INTEGRATED LUMINOSITY >>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Runs through all root ntuple files to get all run numbers
# Sums up integrated luminosity for all used runs (with at least 1 trigger fired)
# Luminosity file taken from: 
# For 2011 data : http://opendata.cern.ch/record/1051
# For 2012 data : http://opendata.cern.ch/record/1052


import ROOT as root
import numpy as np
import os

#file stores luminosity for every run
lumifilename12 = '2012lumi.txt'
lumifilename11 = '2011lumi.txt'
lumifile_url12 = 'http://opendata.cern.ch/record/1052/files/2012lumi.txt'



#path to ntuples
ntuples_dirs = ['ntuples-data/Photon','ntuples-data/DoublePhoton']

# **********************************************************************
# ********** 2011 and 2012 data ****************************************
# **********************************************************************
for directory in ntuples_dirs:
  
  #recorded and delivered luminosity
  lumi_rec = 0
  lumi_del = 0
  #counter
  count_all_files = 0 #all available files
  count_sel_files = 0 #all used files for calc
  count_all_runs = 0 #all available runs
  count_sel_runs = 0 #all used runs for calc
  #create file list
  filelist = os.listdir(directory)
  #set year filename
  filename = ''
  if directory == ntuples_dirs[0]:
    print '\n ********** 2011 DATA ********** '
    filename = lumifilename11
    
  if directory == ntuples_dirs[1]:
    print '\n ********** 2012 DATA ********** '
    filename = lumifilename12

  #create list containing all runs from all root files
  runNums = np.array([],dtype=int)

  # **********************************************************************
  # **************  GET ALL RUN NUMBERS **********************************
  # **********************************************************************
  print('Reading all ntuple files')
  for f in filelist:
    #check for root file
    if f.endswith('.root'):
      count_all_files += 1
      #open root file
      f = root.TFile(os.path.join(directory,f))
      #open tree
      myTree = f.Get("tree;1")
      #this occurs when ROOT file was Zombie
      if type(myTree) != root.TTree:
        continue #next ROOT file
      #increase selected files counter
      count_sel_files += 1
      #check each event
      for event in myTree:
        #get run number and triggers
        run = event.evRunNumber
        trigger = event.Triggers
        #check if run already in run list
        if run not in runNums:
          #at least one trigger needs to be fired 
          if trigger > 0:
            runNums = np.append(runNums,int(run))
            count_all_runs += 1
            
  #delete all double entries ( there shouldnt be any!) 
  runNums = np.unique(runNums)

  # **********************************************************************
  # ************ CALCULATE INTEGRATED LUMINOSITY *************************
  # **********************************************************************
  print('calculating luminosity')
  with open(filename) as lumifile:
    line = lumifile.readline()
    #skip beginning 
    while line[0] != '|':
      line = lumifile.readline()
    line = lumifile.readline()
    while line[0] != '|':
      line = lumifile.readline()
    #reached data entries now
    while line[0] != '+':
      strings = line.split('|')
      #get run number
      run = int(strings[1].split(':')[0])
      if run in runNums:
        #get recorded and delivered lumi
        lumi_rec += float(strings[6])
        lumi_del += float(strings[5])
        #delete run from runNums
        ind = np.where(run == runNums)
        runNums = np.delete(runNums, ind)
        count_sel_runs += 1
      line = lumifile.readline()

  # **********************************************************************
  # ****************** OUTPUT ********************************************
  # **********************************************************************
  print '********** LUMINOSITY **********'
  print 'Recorded lumi: ' + str(lumi_rec*10**(-9)) + '/fb'
  print 'Delivered lumi: ' + str(lumi_del*10**(-9)) + '/fb'
  print '********** RUN STATISTICS **********'
  print 'Number of available runs: ' + str(count_all_runs)
  print 'Number of used runs: ' + str(count_sel_runs)
  print '********** FILE STATISTICS **********'
  print 'Number of all ROOT files: ' + str(count_all_files)
  print 'Number of selected ROOT files: ' + str(count_sel_files)

  
