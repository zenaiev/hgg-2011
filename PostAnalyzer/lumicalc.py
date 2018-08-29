# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>> CALCULATOR FOR INTEGRATED LUMINOSITY >>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Runs through all root ntuple files to get all run numbers
# Sums up integrated luminosity for all used runs (with at least 1 trigger fired)
# Luminosity file taken from: http://opendata.cern.ch/record/1052

import ROOT as root
import numpy as np
import os

#file stores luminosity for every run
lumifile_path = '2012lumi.txt'
lumifile_url = 'http://opendata.cern.ch/record/1052/files/2012lumi.txt'

#recorded and delivered luminosity
lumi_rec = 0
lumi_del = 0

#ntuple counter
count_all_files = 0
count_sel_files = 0

#path to ntuples
ntuples_dir = 'ntuples-data/DoublePhoton'

#get a list of all root file
filelist = os.listdir(ntuples_dir)

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
    f = root.TFile(os.path.join(ntuples_dir,f))
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
          
#delete all double entries ( there shouldnt be any!) 
runNums = np.unique(runNums)

# **********************************************************************
# ************ CALCULATE INTEGRATED LUMINOSITY *************************
# **********************************************************************
print('calculating luminosity')


with open(lumifile_path) as lumifile:
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
    line = lumifile.readline()

# **********************************************************************
# ****************** OUTPUT ********************************************
# **********************************************************************
print '********** LUMINOSITY **********'
print 'Recorded lumi: ' + str(lumi_rec*10**(-9)) + '/fb'
print 'Delivered lumi: ' + str(lumi_del*10**(-9)) + '/fb'
print '********** FILE STATISTICS **********'
print 'Number of all ROOT files: ' + str(count_all_files)
print 'Number of selected ROOT files: ' + str(count_sel_files)

  
