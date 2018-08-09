#!/bin/bash

# compile code (produces two executables)
g++ -g -o hggMakeHist `root-config --cflags --libs` -lMathMore -std=c++11 hggMakeHist.cxx
g++ -g -o hggMakePlots `root-config --cflags --libs` -std=c++11 hggMakePlots.cxx
g++ -g -o pvalPlot `root-config --cflags --libs` -std=c++11 pvalPlot.cxx

# create needed directories if do not exist yet
mkdir -p data mc hist plots
