#!/bin/bash

# compile code (produces two executables)
g++ -o ttbarReco `root-config --cflags --libs` -lMathMore -std=c++11 ttbarReco.cxx
g++ -o ttbarPlots `root-config --cflags --libs` -std=c++11 ttbarPlots.cxx

# create needed directories if do not exist yet
mkdir -p data mc hist plots
