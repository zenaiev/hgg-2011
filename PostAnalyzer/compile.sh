#!/bin/bash

g++ -o ttbarReco `root-config --cflags --libs` -lMathMore -std=c++11 ttbarReco.cxx
g++ -o ttbarPlots `root-config --cflags --libs` -std=c++11 ttbarPlots.cxx
