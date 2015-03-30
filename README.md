# Tprime_MassPlots
Macro Usage:

#Mass plots
root -l

gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so")

.L makeBkgPlotsGeneric_my.C

makeBkgPlotsGeneric_my("where_data_are.root","",false)

Find .eps and .C file here:

https://gfasanel.web.cern.ch/gfasanel/Tprime_analysis/baseline_Tprime_giuseppe_rereco/new/

#Limit plots (1D)
You need a directory with the root files output from combine

python limitPlotter.py -M Asymptotic -p path_limit_root_file

 -r for ratio 

 -7 hadronic category

 -6 leptonic category

 -e expected only

Find .eps and .C here:

https://gfasanel.web.cern.ch/gfasanel/Tprime_analysis/LIMITS_rereco/new/

#2D limits plot: exclusion area
.L MakePlot2D_my.C+

make2D()

It will use border_Tprime_*.txt

Find .eps and .C here:

https://gfasanel.web.cern.ch/gfasanel/Tprime_analysis/2D_EXCL_PLOT/

#2D limits plot: UL in the color scale
root -l -b

.L 2Dplot.cc

make2Dplot(4)

make2Dplot(5)

It will use UL_Tprime*.txt

#Triangle plots:
./runTrianglePlot.py -l './multileptons_v*' -t theory.txt> dump.txt

It will need all the directories multileptons_v* and thery.txt

Find .eps and .C here:

https://gfasanel.web.cern.ch/gfasanel/Tprime_analysis/TRIANGLE_PLOTS/