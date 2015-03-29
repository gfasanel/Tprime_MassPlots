#!/usr/bin/env python
# Original Authors - Nicholas Wardle, Nancy Marinelli, Doug Berry

# Major cleanup from limit-plotter-complete.py
#-------------------------------------------------------------------------
# UserInput
# Plotting in MakeLimitPlot(MG)
from optparse import OptionParser
parser=OptionParser()
parser.add_option("-M","--Method",dest="Method")
parser.add_option("-r","--doRatio",action="store_true")
parser.add_option("-s","--doSmooth",action="store_true")
parser.add_option("-6","--lepcat",action="store_true")
parser.add_option("-7","--hadcat",action="store_true")
parser.add_option("-b","--bayes",dest="bayes")
parser.add_option("-o","--outputLimits",dest="outputLimits")
parser.add_option("-e","--expectedOnly",action="store_true")
parser.add_option("-p","--path",dest="path",default="",type="str")
parser.add_option("-v","--verbose",dest="verbose",action="store_true")
parser.add_option("-a","--append",dest="append",default="",help="Append string to filename")
parser.add_option("","--sideband",dest="sideband",default=False,action="store_true")
parser.add_option("","--addline",action="append",type="str",help="add lines to the plot file.root:color:linestyle:legend entry", default = [])
parser.add_option("","--show",action="store_true")
parser.add_option("","--pval",action="store_true")
parser.add_option("","--addtxt",action="append",type="str", help="Add lines of text under CMS Preliminary",default=[])
parser.add_option("","--square",dest="square",help="Make square plots",action="store_true")
parser.add_option("","--nogrid",dest="nogrid",help="Remove grid from plots",action="store_true")
parser.add_option("","--nobox",dest="nobox",action="store_true",default=False,help="Don't draw box around text")
parser.add_option("-f","--plots_file",dest="plots_file")
parser.add_option("","--ns",action="store_true")

(options,args)=parser.parse_args()

# Standard Imports and calculators
import ROOT
import array,sys,numpy
ROOT.gROOT.ProcessLine(".x tdrstyle.cc")

ROOT.gROOT.SetBatch(True) #come dire root -b 
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPadTickY(0) # I don't want to see ticks on the right side of the plot
ROOT.gStyle.SetPadTickX(0) # I don't want to see ticks on the top side of the plot

#-------------------------------------------------------------------------
# Configuration for the Plotter
bratio=2.28e-03
sigma= []
sigma= [0.571,0.170,0.0569,0.0208,0.00809]

isObsLeTheory = [0,0,0,0,0]
isExpLeTheory = [0,0,0,0,0]
OBS=[0,0]
M=[0,0]
EXPECTED=[0,0,0,0,0]
OBSmasses = []
EXPmasses = []

OBSmassesT = [500,600,700,800,900]
EXPmassesT = [500,600,700,800,900]
epsilon = 10  # USELESS=> make this smaller than your smallest step size

EXP_plot=numpy.arange(500,900,10)#serve per fare le interpolazioni
OBS_plot=numpy.arange(500,900,10)#serve per fare le interpolazioni

for m in OBSmassesT:
        #if "%.1f"%m=="%d.0"%(m+epsilon):continue   # sigh!
    OBSmasses.append(m)
    EXPmasses.append(m)

# Plotting Styles --------------------------------------------------------
OFFSETLOW=0
OFFSETHIGH=0
FONTSIZE=0.034
FILLSTYLE=1001
SMFILLSTYLE=3244
FILLCOLOR_95=ROOT.kYellow
FILLCOLOR_68=ROOT.kGreen
RANGEYABS=[0.000005,0.011] #range giusto
RANGEYRAT=[0.01,30] #era 5.5
RANGEMU = [-4,3.0]
MINPV = 1.0*10E-8
MAXPV = 1.0
Lines = [1.,2.,3.,4.,5.]
MINMH=int(min(EXPmasses))
MAXMH=int(max(EXPmasses))

if(str(options.plots_file)!="None"):
    outputfile = open(options.plots_file, 'w')
if options.show : ROOT.gROOT.SetBatch(False)
if options.addline and not options.pval : sys.exit("Cannot addlines unless running in pvalue")

# ------------------------------------------------------------------------
# SM Signal Normalizer
#if not options.doRatio:
#    ROOT.gROOT.ProcessLine(".L ../libLoopAll.so")
#    signalNormalizer = ROOT.Normalization_8TeV()
ROOT.gROOT.ProcessLine(".L ../libLoopAll.so")
signalNormalizer = ROOT.Normalization_8TeV()
extraString = "SM"
# ------------------------------------------------------------------------
if options.pval:
     EXPmasses=[]
     options.doRatio=True

if options.Method=="MaxLikelihoodFit": 
    options.doRatio = True

if not options.doRatio and options.Method != "Frequentist": 
    ROOT.gROOT.ProcessLine(".L medianCalc.C+g")
    from ROOT import medianCalc
    from ROOT import FrequentistLimits
    

if options.bayes:
  BayesianFile =ROOT.TFile(options.bayes) 
  bayesObs = BayesianFile.Get("observed")
 
Method = options.Method
if not options.path: options.path=Method


EXPName = options.path+"/higgsCombineEXPECTED."+Method
if Method =="MaxLikelihoodFit": EXPName = options.path+"/higgsCombineTest."+Method
#if Method == "Asymptotic" or Method == "AsymptoticNew":  EXPName = options.path+"/higgsCombineTprimeM."+Method  # everyhting contained here
if Method == "Asymptotic" or Method == "AsymptoticNew":  EXPName = options.path+"/higgsCombineTprimeM"  # everyhting contained here
if options.lepcat:
    EXPName = options.path+"/higgsCombineTprime_cat6M"
if options.hadcat:
    EXPName = options.path+"/higgsCombineTprime_cat7M"
print EXPName

if Method == "ProfileLikelihood" or Method=="Asymptotic" or Method=="AsymptoticNew" or Method=="MaxLikelihoodFit":
#  OBSName = options.path+"/higgsCombineTest."+Method
   OBSName = options.path+"/higgsCombineTprimeM"
   if options.lepcat:
       OBSName = options.path+"/higgsCombineTprime_cat6M"
   if options.hadcat:
       OBSName = options.path+"/higgsCombineTprime_cat7M"
if Method == "Bayesian":
  OBSName = options.path+"/higgsCombineOBSERVED.MarkovChainMC"
if Method == "HybridNew":
  OBSName = options.path+"/higgsCombineOBSERVED.HybridNew"


if Method == "HybridNew" or Method == "Asymptotic" or Method == "AsymptoticNew": EXPmasses = OBSmasses[:]
# Done with setip
# --------------------------------------------------------------------------
ROOT.gROOT.ProcessLine( \
   "struct Entry{   \
    double r;       \
   };"
)
from ROOT import Entry
def getOBSERVED(file,entry=0):
  try:
   tree = file.Get("limit")
  except:
   return -1
  br = tree.GetBranch("limit")
  c = Entry()
  br.SetAddress(ROOT.AddressOf(c,'r'))
  tree.GetEntry(entry)
  return c.r

if Method=="HybridNew":
  EXPfiles=[]
  EXPmasses = OBSmasses[:]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon):    # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.quant0.500.root"%m))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.quant0.500.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()
#MIA modifica#  
elif Method=="Asymptotic" or Method=="AsymptoticNew" or Method=="MaxLikelihoodFit":
  EXPfiles=[]
  EXPmasses = OBSmasses[:]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      print EXPNAME+"%d."+Method+".mH125.6.root"
      EXPfiles.append(ROOT.TFile(EXPName+"%d."+Method+".mH125.6.root"%(m+epsilon)))
      print EXPName+"%d."+Method+".mH125.6.root"%(m+epsilon)
      print "MIO cout"
    else:
#      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.root"%m))
      EXPfiles.append(ROOT.TFile(EXPName+"%.0f."%m+Method+".mH125.6.root"))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()
#MIA modifica#
else:
  EXPfiles=[]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.root"%(m+epsilon)))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()

# Get the observed limits - Currently only does up to 1 decimal mass points
OBSfiles = []
if not options.expectedOnly:#ovvero se non giro con l'opzione -e
  for m in OBSmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      OBSfiles.append(ROOT.TFile(OBSName+"%d."+Method+".mH125.6.root"%(m+epsilon)))
    else:
      OBSfiles.append(ROOT.TFile(OBSName+"%.0f."%m+Method+".mH125.6.root"))#questa cosa funziona
    if options.verbose: print "observed MH - ", m, "File - ", OBSfiles[-1].GetName()
#prendere i limiti osservati!
  if Method == "Asymptotic" or Method =="AsymptoticNew" :  obs = [getOBSERVED(O,5) for O in OBSfiles] # observed is last entry in these files
  #print "OSSERVATO =0"
  #print obs[0]
  else: obs = [getOBSERVED(O) for O in OBSfiles]

else:
  obs = [0 for O in OBSmasses]
  OBSfiles = obs[:]

#print "OSSERVATO =0"
#print obs[1]
# -------------------------------------------------------------------------------------------------------------------------------------------
# Set-up the GRAPHS

graph68  = ROOT.TGraphAsymmErrors()
graph95  = ROOT.TGraphAsymmErrors()
graphMed = ROOT.TGraphAsymmErrors()
graphObs = ROOT.TGraphAsymmErrors()
graphOne = ROOT.TGraphAsymmErrors()
graphTheo = ROOT.TGraphAsymmErrors()
dummyGraph= ROOT.TGraphAsymmErrors()

graph68up = ROOT.TGraphErrors()
graph68dn = ROOT.TGraphErrors()
graph95up = ROOT.TGraphErrors()
graph95dn = ROOT.TGraphErrors()
graphmede = ROOT.TGraphErrors()

graphSigma = ROOT.TGraphErrors()

graph68.SetLineColor(1)
graph95.SetLineColor(1)
graph68.SetLineStyle(2)
graph95.SetLineStyle(2)
graph68.SetLineWidth(2)
graph95.SetLineWidth(2)


MG = ROOT.TMultiGraph()

def MakeMlfPlot(MG):
    legend=ROOT.TLegend(0.55,0.19,0.89,0.3)
    legend.SetFillColor(10)
    legend.SetTextFont(42)
    #legend.SetTextSize(FONTSIZE)
    graph68.SetLineStyle(1)
    legend.AddEntry(graph68,"#pm 1#sigma Uncertainty","F")

    if options.square : c = ROOT.TCanvas("c","c",600,600)
    else :c = ROOT.TCanvas("c","c",800,600)

    dhist = ROOT.TH1F("dh","dh",100,MINMH,MAXMH)
    dhist.GetYaxis().SetTitleOffset(1.2)
    dhist.GetXaxis().SetTitleOffset(1.2)
    dhist.GetYaxis().SetTitleSize(0.04)
    dhist.GetXaxis().SetTitleSize(0.04)
    dhist.GetYaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetRangeUser(MINMH,MAXMH)
    dhist.GetYaxis().SetRangeUser(RANGEMU[0],RANGEMU[1])
    dhist.GetXaxis().SetTitle("m_{T} (GeV)")
    dhist.GetYaxis().SetTitle("Best fit #sigma/#sigma_{SM}")
    dhist.Draw("AXIS")

    MG.Draw("L3") 

    # ------------------------------------------------------------------------
    # Additional Lines stored in --addline -----------------------------------
    for lineF in options.addline:

        # Parse the string, should be file.root:color:linestyle:legend entry    
        vals = lineF.split(":")
        ftmp = ROOT.TFile(vals[0])
        grext = ftmp.Get("observed")
        grext.SetLineColor(int(vals[1]))
        grext.SetLineStyle(int(vals[2]))
        grext.SetLineWidth(2)
        legend.AddEntry(grext,vals[3],"L")
        grext.Draw("same")
    # ------------------------------------------------------------------------
        

    c.Update()
    text = ROOT.TLatex()
    text.SetTextColor(ROOT.kRed)
    text.SetTextSize(FONTSIZE)
    text.SetTextFont(42)

    graphOne.Draw("L")
    c.SetGrid(not options.nogrid)
    if not options.nogrid: dhist.Draw("AXIGSAME")

    mytext= ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)
    mytext.SetTextFont(42)
    mytext.SetNDC()
    
    
    box = ROOT.TPave(0.19,0.17,0.42,0.3,2,"NDC")
    box.SetLineColor(1)
    box.SetFillColor(0)
    box.SetShadowColor(0)
    if not options.nobox: box.Draw()
    mytext.DrawLatex(0.2,0.26,"CMS Preliminary")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.2,0.25-(t+1)*0.04,"%s"%(lineT))
    legend.Draw()
    ROOT.gPad.RedrawAxis();
    
    if options.show:raw_input("Looks Ok?")
    c.SaveAs("maxlhplot.pdf")
    c.SaveAs("maxlhplot.png")

#-------------------------------------------------------------------------
def MakePvalPlot(MG):

    legend=ROOT.TLegend(0.55,0.17,0.89,0.35)
    legend.SetFillColor(10)
    legend.SetTextFont(42)
    #legend.SetTextSize(FONTSIZE)
    if not options.expectedOnly: legend.AddEntry(graphObs,"Observed","L")

    if options.square : c = ROOT.TCanvas("c","c",600,600)
    else :c = ROOT.TCanvas("c","c",800,600)

    dhist = ROOT.TH1F("dh","dh",100,MINMH,MAXMH)
    dhist.GetYaxis().SetTitleOffset(1.5)
    dhist.GetXaxis().SetTitleOffset(1.2)
    dhist.GetYaxis().SetTitleSize(0.04)
    dhist.GetXaxis().SetTitleSize(0.04)
    dhist.GetYaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetRangeUser(MINMH,MAXMH)
    dhist.GetYaxis().SetRangeUser(MINPV,MAXPV)
    dhist.GetXaxis().SetTitle("m_{T} (GeV)")
    dhist.GetYaxis().SetTitle("Local p-value")
    dhist.Draw("AXIS")

    MG.Draw("L3")

    # ------------------------------------------------------------------------
    # Additional Lines stored in --addline -----------------------------------
    for lineF in options.addline:

        # Parse the string, should be file.root:color:linestyle:legend entry    
        vals = lineF.split(":")
        ftmp = ROOT.TFile(vals[0])
        grext = ftmp.Get("observed")
        grext.SetLineColor(int(vals[1]))
        grext.SetLineStyle(int(vals[2]))
        grext.SetLineWidth(2)
        legend.AddEntry(grext,vals[3],"L")
        grext.Draw("same")
    # ------------------------------------------------------------------------
        

    c.Update()
    text = ROOT.TLatex()
    text.SetTextColor(ROOT.kRed)
    text.SetTextSize(FONTSIZE)
    text.SetTextFont(42)
    
        
    Vals=[ROOT.RooStats.SignificanceToPValue(L) for L in Lines]
    TLines=[ROOT.TLine(MINMH,V,MAXMH,V) for V in Vals]

    for j,TL in enumerate(TLines):
        TL.SetLineStyle(1)
        TL.SetLineColor(2)
        TL.SetLineWidth(1)
        TL.Draw("same")
        text.DrawLatex(MAXMH+0.2,Vals[j]*0.88,"%d #sigma"%Lines[j])

    c.SetGrid(not options.nogrid)
    c.SetLogy();
    if not options.nogrid: dhist.Draw("AXIGSAME")

    mytext= ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)
    mytext.SetTextFont(42)
    mytext.SetNDC()

    box = ROOT.TPave(0.19,0.17,0.42,0.3,2,"NDC")
    box.SetLineColor(1)
    box.SetFillColor(0)
    box.SetShadowColor(0)
    if not options.nobox: box.Draw()
    mytext.DrawLatex(0.2,0.26,"CMS Preliminary")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.2,0.25-(t+1)*0.04,"%s"%(lineT))
    legend.Draw()
    ROOT.gPad.RedrawAxis();
    
    if options.show:raw_input("Looks Ok?")
    c.SaveAs("pvaluesplot.pdf")
    c.SaveAs("pvaluesplot.png")
#-------------------------------------------------------------------------

def MakeLimitPlot(MG):
    ROOT.gPad.SetLogy()
    #if not options.doRatio: leg=ROOT.TLegend(0.33,0.56,0.89,0.90) #mette la legenda sopra
    #if options.doRatio: leg=ROOT.TLegend(0.53,0.15,0.89,0.25) #mette la legenda sotto
    leg=ROOT.TLegend(0.53,0.20,0.89,0.40)
    leg.SetFillColor(10)

    # Different entries for the different methods
    LegendEntry = ""
    if Method == "ProfileLikelihood": LegendEntry = "PL"
    if Method == "Bayesian": LegendEntry = "Bayesian"
    if Method == "HybridNew": LegendEntry = "CLs"
    if Method == "Asymptotic": LegendEntry = "CLs (Asymptotic)"
    if Method == "AsymptoticNew": LegendEntry = "CLs (Asymptotic)"

    if not options.expectedOnly: leg.AddEntry(graphObs,"Observed","L")
    if options.bayes and not options.expectedOnly: leg.AddEntry(bayesObs,"Observed Bayesian Limit","L")
    leg.SetFillStyle(0);#the legend is trasparent
    leg.AddEntry(graph68,"Expected #pm 1#sigma","FL")
    leg.AddEntry(graph95,"Expected #pm 2#sigma","FL")
    if not options.doRatio: leg.AddEntry(graphOne,"\sigma(T#bar{T})_{theo} #times 2BR(H#rightarrow #gamma #gamma)","L")
    #if not options.doRatio: leg.AddEntry(graphTheo,"Theoretical \sigma(TT) #times 2BR(H#rightarrow #gamma #gamma)","L") #se usi la scala log usa graphTheo

    leg.SetTextFont(42)
    #leg.SetTextSize(FONTSIZE)

    if options.square : C = ROOT.TCanvas("c","c",600,600)
    else: C = ROOT.TCanvas("c","c",700,600)

    C.SetGrid(not options.nogrid)
    dummyHist = ROOT.TH1D("dummy","",1,min(OBSmasses)-OFFSETLOW,max(OBSmasses)+OFFSETHIGH)
    dummyHist.SetTitleSize(0.04,"XY")
    dummyHist.Draw("AXIS")
    MG.Draw("L3")
    if not options.nogrid: dummyHist.Draw("AXIGSAME")
#TITOLI SI MODIFICANO QUI
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetLeftMargin(0.15)# lascia spazio a sinistra
    dummyHist.GetXaxis().SetTitle("m_{T} (GeV)")
    dummyHist.GetXaxis().SetRangeUser(min(OBSmasses)-OFFSETLOW,max(OBSmasses)+OFFSETHIGH)

    if options.doRatio:
     dummyHist.GetYaxis().SetRangeUser(RANGEYRAT[0],RANGEYRAT[1])
     dummyHist.GetYaxis().SetNdivisions(5,int("%d"%(RANGEYRAT[1]-RANGEYRAT[0])),0)
     dummyHist.GetYaxis().SetNdivisions(510)
     dummyHist.GetYaxis().SetTitle("\sigma(H#rightarrow #gamma #gamma)_{95%%CL} / \sigma(H#rightarrow #gamma #gamma)_{%s}"%extraString)
     dummyHist.GetYaxis().SetTitle("\sigma(T#bar{T})_{95%CL} / \sigma(T#bar{T})_{theo}")
    else: 
     dummyHist.GetYaxis().SetRangeUser(RANGEYABS[0],RANGEYABS[1])
     dummyHist.GetYaxis().SetNdivisions(5,int("%d"%(RANGEYABS[1]-RANGEYABS[0])),0)
     dummyHist.GetYaxis().SetTitle("\sigma(T#bar{T}) #times 2BR(H#rightarrow #gamma #gamma)_{95%CL} (pb)")

    dummyHist.GetYaxis().SetTitleOffset(2)
    dummyHist.GetXaxis().SetTitleOffset(1.25)

    MG.SetTitle("")
    mytext = ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)

    mytext.SetNDC()
    mytext.SetTextFont(42)
    mytext.SetTextSize(FONTSIZE)
    
    box = ROOT.TPave(0.19,0.76,0.44,0.89,2,"NDC")
    box.SetLineColor(1)
    box.SetFillColor(0)
    box.SetShadowColor(0)
#    if not options.nobox: box.Draw()
    mytext.DrawLatex(0.2,0.96,"CMS Preliminary #sqrt{s}=8 TeV L=19.7 fb^{-1}")#0.2,0.9
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.2,0.84-(t+1)*(0.04),"%s"%lineT)

    leg.SetBorderSize(0);
    leg.SetTextSize(FONTSIZE);
    leg.Draw()
    ROOT.gPad.RedrawAxis();
    
    dummyHist.GetMaximum();
    line =ROOT.TLine(500,dummyHist.GetMaximum(),900,dummyHist.GetMaximum());
    line.Draw("same");#Otherwise, the top of the frame is hidden below the histos
    if options.show:raw_input("Looks Ok?")

    #Make a bunch of extensions to the plots
    outputname="limit"
    if options.lepcat:
        outputname="limit_lep"
    if options.hadcat:
        outputname="limit_had"
    if options.doSmooth: outputname+="_smooth"
    if options.nogrid: outputname+="_nogrid"    
    outputname+="_"+Method
    if options.doRatio: outputname+="_ratio"
    if options.append!="": outputname+="_"+options.append
    types=[".pdf",".png",".gif",".eps",".ps",".C"]
    if not options.ns:
        for type in types: C.SaveAs("LIMITS_rereco/new/"+outputname+type)

#-------------------------------------------------------------------------


#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#EXPECTED + Bands => QUI DEVO LAVORARE
for i,mass,f in zip(range(len(EXPfiles)),EXPmasses,EXPfiles):
  if options.pval: continue
  sm = 1.
  print mass
  median = array.array('d',[0])
  up68   = array.array('d',[0])
  dn68   = array.array('d',[0])
  up95   = array.array('d',[0])
  dn95   = array.array('d',[0])

#sm sta per sezione d'urto teorica. quindi nel mio caso sigma(TT) per tutti i vari BR
#double Normalization_8TeV::GetXsection(double mass, TString HistName) => mass qui e' quella dell'Higgs, HistName e' il processo
# if not options.doRatio: sm = 2*signalNormalizer.GetBR(125)*signalNormalizer.GetXsection(120,"TprimeM"+"%d"%mass)
  if not options.doRatio: sm = 2*bratio*sigma[i]
  print ""
  print "massa GIUSTO"
  print mass
  print "SEZIONE URTO"
#  if not options.doRatio: print 2*signalNormalizer.GetXsection(120,"TprimeM"+"%d"%mass)
  if not options.doRatio: print 2*sigma[i]
  print "sigma*2BR(H->gg)"
  print sm


  if Method == "Asymptotic" or Method=="AsymptoticNew":   
      median[0] = getOBSERVED(f,2)
      up95[0]   = getOBSERVED(f,4)
      dn95[0]   = getOBSERVED(f,0)
      up68[0]   = getOBSERVED(f,3)
      dn68[0]   = getOBSERVED(f,1)

  elif Method=="MaxLikelihoodFit":
      median[0] = getOBSERVED(f,0)
      up95[0]   = median[0]
      dn95[0]   = median[0]
      up68[0]   = getOBSERVED(f,2)
      dn68[0]   = getOBSERVED(f,1)

  else:
    tree = f.Get("limit")
    medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95)

  graph68.SetPoint(i,float(mass),median[0]*sm)
  graph95.SetPoint(i,float(mass),median[0]*sm)
  graphMed.SetPoint(i,float(mass),median[0]*sm)
  graphOne.SetPoint(i,float(mass),1.*sm)#la teoria
  print "median[0] giusto"
  #print median[0]
  #cambio mio
  EXPECTED[i]=median[0]
  isExpLeTheory[i]=(median[0]<1.);
  print EXPECTED[i]
  print obs[i]
  print "up68"
  print up68[0]
  print "dn68"
  print dn68[0]
  print "up95"
  print up95[0]
  print "dn95"
  print dn95[0]
  if(str(options.plots_file)!="None"):
      outputfile.write(str(mass))
      outputfile.write("  ")
      outputfile.write(str(obs[i]))
      outputfile.write("  ")
      outputfile.write(str(dn95[0]))
      outputfile.write("  ")
      outputfile.write(str(dn68[0]))
      outputfile.write("  ")
      outputfile.write(str(median[0]))
      outputfile.write("  ")
      outputfile.write(str(up68[0]))
      outputfile.write("  ")
      outputfile.write(str(up95[0]))
      outputfile.write("\n")
      
  if Method == "HybridNew":

      up95[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.975.root"))
      dn95[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.025.root"))
      up68[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.840.root"))
      dn68[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.160.root"))

  
  diff95_up = abs(median[0] - up95[0])*sm
  diff95_dn = abs(median[0] - dn95[0])*sm
  diff68_up = abs(median[0] - up68[0])*sm
  diff68_dn = abs(median[0] - dn68[0])*sm
  
  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)
  graphMed.SetPointError(i,0,0,0,0)
  graphOne.SetPointError(i,0,0,0,0)

  newCanvas = ROOT.TCanvas()
  graphOne.Draw("AP")
  graphMed.Draw("PSame")
  graph68.Draw("PSame")
  newCanvas.SaveAs("LIMITS_rereco/Point_Test.pdf")
  
  if options.doSmooth: #*options.doRatio:  # fit the limits
    sm = 1.
#    if not options.doRatio:
#        sm=signalNormalizer.GetBR(125)*signalNormalizer.GetXsection(120,"TprimeM"+"%d"%mass)       

    graphmede.SetPoint(i,float(mass),median[0]*sm)
    graph68up.SetPoint(i,float(mass),up68[0]*sm)
    graph68dn.SetPoint(i,float(mass),dn68[0]*sm)
    graph95up.SetPoint(i,float(mass),up95[0]*sm)
    graph95dn.SetPoint(i,float(mass),dn95[0]*sm)
    graphSigma.SetPoint(i,float(mass),signalNormalizer.GetXsection(120,"TprimeM"+"%d"%mass))

#    graphSigma.SetPointError(i,0,0,0,0)
    graphSigma.SetPoint(i+1,600,0.170)
    graphSigma.SetPoint(i+2,650,0.097)
    graphSigma.SetPoint(i+3,625,0.12)
    graphSigma.SetPoint(i+4,675,0.068)
    graphSigma.SetPoint(i+5,600.001,0.170)
#    graphSigma.SetPoint(i+5,700.001,0.0569)

# Smooth the Bands set -doSmooth
# Since i always fitted to the Absolute, need to see if i want the Ratio instead
if options.doSmooth:
 fitstring = "[0] + [1]*x*x + [2]*x*x*x +[3]*x*x*x*x + [4]*x"
 if options.lepcat:
     fitstring = "[0] +[1]*x*x +[2]*x*x*x*x + [3]*x"
     fitstringlep2 = "(x>=400)*(x<=450)*([4]+[5]*x)"
 fitstringSigma = "(x<=600)*([0] +[1]*x*x +[2]*x*x*x*x+[3]*x*x*x) + (x>=700)*(exp(-[4]*x*x)) +(x>600)*(x<700)*([5]+[6]*x*x+[7]*x*x*x)"
 fitstringSigma2 = "(x>=550)*(x<=600)*([5]+[6]*x*x+[7]*x*x*x)"
#Per far funzionare questo bastardo ho fatto 2 fit e poi quando creo il grafico estraggo secondo una funzione di fit o l'altra a seconda del rande di massa


# medfunc = ROOT.TF1("medfunc",fitstring,109.,150.);
# up68func = ROOT.TF1("up68func",fitstring,109.,150.);
# dn68func = ROOT.TF1("dn68func",fitstring,109.,150.);
# up95func = ROOT.TF1("up95func",fitstring,109.,150.);
# dn95func = ROOT.TF1("dn95func",fitstring,109.,150.);
 medfunc = ROOT.TF1("medfunc",fitstring,399.,902.);
 up68func = ROOT.TF1("up68func",fitstring,399.,902.);
 dn68func = ROOT.TF1("dn68func",fitstring,399.,902.);
 up95func = ROOT.TF1("up95func",fitstring,399.,902.);
 dn95func = ROOT.TF1("dn95func",fitstring,399.,902.);

# medfunclep2 = ROOT.TF1("medfunclep2",fitstringlep2,399.,451.);
# up68funclep2 = ROOT.TF1("up68funclep2",fitstringlep2,399.,451.);
# dn68funclep2 = ROOT.TF1("dn68funclep2",fitstringlep2,399.,451.);
# up95funclep2 = ROOT.TF1("up95funclep2",fitstringlep2,399.,451.);
# dn95funclep2 = ROOT.TF1("dn95funclep2",fitstringlep2,399.,451.);

 Sigmafunc = ROOT.TF1("Sigmafunc",fitstringSigma,399.,902);
 Sigmafunc2 = ROOT.TF1("Sigmafunc2",fitstringSigma2,550.,602);

 graphmede.Fit(medfunc,"R,M,EX0","Q")#va scommentato
 graph68up.Fit(up68func,"R,M,EX0","Q")#va scommentato
 graph68dn.Fit(dn68func,"R,M,EX0","Q")#va scommentato
 graph95up.Fit(up95func,"R,M,EX0","Q")#va scommentato
 graph95dn.Fit(dn95func,"R,M,EX0","Q")#va scommentato

 graphmede.Fit(medfunclep2,"R,M,EX0","Q")#va scommentato
 graph68up.Fit(up68funclep2,"R,M,EX0","Q")#va scommentato
 graph68dn.Fit(dn68funclep2,"R,M,EX0","Q")#va scommentato
 graph95up.Fit(up95funclep2,"R,M,EX0","Q")#va scommentato
 graph95dn.Fit(dn95funclep2,"R,M,EX0","Q")#va scommentato

 graphSigma.Fit(Sigmafunc,"R,M,EX0","Q")
 graphSigma.Fit(Sigmafunc2,"R,M,EX0","Q")


 newCanvas = ROOT.TCanvas()
 graphmede.SetMarkerSize(0.8)
 graphmede.SetMarkerStyle(21)
# graph95up.GetYaxis().SetRangeUser(graph95dn.GetMinimum()-0.2*graph95dn.GetMinimum(),graph95up.GetMaximum()+0.2*graph95up.GetMaximum())
# graph95up.GetYaxis().SetRangeUser(0.0002,0.005)#Se sistemi i range si riesce a far vedere tutto il grafico!
 graph95up.Draw("AP")
 graph68up.Draw("Psame")#senza P non disegna niente il deficiente!
 graphmede.Draw("Psame")#senza P non disegna niente il deficiente!
 graph68dn.Draw("Psame")
 graph95dn.Draw("Psame")
 newCanvas.SaveAs("LIMITS_rereco/smoothTest.pdf")
 graphSigma.Draw("AP")
 newCanvas.SaveAs("LIMITS_rereco/Sigma_Test.pdf")
 for i,mass in zip(range(len(EXPmasses)),EXPmasses):
# for i,mass in zip(range(len(EXP_plot)),EXP_plot):
  sm=1.0
  if not options.doRatio:
#    sm = 2*signalNormalizer.GetBR(125)*signalNormalizer.GetXsection(120,"TprimeM"+"%d"%mass)
    sm = 2*bratio*sigma[i]
#    sm =2*signalNormalizer.GetBR(125)*Sigmafunc.Eval(mass)
    #if(mass>=550)*(mass<=600):
     #   sm =2*signalNormalizer.GetBR(125)*Sigmafunc2.Eval(mass)
  #if options.lepcat:
  #    if (mass<=450)*(mass>=400):
  #          mediansmooth = medfunclep2.Eval(mass)
  #          graphMed.SetPoint(i,mass,mediansmooth*sm)
  #          graph68.SetPoint(i,mass,mediansmooth*sm)
  #          graph95.SetPoint(i,mass,mediansmooth*sm)
  #          graphTheo.SetPoint(i,mass,sm)
  
  #          diff95_up = abs(mediansmooth - up95funclep2.Eval(mass))*sm
  #          diff95_dn = abs(mediansmooth - dn95funclep2.Eval(mass))*sm
  #          diff68_up = abs(mediansmooth - up68funclep2.Eval(mass))*sm
  #          diff68_dn = abs(mediansmooth - dn68funclep2.Eval(mass))*sm
  
  #else:
  print "******************"
  print "******************"
  print options.doSmooth
  mediansmooth = medfunc.Eval(mass)
  graphMed.SetPoint(i,mass,mediansmooth*sm)
  graph68.SetPoint(i,mass,mediansmooth*sm)
  graph95.SetPoint(i,mass,mediansmooth*sm)
  graphTheo.SetPoint(i,mass,sm)
  
  diff95_up = abs(mediansmooth - up95func.Eval(mass))*sm
  diff95_dn = abs(mediansmooth - dn95func.Eval(mass))*sm
  diff68_up = abs(mediansmooth - up68func.Eval(mass))*sm
  diff68_dn = abs(mediansmooth - dn68func.Eval(mass))*sm
  
  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)

# OBSERVED -------- easy as that !
for i,mass in zip(range(len(OBSfiles)),OBSmasses):

    sm = 1.;
    MASSSS=125
    if obs[i] ==-1: continue
#    if not options.doRatio: sm = 2*signalNormalizer.GetBR(125)*signalNormalizer.GetXsection(120,"TprimeM"+"%d"%mass)
    if not options.doRatio: sm = 2*bratio*sigma[i]
    graphObs.SetPoint(i,float(mass),obs[i]*sm)
    print "OSSERVATO"
    print "massa"
    print mass
    print "osservato"
    print obs[i]
    #obs[i] e' l'upper limit osservato
    print "teoria"
    print sm
    print "punto di massa"
    print i
    isObsLeTheory[i]=(obs[i]<1.)
    print "BR check"
    print bratio
    print "Sigma check"
    print sigma[i]
    graphObs.SetPointError(i,0,0,0,0)

print "*****************"
#for j in zip(range(len(OBSfiles))-1):
if not options.expectedOnly:
    for j in 0,1,2,3:
        #print j
        if (isObsLeTheory[j]==1 and isObsLeTheory[j+1]==0):
            M[0]=OBSmasses[j]
            M[1]=OBSmasses[j+1]
            OBS[0]=obs[j]
            OBS[1]=obs[j+1]
            #print "OBS[0] ="
            #print OBS[0]
            #print "OBS[1] ="
            #print OBS[1]
            #print "M[0] ="
            #print M[0]
            #print "M[1] ="
            #print M[1]
    a=0.0
    b=0.0
    x=0.0
    #se la massa esclusa e' maggiore di 500
    if(M[0]!=0):
        a=(OBS[0]-OBS[1])/(M[0]-M[1])
        b=(OBS[0]+OBS[1]-a*(M[0]+M[1]))/2
        x=(1-b)/a
        
    print "(OSSERVATO) Regione esclusa fino a MT"
    print x

if options.expectedOnly:
    for j in 0,1,2,3:
        #print j
        if (isExpLeTheory[j]==1 and isExpLeTheory[j+1]==0):
            M[0]=OBSmasses[j]
            M[1]=OBSmasses[j+1]
            #si chiama OBS, ma come vedi e' Expected
            OBS[0]=EXPECTED[j]
            OBS[1]=EXPECTED[j+1]
            #print "OBS[0] ="
            #print OBS[0]
            #print "OBS[1] ="
            #print OBS[1]
            #print "M[0] ="
            #print M[0]
            #print "M[1] ="
            #print M[1]
    a=0.0
    b=0.0
    x=0.0
    if(M[0]!=0):
        a=(OBS[0]-OBS[1])/(M[0]-M[1])
        b=(OBS[0]+OBS[1]-a*(M[0]+M[1]))/2
        x=(1-b)/a
        
    print "(ATTESO) Regione esclusa fino a MT"
    print x



# Finally setup the graphs and plot
graph95.SetFillColor(FILLCOLOR_95)
graph95.SetFillStyle(FILLSTYLE)
graph68.SetFillColor(FILLCOLOR_68)
graph68.SetFillStyle(FILLSTYLE)
graphMed.SetLineStyle(2)
graphMed.SetLineColor(2)
graphMed.SetMarkerColor(2)
graphMed.SetLineWidth(3)
graphObs.SetLineWidth(2)


if options.bayes:
 bayesObs.SetLineWidth(3)
 bayesObs.SetLineColor(4)
 bayesObs.SetMarkerColor(4)
 bayesObs.SetLineStyle(7)

graphOne.SetLineWidth(3) #Era cosi'
graphOne.SetLineColor(ROOT.kRed)
graphOne.SetMarkerColor(ROOT.kRed)

graphTheo.SetLineWidth(3)
graphTheo.SetLineColor(ROOT.kRed)
graphTheo.SetMarkerColor(ROOT.kRed)
graphObs.SetMarkerStyle(20)
graphObs.SetMarkerSize(1.0) #era 2.0
graphObs.SetLineColor(1)

graphMed.SetLineStyle(2)
graphMed.SetLineColor(ROOT.kBlack)
if not options.pval:MG.Add(graph95)
if not options.pval:MG.Add(graph68)
if not options.pval:MG.Add(graphMed)

if not options.expectedOnly:
  MG.Add(graphObs)
  if options.bayes:
   MG.Add(bayesObs)

if not options.pval: MG.Add(graphOne)
#if not options.pval: MG.Add(graphTheo) usa graphTheo se non usi la scala logaritmica

# Plot -------------------------------------
if options.pval: MakePvalPlot(MG)
elif Method=="MaxLikelihoodFit":  MakeMlfPlot(MG)
else:MakeLimitPlot(MG)
#MakeLimitPlot(MG)
# ------------------------------------------
if options.outputLimits:
  print "Writing Limits To ROOT file --> ",options.outputLimits
  OUTTgraphs = ROOT.TFile(options.outputLimits,"RECREATE")
  graphObs.SetName("observed")
  graphObs.Write()

  if not options.pval:
   graphMed.SetName("median")
   graphMed.Write()
   graph68.SetName("sig1")
   graph68.Write()
   graph95.SetName("sig2")
   graph95.Write()

  OUTTgraphs.Write()

