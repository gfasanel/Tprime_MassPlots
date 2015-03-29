#include <iostream> 
#include <fstream>
#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include "TGraph2D.h"
#include <TH1F.h>
#include <TString.h>
#include <TCut.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TFile.h>
#include "TH2F.h"
#include "TLatex.h"
//#include "TColor.h"

void make2Dplot(int cat){
 int NMAX=10;
 ifstream file;
 //if(cat==4) file.open("limitsResults_COMB.txt");
 //if(cat==5) file.open("limitsResults_COMB_obs.txt");
 if(cat==4) file.open("UL_Tprime_expected.txt");
 if(cat==5) file.open("UL_Tprime_observed.txt");
 int i = 0;

 TGraph2D* gr= new TGraph2D();
 
 while (!file.eof()) {
    	  double mass;
	  double limit;
	  double width;
	  file >> mass >> limit >> width;
	  std::cout<<mass<<" "<<width<<" "<<limit<<std::endl;
	  gr->SetPoint(i, mass, width, limit);
	  i++;
 }
 
 cout<<"Tutti gli ul presi"<<endl;

 TCanvas* c = new TCanvas("c", "c", 1);
 c->cd();
 c->SetLogz();



 //gr->SetMaximum(360);//
 gr->SetMaximum(100);
 //gr->SetMinimum(0.1);//log scale
 gr->GetXaxis()->SetLimits(500,700);
 //gr->GetXaxis()->SetLimits(500,650);

 //gStyle->SetOptTitle(0);
  gr->GetXaxis()->SetTitle("m_{T} [GeV]");
  gr->GetYaxis()->SetTitle("BR(T#rightarrow tH)");
  gr->GetZaxis()->SetTitle("#sigma^{UL} (95%CL)/ #sigma^{theo}");
  gr->GetYaxis()->SetTitleOffset(0.7);
  gr->GetZaxis()->SetTitleOffset(0.7);
  gr->Draw("COLZ");

  double cont[5];
  if(cat==2)for(int i=0; i<5;i++)cont[i]=0.03+i*0.02;
  if( cat==3)for(int i=0; i<5;i++)cont[i]=0.03+i*0.02;
  if(cat==0 )for(int i=0; i<5;i++)cont[i]=0.01+i*0.01;
  if(cat==4 )for(int i=0; i<5;i++)cont[i]=0.005+i*0.005;
  if(cat==1 )for(int i=0; i<5;i++)cont[i]=0.01+i*0.01;


  TList* list1 = (TList*)gr->GetContourList(1.0);
  list1->First();
  //I draw the exclusion border
  if(cat==4){
    ifstream file_expected;
    file_expected.open("border_Tprime_expected.txt");
    int i=0;
    TGraph *contour=new TGraph();
    while (!file_expected.eof()) {
    double mass;
    double br;
    file_expected >> br>> mass;
    std::cout<<mass<<" "<<br/100<<std::endl;
    contour->SetPoint(i, mass, br/100);
    i++;
  }

  contour->SetLineColor(kBlack);  
  contour->SetLineStyle(1);
  contour->SetLineWidth(2);
  contour->Draw("same");
  TLatex *mytext = new TLatex();
  mytext->SetNDC();//without this, it doesn't work properly
  //mytext->SetTextSize(0.034);
  mytext->DrawLatex(0.2,0.9,"CMS Preliminary #sqrt{s}=8 TeV L=19.7 fb^{-1}");
  }else if(cat==5){
    ifstream file_observed;
    file_observed.open("border_Tprime_observed.txt");
    int i=0;
    TGraph *contour=new TGraph();
    while (!file_observed.eof()) {
      double mass;
      double br;
      file_observed >> br>> mass;
      std::cout<<mass<<" "<<br/100<<std::endl;
      contour->SetPoint(i, mass, br/100);
      i++;
    }

  contour->SetLineColor(kBlack);  
  contour->SetLineStyle(1);
  contour->SetLineWidth(2);
  contour->Draw("same");
  TLatex *mytext = new TLatex();
  mytext->SetNDC();//without this, it doesn't work properly
  //mytext->SetTextSize(0.034);
  mytext->DrawLatex(0.2,0.9,"CMS Preliminary #sqrt{s}=8 TeV L=19.7 fb^{-1}");

  }

  string spec;
  if(cat==4){
    spec="expected";
  }else if(cat==5){
    spec="observed";
  }

  c->SaveAs(("2D_EXCL_PLOT/limit_colz_"+spec+".png").c_str());
  c->SaveAs(("2D_EXCL_PLOT/limit_colz_"+spec+".pdf").c_str());
  c->SaveAs(("2D_EXCL_PLOT/limit_colz_"+spec+".eps").c_str());
  c->SaveAs(("2D_EXCL_PLOT/limit_colz_"+spec+".C").c_str());

  //Not Rainbow palette anymmore
  /* c->SaveAs(("2D_EXCL_PLOT/limit_colz_rainbow_"+spec+".png").c_str());
  c->SaveAs(("2D_EXCL_PLOT/limit_colz_rainbow_"+spec+".pdf").c_str());
  c->SaveAs(("2D_EXCL_PLOT/limit_colz_rainbow_"+spec+".eps").c_str());
  c->SaveAs(("2D_EXCL_PLOT/limit_colz_rainbow_"+spec+".C").c_str());
  */

  

}

 
