
#include"TPaveText.h"
#include "TChain.h"
#include "TH1F.h"
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TLatex.h"
#include <iostream>
#include <fstream>
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

using namespace std;
  

TString indirname = "./REDUCED/";
TString savedirname = "./PAS_PLOT";

double lumi = 4860.;

double fracQCD = 0.672;
double fracG = 0.328;




std::string get_sqrtText() {

  char label_sqrt_text[150];
  
  sprintf( label_sqrt_text, "#sqrt{s} = 8 TeV");
  std::string returnString(label_sqrt_text);

  return returnString;

}


TPaveText* get_labelCMS( int legendQuadrant = 0 , std::string year="2012", bool sim=false, std::string run = "ALL") {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for CMS label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.63;
    x2 = 0.8;
    y2 = 0.87;
  } else if( legendQuadrant==2 ) {
    x1 =  0.25;
    y1 = 0.83;
    x2 =  0.42;
    y2 = 0.87;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
    y2 = 0.24;
  } else if( legendQuadrant==0 ) {
    x1 = 0.125;
    y1 = 0.903;
    x2 = 0.4;
    y2 = 0.925;
  }

  
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
 
  std::string leftText;
   
  if(year == "2010")  leftText = "CMS Preliminary 2010, 34 pb^{-1}";
  if (sim)  leftText = "CMS Simulation"; //cwr ->remove 2011
  else {
    if(year == "2011" && run == "RUN2011A")  leftText = "CMS Preliminary RUN2011A 2.034 fb^{-1}";
    if(year == "2011" && run == "RUN2011B")  leftText = "CMS Preliminary 2011, 2.516 fb^{-1}";
    if(year == "2011" && run == "ALL")  leftText = "CMS 4.9 fb^{-1}"; //cwr ->remove 2011
    if(year == "2012" && run == "ALL")  leftText = "CMS Preliminary  #sqrt{s}= 8 TeV L=19.7 fb^{-1}";
    if(year == "none" && run == "ALL")  leftText = "CMS Data"; //cwr ->remove 2011
    if(year == "May2011")leftText = "CMS Preliminary 2011, 858.4 pb^{-1}";

  }
  cmslabel->AddText(leftText.c_str());
  return cmslabel;

}




TPaveText* get_labelSqrt( int legendQuadrant ) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Sqrt label. Using 2." << std::endl;
    legendQuadrant = 2;
  }


  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.78;
    x2 = 0.8;
    y2 = 0.82;
  } else if( legendQuadrant==2 ) {
    x1 = 0.25;
    y1 = 0.78;
    x2 = 0.42;
    y2 = 0.82;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.16;
    x2 = 0.42;
    y2 = 0.2;
  } else if( legendQuadrant==0 ) {
    x1 = 0.65;
    y1 = 0.953;
    x2 = 0.87;
    y2 = 0.975;
  }


  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->AddText("#sqrt{s} = 8 TeV");
  return label_sqrt;

}


void hdensity( TH1F* hout){

   
  for(int i = 1; i<=hout->GetNbinsX(); i++){
   
    double area = hout->GetBinWidth(i);
      
    hout->SetBinContent(i,(double)hout->GetBinContent(i)/area);
    hout->SetBinError(i,(double)hout->GetBinError(i)/area);

    // std::cout << "i: " << i << " j: " << j <<
    //" val: " << hout->GetBinContent(i,j) << " area: " << area <<std::endl; 

  }

}



void make2D(){

  int  kCDF = kYellow-4;
  int  kD0 = kAzure-4;
  int  kCMS = kCyan-7;

  TPaveText* label_cms = get_labelCMS(0,"2012", false, "ALL");
  TPaveText* label_sqrt = get_labelSqrt(0);
  

  //**************************observed cms*******************************//

  TGraph* observed= new TGraph();
  TGraph* observed_contour= new TGraph();
  observed_contour->SetLineWidth(2);

  //observed->SetFillColor(kCMS);//
  observed->SetFillColor(kBlack);//
  //observed->SetFillColor(kOrange);
  observed->SetLineColor(kBlack);
  observed->SetLineWidth(1);
  observed->SetFillStyle(3244);//3013
  //**************************expecteded cms*******************************//
  ifstream file_expected,file_expectedWB0,file_expectedtZ0,file_observed,file_up95,file_up68,file_dn68,file_dn95;
  file_expected.open("border_Tprime_expected.txt");// with BR(T->Wb)=BR(T->tZ)
  file_up95.open("border_Tprime_up95.txt");
  file_up68.open("border_Tprime_up68.txt");
  file_dn68.open("border_Tprime_dn68.txt");
  file_dn95.open("border_Tprime_dn95.txt");
  file_expectedWB0.open("border_Tprime_expectedWB0.txt");
  file_expectedtZ0.open("border_Tprime_expectedtZ0.txt");
  file_observed.open("border_Tprime_observed.txt");

  int i=0;
  int sorter,sorter_end=0;
  TGraph* expected= new TGraph();
  TGraph* expectedWB0= new TGraph();
  TGraph* expectedtZ0= new TGraph();

  TGraph* graph68= new TGraph();
  TGraph* graph95= new TGraph();

  while (!file_expected.eof()) {
    double mass,mass_up95,mass_up68;
    double br;
    file_expected >> br>> mass;
    file_up95 >> br>> mass_up95;
    file_up68 >> br>> mass_up68;

    //std::cout<<mass<<" "<<br/100<<std::endl;
    if(mass==0){mass=500;}
    if(mass_up95==0){mass_up95=500;}
    if(mass_up68==0){mass_up68=500;}

    expected->SetPoint(i, mass, br/100);
    graph68->SetPoint(i, mass_up68, br/100);
    //std::cout<<i<<" "<<mass_up68<<" "<<br/100<<endl;
    graph95->SetPoint(i, mass_up95, br/100);
    sorter=i+1;
    i++;
  }

  //N.B: Qui non va azzerato il contatore i
  while (!file_dn68.eof()) {
    double mass_dn68,mass_dn95;
    double br;
    file_dn68 >> br>> mass_dn68;
    file_dn95 >> br>> mass_dn95;

    if(mass_dn68==0){mass_dn68=500;}
    if(mass_dn95==0){mass_dn95=500;}

    graph68->SetPoint(i, mass_dn68, br/100);
    //std::cout<<i<<" "<<mass_dn68<<" "<<br/100<<endl;
    graph95->SetPoint(i, mass_dn95, br/100);

    sorter_end=i;
    i++;
  }

  //Reordering points, because I want to go around the area and have my sigma-bands
  for(int j=0; j<=(sorter_end-sorter)/2;j++){
    graph68->SwapPoints(sorter+j,sorter_end-j);
    graph95->SwapPoints(sorter+j,sorter_end-j);
  }

  //for(int j=0; j<=13;j++){//Print values
  //  std::cout<<"x  y"<<graph68->GetX()[j]<<" "<<graph68->GetY()[j]<<endl;
  // }


  i=0;
  while (!file_expectedWB0.eof()) {
    double mass;
    double br;
    file_expectedWB0 >> br>> mass;
    //std::cout<<mass<<" "<<br/100<<std::endl;
    if(mass==0){mass=500;}
    expectedWB0->SetPoint(i, mass, br/100);
    i++;
  }

  i=0;
  while (!file_expectedtZ0.eof()) {
    double mass;
    double br;
    file_expectedtZ0 >> br>> mass;
    //std::cout<<mass<<" "<<br/100<<std::endl;
    if(mass==0){mass=500;}
    expectedtZ0->SetPoint(i, mass, br/100);
    i++;
  }

  i=0;
  while (!file_observed.eof()) {
    double mass;
    double br;
    int shift_contour;
    file_observed >> br>> mass;
    //std::cout<<mass<<" "<<br/100<<std::endl;
    if(mass==0){mass=500;shift_contour=i;}
    observed->SetPoint(i, mass, br/100);
    if(mass!=500){
      observed_contour->SetPoint(i-shift_contour,mass,br/100);
    }
    i++;
  }
  observed->SetPoint(i, 0, 1);
  cout<<"Bordi presi"<<endl;

  expected->SetLineWidth(1);
  expected->SetLineStyle(2);

  expectedWB0->SetLineColor(kRed);
  expectedWB0->SetLineWidth(1);
  expectedWB0->SetLineStyle(2);

  expectedtZ0->SetLineColor(kBlue);
  expectedtZ0->SetLineWidth(1);
  expectedtZ0->SetLineStyle(2);

  //BANDS FOR EXPECTED//
  graph95->SetFillColor(kYellow);
  graph95->SetFillStyle(1001);
  graph68->SetFillColor(kGreen);
  graph68->SetFillStyle(1001);

  //PLOTTING//

  TCanvas* c0 = new TCanvas("exclusion limit", "exclusion limit", 1);
  c0->cd();


  TPad * pad2 = new TPad("pad2", "pad2",0.01,0.0001,0.99,0.99);

  pad2->Draw("AH");
  pad2->SetBorderSize(0.);
  pad2->cd();

  
 

  TPad * pad1 = new TPad("pad1", "pad1",0.01,0.002,0.99,0.99);
    
  pad1->Draw();
  pad1->cd();

  TH2F *h2 = new TH2F("h","Axes",100,500,700,10,0.1,1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  h2->GetXaxis()->SetTitle("m_{T} [GeV]");
  h2->GetYaxis()->SetTitle("BR(T#rightarrow tH)");
  h2->GetXaxis()->SetTitleSize(0.048);
  h2->GetYaxis()->SetTitleSize(0.048);

  ////////////
  h2->Draw("");
  graph95->Draw("f");
  graph68->Draw("fsame");
  //expected->Draw("L");//
  expectedWB0->Draw("L");
  expectedtZ0->Draw("L");
  observed->Draw("F");
  observed_contour->Draw("L");
  h2->Draw("AXISSAME");

  TLegend* leg1 = new TLegend(0.581908,0.2286401,0.8847913,0.4023175,NULL,"brNDC");
  leg1->SetFillStyle(0); leg1->SetBorderSize(0); 
  leg1->SetFillColor(0);

  TLegend* leg2 = new TLegend(0.581908,0.4086401,0.8847913,0.5023175,NULL,"brNDC");
  leg2->SetFillStyle(0); leg2->SetBorderSize(0); 
  leg2->SetFillColor(0);
  leg2->AddEntry(observed, "Observed Exclusion", "F");

  TLegend* leg3 = new TLegend(0.581908,0.1286401,0.8847913,0.223175,NULL,"brNDC");
  leg3->SetFillStyle(0); leg3->SetBorderSize(0); 
  leg3->SetFillColor(0);

  TLegend* leg4 = new TLegend(0.72,0.1286401,1.02,0.223175,NULL,"brNDC");
  leg4->SetFillStyle(0); leg4->SetBorderSize(0); 
  leg4->SetFillColor(0);

    
  label_cms->Draw("same");
  //  label_sqrt->Draw("same");
    
  TH1F h;
  h.SetMarkerSize(0);
  h.SetLineColor(kWhite);
    

  //  leg1->AddEntry(expectedWB0, "BR(T#rightarrow Wb)=0; BR(T#rightarrow tZ)= 1-BR(T#rightarrow tH)", "l");
  leg1->AddEntry(expectedWB0, "BR(T#rightarrow Wb) = 0", "l");
  //leg1->AddEntry(expected, "BR(T#rightarrow Wb)=BR(T#rightarrow tZ)= (1-BR(T#rightarrow tH))/2", "l");
  //  leg1->AddEntry(expectedtZ0, "BR(T#rightarrow tZ)=0; BR(T#rightarrow Wb)= 1-BR(T#rightarrow tH)", "l");
  leg1->AddEntry(expectedtZ0, "BR(T#rightarrow tZ) = 0", "l");
  leg3->AddEntry(graph68, "#pm 1 #sigma", "F");
  leg4->AddEntry(graph95, "#pm 2 #sigma", "F");
  leg1->SetHeader("Expected Exclusion");

  TLegendEntry *textobs = (TLegendEntry*)leg2->GetListOfPrimitives()->At(0);
  textobs->SetTextSize(0.038);
  TLegendEntry *header = (TLegendEntry*)leg1->GetListOfPrimitives()->First();
  header->SetTextSize(0.038); 
  TLegendEntry *text1 = (TLegendEntry*)leg1->GetListOfPrimitives()->At(1);
  text1->SetTextSize(0.038);
  //  text1->SetTextFont(22);
  TLegendEntry *text2 = (TLegendEntry*)leg1->GetListOfPrimitives()->At(2);
  text2->SetTextSize(0.038);



  TLegendEntry *header3 = (TLegendEntry*)leg3->GetListOfPrimitives()->First();
  header3->SetTextSize(0.038); 

  TLegendEntry *header4 = (TLegendEntry*)leg4->GetListOfPrimitives()->First();
  header4->SetTextSize(0.038); 


  leg2->Draw("same");
  leg1->Draw("same");
  leg3->Draw("same");
  leg4->Draw("same");

  c0->SaveAs("./2D_EXCL_PLOT/2D_exclusion_final.png");
  c0->SaveAs("./2D_EXCL_PLOT/2D_exclusion_final.pdf");
  c0->SaveAs("./2D_EXCL_PLOT/2D_exclusion_final.eps");
  c0->SaveAs("./2D_EXCL_PLOT/2D_exclusion_final.C");


}



