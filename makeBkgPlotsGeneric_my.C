// Makes partially blinded mass distribution + fit plots for Mass-fac MVA analysis
//Smooth Bands

#ifndef __CINT__
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "TROOT.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

#include <iostream>
#endif 

void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname);

void makeBkgPlotsGeneric_my(std::string filebkg, std::string filesig="", bool blind=true, bool doBands=true, bool baseline=true, bool useBinnedData=false){

	// Globals
	gROOT->SetStyle("Plain");
	gROOT->SetBatch(1);
	gStyle->SetOptStat(0);
	const int ncats = 8;//; 4 inclusive + 2 ttH + 2 Tprime
#
	RooMsgService::instance().setGlobalKillBelow(RooFit::MsgLevel(RooFit::FATAL));

	std::string * labels;
	
	std::string baselinelabels[ncats] = { "Inclusive0"
					  ,"Inclusive1"
					  ,"Inclusive2"
					  ,"Inclusive3"
					  ,"ttH lep"
					  ,"ttH had"
					      ,"T#bar{T}#rightarrowtH#bar{t}H(#rightarrow #gamma #gamma)"//Tprime lep
					      ,"T#bar{T}#rightarrowtH#bar{t}H(#rightarrow #gamma #gamma)"//Tprime had					     
	};
	std::string massfactlabels[9] = { //Not my case: I don't use MVA
		"BDT_{#gamma#gamma} >= 0.91"
		,"0.79  <= BDT_{#gamma#gamma} < 0.91"
		,"0.49 <= BDT_{#gamma#gamma} < 0.79"
		,"-0.05  <= BDT_{#gamma#gamma} < 0.49"
		,"Dijet-tagged class BDT_{VBF} >= 0.985"
		,"Dijet-tagged class 0.93 <= BDT_{VBF} < 0.985"
		,"Muon-tagged class"
		,"Electron-tagged class"
		,"MET-tagged class"
		
	};
	
	if( baseline ) { labels = baselinelabels; }
	else { labels = massfactlabels; }
	
	TFile *fb = TFile::Open(filebkg.c_str());
	TFile *fs = ( filesig.empty() ? fb : TFile::Open(filesig.c_str()) );
	TFile *fsGiuseppe = TFile::Open("../AnalysisScripts/Tprime_allMasses/histo_fasanel.root");	
	TFile *fstth = TFile::Open("../AnalysisScripts/Tprime_allMasses/CMS-HGG_interpolated.root");	

	RooWorkspace *w_bkg  = (RooWorkspace*) fb->Get("cms_hgg_workspace");
	RooRealVar *x = (RooRealVar*) w_bkg->var("CMS_hgg_mass");
	RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");

       	double lumi = intL->getVal()/1000.;


	TLatex *latex = new TLatex();	
	latex->SetTextSize(0.025);
	latex->SetNDC();
	
	TLatex *cmslatex = new TLatex();
	cmslatex->SetTextSize(0.04);
	cmslatex->SetNDC();

	double totalGGHinDIJET = 0;
	double totalVBFinDIJET = 0;
	double totalGGHinINCL = 0;
	double totalVBFinINCL = 0;

	double totalTTHinDIJET = 0;
	double totalWZHinDIJET = 0;
	double totalTTHinINCL = 0;
	double totalWZHinINCL = 0;

	double totalTprimeinINCL = 0;
	double totalTprimeDIJET = 0;

	//	for (int cat=0;cat<ncats;cat++){
	for (int cat=6;cat<ncats;cat++){//Only my Tprime category		
		TCanvas *can = new TCanvas("c","",800,800);
		TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);// Legend is trasparent

		// Get Dataset ->
		RooAbsData *data;
		if (useBinnedData) data =  (RooDataSet*)w_bkg->data(Form("data_mass_cat%d",cat));
		else  data =  (RooDataHist*)w_bkg->data(Form("roohist_data_mass_cat%d",cat));
		data->Print();

		// Background Pdf ->
		/// RooExtendPdf *bkg =  (RooExtendPdf*)w_bkg->pdf(Form("data_pol_model_8TeV_cat%d",cat));
		RooAbsPdf *bkg =  (RooAbsPdf*)w_bkg->pdf(Form("pdf_data_pol_model_8TeV_cat%d",cat));
		bkg->Print();
		bkg->fitTo(*data);
		RooFitResult *r = bkg->fitTo(*data,RooFit::Save(1));
		r->Print();
		
		// Get Signal pdf norms
		std::cout << "Getting Signal Components" << std::endl;
		//TH1F *tthnorm = (TH1F*)fs->Get(Form("th1f_sig_tth_mass_m125.6_cat%d",cat));//dopo che micheli ha rifittato non funziona piu'
		TH1F *tthnorm = (TH1F*)fstth->Get(Form("th1f_sig_tth_mass_m125.6_cat%d",cat));
		TH1F *Tprimenorm = (TH1F*)fsGiuseppe->Get(Form("th1f_sig_TprimeM700_mass_m125.6_cat%d",cat));

		cout<<endl<<endl<<"Tprime Integral all over the range"<<Tprimenorm->Integral()<<endl;		
		cout<<"Integral Tprime mass window"<<Tprimenorm->Integral(Tprimenorm->GetXaxis()->FindBin(123.5),Tprimenorm->GetXaxis()->FindBin(126.5))<<endl;
		cout<<"Integral tth mass window"<<tthnorm->Integral(tthnorm->GetXaxis()->FindBin(123.5),tthnorm->GetXaxis()->FindBin(126.5))<<endl;

		/*if (cat<=3){
		  //totalGGHinINCL+=gghnorm->Integral();
		  //totalVBFinINCL+=vbfnorm->Integral();
		  //totalTTHinINCL+=tthnorm->Integral();
		  //totalWZHinINCL+=wzhnorm->Integral();
		  // totalTprimeinINCL+=Tprimenorm->Integral();
		}else{
		  //totalGGHinDIJET+=gghnorm->Integral();
		  //totalVBFinDIJET+=vbfnorm->Integral();
		  //totalTTHinDIJET+=tthnorm->Integral();
		  //totalWZHinDIJET+=wzhnorm->Integral();
		  totalTprimeDIJET+=Tprimenorm->Integral();
		  }*/
		
		//std::cout << "Rescaling Signal Components" << std::endl;
		//gghnorm->Add(vbfnorm);
		//gghnorm->Add(wzhnorm);
		//gghnorm->Add(tthnorm);


		TH1F *allsig = (TH1F*)Tprimenorm->Clone();
		allsig->Rebin(2);
		allsig->SetLineColor(4);allsig->SetFillColor(38);allsig->SetFillStyle(3001) ;allsig->SetLineWidth(2);
		tthnorm->Rebin(2);
		tthnorm->SetLineColor(kRed);tthnorm->SetFillColor(kRed);tthnorm->SetFillStyle(3001) ;tthnorm->SetLineWidth(2);
		
		/// allsig->SetLineColor(1);
		/// allsig->SetFillColor(38);
		TH1F dumData("d","",80,100,180); dumData.Sumw2();dumData.SetMarkerSize(1.0);dumData.SetMarkerStyle(20);dumData.SetLineWidth(3);
		dumData.Fill(101);
		// TH1F dumSignal("s","",80,100,180); dumSignal.SetLineColor(4);dumSignal.SetFillColor(38);dumSignal.SetFillStyle(3001) ;dumSignal.SetLineWidth(2);
		TH1F dum1Sig("1s","",80,100,180); dum1Sig.SetFillColor(kYellow);dum1Sig.SetFillStyle(1001);
		TH1F dum2Sig("2s","",80,100,180); dum2Sig.SetFillColor(kGreen);dum2Sig.SetFillStyle(1001);
		TH1F dumBkg("b","",80,100,180); dumBkg.SetLineColor(kRed);dumBkg.SetLineWidth(3);
		dumBkg.Draw("P");// dumSignal.Draw("LFsame");
		dumBkg.Draw("Fsame");dum1Sig.Draw("Fsame");dum2Sig.Draw("Lsame");

		// Plot background
		RooPlot *frame = x->frame();

		std::cout << "Plotting Components" << std::endl;
		data->plotOn(frame,RooFit::Binning(80),RooFit::Invisible());
		/// bkg->plotOn(frame,RooFit::VisualizeError(*r,2,1),RooFit::FillColor(kGreen));
		/// bkg->plotOn(frame,RooFit::VisualizeError(*r,1,1),RooFit::FillColor(kYellow));
		bkg->plotOn(frame,RooFit::LineColor(kRed));
		TGraphAsymmErrors *onesigma = 0, *twosigma = 0;
		if( doBands ) {
			onesigma = new TGraphAsymmErrors();
			twosigma = new TGraphAsymmErrors();
			TString name=Form("cat%d",cat);
      doBandsFit(onesigma, twosigma, x, bkg, dynamic_cast<RooCurve*>(frame->getObject(frame->numItems()-1)), 
		 data, frame, name,cat);
		}
		if( blind ) {
		  x->setRange("unblind_up",135,180);//my blind
		  data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_up"));
		  x->setRange("unblind_down",100,115);//my blind
		  data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_down"));
		} else {
		  data->plotOn(frame,RooFit::Binning(80));
		}

		frame->SetTitle("");
		frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
		frame->GetXaxis()->SetNdivisions(5,5,0);
		frame->GetYaxis()->SetTitle("Events / (1 GeV)");
		frame->GetYaxis()->SetTitleOffset(1.2);


		leg->SetTextSize(0.034);//this changes the size of the legend; old 0.024
		leg->AddEntry(&dumData,"Data","PEL");
		leg->AddEntry(&dumBkg,"Bkg Model","L");
		leg->AddEntry(&dum1Sig,"#pm 1#sigma","F");
		leg->AddEntry(&dum2Sig,"#pm 2#sigma","F");
		leg->AddEntry(allsig,"T#bar{T}#rightarrow tH#bar{t}H","F");
		leg->AddEntry(tthnorm,"t#bar{t}H(#rightarrow #gamma#gamma)","F");	       

		
		frame->SetMinimum(0.01);
		frame->SetMaximum(frame->GetMaximum()*1.6); //good if not log scale
		//frame->SetMaximum(frame->GetMaximum()*50);//good for log scale
		frame->Draw();
 		if( doBands ) {
 			twosigma->SetLineColor(kGreen);
 			twosigma->SetFillColor(kGreen);
 			twosigma->SetMarkerColor(kGreen);
 			twosigma->Draw("L3 SAME");     
 			
 			onesigma->SetLineColor(kYellow);
 			onesigma->SetFillColor(kYellow);
 			onesigma->SetMarkerColor(kYellow);
 			onesigma->Draw("L3 SAME");
 			frame->Draw("same");
 		}
		allsig->Draw("samehistF"); //Tprime
		tthnorm->Draw("samehistF");//ttH

		leg->Draw();
		cmslatex->DrawLatex(0.14,0.92,Form("CMS Preliminary #sqrt{s} = 8TeV L = %2.1f fb^{-1}",lumi));
		latex->SetTextSize(0.034);
		latex->DrawLatex(0.15,0.85,labels[cat].c_str());
		if(cat==6){
		latex->DrawLatex(0.15,0.8,"Leptonic Category");
		latex->DrawLatex(0.15,0.75,"M_{T}=700 GeV");
		}else if(cat==7){
		latex->DrawLatex(0.15,0.8,"Hadronic Category");
		latex->DrawLatex(0.15,0.75,"M_{T}=700 GeV");
		}
		
		//do you want log?
		/*
		can->SetLogy();
		can->SaveAs(Form( (baseline ? "baselinecat%d_log.pdf" : "massfacmvacat%d.pdf"),cat));
		can->SaveAs(Form( (baseline ? "baselinecat%d_log.png" : "massfacmvacat%d.png"),cat));
		can->SaveAs(Form( (baseline ? "baselinecat%d_log.eps" : "massfacmvacat%d.eps"),cat));
		can->SaveAs(Form( (baseline ? "baselinecat%d_log.C" : "massfacmvacat%d.C"),cat));		
		*/

		///*
		//can->SaveAs(Form( (baseline ? "baselinecat%d.pdf" : "massfacmvacat%d.pdf"),cat));
		//can->SaveAs(Form( (baseline ? "baselinecat%d.png" : "massfacmvacat%d.png"),cat));
		//can->SaveAs(Form( (baseline ? "baselinecat%d.eps" : "massfacmvacat%d.eps"),cat));
		//can->SaveAs(Form( (baseline ? "baselinecat%d.C" : "massfacmvacat%d.C"),cat));
		//*/
	}

}


using namespace RooFit;


void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname, int cat)
{
	RooRealVar *nlim = new RooRealVar(TString::Format("nlim%s",catname.Data()),"",0.0,0.0,1e+5);

	float ratio9568low=0;
        float ratio9568up=0;
        int counterRatiolow=0;
	int counterRatioup=0;
	double err68=0;
	double err68Up=0;
	double Err68=0;
	double Err68Up=0;

	cout<<"Inside doBandsFit"<<endl;
	cout<<"Number of bins "<<plot->GetXaxis()->GetNbins()<<endl;
	cout<<"Bin @ 125-1.5 "<<plot->GetXaxis()->FindBin(123.5)<<endl;
	cout<<"Bin @ 125+1.5 "<<plot->GetXaxis()->FindBin(125.5)<<endl;
	cout<<"Commented cout for errors 1,2 sigmas"<<endl;

	TGraph *onesigma_up = 0, *onesigma_down = 0;
	onesigma_up = new TGraph();
	onesigma_down = new TGraph();
	double integral_bkg=0;
	double integral_1sup_bkg=0;
	double integral_1sdwn_bkg=0;

	int N=plot->GetXaxis()->GetNbins()+1;
	for (int i=1; i<(plot->GetXaxis()->GetNbins()+1); ++i) {

		double lowedge = plot->GetXaxis()->GetBinLowEdge(i);
		double upedge = plot->GetXaxis()->GetBinUpEdge(i);
		double center = plot->GetXaxis()->GetBinCenter(i);
        
		double nombkg = nomcurve->interpolate(center);

		nlim->setVal(nombkg);
		hmass->setRange("errRange",lowedge,upedge);
		RooAbsPdf *epdf = 0;
		epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        
		RooAbsReal *nll = epdf->createNLL(*datanorm,Extended(),NumCPU(4));
		RooMinimizer minim(*nll);
		minim.setStrategy(0);
		minim.setPrintLevel(-1);
		double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
		double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        
		minim.migrad();
		minim.minos(*nlim);
		
		//This line below shows errors
		//printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
		
		onesigma->SetPoint(i-1,center,nombkg);

		if(i>=30 && i<=32){//I want just the integral below the signal window
		integral_bkg+=nombkg;
		integral_1sup_bkg+=(nombkg + nlim->getErrorHi());
		integral_1sdwn_bkg+=(nombkg + nlim->getErrorLo());//nlim->getErrorLo() fornisce gia' il giusto segno meno
		}

		onesigma_up->SetPoint(i-1,center,nombkg+nlim->getErrorHi());
		onesigma_up->SetPoint(i-1,center,nombkg-nlim->getErrorLo());
		onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());

		Err68+=nlim->getErrorLo();
		Err68Up+=nlim->getErrorHi();

		//Make up for 1 sigma bands (needed only for leptonic)
		err68=-nlim->getErrorLo();
		err68Up=nlim->getErrorHi();

		if(cat==6){
		  if(fabs(nlim->getErrorLo())>0. && fabs(nlim->getErrorHi()) > 0. ){
		    if(nlim->getErrorHi()>0.018 && fabs(nlim->getErrorLo())<0.01425){
		      onesigma->SetPointError(i-1,0.,0.,fabs(nlim->getErrorLo()),nlim->getErrorHi());
		    }else if(nlim->getErrorHi()<0.018 && fabs(nlim->getErrorLo())<0.01425){
		      onesigma->SetPointError(i-1,0.,0.,fabs(nlim->getErrorLo()),0.021);
		      err68Up=0.021;
		    }else if(nlim->getErrorHi()<0.018 && fabs(nlim->getErrorLo())>0.01425){
		      onesigma->SetPointError(i-1,0.,0.,0.01405,0.021);
		      err68=0.01405;
		      err68Up=0.021;
		    }else if(nlim->getErrorHi()>0.018 && fabs(nlim->getErrorLo())>0.01425){
		      onesigma->SetPointError(i-1,0.,0.,0.01405,nlim->getErrorHi());
		      err68=0.01405;
		    }
		  }else if (fabs(nlim->getErrorHi()) > 0.){
		    if(fabs(nlim->getErrorHi())>0.018){
		      onesigma->SetPointError(i-1,0.,0.,0.01405,nlim->getErrorHi());
		      err68=0.01405;
		    }
		    else{
		      onesigma->SetPointError(i-1,0.,0.,0.01405,0.021);
		      err68=0.01405;
		      err68Up=0.021;
		    }
	
		  } else if (fabs(nlim->getErrorLo())>0){
		    if(fabs(nlim->getErrorLo())<0.01425){
		      onesigma->SetPointError(i-1,0.,0.,fabs(nlim->getErrorLo()),0.021);
		      err68Up=0.021;
		    }else{
		      onesigma->SetPointError(i-1,0.,0.,0.01405,0.021);
		      err68=0.01405;
		      err68Up=0.021;
		    }
		  }else {
		    onesigma->SetPointError(i-1,0.,0.,0.01405,0.021);
		    err68=0.01405;
		    err68Up=0.021;
		  }
		}
		//End of make-up for 1 sigma bands

		minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
		// eventually if cl = 0.95 this is the usual 1.92!      
        	minim.migrad();
		minim.minos(*nlim);
		
		twosigma->SetPoint(i-1,center,nombkg);
		twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());      
		ratio9568low+=nlim->getErrorLo()/err68;
		ratio9568up+=nlim->getErrorHi()/err68Up;

		//Make up for 2 sigmas bands
		if(cat==6){
		if(nlim->getErrorLo()==0 && nlim->getErrorHi()==0){
		twosigma->SetPointError(i-1,0.,0.,err68*1.6,err68Up*2.09);
		}else if(nlim->getErrorLo()==0){
		  twosigma->SetPointError(i-1,0.,0.,err68*1.6,nlim->getErrorHi());
		}else if(nlim->getErrorHi()==0){
		twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),err68Up*2.09);
		}
		if(center>130.){
		twosigma->SetPointError(i-1,0.,0.,err68*1.6,err68Up*2.35);
		}
		}else{//Hadronic category
		  twosigma->SetPointError(i-1,0.,0.,err68*1.6,err68Up*2.09);
		}
		//End of make-up for 2 sigmas band
        		
		delete nll;
		delete epdf;
	}
	cout<<"**********"<<endl;
	//cout<<"Mean "<<Err68/N<<endl;
	//cout<<"Mean up "<<Err68Up/N<<endl;
	//cout<<"95/68 up "<<ratio9568up/N<<endl;
	//cout<<"95/68 dwn "<<ratio9568low/N<<endl;

	cout<<"Non-resonant background integral "<<onesigma->Integral()<<endl;
	//cout<<"Non-resonant background integral: one sigma up "<<onesigma_up->Integral()<<endl;
	//cout<<"Non-resonant background integrla: one sigma down "<<onesigma_down->Integral()<<endl;
	//cout<<"Per stampare a schermo tutti i valori di punti e errori del TGraph a 1 sigma, scommenta onesigma->Print()"<<endl;
	//onesigma->Print("V");
	//cout<<"Integral_bkg "<<integral_bkg<<endl;
	//cout<<"Integral_1sup_bkg "<<integral_1sup_bkg<<endl;
	//cout<<"Integral_1sdwn_bkg "<<integral_1sdwn_bkg<<endl;
	//cout<<"Error up "<<(integral_1sup_bkg - integral_bkg)<<endl;
	//cout<<"Error dwn "<<(integral_1sdwn_bkg - integral_bkg)<<endl;

}
