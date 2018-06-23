#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMatrixTSym.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "Fit/FitResult.h"
#include "TLegend.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

Double_t gBinW, gContent, gMean, gSigma;

Double_t gaus_lbg(Double_t *x, Double_t *par)
{
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg;
   if (par[2] == 0) par[2]=1;                
   arg = (x[0] - par[4])/(sqrt2*par[2]);
   Double_t fitval = par[0] + x[0]*par[1]
              + gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg*arg);
   return fitval;
}

void AutoECalibrate(TH1D *hist, Double_t gUpCh, Double_t gUpChErr, string fitopt=""){

	gStyle->SetOptTitle(0);
	TCanvas *c1 = new TCanvas();
	c1->Divide(1,2,0,0,0);
	c1->cd(1);
	const int gate_width = 7;
	const int n = 17;
	hist->Draw();
	hist->GetYaxis()->SetTitle("E_{#gamma} (keV)");
	hist->GetYaxis()->SetTickLength(0);
	hist->GetYaxis()->CenterTitle();
	hist->GetXaxis()->SetRange(-10,gUpCh+600);
	c1->Update();
	Double_t g_ene[n] = {0,39.905,45.602,121.783,244.692,344.276,411.115,443.976,778.903,867.388,964.131,1085.842,1089.767,1112.116,1212.950,1299.124,1408.011};
	Double_t g_err[n] = {0,10e-3,10e-3,4e-3,10e-3,19e-3,3e-3,4e-3,6e-3,6e-3,4e-3,4e-3,14e-3,6e-3,13e-3,9e-3,4e-3};
	if (sizeof(g_ene) != sizeof(g_err)){
		cerr << "Size of arrays do not match!! Check data again!!" << endl;
		return;
	}
	Double_t fit_ene[n];
	Double_t fit_err[n];
	fit_ene[0] = 0; fit_err[0] = 10;
	fit_ene[n-1] = gUpCh; fit_err[n-1] = gUpChErr;
	Double_t grad = g_ene[n-1]/fit_ene[n-1];
	Double_t guess[n-1];
  	gBinW = hist->GetBinWidth(1);

	TF1 fitfunc("gauss_linbg",gaus_lbg, 0, 1, 5);
	fitopt += "RQ";
   	fitopt += "EMS+";
	TFitResultPtr res;
  	fitfunc.SetParName(0,"BgConstant");
  	fitfunc.SetParName(1,"BgSlope   ");
  	fitfunc.SetParName(2,"Sigma     ");
  	fitfunc.SetParName(3,"Content   ");
  	fitfunc.SetParName(4,"Mean      ");

	for(int i = 1; i <= n-1; i++){

		guess[i] = g_ene[i]/grad;

		gContent = hist->Integral(hist->FindBin(guess[i]-gate_width),hist->FindBin(guess[i]+gate_width)); 
  		gMean    = 0.5 * (guess[i]-gate_width + guess[i]+gate_width);  
  		gSigma   = 0.3 * (guess[i]+gate_width - guess[i]-gate_width); 
  		fitfunc.SetParameters(0, 0, gSigma, gContent, gMean); 
  		fitfunc.SetRange(guess[i]-gate_width, guess[i]+gate_width);

		res = hist->Fit("gauss_linbg", fitopt.c_str(), "SAME");

   		fit_ene[i]       = fitfunc.GetParameter(4);
		//cout << "\n Centroid = " << fit_ene[i] << " (" << res->LowerError(4) << ";+" << res->UpperError(4) << ")\n";
		fit_err[i] = (res->LowerError(4)+res->UpperError(4))/2;
		
	}
	
	c1->cd(2);
	TGraphErrors *gr_Ecalib = new TGraphErrors(n,fit_ene,g_ene,fit_err,g_err);

	TF1 quadfit("quadfit","pol2", 0, guess[n-2]);
	quadfit.SetParName(0,"Intercept	:");
  	quadfit.SetParName(1,"Graident	:");
   	quadfit.SetParName(2,"Quadratic	:");
	quadfit.SetParameters(0, grad, 0);
	res = gr_Ecalib->Fit("quadfit", fitopt.c_str(), "SAME");

	Double_t inter = quadfit.GetParameter(0);
	grad = quadfit.GetParameter(1);	
	Double_t quad = quadfit.GetParameter(2);

	Double_t res_ene[n];
	Double_t zero[n] = {0};
	Double_t res_err[n];
	TGraph *gr_uperr = new TGraph(n-1);
	TGraph *gr_loerr = new TGraph(n-1);
	TGraphAsymmErrors *gr_shade = new TGraphAsymmErrors(2*(n-1));
	for(int i = 0; i <= n-1; i++){
		
		res_ene[i] = g_ene[i]-(inter+grad*fit_ene[i]+quad*fit_ene[i]*fit_ene[i]);
		res_err[i] = 0.3;
		gr_uperr->SetPoint( i, fit_ene[i], 10*res->UpperError(0)+res->UpperError(1)*fit_ene[i]+res->UpperError(2)*fit_ene[i]*fit_ene[i]);
		gr_loerr->SetPoint( i, fit_ene[i], 10*res->LowerError(0)+res->LowerError(1)*fit_ene[i]+res->LowerError(2)*fit_ene[i]*fit_ene[i]);
		gr_shade->SetPoint( i, fit_ene[i], 10*res->UpperError(0)+res->UpperError(1)*fit_ene[i]+res->UpperError(2)*fit_ene[i]*fit_ene[i]);
		gr_shade->SetPoint( n-1+i, fit_ene[n-i], 10*res->LowerError(0)+res->LowerError(1)*fit_ene[n-1-i]+res->LowerError(2)*fit_ene[n-i]*fit_ene[n-i]);
	}

	gr_shade->SetPoint(2*(n-1),fit_ene[0], 10*res->LowerError(0));

	//gr_uperr->SetPoint(

	TMultiGraph * mg = new TMultiGraph();
	TGraphErrors *gr_Eresid = new TGraphErrors(n,fit_ene,res_ene,fit_err,res_err);
	gr_Eresid->SetMarkerColor(4);
	gr_Eresid->SetMarkerStyle(24);
	gr_Eresid->SetMarkerSize(0.75);
	TGraph *gr_base = new TGraphErrors(n,fit_ene,zero);
	gr_base->SetLineColor(2);
	gr_base->SetLineWidth(3);
	gr_uperr->SetLineStyle(10);
	gr_uperr->SetLineColorAlpha(1,0.3);
	gr_uperr->SetLineWidth(2);
	gr_loerr->SetLineStyle(10);
	gr_loerr->SetLineColorAlpha(1,0.3);
	gr_loerr->SetLineWidth(2);
	gr_shade->SetFillStyle(1001);;
	gr_shade->SetFillColorAlpha(46,0.5);
	mg->Add(gr_shade,"f");
	mg->Add(gr_base,"C");
	mg->Add(gr_loerr,"C");
	mg->Add(gr_uperr,"C");
	mg->Add(gr_Eresid,"P");
	mg->Draw("A");
	mg->GetXaxis()->SetRange(-10,gUpCh+600);
	mg->GetXaxis()->SetTitle("Channel");
	mg->GetYaxis()->SetTitle("E_{#gamma} - F(ch) (keV)");
	mg->GetYaxis()->SetTickLength(0);
	mg->GetXaxis()->CenterTitle();
	mg->GetYaxis()->CenterTitle();

	c1->Update();

	res->Print();
	cout << "****************************************\n" << endl;

}
