//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//This code calculates the p-value for a given dataset of 2011 and 2012 data
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//directory settings
#include "settings.h"
// C++ library or ROOT header files
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TAxis.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TFitResultPtr.h>
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TMath.h"

double ChiSquared(int intMin, int intMax, TF1* fit_b,TH1D* h_data){
	//set the range for chi2
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	//calculate the chi2
	double chi2 = 0.0;
	for(int i = intMin; i < intMax; i++)
	{
		double intError = fit_b->IntegralError(i,i+1);
		int binNum = (int) i - xmin + 1;
		double data = h_data->GetBinContent(binNum);
		double dataErr = h_data->GetBinError(binNum);
		double fitVal = fit_b->Eval(h_data->GetXaxis()->GetBinCenter(binNum));
		chi2 += (data-fitVal)*(data-fitVal)/(dataErr*dataErr+intError*intError);
	}

	return chi2;
}


//calculate the pValue for given degrees of freedom and chi2
double pValue(int ndf, double chi2)
{
	return TMath::Prob(chi2,ndf);
}

//plot integral error for several ranges
int chi2Evol()
{
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data11 = TFile::Open(TString::Format("%s/data2011_10GeV.root",baseDir.Data()));
	TFile* f_data12 = TFile::Open(TString::Format("%s/data2012_10GeV.root",baseDir.Data()));
	TH1D* h_data11 = (TH1D*) f_data11->Get(TString::Format("h_mgg1"));
	TH1D* h_data = (TH1D*) f_data12->Get(TString::Format("h_mgg1"));
	h_data->Add(h_data11); //h_data holds all events
	//perform background fit
	TF1* fit_b = new TF1("fit_b","pol5",100.0,180.0);
	TFitResultPtr fitRes_b = h_data->Fit(fit_b,"IEMR","same",100.0,180.0);
	//get minimum x value
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	//maximum range to each side
	int r_max = 25;
	double r_list[r_max],chi2_list[r_max];
	//loop over several ranges around 125GeV
	for(int r = 1; r <= r_max;r++)
	{
		r_list[r-1] = r;
		chi2_list[r-1] = ChiSquared(125-r,125+r,fit_b,h_data);
	}

	//create new TGraph
	TCanvas *c1 = new TCanvas("c1","Background integral error",600,600);
	TGraph *gr = new TGraph(r_max,r_list,chi2_list);
	gr->SetTitle("Chi_Squared");
	gr->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
	gr->SetMarkerStyle(7);
	gr->SetMarkerColor(kRed);
	gr->Draw("ap");
	c1->SaveAs("plots/chi2.pdf");
	return 1;
}

int ErrEvol()
{
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data11 = TFile::Open(TString::Format("%s/data2011_10GeV.root",baseDir.Data()));
	TFile* f_data12 = TFile::Open(TString::Format("%s/data2012_10GeV.root",baseDir.Data()));
	TH1D* h_data11 = (TH1D*) f_data11->Get(TString::Format("h_mgg1"));
	TH1D* h_data = (TH1D*) f_data12->Get(TString::Format("h_mgg1"));
	h_data->Add(h_data11); //h_data holds all events
	//perform background fit
	TF1* fit_b = new TF1("fit_b","pol5",100.0,180.0);
	TFitResultPtr fitRes_b = h_data->Fit(fit_b,"IEMR","same",100.0,180.0);
	//get minimum x value
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	//maximum range to each side
	int r_max = 25;
	double r_list[r_max],intErr_list[r_max];
	//loop over several ranges around 125GeV
	for(int r = 1; r <= r_max;r++)
	{
		r_list[r-1] = r;
		double intErr = 0.0;
		for(int i = -r; i < r; i++)
			intErr += fit_b->IntegralError(125+i,125+i+1);
		intErr_list[r] = intErr;
	}

	//create new TGraph
	TCanvas *c1 = new TCanvas("c1","Background integral error",600,600);
	TGraph *gr = new TGraph(r_max,r_list,intErr_list);
	gr->SetTitle("Background error");
	gr->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
	gr->SetMarkerStyle(7);
	gr->SetMarkerColor(kRed);
	gr->Draw("ap");
	c1->SaveAs("plots/intErr.pdf");
	return 1;
}



//***********Main Function*************
int main(int argc, char** argv)
{
	//int ndf = 80-9;
	//double p = pValue(chi2,ndf);
	//("p-val : %f \n", p);
	chi2Evol();
	ErrEvol();
}
