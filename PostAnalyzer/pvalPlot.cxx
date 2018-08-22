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
#include "TMultiGraph.h"

//Calculation of chi2 on a dataset and bkg fit
//params:
//	int intMin	: minimum of range
//  int intMax	: maximum of range
//	TF1* fit_b	: Background fit
//  TH1D* h_data: dataset as histogram
double ChiSquared(int intMin, int intMax, TF1* fit_b,TH1D* h_data){
	//set the range for chi2
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	//calculate the chi2
	double chi2 = 0.0;
	double intError = fit_b->IntegralError(intMin,intMax);
	double dataSum = 0.0;
	double fitSum = fit_b->Integral(intMin,intMax);
	for(int i = intMin; i < intMax; i++)
	{
		int binNum = (int) i - xmin + 1;
		double data = h_data->GetBinContent(binNum);
		dataSum += data;
	}
	//dataSum is already the square of the error
	chi2 = (dataSum - fitSum)*(dataSum - fitSum)  / (dataSum + intError*intError);
	return chi2;
}

//Calculation of chi2 on the mc-signal and bkg fit from data
//params:
//	int intMin	: minimum of range
//  int intMax	: maximum of range
//	TF1* fit_b	: Background fit
//  TH1D* h_mc  : mc-signal as histogram
double ChiSquaredMC(int intMin, int intMax, TF1* fit_b,TH1D* h_mc){
	//set the range for chi2
	double xmin = h_mc->GetXaxis()->GetBinCenter(1)-0.5;
	//calculate the chi2
	double chi2 = 0.0;
	double intError = fit_b->IntegralError(intMin,intMax);
	double dataSum = 0.0;
	double fitSum = fit_b->Integral(intMin,intMax);
	for(int i = intMin; i < intMax; i++)
	{
		int binNum = (int) i - xmin + 1;
		double mc_events = h_mc->GetBinContent(binNum);
		dataSum += mc_events;
	}
	//add bkg to MC signal
	dataSum += fitSum;	
	//dataSum is already the square of the error
	chi2 = (dataSum - fitSum)*(dataSum - fitSum)  / (dataSum + intError*intError);
	return chi2;
}

//calculate the pValue for given degrees of freedom and chi2
double pValue(int ndf, double chi2)
{
	return TMath::Prob(chi2,ndf);
}

//calculate the pvalue for different ranges in around 125 GeV in each class for 2011 data
int pvalClasses11(TString plotName)
{
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data = TFile::Open(TString::Format("%s/data2011_10GeV.root",baseDir.Data()));
	TCanvas* c = new TCanvas("","",600,600);
	c->Divide(2,3,0.01,0.01);
	
	//set titles
	std::vector<TString> vTitle;
  	vTitle.push_back("a) All classes combined");
  	vTitle.push_back("b) Di-jet tagged class");
  	vTitle.push_back("c) Both #gamma in barrel, R_{9}^{min}>0.94");
  	vTitle.push_back("d) Both #gamma in barrel, R_{9}^{min}<0.94");
  	vTitle.push_back("e) One or both #gamma in endcap, R_{9}^{min}>0.94");
  	vTitle.push_back("f) One or both #gamma in endcap, R_{9}^{min}<0.94");

  	//loop over all pads
  	for(int pad = 1; pad <= 6; pad++)
  	{
  		c->cd(pad);
  		//read corresponding class from argument TFile
  		TH1D* h_data = (TH1D*) f_data->Get(TString::Format("h_mgg%d",pad));

  		if(pad != 2)
  		{
  			//fit the background
  			TF1* fit_b = new TF1("fit_b","pol5",100.0,180.0);
  			TFitResultPtr fitRes_b = h_data->Fit(fit_b,"IEMR","same",100.0,180.0);

  			//calculate p-values
  			double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
  			int r_max = 15;
  			double r_list[r_max], pval_list[r_max];
  			for(int r = 1; r <= r_max;r++)
  			{
  				r_list[r-1] = r;
  				double chi2 = ChiSquared(125-r,125+r,fit_b,h_data);
  				pval_list[r-1] = pValue(1,chi2);
  			}
  			//create TGraph for plotting
  			TGraph *gr = new TGraph(r_max,r_list,pval_list);
  			gr->SetTitle(vTitle[pad-1]);
  			gr->GetYaxis()->SetTitle("p-value");
			gr->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
			gr->GetXaxis()->SetTitleSize(0.04);
			gr->GetYaxis()->SetTitleSize(0.04);
			gr->GetYaxis()->SetRangeUser(0,1);
			gr->SetMarkerStyle(20);
			gr->SetMarkerSize(0.35);

			//add a legend for sigma value at width 2
			TLegend* leg = new TLegend(0.06, 0.80, 0.90, 0.93);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0.0);
			printf(" \n Class %d with a p-value at width 2 : %f \n \n", pad, pval_list[1]);
			//TString s = new TString("#sigma_2 = %f", )
			leg->AddEntry((TObject*)NULL, "", "");
			leg->Draw();

			gr->SetMarkerColor(kRed);
			gr->Draw("ap");
  		}

  	}
  	c->SaveAs("plots/" + plotName);
}
//calculate the pvalue for different ranges in around 125 GeV in each class for 2012 data
int pvalClasses12(TString plotName)
{
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data = TFile::Open(TString::Format("%s/data2012_10GeV.root",baseDir.Data()));
	TFile* f_mc = TFile::Open(TString::Format("%s/mcSigReco.root",baseDir.Data()));
	TCanvas* c = new TCanvas("","",600,600);
	c->Divide(2,3,0.01,0.01);

	//set titles
	std::vector<TString> vTitle;
  	vTitle.push_back("a) All classes combined");
  	vTitle.push_back("b) Di-jet tagged class");
  	vTitle.push_back("c) Both #gamma in barrel, R_{9}^{min}>0.94");
  	vTitle.push_back("d) Both #gamma in barrel, R_{9}^{min}<0.94");
  	vTitle.push_back("e) One or both #gamma in endcap, R_{9}^{min}>0.94");
  	vTitle.push_back("f) One or both #gamma in endcap, R_{9}^{min}<0.94");

  	//loop over all pads
  	for(int pad = 1; pad <= 6; pad++)
  	{
  		c->cd(pad);
  		//read corresponding class from argument TFile
  		TH1D* h_data = (TH1D*) f_data->Get(TString::Format("h_mgg%d",pad));
  		TH1D* h_mc = (TH1D*) f_mc->Get(TString::Format("h_mgg%d",pad));

  		if(pad != 2)
  		{
  			//fit the background
  			TF1* fit_b = new TF1("fit_b","pol5",100.0,180.0);
  			TFitResultPtr fitRes_b = h_data->Fit(fit_b,"IEMR","same",100.0,180.0);

  			//calculate p-values
  			double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
  			int r_max = 15;
  			double r_list[r_max], pvalData_list[r_max],pvalMC_list[r_max];
  			for(int r = 1; r <= r_max;r++)
  			{
  				r_list[r-1] = r;
  				double chi2 = ChiSquared(125-r,125+r,fit_b,h_data);
  				double chi2_mc = ChiSquaredMC(125-r,125+r,fit_b,h_mc);
  				pvalData_list[r-1] = pValue(1,chi2);
  				pvalMC_list[r-1] = pValue(1,chi2_mc);
  			}

  			//create TGraph for plotting
  			TMultiGraph *mg = new TMultiGraph();
  			TGraph *gr = new TGraph(r_max,r_list,pvalData_list);
  			TGraph *gr_mc = new TGraph(r_max,r_list,pvalMC_list);



			//set the plot style
			gr->SetMarkerStyle(20);
			gr->SetMarkerSize(0.35);
			gr->SetMarkerColor(kRed);
			mg->Add(gr);
			gr_mc->SetMarkerStyle(20);
			gr_mc->SetMarkerSize(0.35);
			gr_mc->SetMarkerColor(kBlue);
			mg->Add(gr_mc);

			//draw the multigraph
			mg->Draw("ap");
			mg->SetTitle(vTitle[pad-1]);
			mg->GetYaxis()->SetTitle("p-value");
			mg->GetYaxis()->SetTitle("p-value");
			mg->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
			mg->GetXaxis()->SetTitleSize(0.04);
			mg->GetYaxis()->SetTitleSize(0.04);
			mg->GetYaxis()->SetRangeUser(0,1);
			//add a legend for data and mc
			TLegend* leg = new TLegend(0.6, 0.2, 0.9, 0.3);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0.0);
			//printf(" \n Class %d with a p-value at width 2 : %f \n \n", pad, pvalData_list[1]);
			//TString s = new TString("#sigma_2 = %f", )
			leg->AddEntry(gr, "Data 2012","pe");
			leg->AddEntry(gr_mc, "MC 2012","pe");
			leg->Draw();
  		}

  	}
  	c->SaveAs("plots/" + plotName);
}

int pval12(TString plotName)
{
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data = TFile::Open(TString::Format("%s/data2012_10GeV.root",baseDir.Data()));
	TFile* f_mc = TFile::Open(TString::Format("%s/mcSigReco.root",baseDir.Data()));
	TCanvas* c = new TCanvas("","",600,600);

	//read corresponding class from argument TFile
	TH1D* h_data = (TH1D*) f_data->Get(TString::Format("h_mgg%d",1));
	TH1D* h_mc = (TH1D*) f_mc->Get(TString::Format("h_mgg%d",1));

	//fit the background
	TF1* fit_b = new TF1("fit_b","pol5",100.0,180.0);
	TFitResultPtr fitRes_b = h_data->Fit(fit_b,"IEMR","same",100.0,180.0);

	//calculate p-values
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	int r_max = 15;
	double r_list[r_max], pvalData_list[r_max],pvalMC_list[r_max];
	for(int r = 1; r <= r_max;r++)
	{
		r_list[r-1] = r;
		double chi2 = ChiSquared(125-r,125+r,fit_b,h_data);
		double chi2_mc = ChiSquaredMC(125-r,125+r,fit_b,h_mc);
		pvalData_list[r-1] = pValue(1,chi2);
		pvalMC_list[r-1] = pValue(1,chi2_mc);
	}

	//create TGraph for plotting
	TMultiGraph *mg = new TMultiGraph();
	TGraph *gr = new TGraph(r_max,r_list,pvalData_list);
	TGraph *gr_mc = new TGraph(r_max,r_list,pvalMC_list);

	printf("width : %f \n", r_list[1]);
	printf("pvalMC : %f \n", pvalMC_list[1]);
	printf("pvalData : %f \n", pvalData_list[1]);

	//set the plot style
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(kRed);
	mg->Add(gr);
	gr_mc->SetMarkerStyle(20);
	gr_mc->SetMarkerSize(0.6);
	gr_mc->SetMarkerColor(kBlue);
	mg->Add(gr_mc);

	//draw the multigraph
	mg->Draw("ap");
	mg->SetTitle("2012 Open data - combined classes");
	mg->GetYaxis()->SetTitle("p-value");
	mg->GetYaxis()->SetTitle("p-value");
	mg->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
	mg->GetXaxis()->SetTitleSize(0.04);
	mg->GetYaxis()->SetTitleSize(0.04);
	mg->GetYaxis()->SetRangeUser(0,1);
	//add a legend for data and mc
	TLegend* leg = new TLegend(0.6, 0.15, 0.9, 0.25);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0.0);
	leg->AddEntry(gr, "Data 2012","pe");
	leg->AddEntry(gr_mc, "MC 2012","pe");
	leg->Draw();



  	c->SaveAs("plots/" + plotName);
}
//***********Main Function*************
int main(int argc, char** argv)
{
	//Classes plots
	pvalClasses11("pval_2011_classes.pdf");
	pvalClasses12("pval_2012_classes.pdf");
	//combined plots for 2012
	pval12("pval_2012.eps");


}
