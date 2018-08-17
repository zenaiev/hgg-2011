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

//plot integral error for several ranges
int chi2EvolData(TH1D* h_data, TF1* fit_b, TString plotName)
{
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
	c1->SaveAs("plots/" + plotName);
	return 1;
}

int chi2EvolMC(TH1D* h_mc, TF1* fit_b, TString plotName)
{
	//get minimum value of histogram
	double xmin = h_mc->GetXaxis()->GetBinCenter(1)-0.5;
	//maximum range to each side
	int r_max = 25;
	double r_list[r_max], chi2_list[r_max];
	for(int r = 1; r <= r_max; r++)
	{
		r_list[r-1] = r;
		chi2_list[r-1] = ChiSquaredMC(125-r,125+r,fit_b,h_mc);
	}
	TCanvas *c1 = new TCanvas("c1","Background integral error",600,600);
	TGraph *gr = new TGraph(r_max,r_list,chi2_list);
	gr->SetTitle("Chi_Squared MC");
	gr->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
	gr->SetMarkerStyle(7);
	gr->SetMarkerColor(kRed);
	gr->Draw("ap");
	c1->SaveAs("plots/" + plotName);
	return 1;


}

int pvalEvolData(TH1D* h_data, TF1* fit_b, TString plotName)
{

	//get minimum x value
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	//maximum range to each side
	int r_max = 25;
	double r_list[r_max],pval_list[r_max];
	//loop over several ranges around 125GeV
	for(int r = 1; r <= r_max;r++)
	{
		r_list[r-1] = r;
		double chi2 = ChiSquared(125-r,125+r,fit_b,h_data);
		//degrees of freedom is 1 because all bins are summed/ integrated bkg
		int ndf =1;
		pval_list[r-1] = pValue(ndf,chi2);
	}

	//create new TGraph
	TCanvas *c1 = new TCanvas("c1","Background integral error",600,600);
	TGraph *gr = new TGraph(r_max,r_list,pval_list);
	gr->SetTitle("p-value");
	gr->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
	gr->SetMarkerStyle(7);
	gr->SetMarkerColor(kRed);
	gr->Draw("ap");
	c1->SaveAs("plots/" + plotName);
	return 1;
}

int pvalEvolMC(TH1D* h_mc, TF1* fit_b, TString plotName)
{
	//get minimum value of histogram
	double xmin = h_mc->GetXaxis()->GetBinCenter(1)-0.5;
	//maximum range to each side
	int r_max = 25;
	double r_list[r_max], pval_list[r_max];
	for(int r = 1; r <= r_max; r++)
	{
		r_list[r-1] = r;
		double chi2 = ChiSquaredMC(125-r,125+r,fit_b,h_mc);
		pval_list[r-1] = pValue(1,chi2);
	}
	TCanvas *c1 = new TCanvas("c1","Background integral error",600,600);
	TGraph *gr = new TGraph(r_max,r_list,pval_list);
	gr->SetTitle("p-value MC");
	gr->GetXaxis()->SetTitle("Width of range (around 125 GeV)");
	gr->SetMarkerStyle(7);
	gr->SetMarkerColor(kRed);
	gr->Draw("ap");
	c1->SaveAs("plots/" + plotName);
	return 1;
	

}

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
	//Read 2011 and 2012 data
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data11 = TFile::Open(TString::Format("%s/data2011_10GeV.root",baseDir.Data()));
	TFile* f_data12 = TFile::Open(TString::Format("%s/data2012_10GeV.root",baseDir.Data()));
	TH1D* h_data11 = (TH1D*) f_data11->Get(TString::Format("h_mgg1"));
	TH1D* h_data12 = (TH1D*) f_data12->Get(TString::Format("h_mgg1"));
	//print number of entries for each year
	printf("Entries 2011: %f \n", h_data11->GetEntries());
	printf("Entries 2012: %f \n", h_data12->GetEntries());
	//combine data
	TH1D* h_data = (TH1D*) h_data11->Clone();
	h_data->Add(h_data12);
	printf("Entries combined: %f \n", h_data->GetEntries());

	//Classes plots
	pvalClasses11("pval_2011_classes.pdf");
	pvalClasses12("pval_2012_classes.pdf");
	pval12("pval_2012.eps");


}
