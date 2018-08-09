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

double ChiSquared(int intMin, int intMax, double fitMin = 100.0, double fitMax = 180.0){
	//import data
	printf("intmin : %d \n", intMin);
	printf("intMax : %d \n", intMax);
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data11 = TFile::Open(TString::Format("%s/data2011_10GeV.root",baseDir.Data()));
	TFile* f_data12 = TFile::Open(TString::Format("%s/data2012_10GeV.root",baseDir.Data()));
	TH1D* h_data11 = (TH1D*) f_data11->Get(TString::Format("h_mgg1"));
	TH1D* h_data = (TH1D*) f_data12->Get(TString::Format("h_mgg1"));
	h_data->Add(h_data11); //h_data holds all events
	//perform background fit
	TF1* fit_b = new TF1("fit_b","pol5",fitMin,fitMax);
	TFitResultPtr fitRes = h_data->Fit(fit_b,"IEMR","same",fitMin,fitMax);

	//set the range for chi2
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	
	//calculate the chi2
	double chi2 = 0.0;
	printf("Fit done \n");
	for(int i = intMin; i < intMax; i++)
	{
		double intError = fit_b->IntegralError(i-xmin,i-xmin+1); //slow!
		int binNum = (int) i - xmin + 1;
		double data = h_data->GetBinContent(binNum);
		double dataErr = h_data->GetBinError(binNum);
		double fitVal = fit_b->Eval(h_data->GetXaxis()->GetBinCenter(binNum));
		chi2 += (data-fitVal)*(data-fitVal)/(dataErr*dataErr+intError*intError);
	}

	return chi2;
}

double ChiSquaredFits(double m_H,int intMin,int intMax,double fitMin = 100.0, double fitMax = 180.0)
{
	TString baseDir = gHistDir; //from settings.h
	TFile* f_data11 = TFile::Open(TString::Format("%s/data2011_10GeV.root",baseDir.Data()));
	TFile* f_data12 = TFile::Open(TString::Format("%s/data2012_10GeV.root",baseDir.Data()));
	TH1D* h_data11 = (TH1D*) f_data11->Get(TString::Format("h_mgg1"));
	TH1D* h_data = (TH1D*) f_data12->Get(TString::Format("h_mgg1"));
	h_data->Add(h_data11); //h_data holds all events
	//perform background fit
	TF1* fit_b = new TF1("fit_b","pol5",fitMin,fitMax);
	TFitResultPtr fitRes_b = h_data->Fit(fit_b,"IEMR","same",fitMin,fitMax);
	//perform background + signal fit
	TF1* fit_sb = new TF1("fit_sb","gaus(0)+pol5(3)",fitMin,fitMax);
	//set start parameters
	fit_sb->SetParameter(0, 0.0);
	fit_sb->SetParameter(1, m_H);
	fit_sb->SetParameter(2, 2.0);
	fit_sb->FixParameter(1, m_H); // fix the Higgs mass
	//set initial parameters for background
	for(int p = 3; p <= 9; p++)
		fit_sb->SetParameter(p,fit_b->GetParameter(p-3));
	//fit
	TFitResultPtr fitRes_sb = h_data->Fit(fit_sb,"IEMR","same",fitMin,fitMax);
	double xmin = h_data->GetXaxis()->GetBinCenter(1)-0.5;
	//calculate chi2
	double chi2 = 0.0;
	for(int i = intMin; i < intMax; i++)
	{
		int binNum = (int) i - xmin +1;
		double bVal = fit_b->Eval(h_data->GetXaxis()->GetBinCenter(binNum));;
		double sbVal = fit_sb->Eval(h_data->GetXaxis()->GetBinCenter(binNum));;
		chi2 += (sbVal - bVal)*(sbVal - bVal)/bVal;
	}
	return chi2;

}

//calculate the pValue for given degrees of freedom and chi2
double pValue(int ndf, double chi2)
{
	return TMath::Prob(chi2,ndf);
}



//***********Main Function*************
int main(int argc, char** argv)
{
	printf("p-value plot \n");
	double chi2 = ChiSquared(124,126);
	printf("chi2 : %f \n", chi2);
	int ndf = 80-9;
	double p = pValue(chi2,ndf);
	printf("p-val : %f \n", p);
}
