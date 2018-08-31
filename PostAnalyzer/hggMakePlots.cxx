// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This code processes ROOT histograms for hgg analysis,
// (produced by hggMakeHist.cxx), and makes final plots and numbers
// Run: ./hggMakePlots
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// additional files from this analysis 
#include "plots.h"
#include "settings.h"
// C++ library or ROOT header files
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"

//Plots 2011 Histograms for each class and all classes combined in 1 GeV bins
//can be compared to plot from: CMS-HIG-11-033
int Plot2011()
{
  // set user style
  Style();

  // directory with input histograms
  TString baseDir = gHistDir;
  // directory for output plots (must exist)
  TString plotDir = gPlotsDir;

  // 2012
  TCanvas* c = new TCanvas("", "", 600, 600);
  c->Divide(2, 3, 0.01, 0.01);

  
  //Titles for all pads
  std::vector<TString> vTitle;
  vTitle.push_back("a) All classes combined");
  vTitle.push_back("b) Di-jet tagged class");
  vTitle.push_back("c) Both #gamma in barrel, R_{9}^{min}>0.94");
  vTitle.push_back("d) Both #gamma in barrel, R_{9}^{min}<0.94");
  vTitle.push_back("e) One or both #gamma in endcap, R_{9}^{min}>0.94");
  vTitle.push_back("f) One or both #gamma in endcap, R_{9}^{min}<0.94");

  for(int pad = 1; pad <= 6; pad++)
  {
    //Load corresponding data
    c->cd(pad);
    TFile* f_data = TFile::Open(TString::Format("%s/data2011_10GeV.root", baseDir.Data()));
    TH1D* h_data = (TH1D*)f_data->Get(TString::Format("h_mgg2011_1_%d", pad));
    //set style of plot
    h_data->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    h_data->GetXaxis()->SetTitleOffset(3.0);
    h_data->GetYaxis()->SetTitle("Events / (1 GeV)");
    h_data->GetYaxis()->SetTitleOffset(3.0);
    h_data->GetYaxis()->SetNdivisions(505);
    SetHistoAxisFonts(h_data);
    ScaleHistoFonts(h_data, 0.5);
    h_data->SetMarkerColor(1);
    h_data->SetLineColor(1);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.35);
    h_data->SetTitle("");
    h_data->Draw();

    if(pad != 2)
    {
      double fitMin = 100.0;
      double fitMax = 180.0;
      // background only gives:
      // pol5 gives chi2/dif = 260/74
      // pol7 gives chi2/dof = 93/72
      //pol5 fit of data (background only)
      TF1* fit_b = new TF1("fit_b","pol5", fitMin, fitMax);
      h_data->Fit(fit_b, "IEMR", "same", fitMin, fitMax);
      //pol5 + gaus fit (signal + background)
      TF1* fit_sb = new TF1("fit_sb","gaus(0)+pol5(3)", fitMin, fitMax);
      //Fixing parameters for Gaus fit
      fit_sb->SetParameter(0, 0.0);
      fit_sb->SetParameter(1, 125.0);
      fit_sb->SetParameter(2, 2.0);
      fit_sb->FixParameter(0, 0.0);
      fit_sb->FixParameter(1, 125.0);
      fit_sb->FixParameter(2, 2.0);
      //Set parameters for pol5 fit from first fit
      for(int p = 3; p <= 8; p++)
        fit_sb->SetParameter(p, fit_b->GetParameter(p - 3));
      h_data->Fit(fit_sb, "I0", "", fitMin, fitMax);
      TF1* fit_sb_save = h_data->GetFunction("fit_sb");
      fit_sb_save = new TF1(*fit_sb_save);

      // create a TGraphErrors to hold the confidence intervals
      const int ngr = 100;
      TGraphErrors* gunc68 = new TGraphErrors(ngr);
      TGraphErrors* gunc95 = new TGraphErrors(ngr);
      for (int i = 0; i < ngr; i++)
      {
        double xmin = h_data->GetXaxis()->GetBinCenter(1);
        double xmax = h_data->GetXaxis()->GetBinCenter(h_data->GetNbinsX());
        double x = (i + 0.5) / ngr * (xmax - xmin) + xmin;
        double y = fit_sb->Eval(x);
        gunc68->SetPoint(i, x, y);
        gunc95->SetPoint(i, x, y);
      }
      /*Compute the confidence intervals at the x points of the created graph*/
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gunc68, 0.68);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gunc95);
      gunc68->SetFillColor(kYellow);
      gunc68->SetLineColor(kYellow);
      gunc95->SetFillColor(kGreen);
      gunc95->SetLineColor(kGreen);
      gunc95->Draw("4");
      gunc68->Draw("4");
      fit_sb_save->SetLineColor(kRed);
      fit_sb_save->SetLineWidth(0.75);
      fit_sb_save->Draw("same");

      h_data->Draw("same");

      if(pad == 1)
      {
        TLegend* leg = new TLegend(0.50, 0.45, 0.90, 0.80);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.0);
        leg->AddEntry(h_data, "Data 2011", "pe");
        leg->AddEntry(fit_sb, "Bkg Model", "l");
        leg->AddEntry(gunc68, "#pm 1#sigma", "f");
        leg->AddEntry(gunc95, "#pm 2#sigma", "f");
        leg->Draw();
      }
    }

    TLegend* leg = new TLegend(0.06, 0.80, 0.90, 0.93);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0.0);
    leg->AddEntry((TObject*)NULL, vTitle[pad - 1], "");
    leg->Draw();

    if(pad == 2)
    {
      TLegend* leg = new TLegend(0.10, 0.30, 0.90, 0.7);
      leg->SetBorderSize(0.0);
      leg->SetFillStyle(0);
      leg->SetTextColor(kGray);
      leg->AddEntry((TObject*)(NULL), "not attempted", "");
      leg->Draw();
    }
  }

  c->SaveAs(plotDir + "/hgg-2011.pdf");
  c->SaveAs(plotDir + "/hgg-2011.eps");
}

//Plot 2012 Histograms for each class in 1.5 GeV bins
// can be compared to internal analysis note: CMS AN-2013/008
int Plot2012()
{
  // set user style
  Style();

  // directory with input histograms
  TString baseDir = gHistDir;
  // directory for output plots (must exist)
  TString plotDir = gPlotsDir;

  //Canvas
  TCanvas* c = new TCanvas("", "", 600, 600);
  c->Divide(2,2);
  
  //Titles for all pads : event classes
  std::vector<TString> vTitle;
  vTitle.push_back("a) Both #gamma in barrel, R_{9}^{min}>0.94");
  vTitle.push_back("b) Both #gamma in barrel, R_{9}^{min}<0.94");
  vTitle.push_back("c) One or both #gamma in endcap, R_{9}^{min}>0.94");
  vTitle.push_back("d) One or both #gamma in endcap, R_{9}^{min}<0.94");
  
  //loop over all pads/eventclasses
  for(int pad = 1; pad <= 4; pad++)
  {
    //change to pad
    c->cd(pad);
    TPad* pad1 = new TPad("", "", 0.00, 0.40, 1.00, 1.00);
    pad1->SetMargin(0.11, 0.03, 0.005, 0.07);
    
    //set maximum value yLabel
    double yMax = 2000;
    if(pad == 1) yMax = 1000;
    if(pad == 2) yMax = 3500;
    if(pad == 3) yMax = 1000;
    if(pad == 4) yMax = 4500;
    
    
    TH2F* hr = new TH2F("", "", 1.0, 100.0, 180.0,1,0.01,yMax);
    hr->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    hr->GetXaxis()->SetTitleOffset(3.0);
    hr->GetYaxis()->SetTitle("Events / 1.5 GeV");
    hr->GetYaxis()->SetTitleOffset(3.0);
    SetHistoAxisFonts(hr);
    ScaleHistoFonts(hr, 0.5);
    pad1->cd();
    hr->Draw();

    TFile* f_data = TFile::Open(TString::Format("%s/data2012_15GeV.root", baseDir.Data()));
    TH1D* h_data = (TH1D*)f_data->Get(TString::Format("h_mgg2012_2_%d",pad+2));

    TFile* f_mc = TFile::Open(TString::Format("%s/mcSigReco.root", baseDir.Data()));
    TH1D* h_mc = (TH1D*)f_mc->Get(TString::Format("h_mgg%d",pad+2));

    h_data->SetLineColor(1);
    h_data->SetMarkerColor(1);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.3);
    h_data->SetTitle("");
    h_data->GetXaxis()->SetTitle("");
    h_data->GetYaxis()->SetTitle("");

    // fit
    double fitMin = 100.0;
    double fitMax = 180.0;
    // background only gives:
    // pol5 gives chi2/dif = 260/74
    // pol7 gives chi2/dof = 93/72
    //pol2 fit (background only)
    TF1* fit_b = new TF1("fit_b","pol5", fitMin, fitMax);
    h_data->Fit(fit_b, "IEMR", "same", fitMin, fitMax);
    //pol2 + gaus fit (signal + background)
    TF1* fit_sb = new TF1("fit_sb","gaus(0)+pol5(3)", fitMin, fitMax);
    fit_sb->SetParameter(0, 0.0);
    fit_sb->SetParameter(1, 125.0);
    fit_sb->SetParameter(2, 2.0);
    //Fix only mass parameter
    //fit_sb->FixParameter(0, 0.0);
    fit_sb->FixParameter(1, 125.0);
    //fit_sb->FixParameter(2, 2.0);
    //Set other parameters from background
    for(int p = 3; p <= 9; p++)
      fit_sb->SetParameter(p, fit_b->GetParameter(p - 3));
    //h->Fit(fit, "IEMR", "same");
    h_data->Fit(fit_sb, "I", "same", fitMin, fitMax);
    TF1* fit_sb_save = h_data->GetFunction("fit_sb");
    fit_sb_save = new TF1(*fit_sb_save);
    h_data->Fit(fit_b, "I0", "", fitMin, fitMax); // Why again?

    TF1* fit_s = new TF1("fit_sb","gaus", fitMin, fitMax); //should be TF1("fit_s",..) ?
    for(int p = 0; p <= 2; p++)
      fit_s->SetParameter(p, fit_sb->GetParameter(p));

    h_mc->SetLineColor(kBlue);
    h_mc->Draw("h0 same");
    //h_mc->Print("all");
    fit_sb_save->SetLineColor(kRed);
    fit_sb_save->Draw("same");
    h_data->Draw("e0 same");

    hr->Draw("axis same");

    //Second pad
    TPad* pad2 = new TPad("", "", 0.00, 0.00, 1.00, 0.40);
    pad2->SetMargin(0.11, 0.03, 0.20, 0.005);
    TH2F* hrr = new TH2F("", "", 1, 100.0, 180.0, 1, -49.99, 49.99);
    hrr->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hrr->GetXaxis()->SetTitleOffset(5.0); //not working!! needs a fix somewhere
    hrr->GetYaxis()->SetTitle("Events / 1.5 GeV");
    hrr->GetYaxis()->SetTitleOffset(3.0);
    hrr->GetYaxis()->SetNdivisions(405);
    SetHistoAxisFonts(hrr);
    ScaleHistoFonts(hrr, 0.5);

    pad2->cd();
    hrr->Draw();

    //errors
    TGraphErrors* gunc68 = NULL;
    TGraphErrors* gunc95 = NULL;
    TH1D* hsig = CalculateRatio(h_data, fit_b, gunc68, gunc95);
    gunc68->SetFillColor(kGreen);
    gunc68->SetLineStyle(2);
    gunc68->SetLineColor(kRed);
    gunc95->SetFillColor(kYellow);
    gunc95->SetLineStyle(2);
    gunc95->SetLineColor(kRed);
    gunc95->Draw("4");
    gunc68->Draw("4");
    //set zero bkg fit for second pad
    TF1* fit_b0 = new TF1("fit_b0","pol5", fitMin, fitMax);
    for(int p = 3; p <= 9; p++)
      fit_b0->SetParameter(p - 3, 0.0);
    fit_b0->SetLineColor(kRed);
    fit_b0->SetLineStyle(2);
    fit_b0->Draw("same");
    fit_s->SetLineColor(kRed);
    fit_s->Draw("same");
    h_mc->Draw("h0 same");
    hsig->Draw("e0 same");
    hrr->Draw("axis same");

    pad1->cd();
    
    //legend 
    TLegend* leg = new TLegend(0.65, 0.50, 0.88, 0.88);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0.0);
    leg->AddEntry(h_data, "Data 2012", "pe");
    leg->AddEntry(fit_sb_save, "S+B fits (sum)", "l");
    leg->AddEntry(fit_b0, "B component", "l");
    leg->AddEntry(gunc68, "#pm 1#sigma", "lf");
    leg->AddEntry(gunc95, "#pm 2#sigma", "lf");
    leg->AddEntry(h_mc, "SM prediction", "l");
    leg->Draw();
    
    TLegend* legTitle = new TLegend(0.06, 0.80, 0.6, 0.93);
    legTitle->SetFillStyle(0);
    legTitle->SetBorderSize(0.0);
    legTitle->AddEntry((TObject*)NULL, vTitle[pad - 1], "");
    legTitle->Draw();

    c->cd(pad);
    pad1->Draw();
    pad2->Draw();
    
    //Text for bg component deletion
    TPaveText* pt = new TPaveText(0.55, 0.30, 0.92, 0.38);
    pt->SetBorderSize(0.0);
    pt->SetFillStyle(0);
    pt->AddText("B component subtracted");
    pt->Draw();
  }

  c->SaveAs(plotDir + "/hgg-2012.pdf");
  c->SaveAs(plotDir + "/hgg-2012.eps");
}

//plot 2011 + 2012 data for all classes combined
//can be compared to: CMS-HIG-13-001
int PlotCombFit()
{
  // set user style
  Style();

  // directory with input histograms
  TString baseDir = gHistDir;
  // directory for output plots (must exist)
  TString plotDir = gPlotsDir;

  // 2012
  TCanvas* c = new TCanvas("", "", 600, 600);
  //c->SetMargin(0.0, 0.0, 0.0, 0.0);

  TPad* pad1 = new TPad("", "", 0.00, 0.40, 1.00, 1.00);
  pad1->SetMargin(0.11, 0.03, 0.005, 0.07);
  TH2F* hr = new TH2F("", "", 1, 100.0, 180.0, 1, 0.01, 7500.0);
  hr->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  hr->GetYaxis()->SetTitle("Events / GeV");
  //hr->GetYaxis()->SetTitleOffset(0.0);
  SetHistoAxisFonts(hr);
  //ScaleHistoFonts(hr, 1.5);
  pad1->cd();
  hr->Draw();

  TFile* f_data12 = TFile::Open(TString::Format("%s/data2012_10GeV.root", baseDir.Data()));
  TFile* f_data11 = TFile::Open(TString::Format("%s/data2011_10GeV.root", baseDir.Data()));
  
  TH1D* h_data12 = (TH1D*)f_data12->Get(TString::Format("h_mgg2012_1_1"));
  TH1D* h_data = (TH1D*)f_data11->Get(TString::Format("h_mgg2011_1_1"));
  h_data->Add(h_data12);
  TFile* f_mc = TFile::Open(TString::Format("%s/mcSigReco.root", baseDir.Data()));
  TH1D* h_mc = (TH1D*)f_mc->Get(TString::Format("h_mgg1"));

  h_data->SetLineColor(1);
  h_data->SetMarkerColor(1);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.6);
  h_data->SetTitle("");
  h_data->GetXaxis()->SetTitle("");
  h_data->GetYaxis()->SetTitle("");

  // fit
  double fitMin = 100.0;
  double fitMax = 180.0;
  // background only gives:
  // pol5 gives chi2/dif = 260/74
  // pol7 gives chi2/dof = 93/72
  TF1* fit_b = new TF1("fit_b","pol5", fitMin, fitMax);
  h_data->Fit(fit_b, "IEMR", "same", fitMin, fitMax);
  TF1* fit_sb = new TF1("fit_sb","gaus(0)+pol5(3)", fitMin, fitMax);
  fit_sb->SetParameter(0, 0.0);
  fit_sb->SetParameter(1, 125.0);
  fit_sb->SetParameter(2, 2.0);
  //fit_sb->FixParameter(0, 0.0);
  fit_sb->FixParameter(1, 125.0);
  //fit_sb->FixParameter(2, 2.0);
  for(int p = 3; p <= 9; p++)
    fit_sb->SetParameter(p, fit_b->GetParameter(p - 3));
  //h->Fit(fit, "IEMR", "same");
  h_data->Fit(fit_sb, "I", "same", fitMin, fitMax);
  TF1* fit_sb_save = h_data->GetFunction("fit_sb");
  fit_sb_save = new TF1(*fit_sb_save);
  h_data->Fit(fit_b, "I0", "", fitMin, fitMax);

  TF1* fit_s = new TF1("fit_sb","gaus", fitMin, fitMax);
  for(int p = 0; p <= 2; p++)
    fit_s->SetParameter(p, fit_sb->GetParameter(p));

  h_mc->SetLineColor(kBlue);
  h_mc->Draw("h0 same");
  //h_mc->Print("all");
  fit_sb_save->SetLineColor(kRed);
  fit_sb_save->Draw("same");
  h_data->Draw("e0 same");

  hr->Draw("axis same");

  // ratio
  TPad* pad2 = new TPad("", "", 0.00, 0.00, 1.00, 0.40);
  pad2->SetMargin(0.11, 0.03, 0.20, 0.005);
  TH2F* hrr = new TH2F("", "", 1, 100.0, 180.0, 1, -199.99, 199.99);
  hrr->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  hrr->GetYaxis()->SetTitle("Events / GeV");
  hrr->GetXaxis()->SetTitleOffset(2.5);
  hrr->GetYaxis()->SetTitleOffset(1.5);
  hrr->GetYaxis()->SetNdivisions(405);
  SetHistoAxisFonts(hrr);
  pad2->cd();
  hrr->Draw();
  TGraphErrors* gunc68 = NULL;
  TGraphErrors* gunc95 = NULL;
  TH1D* hsig = CalculateRatio(h_data, fit_b, gunc68, gunc95);
  //gunc95->Print("all");
  gunc68->SetFillColor(kGreen);
  gunc68->SetLineStyle(2);
  gunc68->SetLineColor(kRed);
  gunc95->SetFillColor(kYellow);
  gunc95->SetLineStyle(2);
  gunc95->SetLineColor(kRed);
  gunc95->Draw("4");
  gunc68->Draw("4");
  TF1* fit_b0 = new TF1("fit_b0","pol2", fitMin, fitMax);
  for(int p = 3; p <= 5; p++)
    fit_b0->SetParameter(p - 3, 0.0);
  fit_b0->SetLineColor(kRed);
  fit_b0->SetLineStyle(2);
  fit_b0->Draw("same");
  fit_s->SetLineColor(kRed);
  fit_s->Draw("same");
  h_mc->Draw("h0 same");
  hsig->Draw("e0 same");
  hrr->Draw("axis same");

  pad1->cd();
  TLegend* leg = new TLegend(0.50, 0.40, 0.88, 0.88);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0.0);
  leg->AddEntry(h_data, "Data 2011 + Data 2012", "pe");
  leg->AddEntry(fit_sb_save, "S+B fits (sum)", "l");
  leg->AddEntry(fit_b0, "B component", "l");
  leg->AddEntry(gunc68, "#pm 1#sigma", "lf");
  leg->AddEntry(gunc95, "#pm 2#sigma", "lf");
  leg->AddEntry(h_mc, "SM prediction", "l");
  leg->Draw();

  c->cd();
  pad1->Draw();
  pad2->Draw();
  TPaveText* pt = new TPaveText(0.55, 0.30, 0.92, 0.38);
  pt->SetBorderSize(0.0);
  pt->SetFillStyle(0);
  pt->AddText("B component subtracted");
  pt->Draw();
  TPaveText* ptlum = new TPaveText(0.15,0.8,0.5,1);
  ptlum->SetBorderSize(0.0);
  ptlum->SetFillStyle(0);
  ptlum->AddText("11.58 fb^{-1} (8 TeV) + 2.33 fb^{-1} (7 TeV)");
  ptlum->Draw();
  c->SaveAs(plotDir + "/hgg-comb.pdf");
  c->SaveAs(plotDir + "/hgg-comb.eps");
}

//plot 2011 + 2012 data for each class
// no comparison plot existing
int PlotCombFitClasses()
{
  // set user style
  Style();

  // directory with input histograms
  TString baseDir = gHistDir;
  // directory for output plots (must exist)
  TString plotDir = gPlotsDir;

  //Canvas
  TCanvas* c = new TCanvas("", "", 600, 600);
  c->Divide(2,2);
  
  //Titles for all pads : event classes
  std::vector<TString> vTitle;
  vTitle.push_back("a) Both #gamma in barrel, R_{9}^{min}>0.94");
  vTitle.push_back("b) Both #gamma in barrel, R_{9}^{min}<0.94");
  vTitle.push_back("c) One or both #gamma in endcap, R_{9}^{min}>0.94");
  vTitle.push_back("d) One or both #gamma in endcap, R_{9}^{min}<0.94");
  
  //loop over all pads/eventclasses
  for(int pad = 1; pad <= 4; pad++)
  {
    //change to pad
    c->cd(pad);
    TPad* pad1 = new TPad("", "", 0.00, 0.40, 1.00, 1.00);
    pad1->SetMargin(0.11, 0.03, 0.005, 0.07);
    
    //set maximum value yLabel
    double yMax = 2000;
    if(pad == 1) yMax = 1000;
    if(pad == 2) yMax = 3500;
    if(pad == 3) yMax = 1000;
    if(pad == 4) yMax = 4500;
    
    
    TH2F* hr = new TH2F("", "", 1.0, 100.0, 180.0,1,0.01,yMax);
    hr->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    hr->GetXaxis()->SetTitleOffset(3.0);
    hr->GetYaxis()->SetTitle("Events / GeV");
    hr->GetYaxis()->SetTitleOffset(3.0);
    SetHistoAxisFonts(hr);
    ScaleHistoFonts(hr, 0.5);
    pad1->cd();
    hr->Draw();

    TFile* f_data12 = TFile::Open(TString::Format("%s/data2012_10GeV.root", baseDir.Data()));
    TFile* f_data11 = TFile::Open(TString::Format("%s/data2011_10GeV.root", baseDir.Data()));
    TH1D* h_data12 = (TH1D*)f_data12->Get(TString::Format("h_mgg2012_1_%d",pad+2));
    TH1D* h_data = (TH1D*) f_data11->Get(TString::Format("h_mgg2011_1_%d", pad+2));
    h_data->Add(h_data12);

    TFile* f_mc = TFile::Open(TString::Format("%s/mcSigReco.root", baseDir.Data()));
    TH1D* h_mc = (TH1D*)f_mc->Get(TString::Format("h_mgg%d",pad+2));

    h_data->SetLineColor(1);
    h_data->SetMarkerColor(1);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.3);
    h_data->SetTitle("");
    h_data->GetXaxis()->SetTitle("");
    h_data->GetYaxis()->SetTitle("");

    // fit
    double fitMin = 100.0;
    double fitMax = 180.0;
    // background only gives:
    // pol5 gives chi2/dif = 260/74
    // pol7 gives chi2/dof = 93/72
    //pol2 fit (background only)
    TF1* fit_b = new TF1("fit_b","pol5", fitMin, fitMax);
    h_data->Fit(fit_b, "IEMR", "same", fitMin, fitMax);
    //pol2 + gaus fit (signal + background)
    TF1* fit_sb = new TF1("fit_sb","gaus(0)+pol5(3)", fitMin, fitMax);
    fit_sb->SetParameter(0, 0.0);
    fit_sb->SetParameter(1, 125.0);
    fit_sb->SetParameter(2, 2.0);
    //Fix only mass parameter
    //fit_sb->FixParameter(0, 0.0);
    fit_sb->FixParameter(1, 125.0);
    //fit_sb->FixParameter(2, 2.0);
    //Set other parameters from background
    for(int p = 3; p <= 9; p++)
      fit_sb->SetParameter(p, fit_b->GetParameter(p - 3));
    //h->Fit(fit, "IEMR", "same");
    h_data->Fit(fit_sb, "I", "same", fitMin, fitMax);
    TF1* fit_sb_save = h_data->GetFunction("fit_sb");
    fit_sb_save = new TF1(*fit_sb_save);
    h_data->Fit(fit_b, "I0", "", fitMin, fitMax);

    TF1* fit_s = new TF1("fit_sb","gaus", fitMin, fitMax);
    for(int p = 0; p <= 2; p++)
      fit_s->SetParameter(p, fit_sb->GetParameter(p));

    h_mc->SetLineColor(kBlue);
    h_mc->Draw("h0 same");
    //h_mc->Print("all");
    fit_sb_save->SetLineColor(kRed);
    fit_sb_save->Draw("same");
    h_data->Draw("e0 same");

    hr->Draw("axis same");
    //Second pad
    // ratio
    TPad* pad2 = new TPad("", "", 0.00, 0.00, 1.00, 0.40);
    pad2->SetMargin(0.11, 0.03, 0.20, 0.005);
    TH2F* hrr = new TH2F("", "", 1, 100.0, 180.0, 1, -49.99, 49.99);
    hrr->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    hrr->GetXaxis()->SetTitleOffset(5.0);
    hrr->GetYaxis()->SetTitle("Events / GeV");
    hrr->GetYaxis()->SetTitleOffset(3.0);
    hrr->GetYaxis()->SetNdivisions(405);
    
    SetHistoAxisFonts(hrr);
    ScaleHistoFonts(hrr, 0.5);
    pad2->cd();
    hrr->Draw();
    TGraphErrors* gunc68 = NULL;
    TGraphErrors* gunc95 = NULL;
    TH1D* hsig = CalculateRatio(h_data, fit_b, gunc68, gunc95);
    //gunc95->Print("all");
    gunc68->SetFillColor(kGreen);
    gunc68->SetLineStyle(2);
    gunc68->SetLineColor(kRed);
    gunc95->SetFillColor(kYellow);
    gunc95->SetLineStyle(2);
    gunc95->SetLineColor(kRed);
    gunc95->Draw("4");
    gunc68->Draw("4");
    TF1* fit_b0 = new TF1("fit_b0","pol5", fitMin, fitMax);
    for(int p = 3; p <= 9; p++)
      fit_b0->SetParameter(p - 3, 0.0);
    fit_b0->SetLineColor(kRed);
    fit_b0->SetLineStyle(2);
    fit_b0->Draw("same");
    fit_s->SetLineColor(kRed);
    fit_s->Draw("same");
    h_mc->Draw("h0 same");
    hsig->Draw("e0 same");
    hrr->Draw("axis same");

    pad1->cd();
    
    //add errors 
    
    //legend 
    TLegend* leg = new TLegend(0.65, 0.50, 0.88, 0.88);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0.0);
    leg->AddEntry(h_data, "Data 2011 + Data 2012", "pe");
    leg->AddEntry(fit_sb_save, "S+B fits (sum)", "l");
    leg->AddEntry(fit_b0, "B component", "l");
    leg->AddEntry(gunc68, "#pm 1#sigma", "lf");
    leg->AddEntry(gunc95, "#pm 2#sigma", "lf");
    leg->AddEntry(h_mc, "SM prediction", "l");
    leg->Draw();
    
    TLegend* legTitle = new TLegend(0.06, 0.80, 0.6, 0.93);
    legTitle->SetFillStyle(0);
    legTitle->SetBorderSize(0.0);
    legTitle->AddEntry((TObject*)NULL, vTitle[pad - 1], "");
    legTitle->Draw();

    c->cd(pad);
    pad1->Draw();
    pad2->Draw();
    
    //Add text for bg component deletion
    TPaveText* pt = new TPaveText(0.55, 0.30, 0.92, 0.38);
    pt->SetBorderSize(0.0);
    pt->SetFillStyle(0);
    pt->AddText("B component subtracted");
    pt->Draw();
  }

  c->SaveAs(plotDir + "/hgg-comb-classes.pdf");
  c->SaveAs(plotDir + "/hgg-comb-classes.eps");
}


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>> Main function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char** argv)
{
  //plot all
  Plot2011();
  Plot2012();
  PlotCombFit();
  PlotCombFitClasses();
  return 0;
}
