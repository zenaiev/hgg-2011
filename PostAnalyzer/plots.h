#include <TGaxis.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include "TVirtualFitter.h"
#include <TStyle.h>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>> Prepare plot style >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Modify as you want, if needed consult
// https://root.cern.ch/doc/master/classTStyle.html
//
void Style()
{
  gStyle->SetOptStat(000000000);
  gStyle->SetTitle(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetPadGridX(1);
  //gStyle->SetPadGridY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetEndErrorSize(5);
  TGaxis::SetMaxDigits(3);
  gStyle->SetErrorX(0.0);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.08);
  //gStyle->SetNdivisions(206, "xyz");
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void ScaleAxisFonts(TAxis* axis, const double scale)
{
  axis->SetTitleSize(scale * axis->GetTitleSize());
  axis->SetLabelSize(scale * axis->GetLabelSize());
}

void ScaleHistoFonts(TH1* h, const double scale)
{
  TAxis* xaxis = h->GetXaxis();
  ScaleAxisFonts(xaxis, scale);
  TAxis* yaxis = h->GetYaxis();
  ScaleAxisFonts(yaxis, scale);
}

void SetAxisFonts(TAxis* axis, const int font = 63, const int size = 20 )
{
  axis->SetTitleFont(font);
  axis->SetTitleSize(size);
  axis->SetLabelFont(font);
  axis->SetLabelSize(size);
}

void SetHistoAxisFonts(TH1* h, const int font = 63)
{
  TAxis* xaxis = h->GetXaxis();
  SetAxisFonts(xaxis, font);
  TAxis* yaxis = h->GetYaxis();
  SetAxisFonts(yaxis, font);
}

TH1D* CalculateRatio(const TH1D* hsum, const TF1* fb, TGraphErrors*& gunc68, TGraphErrors*& gunc95)
{
  // subtracted background function (0 +- unc.)
  //fb0 = new TF1(*fb);
  // create a TGraphErrors to hold the confidence intervals
  const int ngr = 60;
  gunc68 = new TGraphErrors(ngr);
  gunc95 = new TGraphErrors(ngr);
  //gunc->SetTitle("Fitted line with .95 conf. band");
  for (int i = 0; i < ngr; i++)
  {
    double xmin = hsum->GetXaxis()->GetBinCenter(1);
    double xmax = hsum->GetXaxis()->GetBinCenter(hsum->GetNbinsX());
    double x = (i + 0.5) / ngr * (xmax - xmin) + xmin;
    gunc68->SetPoint(i, x, 0);
    gunc95->SetPoint(i, x, 0);
  }
  /*Compute the confidence intervals at the x points of the created graph*/
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gunc68, 0.68);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gunc95);
  //Now the graph contains function values as its y-coordinates
  //and confidence intervals as the errors on these coordinates
  //Draw the graph, the function and the confidence intervals

  TH1D* hsig = new TH1D(*hsum);
  for(int b = 1; b <= hsum->GetNbinsX(); b++)
  {
    double x = hsum->GetXaxis()->GetBinCenter(b);
    double bg = fb->Eval(x);
    double sig = hsum->GetBinContent(b) - bg;
    // remove points outside of fit range
    if(x < fb->GetXmin() || x > fb->GetXmax())
      sig = 9e9;
    hsig->SetBinContent(b, sig);
  }

  for (int i = 0; i < ngr; i++)
  {
    double x = gunc95->GetX()[i];
    double bg = fb->Eval(x);
    gunc68->GetY()[i] -= bg;
    gunc95->GetY()[i] -= bg;
    //gunc->SetPointError(i, 0.0, gunc->GetErrorY(i) - bg);
    // remove points outside of fit range
    if(x < fb->GetXmin() || x > fb->GetXmax())
    {
      //gunc68->RemovePoint(i);
      //gunc95->RemovePoint(i);
      //gunc68->GetX()[i] = 0;
      //gunc95->GetX()[i] = 0;
      gunc68->SetPointError(i, 0.0, 0.0);
      gunc95->SetPointError(i, 0.0, 0.0);
    }
  }

  return hsig;
}
