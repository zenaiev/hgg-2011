#include "settings.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

TGraphAsymmErrors* HtoGragh(TH1* h, double xpos = 0.5, TH1* hline = NULL)
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors;
  for(int b = 0; b < h->GetNbinsX(); b++)
  {
    double xmin = h->GetXaxis()->GetBinLowEdge(b + 1);
    double xmax = h->GetXaxis()->GetBinUpEdge(b + 1);
    double bw = xmax - xmin;
    double x = xmin + xpos * bw;
    double y = h->GetBinContent(b + 1);
    double yer = h->GetBinError(b + 1);
    int point = g->GetN();
    g->SetPoint(point, x, y);
    double xerl = 0.0;
    double xerh = 0.0;
    g->SetPointError(point, 0.0, 0.0, yer, yer);
  }
  return g;
}

class ZPlotCSInput
{
  public:
    // directories
    TString baseDir;
    TString plotDir;
    // MC backgrounds
    std::vector<TString> VecMCBackgr;
    // variables
    std::vector<TH2F*> VecHR;
    std::vector<TString> VecVar;
    // channels
    std::vector<TString> VecTitle;
    std::vector<int> VecColor;
    std::vector<int> VecStyle;
    // absolute or normalised
    bool Norm;
    // reference from paper
    bool Paper;
};

// paper cs ptt
void GetPaperCS(const TString& var, TGraphAsymmErrors* gstat, TGraphAsymmErrors* gsyst)
{
  double *sig, *stat, *syst;
  if(var == "ptt")
  {
    sig = new double[5] {5.10e-3,6.26e-3,2.96e-3,0.70e-3,0.12e-3};
    stat = new double[5] {2.2    ,2.6    ,2.6    ,3.5    ,7.5    };
    syst = new double[5] {5.6    ,3.9    ,4.9    ,6.2    ,5.4    };
  }
  if(var == "yt")
  {
    sig = new double[8] {0.091,0.255,0.302,0.351,0.371,0.306,0.241,0.090};
    stat = new double[8] {4.1,3.1,3.3,3.2,3.2,3.4,3.3,4.0};
    syst = new double[8] {6.8,4.9,4.0,3.8,3.8,4.0,4.9,6.8};
  }
  if(var == "ytt")
  {
    sig = new double[6] {0.03,0.219,0.418,0.393,0.218,0.04};
    stat = new double[6] {13.7,3.1,2.4,2.5,3.2,10.4};
    syst = new double[6] {14.7,4.4,3.6,3.6,4.4,14.7};
  }
  if(var == "pttt")
  {
    sig = new double[4] {1.6e-2,0.97e-2,0.32e-2,0.05e-2};
    stat = new double[4] {2.9,2.1,2.5,3.7};
    syst = new double[4] {24.9,10.7,13.2,6.9};
  }
  if(var == "mtt")
  {
    sig = new double[8] {0.0,5.26e-3,4.58e-3,2.46e-3,1.07e-3,0.39e-3,0.08e-3,0.01e-3};
    stat = new double[8] {0.0,5.4,3.8,4.9,6.1,6.9,13.3,22.4};
    syst = new double[8] {0.0,10.4,4.1,7.6,3.9,11.4,27.0,43.6};
  }
  for(int p = 0; p < gstat->GetN(); p++)
  {
    gstat->SetPoint(p, (gstat->GetX())[p], sig[p]);
    gsyst->SetPoint(p, (gsyst->GetX())[p], sig[p]);
    double tot = TMath::Sqrt(stat[p] * stat[p] + syst[p] * syst[p]);
    gstat->SetPointError(p, 0.0, 0.0, (gstat->GetY())[p] * stat[p] / 100.0, (gstat->GetY())[p] * stat[p] / 100.0);
    gsyst->SetPointError(p, 0.0, 0.0, (gsyst->GetY())[p] * syst[p] / 100.0, (gsyst->GetY())[p] * tot / 100.0);
  }
}


void PlotCS(const ZPlotCSInput& in)
{
  TCanvas* c_cs;
  //c_cs = new TCanvas(TString::Format("ccs"), "", 800, 1200);
  //c_cs->Divide(2,3, 0.0001);
  c_cs = new TCanvas(TString::Format("ccs"), "", 1200, 800);
  c_cs->Divide(3,2, 0.0001);
  for(int v = 0; v < 5; v++)
  {
    c_cs->cd(v + 1);
    if(v == 4)
      gPad->SetLogy();
    in.VecHR[v]->Draw();
    TString var = in.VecVar[v];
    TH1D *hcombsig, *hcombreco, *hcombgen;
    double nsig = 0;
    double nreco = 0;
    double ngen = 0;
    TLegend* leg = new TLegend(0.46, 0.68, 0.92, 0.90);
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (int ch = 1; ch < 4; ch++)
    {
      // MC backgr
      TH1D* hbackgr;
      for(int mc = 0; mc < in.VecMCBackgr.size(); mc++)
      {
        TFile* f = TFile::Open(TString::Format("%s/mc%sReco-c%d.root", in.baseDir.Data(), in.VecMCBackgr[mc].Data(), ch));        
        TH1D* h = (TH1D*)f->Get(TString::Format("h_%s_cs", var.Data()));
        if(mc == 0)
          hbackgr = new TH1D(*h);
        else
          hbackgr->Add(h);
      }
      //hbackgr->Print("all");
      // data
      TFile* f = TFile::Open(TString::Format("%s/data-c%d.root", in.baseDir.Data(), ch));
      TH1D* hsig = new TH1D(*(TH1D*)f->Get(TString::Format("h_%s_cs", var.Data())));
      //hsig->Print("all");
      hsig->Add(hbackgr, -1.0);
      nsig += hsig->Integral(0, hsig->GetNbinsX());
      if(ch == 1)
        hcombsig = new TH1D(*hsig);
      else
        hcombsig->Add(hsig);
      // MC reco
      f = TFile::Open(TString::Format("%s/mcSigReco-c%d.root", in.baseDir.Data(), ch));
      TH1D* hacc = new TH1D(*(TH1D*)f->Get(TString::Format("h_%s_cs", var.Data())));
      nreco += hacc->Integral(0, hacc->GetNbinsX());
      if(ch == 1)
        hcombreco = new TH1D(*hacc);
      else
        hcombreco->Add(hacc);
      // MC gen
      f = TFile::Open(TString::Format("%s/mcSigGen-c%d.root", in.baseDir.Data(), ch));
      TH1D* hgen = (TH1D*)f->Get(TString::Format("h_%s_cs", var.Data()));
      ngen += hgen->Integral(0, hgen->GetNbinsX());
      if(ch == 1)
        hcombgen = new TH1D(*hgen);
      else
        hcombgen->Add(hgen);
      // acceptance
      hacc->Divide(hgen);
      // cross section
      TH1D* hcs = hsig;
      hcs->Divide(hacc);
      double cs, cserr;
      cs = hcs->IntegralAndError(0, hcs->GetNbinsX() + 1, cserr);
      double br = 0.0115 * 2500.0;
      if(ch == 3)
        br *= 2.0;
      printf("CS c%d %s:  %.1f +- %.1f\n", ch, var.Data(), cs / br, cserr / br);
      if(in.Norm)
        hcs->Scale(1.0 / hcs->Integral(), "width");
      //hcs->Print("all");
      //return 1;
      // draw
      //TGraphAsymmErrors* gcs = HtoGragh(hcs, (ch == 1) ? 0.2 : ((ch == 2) ? 0.35 : 0.65));
      TGraphAsymmErrors* gcs = HtoGragh(hcs, 0.10 + 0.1 * ch);
      gcs->SetLineColor(in.VecColor[ch]);
      gcs->SetMarkerStyle(in.VecStyle[ch]);
      gcs->SetMarkerColor(in.VecColor[ch]);
      gcs->SetMarkerSize(0.9);
      leg->AddEntry(gcs, in.VecTitle[ch], "p");
      gcs->Draw("pz0");
    }
    // combined
    TH1D* hcombacc = hcombreco;
    hcombacc->Divide(hcombgen);
    TH1D* hcombcs = hcombsig;
    hcombcs->Divide(hcombacc);
    double cs, cserr;
    cs = hcombcs->IntegralAndError(0, hcombcs->GetNbinsX() + 1, cserr);
    double br = 0.046 * 2500.0;
    //printf("CS dilepton %s:  %.1f +- %.1f (%.1f)\n", var.Data(), cs / br, cserr / br, nsig / nreco * ngen / br);
    printf("CS dilepton %s:  %.1f +- %.1f\n", var.Data(), cs / br, cserr / br, nsig / nreco * ngen / br);
    if(in.Norm)
      hcombcs->Scale(1.0 / hcombcs->Integral(), "width");
    TGraphAsymmErrors* gcombcs = HtoGragh(hcombcs);
    gcombcs->SetLineColor(in.VecColor[0]);
    gcombcs->SetMarkerStyle(in.VecStyle[0]);
    gcombcs->SetMarkerColor(in.VecColor[0]);
    gcombcs->SetMarkerSize(0.9);
    // plot data
    leg->AddEntry(gcombcs, in.VecTitle[0], "pel");
    gcombcs->Draw("pz0");
    hcombcs->SetLineColor(1);
    hcombcs->Draw("hist same");
    hcombcs->Print("all");
    // plot paper results
    if(in.Paper)
    {
      TGraphAsymmErrors* gstat = HtoGragh(hcombcs, 0.7);
      gstat->SetLineColor(1);
      gstat->SetMarkerStyle(24);
      gstat->SetMarkerColor(1);
      gstat->SetMarkerSize(0.9);
      TGraphAsymmErrors* gsyst = HtoGragh(hcombcs, 0.7);
      gsyst->SetLineColor(1);
      gsyst->SetMarkerStyle(24);
      gsyst->SetMarkerColor(1);
      gsyst->SetMarkerSize(0.9);
      GetPaperCS(var, gstat, gsyst);
      gstat->Draw("p0");
      gsyst->Draw("pz0");
      //gstat->Print("all");
      leg->AddEntry(gstat, "CMS-TOP-11-013", "pe");
    }
    leg->Draw();
    in.VecHR[v]->Draw("axis same");
  }
  // save plot
  TString name = TString::Format("%s/cs%s", in.plotDir.Data(), (in.Norm) ? "_norm" : "");
  c_cs->SaveAs(name + ".eps");
  c_cs->SaveAs(name + ".pdf");
}
