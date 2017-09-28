#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TStyle.h"
#include "TArrow.h"

#include "include/getLinBins.h"

void prettyCanv(TCanvas* canv_p)
{
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetBottomMargin(canv_p->GetLeftMargin());

  return;
}


int solveFixedPoint()
{
  gStyle->SetOptStat(0);

  const double vx = 0.1;

  const double kx = 4.;
  const double ky = 2.;

  //  const double gammax = 10;

  TFile* outFile_p = new TFile("outFile.root", "RECREATE");
  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  TH1F* dummyHist_p = new TH1F("dummyHist_h", ";;", 3, 0, 4);
  dummyHist_p->SetMaximum(4);
  dummyHist_p->SetMinimum(0);
  dummyHist_p->DrawCopy("");

  const int nLinBins = 100000;
  const Float_t low = .1;
  const Float_t hi = 4.;
  Double_t linBins[nLinBins+1];
  getLinBins(low, hi, nLinBins, linBins);

  TH1F* yNull_p = new TH1F("yNull_h", ";x;y", nLinBins, linBins);
  TH1F* xNull_p = new TH1F("xNull_h", ";x;y", nLinBins, linBins);

  for(int i = 0; i < xNull_p->GetNbinsX(); ++i){
    double xVal = yNull_p->GetBinCenter(i+1);

    yNull_p->SetBinContent(i+1, ky*xVal);
    yNull_p->SetBinError(i+1, 0);

    xNull_p->SetBinContent(i+1, kx*xVal*xVal/((xVal - vx)*(1 + xVal*xVal)) - 1);
    xNull_p->SetBinError(i+1, 0);
  }
  
  yNull_p->SetMarkerSize(0.5);
  yNull_p->SetMarkerStyle(20);
  yNull_p->SetMarkerColor(kBlue);

  xNull_p->SetMarkerSize(0.5);
  xNull_p->SetMarkerStyle(20);
  xNull_p->SetMarkerColor(kRed);

  yNull_p->Draw("P SAME");
  xNull_p->Draw("P SAME");

  outFile_p->cd();
  canv_p->Write("", TObject::kOverwrite);

  canv_p->SaveAs("canvFixedPoint.pdf");
  delete canv_p;
  delete yNull_p;
  delete xNull_p;
  delete dummyHist_p;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal = solveFixedPoint();
  return retVal;
}
