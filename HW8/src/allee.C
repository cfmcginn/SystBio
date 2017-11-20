#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"

#include "include/getLogBins.h"
#include "include/kirchnerPalette.h"

int allee()
{
  kirchnerPalette col;

  double aParam = 5.;
  double bParam = 1.;

  const Int_t nScan = 10000;
  const double rLow = 0.001;
  const double rHi = 3*aParam*bParam*bParam;
  Double_t rBins[nScan+1];
  getLogBins(rLow, rHi, nScan, rBins);

  TGraph* graph_p = new TGraph();
  graph_p->SetMarkerStyle(20);
  graph_p->SetMarkerSize(.5);
  graph_p->SetMarkerColor(col.getColor(0));
  graph_p->SetLineColor(col.getColor(0));

  double minVal = 0;
  double maxVal = 0;

  for(Int_t i = 0; i < nScan; ++i){
    double nVal1 = bParam + TMath::Sqrt((rBins[i]/aParam));
    double nVal2 = bParam - TMath::Sqrt((rBins[i]/aParam));

    if(nVal1 < minVal) minVal = nVal1;
    if(nVal2 < minVal) minVal = nVal2;

    if(nVal1 > maxVal) maxVal = nVal1;
    if(nVal2 > maxVal) maxVal = nVal2;

    graph_p->SetPoint(2*i, rBins[i], nVal1);
    graph_p->SetPoint(2*i + 1, rBins[i], nVal2);
  }

  TH1F* dummy_p = new TH1F("dummy_h", "", nScan, rBins);
  dummy_p->SetMaximum(maxVal*1.1);
  dummy_p->SetMinimum(minVal*1.1);
  
  TCanvas* canv_p = new TCanvas("allee_c", "allee_c", 500, 500);
  dummy_p->DrawCopy();
  graph_p->Draw("P");
  gPad->SetLogx();

  canv_p->SaveAs("pdfDir/allee.pdf");

  delete canv_p;
  delete dummy_p;
  delete graph_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += allee();
  return retVal;
}
