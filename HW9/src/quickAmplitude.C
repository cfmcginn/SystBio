#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TStyle.h"

#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int quickAmplitude()
{
  kirchnerPalette col;

  TH1F* dummy_p = new TH1F("dummy_h", ";N;Amplitude", 10, 100, 50000);
  dummy_p->SetMaximum(.5);
  dummy_p->SetMinimum(.005);
  dummy_p->GetXaxis()->CenterTitle();
  dummy_p->GetYaxis()->CenterTitle();

  TGraph* graph_p = new TGraph();
  graph_p->SetPoint(0, 300, .28);
  graph_p->SetPoint(1, 3000, .08);
  graph_p->SetPoint(2, 30000, .02);

  graph_p->SetMarkerColor(col.getColor(0));
  graph_p->SetLineColor(col.getColor(0));
  graph_p->SetMarkerStyle(20);
  graph_p->SetMarkerSize(.8);

  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);

  dummy_p->DrawCopy();
  graph_p->Draw("P");

  gStyle->SetOptStat(0);

  gPad->SetLogx();
  gPad->SetLogy();

  canv_p->SaveAs("pdfDir/amplitude.pdf");

  delete canv_p;

  delete graph_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += quickAmplitude();
  return retVal;
}
