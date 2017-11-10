#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TH1F.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TLatex.h"

#include "include/getLinBins.h"
#include "include/plotUtilities.h"

int nullClines()
{
  const double beta = 1.;
  //  const double kappa = 1.;

  const int nCBins = 200;
  const double cMin = 0;
  const double cMax = 2;
  Double_t cBins[nCBins+1];
  getLinBins(cMin, cMax, nCBins, cBins);

  const int nIBins = 200;
  const double iMin = 0;
  const double iMax = 2;
  Double_t iBins[nIBins+1];
  getLinBins(iMin, iMax, nIBins, iBins);

  TH1F* hist_p = new TH1F("hist_h", ";C;I", 10000, cMin, cMax);
  hist_p->SetMinimum(iMin);
  hist_p->SetMaximum(iMax);

  TH1F* histNullI_p = new TH1F("histNullI_h", ";C;I", 10000, cMin, cMax);
  histNullI_p->SetMinimum(iMin);
  histNullI_p->SetMaximum(iMax);

  for(Int_t i = 0; i < histNullI_p->GetNbinsX(); ++i){
    histNullI_p->SetBinContent(i+1, histNullI_p->GetBinCenter(i+1));
    histNullI_p->SetBinError(i+1, 0);
  }


  TH1F* histNullC_p = new TH1F("histNullC_h", ";C;I", 10000, cMin, cMax);
  histNullC_p->SetMinimum(iMin);
  histNullC_p->SetMaximum(iMax);

  for(Int_t i = 0; i < histNullC_p->GetNbinsX(); ++i){
    Double_t val = histNullC_p->GetBinCenter(i+1);
    if(val > 1) continue;
    val = TMath::Sqrt((val - val*val)/beta);
    if(val < 0) continue;

    histNullC_p->SetBinContent(i+1, val);
    histNullC_p->SetBinError(i+1, 0);
  }

  TArrow* line_p = new TArrow();
  line_p->SetLineStyle(1);


  TCanvas* canvNullI_p = new TCanvas("canvNullI_c", "canvNullI_c", 500, 500);
  prettyCanv(canvNullI_p);
  prettyTH1(hist_p, .6, 20, 1);
  gStyle->SetOptStat(0);
  hist_p->DrawCopy();
  histNullI_p->DrawCopy("SAME");


  for(int c = 1; c < 10; ++c){
    double xPos = cMax/10.*c;

    for(int i = 1; i < 10; ++i){
      double yPos = iMax/10.*i;
      
      if(xPos < yPos) line_p->DrawArrow(xPos-.02, yPos, xPos+.02, yPos, 0.01, ">");
      else line_p->DrawArrow(xPos-.02, yPos, xPos+.02, yPos, 0.01, "<");

    }
  }


  TCanvas* canvNullC_p = new TCanvas("canvNullC_c", "canvNullC_c", 500, 500);
  prettyCanv(canvNullC_p);
  prettyTH1(histNullC_p, .6, 20, 1);
  gStyle->SetOptStat(0);
  hist_p->DrawCopy();
  histNullC_p->DrawCopy("SAME");

  for(int c = 1; c < 10; ++c){
    double xPos = cMax/10.*c;
    double val = TMath::Sqrt((xPos - xPos*xPos)/beta);

    for(int i = 1; i < 10; ++i){
      double yPos = iMax/10.*i;
      
      if(val > yPos) line_p->DrawArrow(xPos, yPos-.02, xPos, yPos+.02, 0.01, ">");
      else line_p->DrawArrow(xPos, yPos-.02, xPos, yPos+.02, 0.01, "<");

    }
  }

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(20);
  label_p->SetNDC();

  TCanvas* canvOsc_p = new TCanvas("canvOsc_c", "canvOsc_c", 500, 500);
  prettyCanv(canvOsc_p);
  hist_p->DrawCopy();
  histNullI_p->DrawCopy("SAME");
  histNullC_p->DrawCopy("SAME");
  gStyle->SetOptStat(0);

  canvOsc_p->cd();
  canvOsc_p->SaveAs("pdfDir/nullOsc.pdf");

  canvNullI_p->cd();
  label_p->DrawLatex(.25, .95, "di/dt = 0");
  canvNullI_p->SaveAs("pdfDir/nullClineI.pdf");


  canvNullC_p->cd();
  label_p->DrawLatex(.25, .95, "dc/dt = 0");
  canvNullC_p->SaveAs("pdfDir/nullClineC.pdf");


  delete hist_p;
  delete histNullI_p;
  delete histNullC_p;
  delete canvOsc_p;
  delete canvNullI_p;
  delete canvNullC_p;

  delete label_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += nullClines();
  return retVal;
}
