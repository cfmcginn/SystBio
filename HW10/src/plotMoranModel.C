#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"

int plotMoranModel(const std::string inFileName, const int nPop)
{
  if(nPop != 10 && nPop != 500) return 1;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  TH1F* hist_UE_p = (TH1F*)inFile_p->Get(("sGen_nSim10000_NPop" + std::to_string(nPop) + "_mutRate0p000000_h").c_str());
  TH1F* hist1_p = (TH1F*)inFile_p->Get(("sGen_Win_nSim10000_NPop" + std::to_string(nPop) + "_mutRate0p000000_h").c_str());
  TH1F* hist2_p = (TH1F*)inFile_p->Get(("sGen_Win_nSim10000_NPop" + std::to_string(nPop) + "_mutRate0p000020_h").c_str());
  TH1F* hist3_p = (TH1F*)inFile_p->Get(("sGen_Win_nSim10000_NPop" + std::to_string(nPop) + "_mutRate0p010000_h").c_str());

  kirchnerPalette col;

  hist_UE_p->GetXaxis()->CenterTitle();
  hist_UE_p->GetYaxis()->CenterTitle();

  hist_UE_p->GetXaxis()->SetTitle("S of Winning Mutants");

  hist_UE_p->SetMarkerStyle(20);
  hist1_p->SetMarkerStyle(20);
  hist2_p->SetMarkerStyle(20);
  hist3_p->SetMarkerStyle(20);

  hist_UE_p->SetMarkerSize(1);
  hist1_p->SetMarkerSize(1);
  hist2_p->SetMarkerSize(1);
  hist3_p->SetMarkerSize(1);

  hist_UE_p->SetMarkerColor(col.getColor(0));
  hist1_p->SetMarkerColor(col.getColor(2));
  hist2_p->SetMarkerColor(col.getColor(3));
  hist3_p->SetMarkerColor(col.getColor(4));

  hist_UE_p->SetLineColor(col.getColor(0));
  hist1_p->SetLineColor(col.getColor(2));
  hist2_p->SetLineColor(col.getColor(3));
  hist3_p->SetLineColor(col.getColor(4));

  gStyle->SetOptStat(0);

  hist_UE_p->SetMaximum(.12);
  hist_UE_p->SetMinimum(0.);

  hist_UE_p->DrawCopy("E1 P");
  hist1_p->DrawCopy("E1 P SAME");
  hist2_p->DrawCopy("E1 P SAME");
  hist3_p->DrawCopy("E1 P SAME");

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);

  label_p->DrawLatex(.15, .95, ("N_{Pop}=" + std::to_string(nPop) + "; Simulate 10000 times").c_str());

  TLegend* leg_p = new TLegend(0.5, 0.9, 0.5, 0.9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);

  leg_p->AddEntry(hist_UE_p, "Sampling Distrib.(s_{0}=0.02)", "P L");
  leg_p->AddEntry(hist1_p, "#mu=0.0", "P L");
  leg_p->AddEntry(hist2_p, "#mu=0.00002", "P L");
  leg_p->AddEntry(hist3_p, "#mu=0.01", "P L");

  leg_p->Draw("SAME");

  canv_p->SaveAs(("pdfDir/plotMoran_NPop" + std::to_string(nPop) + ".pdf").c_str());

  delete leg_p;
  delete label_p;
  delete canv_p;
  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./plotMoranModel.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += plotMoranModel(argv[1], std::stoi(argv[2]));
  return retVal;
}
