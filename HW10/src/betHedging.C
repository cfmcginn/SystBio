#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

#include "include/doGlobalDebug.h"
#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int betHedging(const int n0 = 100)
{
  TRandom3* randGen_p = new TRandom3(0);

  const double lLow = 0.7;
  const double lHi = 1.35;

  const int nTimeSteps = 300;
  const int nTrials = 100000;

  TFile* outFile_p = new TFile("output/betHedging.root", "RECREATE");
  TH1F* simHist_h = new TH1F("simHist_h", ";;", 100, 0, 2.*n0);
  TH1F* simHist_n0_h = new TH1F("simHist_n0_h", ";n(t=300)/n_{0};Counts", 100, 0, 2.);
  simHist_n0_h->GetXaxis()->CenterTitle();
  simHist_n0_h->GetYaxis()->CenterTitle();

  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  canv_p->cd();
  TH1F* dummyHist_h = new TH1F("dummyHist_h", ";Time (Generations);n(t)/n_{0}", 300, 0, 300);
  dummyHist_h->GetXaxis()->CenterTitle();
  dummyHist_h->GetYaxis()->CenterTitle();
  dummyHist_h->DrawCopy();
  dummyHist_h->SetMaximum(2.);
  dummyHist_h->SetMinimum(0.);
  dummyHist_h->DrawCopy();
  
  std::vector<double> vals;
  kirchnerPalette col;

  for(int i = 0; i < nTrials; ++i){
    canv_p->cd();      

    TGraph* graph_p = new TGraph();

    double count = n0;
    for(int j = 0; j < nTimeSteps; ++j){
      graph_p->SetPoint(j, j, count/(double)n0);

      double draw1 = randGen_p->Uniform(0,1);
      if(draw1 < .5) count *= lLow;
      else count *= lHi;
    }

    graph_p->SetPoint(nTimeSteps, nTimeSteps, count/(double)n0);
    if(i < 6){
      graph_p->SetMarkerColor(col.getColor(i));
      graph_p->SetLineColor(col.getColor(i));
      graph_p->SetMarkerStyle(20);
      graph_p->SetMarkerSize(.5);
      canv_p->cd();
      graph_p->Draw("P");
    }

    vals.push_back(count);
    simHist_h->Fill(count);
    simHist_n0_h->Fill(count/(double)n0);
  }

  std::sort(std::begin(vals), std::end(vals));

  double median = (vals.at(vals.size()/2 - 1) + vals.at(vals.size()/2))/2.;
  double mean = simHist_h->GetMean();
  
  std::cout << "n0: " << n0 << std::endl;
  std::cout << "Median: " << median << ", " << median/(double)n0 << std::endl;
  std::cout << "Mean: " << mean << ", " << mean/(double)n0 << std::endl;

  outFile_p->cd();
  simHist_h->Write("", TObject::kOverwrite);
  delete simHist_h;

  canv_p->cd();
  gStyle->SetOptStat(0);
  canv_p->SaveAs("pdfDir/timeSeriesBetHedge.pdf");
  delete canv_p;

  canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  canv_p->cd();

  simHist_n0_h->DrawCopy("HIST");
  gStyle->SetOptStat(0);

  canv_p->SaveAs("pdfDir/betHedge300Hist.pdf");
  delete canv_p;

  simHist_n0_h->Write("", TObject::kOverwrite);
  delete simHist_n0_h;


  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 2){
    std::cout << "Usage: ./betHedging.exe <n0>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 1) retVal += betHedging();
  else if(argc == 2) retVal += betHedging(std::stoi(argv[1]));
  return retVal;
}
