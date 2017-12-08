#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

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

  std::vector< std::vector<double> > snapShots;
  snapShots.reserve(nTimeSteps+1);
  for(Int_t i = 0; i < nTimeSteps+1; ++i){
    std::vector<double> temp;
    snapShots.push_back(temp);
  }


  TFile* outFile_p = new TFile("output/betHedging.root", "RECREATE");
  TH1F* simHist_h = new TH1F("simHist_h", ";;", 100, 0, 2.*n0);
  TH1F* simHist_n0_h = new TH1F("simHist_n0_h", ";n(t=300)/n_{0};Counts", 100, 0, 2.);
  simHist_n0_h->GetXaxis()->CenterTitle();
  simHist_n0_h->GetYaxis()->CenterTitle();

  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
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

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TLegend* leg_p = new TLegend(0.6, 0.9, 0.6, 0.9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);

  const Int_t nGraph = 6;
  TGraph* graph_p[nGraph];
  TH1F* dummyH_p[nGraph];

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(int i = 0; i < nTrials; ++i){
    canv_p->cd();      

    if(i < nGraph) graph_p[i] = new TGraph();

    double count = n0;
    for(int j = 0; j < nTimeSteps; ++j){
      snapShots.at(j).push_back(count/(double)n0);

      if(i < nGraph) graph_p[i]->SetPoint(j, j, count/(double)n0);

      double draw1 = randGen_p->Uniform(0,1);
      if(draw1 < .5) count *= lLow;
      else count *= lHi;
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    snapShots.at(nTimeSteps).push_back(count/(double)n0);

    if(i < nGraph){
      graph_p[i]->SetPoint(nTimeSteps, nTimeSteps, count/(double)n0);
      graph_p[i]->SetMarkerColor(col.getColor(i));
      graph_p[i]->SetLineColor(col.getColor(i));
      graph_p[i]->SetMarkerStyle(20);
      graph_p[i]->SetMarkerSize(.5);

      dummyH_p[i] = new TH1F(("dummy_" + std::to_string(i)).c_str(), "", 10, 0, 10);
      dummyH_p[i]->SetMarkerColor(col.getColor(i));
      dummyH_p[i]->SetLineColor(col.getColor(i));
      dummyH_p[i]->SetMarkerStyle(20);
      dummyH_p[i]->SetMarkerSize(.5);

      leg_p->AddEntry(dummyH_p[i], ("Example " + std::to_string(i)).c_str(), "P L");

      canv_p->cd();
      graph_p[i]->Draw("P");
    }

    vals.push_back(count);
    simHist_h->Fill(count);
    simHist_n0_h->Fill(count/(double)n0);
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::sort(std::begin(vals), std::end(vals));

  double median = (vals.at(vals.size()/2 - 1) + vals.at(vals.size()/2))/2.;
  double mean = simHist_h->GetMean();

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::cout << "n0: " << n0 << std::endl;
  std::cout << "Median: " << median << ", " << median/(double)n0 << std::endl;
  std::cout << "Mean: " << mean << ", " << mean/(double)n0 << std::endl;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();
  simHist_h->Write("", TObject::kOverwrite);
  delete simHist_h;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  canv_p->cd();
  gStyle->SetOptStat(0);
  leg_p->Draw("SAME");
  canv_p->SaveAs("pdfDir/timeSeriesBetHedge.pdf");
  delete canv_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  canv_p->cd();

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  simHist_n0_h->DrawCopy("HIST");
  gStyle->SetOptStat(0);

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  canv_p->SaveAs("pdfDir/betHedge300Hist.pdf");
  delete canv_p;

  canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  TH1F* dummy2_p = new TH1F("dummy2_h", ";Time;Mean or Median n(t)", 10, 0, 300);
  dummy2_p->SetMaximum(2.);
  dummy2_p->SetMinimum(0.0);
  dummy2_p->GetXaxis()->CenterTitle();
  dummy2_p->GetYaxis()->CenterTitle();
  dummy2_p->Draw();

  TGraph* graphMean_p = new TGraph();
  TGraph* graphMedian_p = new TGraph();

  TH1F* histMean_p = new TH1F("histMean", "", 10, 0, 10);
  TH1F* histMedian_p = new TH1F("histMedian", "", 10, 0, 10);

  histMean_p->SetMarkerStyle(20);
  histMean_p->SetMarkerSize(.6);
  histMean_p->SetMarkerColor(col.getColor(0));
  histMean_p->SetLineColor(col.getColor(0));

  histMedian_p->SetMarkerStyle(20);
  histMedian_p->SetMarkerSize(.6);
  histMedian_p->SetMarkerColor(col.getColor(2));
  histMedian_p->SetLineColor(col.getColor(2));

  graphMean_p->SetMarkerStyle(20);
  graphMean_p->SetMarkerSize(.6);
  graphMean_p->SetMarkerColor(col.getColor(0));
  graphMean_p->SetLineColor(col.getColor(0));

  graphMedian_p->SetMarkerStyle(20);
  graphMedian_p->SetMarkerSize(.6);
  graphMedian_p->SetMarkerColor(col.getColor(2));
  graphMedian_p->SetLineColor(col.getColor(2));

  double max = 0;
  double min = 1;

  for(unsigned int i = 0; i < nTimeSteps+1; ++i){
    std::vector<double> temp = snapShots.at(i);
    std::sort(std::begin(temp), std::end(temp));

    double mean = 0.0;
    for(unsigned int j = 0; j < temp.size(); ++j){
      mean += temp.at(j);
    }
    mean /= (double)temp.size();

    median = temp.at(temp.size()/2-1) + temp.at(temp.size()/2);
    median /= 2.;

    if(mean > max) max = mean;
    std::cout << mean << ", ";
    //    std::cout << temp.at(temp.size()-1) << ", ";

    graphMean_p->SetPoint(i, i, mean);
    graphMedian_p->SetPoint(i, i, median);

    if(median < min) min = median;
  }
  std::cout << std::endl;
  dummy2_p->SetMaximum(max);
  dummy2_p->SetMinimum(min);
  gPad->Modified();
  gPad->SetLogy();

  TLegend* leg2_p = new TLegend(0.6, 0.5, 0.9, 0.7);
  leg2_p->SetBorderSize(0);
  leg2_p->SetFillStyle(0);
  leg2_p->SetFillColor(0);
  leg2_p->AddEntry(histMean_p, "Mean", "P L");
  leg2_p->AddEntry(histMedian_p, "Median", "P L");


  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  label_p->DrawLatex(.2, .95, "Simulate 100000 times");
  leg2_p->Draw("SAME");
  graphMedian_p->Draw("P");
  graphMean_p->Draw("P");

  canv_p->SaveAs("pdfDir/betHedgeMeanMedian.pdf");

  delete canv_p;

  delete graphMean_p;
  delete graphMedian_p;
  delete dummy2_p;

  simHist_n0_h->Write("", TObject::kOverwrite);
  delete simHist_n0_h;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

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
