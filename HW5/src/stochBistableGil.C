#include <iostream>

#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TStyle.h"

#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int stochBistableGil()
{
  TRandom3* randGen_p = new TRandom3(0);

  kirchnerPalette col;

  //  double gammaVal = 1.;
  double k1k2_1 = .0001;
  double v0 = 12.5;
  double v1 = 200.;

  double maxTime = 1000000;
  
  std::vector<double> timeVals;
  std::vector<int> counts;
  
  timeVals.push_back(0);
  counts.push_back(0);

  while(timeVals.at(timeVals.size()-1) < maxTime){
    int count = counts.at(counts.size()-1);
    double rateProd = (v0 + v1*k1k2_1*count*count)/(1 + k1k2_1*count*count);
    double rateDeg = count;

    double timeProd = 1./rateProd*TMath::Log(1./randGen_p->Uniform(0.0000001, 1));
    double timeDeg = 100000000;
    if(count != 0) timeDeg = 1./rateDeg*TMath::Log(1./randGen_p->Uniform(0.0000001, 1));
    
    if(timeProd < timeDeg){
      timeVals.push_back(timeVals.at(timeVals.size()-1) + timeProd);
      counts.push_back(counts.at(counts.size()-1) + 1);
    }
    else{
      timeVals.push_back(timeVals.at(timeVals.size()-1) + timeDeg);
      counts.push_back(counts.at(counts.size()-1) - 1);
    }
  }

  TCanvas* canvBistableGil_c = new TCanvas("canvBistableGil_c", "canvBistableGil_c", 800, 500);
  prettyCanv(canvBistableGil_c);

  TGraph* graphBistableGil_p = new TGraph();
  TH1F* dummy_h = new TH1F("dummy_h", ";Time (s);Counts (Bistable system, Problem 4)", 10, 0, 5100);
  prettyTH1(dummy_h, 10., 10, 10);

  graphBistableGil_p->SetMarkerColor(col.getColor(0));
  graphBistableGil_p->SetMarkerSize(.2);
  graphBistableGil_p->SetMarkerStyle(20);

  double lessThan50Sum = 0;
  double greaterThan50Sum = 0;
  double max = 0;

  const int points = 1000;
  for(unsigned int i = 0; i < timeVals.size(); ++i){
    if(i%points == 0) graphBistableGil_p->SetPoint(i/points, timeVals.at(i), counts.at(i));

    if(counts.at(i) > max) max = counts.at(i);

    double deltaTime = timeVals.at(i);
    if(i>0) deltaTime -= timeVals.at(i-1);

    if(counts.at(i) < 50) lessThan50Sum += deltaTime;
    else greaterThan50Sum += deltaTime;
  }

  std::cout << lessThan50Sum << ", " << greaterThan50Sum << std::endl;

  gStyle->SetOptStat(0);
  canvBistableGil_c->cd();
  dummy_h->SetMaximum(max);
  dummy_h->DrawCopy();
  graphBistableGil_p->Draw("P");
  canvBistableGil_c->SaveAs("pdfDir/canvBistableGil.pdf");

  delete graphBistableGil_p;
  delete canvBistableGil_c;

  delete randGen_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += stochBistableGil();
  return retVal;
}
