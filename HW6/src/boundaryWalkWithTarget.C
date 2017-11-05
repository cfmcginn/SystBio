#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>

#include "TFile.h"
#include "TH2I.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"

#include "include/plotUtilities.h"

int boundaryWalkWithTarget(const int nDPos, const int nSize)
{
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  label_p->SetNDC();

  const int nSimul = 10000;
  TRandom3* randGen_p = new TRandom3(0);

  TFile* outFile_p = new TFile("output/boundaryWalkWithTarget.root", "RECREATE");
  TH1F* stepsToTarget_p = new TH1F("stepsToTarget_h", ";Steps (N);Counts", 100, -0.5, 99999.5);
  stepsToTarget_p->Sumw2();

  std::vector<double> meanMedian;

  //  const int printInterval = 1000;
  
  double nCounts = 0;

  for(int simIter = 0; simIter < nSimul; ++simIter){
    //    if(simIter%printInterval == 0) std::cout << "Sim " << simIter << "/" << nSimul << std::endl;

    int posX = 0;
    int posY = 0;

    int targetX = 0;
    double targetXPull = randGen_p->Uniform(-nDPos - .5, nDPos + .5);
    for(Int_t i = 0; i < 41; ++i){
      if(targetXPull >= -nDPos - .5 + i && targetXPull <= -nDPos - .5 + (i+1)){
	targetX = -nDPos + i;
	break;
      }
    }

    bool targetYIsPos = true;
    if(randGen_p->Uniform(0,1) < .5) targetYIsPos = false;
    int targetY = nDPos - TMath::Abs(targetX);
    if(!targetYIsPos) targetY *= -1;

    //    if(simIter%printInterval == 0) std::cout << "targetX,Y: " << targetX << ", " << targetY << std::endl;

    int stepCount = 0;

    while(posX != targetX || posY != targetY){
      bool isX = true;
      if(randGen_p->Uniform(0,1) > .5) isX = false;
      bool isPos = true;
      if(randGen_p->Uniform(0,1) > .5) isPos = false;

      if(isX && isPos) ++posX;
      else if(isX && !isPos) --posX;
      else if(!isX && isPos) ++posY;
      else if(!isX && !isPos) --posY;
      
      if(posX == nSize/2 + 1) posX = -(nSize/2);
      else if(posX == -nSize/2 - 1) posX = (nSize/2);
      else if(posY == nSize/2 + 1) posY = -(nSize/2);
      else if(posY == -nSize/2 - 1) posY = (nSize/2);

      ++stepCount;
    }

    stepsToTarget_p->Fill(stepCount);
    meanMedian.push_back(stepCount);
    ++nCounts;
  }

  std::sort(meanMedian.begin(), meanMedian.end());
  double mean = std::accumulate(meanMedian.begin(), meanMedian.end(), 0.0)/(double)meanMedian.size();
  double median = (meanMedian.at(meanMedian.size()/2-1) + meanMedian.at(meanMedian.size()/2))/2.;
  double stddev = std::inner_product(meanMedian.begin(), meanMedian.end(), meanMedian.begin(), 0.0);
  stddev = TMath::Sqrt(stddev/(double)meanMedian.size() - mean*mean);

  std::cout << "Mean, median, stddev: " << mean << ", " << median << ", " << stddev << std::endl;

  outFile_p->cd();

  stepsToTarget_p->Write("", TObject::kOverwrite);

  TCanvas* canv_p = new TCanvas("stepsToTarget_c", "stepsToTarget_c", 500, 500);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);

  canv_p->SetBottomMargin(1.3*canv_p->GetBottomMargin());
  canv_p->SetLeftMargin(canv_p->GetBottomMargin());

  stepsToTarget_p->GetXaxis()->CenterTitle();
  stepsToTarget_p->GetYaxis()->CenterTitle();
  stepsToTarget_p->Scale(1./nCounts);

  stepsToTarget_p->GetXaxis()->SetNdivisions(505);
  stepsToTarget_p->DrawCopy("HIST E1");
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  label_p->DrawLatex(.65, .94, "2D Diffusion");
  label_p->DrawLatex(.65, .88, ("D="+ std::to_string(nDPos) + "; N=" + std::to_string(nSize)).c_str());
  label_p->DrawLatex(.65, .82, ("Mean=" + prettyString(stepsToTarget_p->GetMean(), 1, false)).c_str());
  label_p->DrawLatex(.65, .76, ("StdDev=" + prettyString(stepsToTarget_p->GetStdDev(), 1, false)).c_str());
  label_p->DrawLatex(.65, .7, ("Median=" + prettyString(median, 1, false)).c_str());

  canv_p->SaveAs(("pdfDir/stepsToTarget_nSize" + std::to_string(nSize) + "_D" + std::to_string(nDPos) + ".pdf").c_str());

  delete canv_p;
  delete stepsToTarget_p;

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./boundaryWalkWithTarget.exe <inputD> <inputL>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += boundaryWalkWithTarget(std::stoi(argv[1]), std::stoi(argv[2]));
  return retVal;
}
