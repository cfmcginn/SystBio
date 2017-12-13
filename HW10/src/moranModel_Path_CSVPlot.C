#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"

#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"

int moranModel_Path_CSVPlot(const std::string inFileName)
{
  std::vector<double> muteVals;
  std::vector<double> probValsA;
  std::vector<double> probValsAErr;
  std::vector<double> probValsB;
  std::vector<double> probValsBErr;

  std::vector<std::string> totVals;
  std::ifstream file(inFileName.c_str());
  std::string tempStr;
  unsigned int counter = 0;
  
  while(std::getline(file, tempStr)){
    if(tempStr.size() == 0) continue;
    if(counter==0){++counter; continue;}

    totVals.push_back(tempStr);
    double mute = std::stof(tempStr.substr(0, tempStr.find(",")));
    muteVals.push_back(mute);
    tempStr.replace(0,tempStr.find(",")+1, "");

    mute = std::stof(tempStr.substr(0, tempStr.find(",")));
    probValsA.push_back(mute);
    tempStr.replace(0,tempStr.find(",")+1, "");

    mute = std::stof(tempStr.substr(0, tempStr.find(",")));
    probValsAErr.push_back(mute);
    tempStr.replace(0,tempStr.find(",")+1, "");

    mute = std::stof(tempStr.substr(0, tempStr.find(",")));
    probValsB.push_back(mute);
    tempStr.replace(0,tempStr.find(",")+1, "");

    mute = std::stof(tempStr.substr(0, tempStr.find(",")));
    probValsBErr.push_back(mute);
    tempStr.replace(0,tempStr.find(",")+1, "");
    

    ++counter;
  }
  file.close();
  
  const Int_t nMuteBins = muteVals.size();
  Double_t muteBins[nMuteBins+1];
  for(Int_t muI = 1; muI < nMuteBins; ++muI){
    muteBins[muI] = (muteVals[muI] + muteVals[muI-1])/2.;
  }
  muteBins[0] = muteBins[1]/2.;
  muteBins[nMuteBins] = muteBins[nMuteBins-1]*2.;

  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  kirchnerPalette col;

  TH1F* hist_p = new TH1F("hist_h", ";Mutation Rate;Path Fraction", nMuteBins, muteBins);
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();
  hist_p->SetMaximum(1.1);
  hist_p->SetMinimum(0.0);
  hist_p->DrawCopy();

  gStyle->SetOptStat(0);
  gPad->SetLogx();

  TH1F* histA_p = new TH1F("histA_h", ";Mutation Rate;Path Fraction", nMuteBins, muteBins);
  histA_p->SetMarkerSize(0.0);
  histA_p->SetMarkerStyle(20);
  histA_p->SetFillColor(col.getColor(0));
  histA_p->SetLineColor(1);
  histA_p->SetLineWidth(2);

  TH1F* histB_p = new TH1F("histB_h", ";Mutation Rate;Path Fraction", nMuteBins, muteBins);
  histB_p->SetMarkerSize(0.0);
  histB_p->SetMarkerStyle(20);
  histB_p->SetFillColor(col.getColor(2));
  histB_p->SetLineColor(1);
  histB_p->SetLineWidth(2);

  for(Int_t muI = 0; muI < nMuteBins; ++muI){
    Double_t totVal = probValsA.at(muI);
    histA_p->SetBinContent(muI+1, totVal);
    histA_p->SetBinError(muI+1, probValsAErr.at(muI));

    totVal += probValsB.at(muI);
    histB_p->SetBinContent(muI+1, totVal);
    histB_p->SetBinError(muI+1, probValsBErr.at(muI));
  }
  
  histB_p->DrawCopy("HIST SAME");
  histA_p->DrawCopy("HIST E1 SAME");

  gPad->RedrawAxis();
  gPad->Modified();

  canv_p->SaveAs("pdfDir/tempCSVPlot.pdf");
  delete canv_p;

  delete hist_p;
  delete histA_p;
  delete histB_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./moranModel_Path_CSVPlot.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += moranModel_Path_CSVPlot(argv[1]);
  return retVal;
}
