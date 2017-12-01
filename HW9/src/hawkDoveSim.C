#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TLatex.h"

#include "include/getLinBins.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"

double getMeanPorF(const int arrSize, Double_t arrVal[], Double_t arrWeight[])
{
  double arrMean = 0.0;
  double arrDenom = 0.0;
  for(Int_t arrI = 0; arrI < arrSize; ++arrI){
    arrMean += arrVal[arrI]*arrWeight[arrI];
    arrDenom += arrWeight[arrI];
  }
  
  return (arrMean/arrDenom);
}

int hawkDoveSim(const double muParam = 0.001)
{
  TDatime* date = new TDatime();
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  label_p->SetNDC();

  TFile* outFile_p = new TFile("output/hawkDoveSim.root", "RECREATE");

  const int nProbBins = 100;
  Double_t probBins[nProbBins+1];
  getLinBins(0., 1., nProbBins, probBins);

  TH1F* probability_p = new TH1F("prob_h", ";p_{Hawk};Prob(p_{Hawk})", nProbBins, probBins);
  TH1F* dpdt_p = new TH1F("dpdt_h", ";p_{Hawk};dProb(p_{Hawk})/dt", nProbBins, probBins);

  const float bParam = 1.;
  const float cParam = 4.;

  const int nPop = 10001;
  const float pLow = 0.0;
  const float pHi = 1.0;
  Double_t pVals[nPop];
  Double_t dpdtVals[nPop];
  Double_t pop[nPop];
  getLinBins(pLow, pHi, nPop-1, pVals);
  for(int i = 0; i < nPop; ++i){pop[i] = 1;}

  const int nTimeSteps = 1000000;
  const double timeStep = 0.05;

  for(Int_t timeI = 0; timeI < nTimeSteps; ++timeI){
    if(timeI%1000 == 0) std::cout << "nTimeStep: " << timeI << "/" << nTimeSteps << std::endl;

    Double_t pMean = getMeanPorF(nPop, pVals, pop);

    Double_t pay[nPop];
    for(Int_t pI = 0; pI < nPop; ++pI){
      pay[pI] = (bParam - cParam)*pVals[pI]*pMean/2.;
      pay[pI] += bParam*pVals[pI]*(1-pMean);
      pay[pI] += bParam*(1-pVals[pI])*(1-pMean)/2.;
    }

    Double_t payMean = getMeanPorF(nPop, pay, pop);

    for(Int_t pI = 0; pI < nPop; ++pI){
      Double_t deltaP = (pay[pI] - payMean)*pop[pI] - muParam*pop[pI] + muParam;
      dpdtVals[pI] = deltaP;
      pop[pI] += deltaP*timeStep;
    }
  }

  Double_t pMean = getMeanPorF(nPop, pVals, pop);

  Double_t totPop = 0;
  for(Int_t i = 0; i < nPop; ++i){
    totPop += pop[i];
    probability_p->Fill(pVals[i], pop[i]);
    dpdt_p->Fill(pVals[i], dpdtVals[i]);
  }

  outFile_p->cd();
  probability_p->Write("", TObject::kOverwrite);
  dpdt_p->Write("", TObject::kOverwrite);
  
  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  probability_p->GetXaxis()->CenterTitle();
  probability_p->GetYaxis()->CenterTitle();
  probability_p->Scale(1./probability_p->Integral());
  probability_p->SetMaximum(.1);
  probability_p->SetMinimum(0.);
  probability_p->DrawCopy("HIST");
  label_p->DrawLatex(.25, .95, ("#mu=" + prettyString(muParam, 4, false)).c_str());
  label_p->DrawLatex(.45, .95, ("Time=" + std::to_string(int(nTimeSteps*timeStep)) + "(s)").c_str());
  label_p->DrawLatex(.75, .95, ("#Deltat=" + prettyString(timeStep, 2, false) + "(s)").c_str());
  label_p->DrawLatex(.75, .55, ("Mean=" + prettyString(probability_p->GetMean(),4,false)).c_str());
  label_p->DrawLatex(.75, .48, ("b=" + std::to_string(int(bParam))).c_str());
  label_p->DrawLatex(.75, .41, ("c=" + std::to_string(int(cParam))).c_str());

  gStyle->SetOptStat(0);


  canv_p->SaveAs(("pdfDir/hawkDoveProb_Mu" + prettyString(muParam, 4, true) + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());
  delete canv_p;

  canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  dpdt_p->GetXaxis()->CenterTitle();
  dpdt_p->GetYaxis()->CenterTitle();
  dpdt_p->Scale(1./probability_p->Integral());
  dpdt_p->DrawCopy("HIST");
  label_p->DrawLatex(.25, .95, ("#mu=" + prettyString(muParam, 4, false)).c_str());
  label_p->DrawLatex(.45, .95, ("Time=" + std::to_string(int(nTimeSteps*timeStep)) + "(s)").c_str());
  label_p->DrawLatex(.75, .95, ("#Deltat=" + prettyString(timeStep, 2, false) + "(s)").c_str());
  gStyle->SetOptStat(0);
  canv_p->SaveAs(("pdfDir/hawkDoveProb_DPDT_Mu" + prettyString(muParam, 4, true) + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());
  delete canv_p;


  delete probability_p;
  delete dpdt_p;
  outFile_p->Close();
  delete outFile_p;

  std::cout << totPop << ", " << pMean << std::endl;

  delete date;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 2){
    std::cout << "Usage: ./hawkDoveSime.exe <muParam-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 1) retVal += hawkDoveSim();
  else if(argc == 2) retVal += hawkDoveSim(std::stof(argv[1]));
  return retVal;
}
