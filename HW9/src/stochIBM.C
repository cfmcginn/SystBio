#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TLegend.h"

#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"

int stochIBM(const double ntot = 3000)
{
  TRandom3* randGen_p = new TRandom3(0);

  const double d1Param = 0.1;
  const double d2Param = 0.05;
  const double bParam = 0.1;
  const double p1Param = 0.25;
  const double p2Param = 0.05;
  /*
  const double d1Param = 0.1;
  const double d2Param = 0.05;
  const double bParam = 0.1;
  const double p1Param = 0.25;
  const double p2Param = 0.05;
  */
  const double rParam = 2*bParam - d2Param;
  const double kParam = 1. - d2Param/(2.*bParam);
  const double pParam = (p1Param + p2Param + bParam);

  double stablef1 = rParam/(2.*pParam)*(1. - d1Param/(2.*p1Param*kParam));
  double stablef2 = d1Param/(2.*p1Param);

  double tauParam = -2.*bParam*stablef2;
  double deltaParam = 4*pParam*p1Param*stablef1*stablef2;

  const double f10 = 0.5;
  const double f20 = 0.5;

  const double n10 = ntot/2.;
  const double n20 = ntot/2.;
  
  const int nTimeSteps = 1000000;
  double deltaTimeStep = .001;

  double currf1 = f10;
  double currf2 = f20;
  double currTime = 0.;

  double currN1 = n10;
  double currN2 = n20;
  double currTimeGill = 0.;

  TGraph* graphf1_p = new TGraph();
  TGraph* graphf2_p = new TGraph();

  graphf1_p->SetPoint(0, currTime, currf1);
  graphf2_p->SetPoint(0, currTime, currf2);

  TGraph* graphN1_p = new TGraph();
  TGraph* graphN2_p = new TGraph();

  std::vector<double> absDeltaN1;
  std::vector<double> absDeltaN2;

  graphN1_p->SetPoint(0, currTime, currN1/ntot);
  graphN2_p->SetPoint(0, currTime, currN2/ntot);

  double nextThreshDet = .01;
  int pointDet = 1;
  for(Int_t i = 0; i < nTimeSteps; ++i){
    if(i%1000 == 0) std::cout << "Step " << i << "/" << nTimeSteps << std::endl;

    double df1dt = 2*p1Param*currf1*currf2 - d1Param*currf1;
    double df2dt = rParam*currf2*(1. - currf2/kParam) - 2*pParam*currf1*currf2;

    currf1 += df1dt*deltaTimeStep;
    currf2 += df2dt*deltaTimeStep;
    currTime += deltaTimeStep;

    if(currTime > nextThreshDet){
      graphf1_p->SetPoint(pointDet, currTime, currf1);
      graphf2_p->SetPoint(pointDet, currTime, currf2);
      ++pointDet;
      nextThreshDet += .01;
    }

  }

  int gilPoint = 1;
  double nextThresh = .5;
  while(currTimeGill < currTime){
    double time1 = randGen_p->Exp(1./(d1Param*currN1));
    double time2 = randGen_p->Exp(1./(2.*bParam*currN2*(ntot - currN1 - currN2)/ntot));
    double time3 = randGen_p->Exp(1./(2.*p2Param*currN1*currN2/ntot + d2Param*currN2));
    double time4 = randGen_p->Exp(1./(2*p1Param*currN1*currN2/ntot));

    if(time1 < time2 && time1 < time3 && time1 < time4 && currN1 != 1){
      --currN1;
      currTimeGill += time1;
    }
    else if(time2 < time3 && time2 < time4){
      ++currN2;
      currTimeGill += time2;
    }
    else if(time3 < time4 && currN2 != 1){
      --currN2;
      currTimeGill += time3;
    }
    else if(currN2 != 1){
      ++currN1;
      --currN2;
      currTimeGill += time4;
    }
    else{
      ++currN2;
      currTimeGill += time2;
    }
    
    if(currTimeGill > nextThresh){
      graphN1_p->SetPoint(gilPoint, currTimeGill, currN1/ntot);
      graphN2_p->SetPoint(gilPoint, currTimeGill, currN2/ntot);
      ++gilPoint;
      nextThresh += .5;
    }

    if(currTimeGill > 300){
      absDeltaN1.push_back(TMath::Abs(currN1/ntot - stablef1));
      absDeltaN2.push_back(TMath::Abs(currN2/ntot - stablef2));
    }
  }


  std::cout << absDeltaN1.size() << ", " << absDeltaN2.size() << std::endl;

  std::vector<double> max10N1;

  for(int i = 0; i < 10; ++i){
    double max = -1;
    int maxPos = -1;
    for(unsigned int j = 0; j < absDeltaN1.size(); ++j){
      if(absDeltaN1.at(j) > max){
	maxPos = j;
	max = absDeltaN1.at(j);
      }
    }

    max10N1.push_back(max);

    for(int l = 0; l < 500; ++l){
      absDeltaN1.erase(absDeltaN1.begin()+maxPos+l+1);
    }

    for(int l = 0; l < 500; ++l){
      absDeltaN1.erase(absDeltaN1.begin()+maxPos);
      --maxPos;
    }
  }

  for(unsigned int i = 0; i < max10N1.size(); ++i){
    std::cout << max10N1.at(i) << std::endl;
  }

  std::cout << stablef1 << ", " << stablef2 << std::endl;

  kirchnerPalette col;
  graphf1_p->SetMarkerColor(col.getColor(0));
  graphf2_p->SetMarkerColor(col.getColor(2));
  graphf1_p->SetMarkerStyle(20);
  graphf2_p->SetMarkerStyle(20);
  graphf1_p->SetMarkerSize(.8);
  graphf2_p->SetMarkerSize(.8);
  graphf1_p->SetLineColor(col.getColor(0));
  graphf2_p->SetLineColor(col.getColor(2));

  graphN1_p->SetMarkerColor(col.getColor(0));
  graphN2_p->SetMarkerColor(col.getColor(2));
  graphN1_p->SetMarkerStyle(20);
  graphN2_p->SetMarkerStyle(20);
  graphN1_p->SetMarkerSize(.2);
  graphN2_p->SetMarkerSize(.2);
  graphN1_p->SetLineColor(col.getColor(0));
  graphN2_p->SetLineColor(col.getColor(2));

  TH1F* dummyHist_p = new TH1F("dummyHist_h", ";Time;f1,f2", 10, 0, currTime);
  dummyHist_p->SetMaximum(0.7);
  dummyHist_p->SetMinimum(0.0);
  dummyHist_p->GetXaxis()->CenterTitle();
  dummyHist_p->GetYaxis()->CenterTitle();

  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  dummyHist_p->DrawCopy();
  graphf1_p->Draw("SAME");
  graphf2_p->Draw("SAME");

  graphN1_p->Draw("SAME P");
  graphN2_p->Draw("SAME P");

  gStyle->SetOptStat(0);

  TDatime* date = new TDatime();
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  label_p->SetNDC();
  label_p->DrawLatex(.2, .95, ("N=" + std::to_string(int(ntot))).c_str());

  label_p->DrawLatex(.75, .9, ("d1=" + prettyString(d1Param, 2, false)).c_str());
  label_p->DrawLatex(.75, .85, ("d2=" + prettyString(d2Param, 2, false)).c_str());
  label_p->DrawLatex(.75, .8, ("b=" + prettyString(bParam, 2, false)).c_str());
  label_p->DrawLatex(.75, .75, ("p1=" + prettyString(p1Param, 2, false)).c_str());
  label_p->DrawLatex(.75, .7, ("p2=" + prettyString(p2Param, 2, false)).c_str());

  TLegend* leg_p = new TLegend(.5, .9, .5, .9);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(18);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);

  TH1F* dumm1 = new TH1F("dumm1", "", 10, 0, 1);
  dumm1->SetMarkerColor(col.getColor(0));
  dumm1->SetMarkerStyle(20);
  dumm1->SetMarkerSize(.4);
  dumm1->SetLineColor(col.getColor(0));

  TH1F* dumm2 = new TH1F("dumm2", "", 10, 0, 1);
  dumm2->SetMarkerColor(col.getColor(2));
  dumm2->SetMarkerStyle(20);
  dumm2->SetMarkerSize(.4);
  dumm2->SetLineColor(col.getColor(2));

  leg_p->AddEntry(dumm1, "Predator", "L P");
  leg_p->AddEntry(dumm2, "Prey", "L P");

  leg_p->Draw("SAME");

  canv_p->SaveAs(("pdfDir/stochIBM_NTot" + std::to_string(int(ntot)) + "_p1" + prettyString(p1Param,2,true) + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());

  delete canv_p;

  delete graphf1_p;
  delete graphf2_p;

  delete graphN1_p;
  delete graphN2_p;

  delete date;

  delete randGen_p;

  std::cout << "Tau^2, 4delta^2: " << tauParam*tauParam << ", " << 4*deltaParam << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 1 && argc !=2){
    std::cout << "Usage: ./stochIBM.exe <ntot-optionale>" << std::endl;
    return 0;
  }
  
  int retVal = 0;
  if(argc == 1) retVal += stochIBM();
  else if(argc == 2) retVal += stochIBM(std::stof(argv[1]));
  return 0;
}
