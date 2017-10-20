#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1F.h"
#include "TStyle.h"

#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int stochPotentialProb()
{
  kirchnerPalette col;

  double gammaVal = 1.;
  double k1k2 = .0001;
  double v0 = 12.5;
  double v1 = 200.;

  double intervalEmu = .00001;
  double maxX = 200;
  
  std::vector<double> xVals;
  std::vector<double> intVals;
  std::vector<double> yValsE;
  std::vector<double> yValsDenom;
  std::vector<double> yValsNum;
  std::vector<double> yValsProb;
  
  xVals.push_back(0);
  intVals.push_back(v0*intervalEmu);
  yValsE.push_back(TMath::Exp(-v0*intervalEmu));
  yValsDenom.push_back(v0*intervalEmu);
  yValsNum.push_back(v0*intervalEmu);
  yValsProb.push_back(yValsE.at(yValsE.size()-1)/yValsDenom.at(yValsDenom.size()-1));

  while(xVals.at(xVals.size()-1) < maxX){
    double xVal = xVals.at(xVals.size()-1) + intervalEmu;
    double prodVal = (v0 + v1*k1k2*xVal*xVal)/(1 + k1k2*xVal*xVal);
    double degVal = gammaVal*xVal;

    double valInt = -2*intervalEmu*(prodVal - degVal)/(prodVal + degVal) + intVals.at(intVals.size()-1);
    double valE = TMath::Exp(-valInt);
    double valDenom = prodVal + degVal;
    double valNum = prodVal - degVal;

    xVals.push_back(xVal);
    intVals.push_back(valInt);
    yValsE.push_back(valE);
    yValsDenom.push_back(valDenom);
    yValsNum.push_back(valNum);
    yValsProb.push_back(yValsE.at(yValsE.size()-1)/yValsDenom.at(yValsDenom.size()-1));
  }

  TCanvas* canvPotentialE_c = new TCanvas("canvPotentialE_c", "canvPotentialE_c", 500, 500);
  prettyCanv(canvPotentialE_c);
  TH1F* dummyE_p = new TH1F("dummyE_h", ";Counts;Potential #phi", 200, 0, 200);

  TCanvas* canvPotentialDenom_c = new TCanvas("canvPotentialDenom_c", "canvPotentialDenom_c", 500, 500);
  prettyCanv(canvPotentialDenom_c);
  TH1F* dummyDenom_p = new TH1F("dummyDenom_h", ";Counts;f(n)+g(n)", 200, 0, 200);

  TCanvas* canvPotentialNum_c = new TCanvas("canvPotentialNum_c", "canvPotentialNum_c", 500, 500);
  prettyCanv(canvPotentialNum_c);
  TH1F* dummyNum_p = new TH1F("dummyNum_h", ";Counts;f(n)-g(n)", 200, 0, 200);

  TCanvas* canvPotentialProb_c = new TCanvas("canvPotentialProb_c", "canvPotentialProb_c", 500, 500);
  prettyCanv(canvPotentialProb_c);
  TH1F* dummyProb_p = new TH1F("dummyProb_h", ";Counts;Probability (Unnormalized)", 200, 0, 200);

  TCanvas* canvPotentialProb2_c = new TCanvas("canvPotentialProb2_c", "canvPotentialProb2_c", 500, 500);
  prettyCanv(canvPotentialProb2_c);
  TH1F* dummyProb2_p = new TH1F("dummyProb2_h", ";Counts;Probability p(n) = #frac{A}{f+g} e^{-#phi}", 200, 0, 200);

  prettyTH1(dummyE_p, 10, 10, 10);
  prettyTH1(dummyDenom_p, 10, 10, 10);
  prettyTH1(dummyNum_p, 10, 10, 10);
  prettyTH1(dummyProb_p, 10, 10, 10);
  prettyTH1(dummyProb2_p, 10, 10, 10);

  dummyE_p->GetXaxis()->SetTitleFont(43);
  dummyE_p->GetYaxis()->SetTitleFont(43);
  dummyE_p->GetXaxis()->SetTitleSize(26);
  dummyE_p->GetYaxis()->SetTitleSize(26);


  dummyDenom_p->GetXaxis()->SetTitleFont(43);
  dummyDenom_p->GetYaxis()->SetTitleFont(43);
  dummyDenom_p->GetXaxis()->SetTitleSize(26);
  dummyDenom_p->GetYaxis()->SetTitleSize(26);


  dummyNum_p->GetXaxis()->SetTitleFont(43);
  dummyNum_p->GetYaxis()->SetTitleFont(43);
  dummyNum_p->GetXaxis()->SetTitleSize(26);
  dummyNum_p->GetYaxis()->SetTitleSize(26);

  dummyProb_p->GetXaxis()->SetTitleFont(43);
  dummyProb_p->GetYaxis()->SetTitleFont(43);
  dummyProb_p->GetXaxis()->SetTitleSize(26);
  dummyProb_p->GetYaxis()->SetTitleSize(26);

  dummyProb2_p->GetXaxis()->SetTitleFont(43);
  dummyProb2_p->GetYaxis()->SetTitleFont(43);
  dummyProb2_p->GetXaxis()->SetTitleSize(26);
  dummyProb2_p->GetYaxis()->SetTitleSize(18);

  
  TGraph* graphPotentialE_p = new TGraph();
  TGraph* graphPotentialDenom_p = new TGraph();
  TGraph* graphPotentialNum_p = new TGraph();
  TGraph* graphPotentialProb_p = new TGraph();
  TGraph* graphPotentialProb2_p = new TGraph();
  
  graphPotentialE_p->SetMarkerSize(.4);
  graphPotentialE_p->SetMarkerStyle(20);
  graphPotentialE_p->SetMarkerColor(col.getColor(0));

  graphPotentialDenom_p->SetMarkerSize(.4);
  graphPotentialDenom_p->SetMarkerStyle(20);
  graphPotentialDenom_p->SetMarkerColor(col.getColor(0));

  graphPotentialNum_p->SetMarkerSize(.4);
  graphPotentialNum_p->SetMarkerStyle(20);
  graphPotentialNum_p->SetMarkerColor(col.getColor(0));

  graphPotentialProb_p->SetMarkerSize(.4);
  graphPotentialProb_p->SetMarkerStyle(20);
  graphPotentialProb_p->SetMarkerColor(col.getColor(0));

  graphPotentialProb2_p->SetMarkerSize(.4);
  graphPotentialProb2_p->SetMarkerStyle(20);
  graphPotentialProb2_p->SetMarkerColor(col.getColor(0));

  const int points = 100000;

  double maxE = 0;
  double maxDenom = 0;
  double maxNum = 0;
  double maxProb = 0;
  double maxProb2 = 0;

  double minE = 10000;
  double minDenom = 10000;
  double minNum = 10000;
  double minProb = 10000;
  double minProb2 = 10000;

  double probSum = 0;
  double probSumLess50 = 0;
  double probSumGreater50 = 0;
  for(unsigned int i = 0; i < xVals.size(); ++i){
    if(i%points == 0 && i != 0) graphPotentialE_p->SetPoint((i-1)/points, xVals.at(i), yValsE.at(i));
    if(i%points == 0 && i != 0) graphPotentialDenom_p->SetPoint((i-1)/points, xVals.at(i), yValsDenom.at(i));
    if(i%points == 0 && i != 0) graphPotentialNum_p->SetPoint((i-1)/points, xVals.at(i), yValsNum.at(i));
    if(i%points == 0 && i != 0) graphPotentialProb_p->SetPoint((i-1)/points, xVals.at(i), yValsProb.at(i));

    if(maxE < yValsE.at(i)) maxE = yValsE.at(i);
    if(minE > yValsE.at(i)) minE = yValsE.at(i);

    if(maxDenom < yValsDenom.at(i)) maxDenom = yValsDenom.at(i);
    if(minDenom > yValsDenom.at(i)) minDenom = yValsDenom.at(i);

    if(maxNum < yValsNum.at(i)) maxNum = yValsNum.at(i);
    if(minNum > yValsNum.at(i)) minNum = yValsNum.at(i);

    if(maxProb < yValsProb.at(i)) maxProb = yValsProb.at(i);
    if(minProb > yValsProb.at(i)) minProb = yValsProb.at(i);

    probSum += yValsProb.at(i);
    if(xVals.at(i) < 50) probSumLess50 += yValsProb.at(i);
    else probSumGreater50 += yValsProb.at(i);
  }
  
  
  int npoints = 0;
  double sumNorm = 0;
  for(unsigned int i = 0; i < xVals.size(); ++i){
    if(i%points == 0 && i != 0){
      graphPotentialProb2_p->SetPoint((i-1)/points, xVals.at(i), points*yValsProb.at(i)/probSum);
      sumNorm += points*yValsProb.at(i)/probSum;
      ++npoints;
    }

    if(i == 0) continue;
    if(maxProb2 < yValsProb.at(i)*points/probSum) maxProb2 = yValsProb.at(i)*points/probSum;
    if(minProb2 > yValsProb.at(i)*points/probSum) minProb2 = yValsProb.at(i)*points/probSum;
  }
  std::cout << "ProbSum: " << probSum << ", " << probSumLess50 << ", " << probSumGreater50 << ", " << npoints << std::endl;

  std::cout << "Sum norm: " << sumNorm << std::endl;

  maxE *= 5;
  minE /= 2;

  maxDenom *= 5;
  minDenom /= 2;

  maxNum *= 5;
  minNum /= 2;


  std::cout << "Maxprob: " << maxProb << ", " << minProb << std::endl;
  maxProb *= 5;
  minProb /= 2;

  maxProb2 *= 1.2;
  minProb2 /= 2;

  canvPotentialE_c->cd();
  dummyE_p->SetMaximum(maxE);
  dummyE_p->SetMinimum(minE);
  gStyle->SetOptStat(0);
  dummyE_p->DrawCopy();
  gPad->SetLogx();
  gPad->SetLogy();
  graphPotentialE_p->Draw("P");
  canvPotentialE_c->SaveAs("pdfDir/canvPotential_E.pdf");
  delete graphPotentialE_p;
  delete canvPotentialE_c;

  canvPotentialDenom_c->cd();
  dummyDenom_p->SetMaximum(maxDenom);
  dummyDenom_p->SetMinimum(minDenom);
  dummyDenom_p->DrawCopy();
  gPad->SetLogx();
  gPad->SetLogy();
  graphPotentialDenom_p->Draw("P");
  canvPotentialDenom_c->SaveAs("pdfDir/canvPotential_Denom.pdf");
  delete graphPotentialDenom_p;
  delete canvPotentialDenom_c;

  canvPotentialNum_c->cd();
  dummyNum_p->SetMaximum(maxNum);
  dummyNum_p->SetMinimum(minNum);
  dummyNum_p->DrawCopy();
  graphPotentialNum_p->Draw("P");
  canvPotentialNum_c->SaveAs("pdfDir/canvPotential_Num.pdf");
  delete graphPotentialNum_p;
  delete canvPotentialNum_c;

  canvPotentialProb_c->cd();
  dummyProb_p->SetMaximum(maxProb);
  dummyProb_p->SetMinimum(minProb);
  dummyProb_p->DrawCopy();
  gPad->SetLogx();
  gPad->SetLogy();
  graphPotentialProb_p->Draw("P");
  canvPotentialProb_c->SaveAs("pdfDir/canvPotential_Prob.pdf");
  delete graphPotentialProb_p;
  delete canvPotentialProb_c;

  canvPotentialProb2_c->cd();
  dummyProb2_p->SetMaximum(maxProb2);
  dummyProb2_p->SetMinimum(minProb2);
  dummyProb2_p->DrawCopy();
  graphPotentialProb2_p->Draw("P");
  canvPotentialProb2_c->SaveAs("pdfDir/canvPotential_Prob2.pdf");
  delete graphPotentialProb2_p;
  delete canvPotentialProb2_c;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += stochPotentialProb();
  return retVal;
}
