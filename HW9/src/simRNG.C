#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int simRNG()
{
  kirchnerPalette col;

  const Int_t nSim = 100000;

  const Int_t nRandNorm = 75;
  const Double_t meanNorm = 0.0;
  const Double_t sigmaNorm = 1.35;

  const Int_t nRandExp = 400;
  const Double_t meanExp = 0.5;
  
  std::vector<double> normVals;
  normVals.reserve(nSim);

  std::vector<double> expVals;
  expVals.reserve(nSim);

  TRandom3* randGen_p = new TRandom3(0);
  
  TFile* outFile_p = new TFile("output/simRNG.root", "RECREATE");

  for(Int_t i = 0; i < nSim; ++i){
    double normMax= -10000000;
    double expMax= -10000000;

    for(Int_t j = 0; j < nRandNorm; ++j){
      double temp = randGen_p->Gaus(meanNorm, sigmaNorm);
      if(temp > normMax) normMax= temp;
    }

    for(Int_t j = 0; j < nRandExp; ++j){
      double temp = randGen_p->Exp(meanExp);
      if(temp > expMax) expMax= temp;
    }

    normVals.push_back(normMax);
    expVals.push_back(expMax);
  }
  
  std::sort(std::begin(normVals), std::end(normVals));
  std::sort(std::begin(expVals), std::end(expVals));

  std::cout << "Norm low,high: " << normVals.at(0) << ", " << normVals.at(normVals.size()-1) << std::endl;
  std::cout << "Exp low,high: " << expVals.at(0) << ", " << expVals.at(expVals.size()-1) << std::endl;


  Double_t normMin = normVals.at(0);
  Double_t normMax = normVals.at(normVals.size()-1);

  Double_t expMin = expVals.at(0);
  Double_t expMax = expVals.at(expVals.size()-1);

  Double_t globalMin = TMath::Min(normMin, expMin);
  Double_t globalMax = TMath::Max(normMax, expMax);


  Double_t globalInt = globalMax - globalMin;
  globalMin -= globalInt/10.;
  globalMax += globalInt/10.;

  const int nGlobalBins = 200;

  TH1F* histNorm_p = new TH1F("histNorm_h", ";S_{n} = Max_{n}(s_{i});Probability", nGlobalBins, globalMin, globalMax);
  TH1F* histExp_p = new TH1F("histExp_h", ";S_{n} = Max_{n}(s_{i});Probability", nGlobalBins, globalMin, globalMax);
  TH1F* histNorm_Check_p = new TH1F("histNorm_Check_h", ";;", nGlobalBins, globalMin, globalMax);
  TH1F* histExp_Check_p = new TH1F("histExp_Check_h", ";;", nGlobalBins, globalMin, globalMax);
  TH1F* histNorm_Cum_p = new TH1F("histNorm_Cum_h", ";S_{n} = Max_{n}(s_{i});Cumulative Probability", nGlobalBins, globalMin, globalMax);
  TH1F* histExp_Cum_p = new TH1F("histExp_Cum_h", ";S_{n} = Max_{n}(s_{i});Cumulative Probability", nGlobalBins, globalMin, globalMax);
  TH1F* histNorm_CumCheck_p = new TH1F("histNorm_CumCheck_h", ";;", nGlobalBins, globalMin, globalMax);
  TH1F* histExp_CumCheck_p = new TH1F("histExp_CumCheck_h", ";;", nGlobalBins, globalMin, globalMax);

  histNorm_p->Sumw2();
  histExp_p->Sumw2();
  histNorm_Cum_p->Sumw2();
  histExp_Cum_p->Sumw2();
  histNorm_CumCheck_p->Sumw2();
  histExp_CumCheck_p->Sumw2();

  TF1* gausCheck_p = new TF1("gausCheck_p", "gaus", -1000, 1000);
  gausCheck_p->SetParameter(0, 1);
  gausCheck_p->SetParameter(1, meanNorm);
  gausCheck_p->SetParameter(2, sigmaNorm);

  TF1* expCheck_p = new TF1("expCheck_p", "expo", -1000, 1000);
  expCheck_p->SetParameter(0, 1);
  expCheck_p->SetParameter(1, -1./meanExp);

  expCheck_p->Print("ALL");

  for(Int_t i = 0; i < histNorm_CumCheck_p->GetNbinsX(); ++i){
    histNorm_CumCheck_p->SetBinContent(i+1, gausCheck_p->Eval(histNorm_CumCheck_p->GetBinCenter(i+1)));
    histExp_CumCheck_p->SetBinContent(i+1, expCheck_p->Eval(histExp_CumCheck_p->GetBinCenter(i+1)));
  }

  double lostValNorm = 0;
  double lostValExp = 0;

  for(Int_t i = 0; i < nGlobalBins; ++i){
    Double_t binCenter = histNorm_CumCheck_p->GetBinLowEdge(1) - histNorm_CumCheck_p->GetBinWidth(i)*(i+.5);

    lostValNorm += gausCheck_p->Eval(binCenter);
    if(binCenter > 0) lostValExp += expCheck_p->Eval(binCenter);
  }

  std::cout << lostValNorm << std::endl;
  std::cout << lostValExp << std::endl;

  double cumLostValNorm = lostValNorm/(histNorm_CumCheck_p->Integral() + lostValNorm);
  double cumLostValExp = lostValExp/(histExp_CumCheck_p->Integral() + lostValExp);

  histNorm_CumCheck_p->Scale(1./(histNorm_CumCheck_p->Integral() + lostValNorm));
  histExp_CumCheck_p->Scale(1./(histExp_CumCheck_p->Integral() + lostValExp));

  histNorm_CumCheck_p->SetBinContent(1, histNorm_CumCheck_p->GetBinContent(1) + cumLostValNorm);
  histExp_CumCheck_p->SetBinContent(1, histExp_CumCheck_p->GetBinContent(1) + cumLostValExp);

  for(Int_t i = 1; i < histNorm_CumCheck_p->GetNbinsX(); ++i){
    histNorm_CumCheck_p->SetBinContent(i+1, histNorm_CumCheck_p->GetBinContent(i) + histNorm_CumCheck_p->GetBinContent(i+1));
    histExp_CumCheck_p->SetBinContent(i+1, histExp_CumCheck_p->GetBinContent(i) + histExp_CumCheck_p->GetBinContent(i+1));
  }

  for(Int_t i = 0; i < histNorm_CumCheck_p->GetNbinsX(); ++i){
    histNorm_Check_p->SetBinContent(i+1, nRandNorm*TMath::Power(histNorm_CumCheck_p->GetBinContent(i+1), nRandNorm-1)*gausCheck_p->Eval(histNorm_Check_p->GetBinCenter(i+1)));
    histExp_Check_p->SetBinContent(i+1, nRandExp*TMath::Power(histExp_CumCheck_p->GetBinContent(i+1), nRandExp-1)*expCheck_p->Eval(histExp_Check_p->GetBinCenter(i+1)));

    histNorm_CumCheck_p->SetBinContent(i+1, TMath::Power(histNorm_CumCheck_p->GetBinContent(i+1), nRandNorm));
    histExp_CumCheck_p->SetBinContent(i+1, TMath::Power(histExp_CumCheck_p->GetBinContent(i+1), nRandExp));
  }

  histNorm_Check_p->Scale(1./histNorm_Check_p->Integral());
  histExp_Check_p->Scale(1./histExp_Check_p->Integral());

  for(unsigned int i = 0; i < normVals.size(); ++i){
    histNorm_p->Fill(normVals.at(i));
  }

  for(unsigned int i = 0; i < expVals.size(); ++i){
    histExp_p->Fill(expVals.at(i));
  }

  histNorm_p->Scale(1./histNorm_p->Integral());
  histExp_p->Scale(1./histExp_p->Integral());

  histNorm_Cum_p->SetBinContent(1, histNorm_p->GetBinContent(1));
  histExp_Cum_p->SetBinContent(1, histExp_p->GetBinContent(1));

  for(Int_t bI = 1; bI < histNorm_p->GetNbinsX(); ++bI){
    histNorm_Cum_p->SetBinContent(bI+1, histNorm_Cum_p->GetBinContent(bI) + histNorm_p->GetBinContent(bI+1));
    histExp_Cum_p->SetBinContent(bI+1, histExp_Cum_p->GetBinContent(bI) + histExp_p->GetBinContent(bI+1));
  }

  outFile_p->cd();

  histNorm_p->Write("", TObject::kOverwrite);
  histExp_p->Write("", TObject::kOverwrite);

  histNorm_Check_p->Write("", TObject::kOverwrite);
  histExp_Check_p->Write("", TObject::kOverwrite);

  histNorm_Cum_p->Write("", TObject::kOverwrite);
  histExp_Cum_p->Write("", TObject::kOverwrite);

  histNorm_CumCheck_p->Write("", TObject::kOverwrite);
  histExp_CumCheck_p->Write("", TObject::kOverwrite);

  TDatime* date = new TDatime();
  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  histNorm_p->GetXaxis()->CenterTitle();
  histNorm_p->GetYaxis()->CenterTitle();
  histNorm_p->SetMarkerColor(col.getColor(0));
  histNorm_p->SetLineColor(col.getColor(0));
  histNorm_p->SetMarkerStyle(20);
  histNorm_p->SetMarkerSize(.5);

  histExp_p->GetXaxis()->CenterTitle();
  histExp_p->GetYaxis()->CenterTitle();
  histExp_p->SetMarkerColor(col.getColor(2));
  histExp_p->SetLineColor(col.getColor(2));
  histExp_p->SetMarkerStyle(20);
  histExp_p->SetMarkerSize(.5);

  histNorm_p->SetMinimum(0.0);
  histNorm_p->SetMaximum(TMath::Max(histNorm_p->GetMaximum(), histExp_p->GetMaximum())*1.1);

  histNorm_p->DrawCopy("E1 P");
  histExp_p->DrawCopy("E1 P SAME");

  gStyle->SetOptStat(0);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);

  TLegend* leg_p = new TLegend(0.5, .3, 0.95, .7);
  leg_p->SetFillStyle(0);
  leg_p->SetBorderSize(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(18);
  leg_p->AddEntry(histNorm_p, ("Normal (#mu=" + prettyString(meanNorm, 1, false) + ";#sigma=" + prettyString(sigmaNorm, 2, false) + ")").c_str(), "P L");
  leg_p->AddEntry(histExp_p, ("Exponential (#tau=" + prettyString(meanExp,2,false) + ")").c_str() , "P L");

  leg_p->Draw("SAME");

  label_p->DrawLatex(.25, .95, "Probability Density");

  canv_p->SaveAs(("pdfDir/simulate_" + std::to_string(date->GetDate()) + ".pdf").c_str());
  delete canv_p;

  canv_p = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_p);
  histNorm_Cum_p->GetXaxis()->CenterTitle();
  histNorm_Cum_p->GetYaxis()->CenterTitle();
  histNorm_Cum_p->SetMarkerColor(col.getColor(0));
  histNorm_Cum_p->SetLineColor(col.getColor(0));
  histNorm_Cum_p->SetMarkerStyle(20);
  histNorm_Cum_p->SetMarkerSize(.5);

  histExp_Cum_p->GetXaxis()->CenterTitle();
  histExp_Cum_p->GetYaxis()->CenterTitle();
  histExp_Cum_p->SetMarkerColor(col.getColor(2));
  histExp_Cum_p->SetLineColor(col.getColor(2));
  histExp_Cum_p->SetMarkerStyle(20);
  histExp_Cum_p->SetMarkerSize(.5);

  histNorm_Cum_p->SetMinimum(0.0);
  histNorm_Cum_p->SetMaximum(TMath::Max(histNorm_Cum_p->GetMaximum(), histExp_Cum_p->GetMaximum())*1.1);

  histNorm_Cum_p->DrawCopy("E1 P");
  histExp_Cum_p->DrawCopy("E1 P SAME");

  gStyle->SetOptStat(0);

  leg_p->Draw("SAME");
  label_p->DrawLatex(.25, .95, "Cumulative Probability ");

  canv_p->SaveAs(("pdfDir/simulate_Cumulative_" + std::to_string(date->GetDate()) + ".pdf").c_str());
  delete canv_p;



  delete histNorm_p;
  delete histExp_p;

  delete histNorm_Check_p;
  delete histExp_Check_p;

  delete histNorm_Cum_p;
  delete histExp_Cum_p;

  delete histNorm_CumCheck_p;
  delete histExp_CumCheck_p;

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += simRNG();
  return retVal;
}
