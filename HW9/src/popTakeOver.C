#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"

#include "include/getLinBins.h"

int popTakeOver()
{
  const int nSize1 = 10;
  const int nSize2 = 1000;

  const Int_t nSBins = 10000;
  float sLow = 0.;
  float sHi = 0.01;
  Double_t sBins[nSBins+1];
  getLinBins(sLow, sHi, nSBins, sBins);

  TFile* outFile_p = new TFile("output/popTakeOver.root", "RECREATE");
  TH1F* uniform1_p = new TH1F("uniform1_h", ";;", nSBins, sBins);
  TH1F* uniform2_p = new TH1F("uniform2_h", ";;", nSBins, sBins);

  TH1F* expo1_p = new TH1F("expo1_h", ";;", nSBins, sBins);
  TH1F* expo2_p = new TH1F("expo2_h", ";;", nSBins, sBins);

  TF1* f1_p = new TF1("f1_p", "expo", sLow, sHi);
  f1_p->SetParameter(0, 1);
  f1_p->SetParameter(1, -1./0.01);

  f1_p->Print("ALL");

  for(Int_t i = 0; i < uniform1_p->GetNbinsX(); ++i){
    double cent = 1. + uniform1_p->GetBinCenter(i+1);
    double val1 = (1. - 1./cent)/(1. - 1./TMath::Power(cent, nSize1));
    double val2 = (1. - 1./cent)/(1. - 1./TMath::Power(cent, nSize2));

    uniform1_p->SetBinContent(i+1, val1);
    uniform2_p->SetBinContent(i+1, val2);

    expo1_p->SetBinContent(i+1, val1*f1_p->Eval(cent-1.));
    expo2_p->SetBinContent(i+1, val2*f1_p->Eval(cent-1.));
  }

  outFile_p->cd();

  uniform1_p->Write("", TObject::kOverwrite);
  uniform2_p->Write("", TObject::kOverwrite);

  expo1_p->Write("", TObject::kOverwrite);
  expo2_p->Write("", TObject::kOverwrite);

  delete uniform1_p;
  delete uniform2_p;

  delete expo1_p;
  delete expo2_p;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += popTakeOver();
  return retVal;
}
