#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1F.h"

#include "include/kirchnerPalette.h"

void simulData()
{
  kirchnerPalette col;

  const double startE = .01;
  const double startS = 1.;
  const double startES = 0.;
  const double startP = 0.;
  
  const double rateK1P = 1;
  const double rateK1M = .1;
  const double rateK2 = .001;

  const double interval5 = .01;
  const double interval100k = 1;
  const unsigned long long steps5 = 5./interval5;
  const unsigned long long steps100k = 100000./interval100k;

  TFile* outFile_p = new TFile("simulDataSystBio.root", "RECREATE");
  TH1F* eHist5_p = new TH1F("eHist5_h", ";;", steps5, 0, 5);
  eHist5_p->Sumw2();
  eHist5_p->SetMarkerColor(col.getColor(0));

  TH1F* sHist5_p = new TH1F("sHist5_h", ";;", steps5, 0, 5);
  sHist5_p->Sumw2();
  sHist5_p->SetMarkerColor(col.getColor(1));

  TH1F* esHist5_p = new TH1F("esHist5_h", ";;", steps5, 0, 5);
  esHist5_p->Sumw2();
  esHist5_p->SetMarkerColor(col.getColor(2));

  TH1F* pHist5_p = new TH1F("pHist5_h", ";;", steps5, 0, 5);
  pHist5_p->Sumw2();
  pHist5_p->SetMarkerColor(col.getColor(3));


  TH1F* eHist100k_p = new TH1F("eHist100k_h", ";;", steps100k, 0, 100000);
  eHist100k_p->Sumw2();
  eHist100k_p->SetMarkerColor(col.getColor(0));

  TH1F* sHist100k_p = new TH1F("sHist100k_h", ";;", steps100k, 0, 100000);
  sHist100k_p->Sumw2();
  sHist100k_p->SetMarkerColor(col.getColor(1));

  TH1F* esHist100k_p = new TH1F("esHist100k_h", ";;", steps100k, 0, 100000);
  esHist100k_p->Sumw2();
  esHist100k_p->SetMarkerColor(col.getColor(2));

  TH1F* pHist100k_p = new TH1F("pHist100k_h", ";;", steps100k, 0, 100000);
  pHist100k_p->Sumw2();
  pHist100k_p->SetMarkerColor(col.getColor(3));

  std::vector<double> pRateVsS;
  std::vector<double> pRateVsS_SVal;
  pRateVsS.reserve(100000);
  pRateVsS_SVal.reserve(100000);

  double currE = startE;
  double currS = startS;
  double currES = startES;
  double currP = startP;
  //  double currPRate

  int iter = 0;
  eHist5_p->SetBinContent(iter+1, currE);
  eHist5_p->SetBinError(iter+1, 0);

  sHist5_p->SetBinContent(iter+1, currS);
  sHist5_p->SetBinError(iter+1, 0);

  esHist5_p->SetBinContent(iter+1, currES);
  esHist5_p->SetBinError(iter+1, 0);

  pHist5_p->SetBinContent(iter+1, currP);
  pHist5_p->SetBinError(iter+1, 0);

  eHist100k_p->SetBinContent(iter+1, currE);
  eHist100k_p->SetBinError(iter+1, 0);

  sHist100k_p->SetBinContent(iter+1, currS);
  sHist100k_p->SetBinError(iter+1, 0);

  esHist100k_p->SetBinContent(iter+1, currES);
  esHist100k_p->SetBinError(iter+1, 0);

  pHist100k_p->SetBinContent(iter+1, currP);
  pHist100k_p->SetBinError(iter+1, 0);

  ++iter;

  for(unsigned long long i = 0; i < steps5; ++i){
    if(i > (unsigned long long)eHist5_p->GetNbinsX()) break;

    double deltaE = interval5*((rateK1M + rateK2)*currES - rateK1P*currE*currS);
    currE += deltaE;

    double deltaS = interval5*((rateK1M)*currES - rateK1P*currE*currS);
    currS += deltaS;

    double deltaES = interval5*((rateK1P)*currE*currS - (rateK1M + rateK2)*currES);
    currES += deltaES;

    double deltaP = interval5*((rateK2)*currES);
    currP += deltaP;
    
    eHist5_p->SetBinContent(iter+1, currE);
    eHist5_p->SetBinError(iter+1, 0);

    sHist5_p->SetBinContent(iter+1, currS);
    sHist5_p->SetBinError(iter+1, 0);

    esHist5_p->SetBinContent(iter+1, currES);
    esHist5_p->SetBinError(iter+1, 0);

    pHist5_p->SetBinContent(iter+1, currP);
    pHist5_p->SetBinError(iter+1, 0);
    
    ++iter;
  }

  currE = startE;
  currS = startS;
  currES = startES;
  currP = startP;

  pRateVsS.push_back(rateK2*currES);
  pRateVsS_SVal.push_back(currS);

  iter = 1;
  for(unsigned long long i = 0; i < steps100k; ++i){
    if(i > (unsigned long long)eHist100k_p->GetNbinsX()) break;
    double deltaE = interval100k*((rateK1M + rateK2)*currES - rateK1P*currE*currS);
    currE += deltaE;

    double deltaS = interval100k*((rateK1M)*currES - rateK1P*currE*currS);
    currS += deltaS;

    double deltaES = interval100k*((rateK1P)*currE*currS - (rateK1M + rateK2)*currES);
    currES += deltaES;

    double deltaP = interval100k*((rateK2)*currES);
    currP += deltaP;
    
    eHist100k_p->SetBinContent(iter+1, currE);
    eHist100k_p->SetBinError(iter+1, 0);

    sHist100k_p->SetBinContent(iter+1, currS);
    sHist100k_p->SetBinError(iter+1, 0);

    esHist100k_p->SetBinContent(iter+1, currES);
    esHist100k_p->SetBinError(iter+1, 0);

    pHist100k_p->SetBinContent(iter+1, currP);
    pHist100k_p->SetBinError(iter+1, 0);

    pRateVsS.push_back(rateK2*currES);
    pRateVsS_SVal.push_back(currS);
    
    ++iter;
  }


  const int th1N = pRateVsS.size();
  Double_t binsInS[th1N+1];

  for(unsigned int i = 1; i < pRateVsS.size(); ++i){
    binsInS[i] = (pRateVsS_SVal.at(i-1) + pRateVsS_SVal.at(i))/2;
  }

  binsInS[0] = binsInS[1] - (binsInS[2]-binsInS[1]);
  binsInS[th1N] = binsInS[th1N-1] + (binsInS[th1N-1]-binsInS[th1N-2]);
  std::cout << th1N << std::endl;
  std::cout << binsInS[0] << ", " << binsInS[1] << ", " << binsInS[th1N-1] << ", " << binsInS[th1N] << std::endl;
  //  binsInS[0] = 
  

  outFile_p->cd();


  std::cout << __LINE__ << std::endl;

  eHist5_p->Write("", TObject::kOverwrite);
  std::cout << __LINE__ << std::endl;
  delete eHist5_p;

  std::cout << __LINE__ << std::endl;

  sHist5_p->Write("", TObject::kOverwrite);
  delete sHist5_p;

  esHist5_p->Write("", TObject::kOverwrite);
  delete esHist5_p;

  pHist5_p->Write("", TObject::kOverwrite);
  delete pHist5_p;

  eHist100k_p->Write("", TObject::kOverwrite);
  delete eHist100k_p;

  sHist100k_p->Write("", TObject::kOverwrite);
  delete sHist100k_p;

  esHist100k_p->Write("", TObject::kOverwrite);
  delete esHist100k_p;

  pHist100k_p->Write("", TObject::kOverwrite);
  delete pHist100k_p;

  std::cout << __LINE__ << std::endl;
  
  outFile_p->Close();
  delete outFile_p;

  return;
}
