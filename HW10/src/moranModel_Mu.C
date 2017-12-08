#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"

int moranModel(const int nPopSize = 10)
{
  TRandom3* randGen_p = new TRandom3(0);
  const double s0 = 0.02;
  const unsigned int nSim = 50000;

  TFile* outFile_p = new TFile("output/moranModel.root", "RECREATE");
  TH1F* sGen_p = new TH1F("sGen_h", ";s_{Gen.};Counts", 200, 0, .2);
  TH1F* sGen_Win_p = new TH1F("sGen_Win_h", ";s_{Gen. Win};Counts", 200, 0, .2);

  sGen_p->Sumw2();
  sGen_Win_p->Sumw2();

  std::vector<double> sGen_Win_;

  while(sGen_Win_.size() < nSim){
    const double sGen = randGen_p->Exp(s0);
    sGen_p->Fill(sGen);

    double startPop = nPopSize-1.;
    double mutPop = 1.;

    while(mutPop != 0 && mutPop != nPopSize){
      double rel = 1. + sGen;
      double prob1 = rel*mutPop/(rel*mutPop + startPop);
      double draw1 = randGen_p->Uniform(0,1);
      double draw2 = randGen_p->Uniform(0,1);

      if(draw1 < (double(startPop))/(double(nPopSize))) --startPop;
      else --mutPop;

      if(draw2 < prob1) ++mutPop;
      else ++startPop;
    }

    if(mutPop == nPopSize){
      if(sGen_Win_.size()%100 == 0) std::cout << sGen_Win_.size() << std::endl;
      sGen_Win_.push_back(sGen);
    }
  }

  for(unsigned int i = 0; i < sGen_Win_.size(); ++i){
    sGen_Win_p->Fill(sGen_Win_.at(i));
  }
  
  outFile_p->cd();
  
  sGen_p->Scale(1./sGen_p->Integral());
  sGen_Win_p->Scale(1./sGen_Win_p->Integral());

  sGen_p->Write("", TObject::kOverwrite);
  delete sGen_p;

  sGen_Win_p->Write("", TObject::kOverwrite);
  delete sGen_Win_p;

  outFile_p->Close();
  delete outFile_p;
  
  delete randGen_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 2){
    std::cout << "Usage: ./moranModel.exe <nPopSize-optional>" << std::endl;
  }

  int retVal = 0;
  if(argc == 1) retVal += moranModel();
  else retVal += moranModel(std::stoi(argv[1]));
  return retVal;
}
