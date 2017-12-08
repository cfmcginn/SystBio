#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include "include/doGlobalDebug.h"
#include "include/plotUtilities.h"

int moranModel(const int nPopSize=10, const double mutRate=0.0)
{
  TRandom3* randGen_p = new TRandom3(0);
  const double s0 = 0.02;

  TFile* outFile_p = new TFile("output/moranModel.root", "UPDATE");
  TH1F* sGen_p = new TH1F(("sGen_NPop" + std::to_string(nPopSize) + "_mutRate" + prettyString(mutRate, 6, true) + "_h").c_str(), ";s_{Gen.};Counts", 100, 0, s0*10.);
  TH1F* sGen_Win_p = new TH1F(("sGen_Win_NPop" + std::to_string(nPopSize) + "_mutRate" + prettyString(mutRate, 6, true) + "_h").c_str(), ";s_{Gen. Win};Counts", 100, 0, s0*10.);
  sGen_p->Sumw2();
  sGen_Win_p->Sumw2();

  const unsigned int nSim = 5000;
  std::vector<double> sGen_Win_;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  while(nSim > sGen_Win_.size()){
    std::vector<int> pop;
    std::vector<double> relFit;

    pop.push_back(nPopSize-1);
    pop.push_back(1);
    relFit.push_back(1.);
    double sGenStart = randGen_p->Exp(s0);
    sGen_p->Fill(sGenStart);
    relFit.push_back(1. + sGenStart);

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    int pops = 1;
    while(pops != nPopSize){
      double denom = 0.;
      for(unsigned int jI = 0; jI < pop.size(); ++jI){
	denom += pop.at(jI)*relFit.at(jI);
      }

      std::vector<double> probRep;
      std::vector<double> probDie;

      for(unsigned int jI = 0; jI < pop.size(); ++jI){
	probRep.push_back(pop.at(jI)*relFit.at(jI)/denom);
	probDie.push_back(double(pop.at(jI))/(double)nPopSize);
      }
     
      double draw1 = randGen_p->Uniform(0,1);
      double draw2 = randGen_p->Uniform(0,1);
      double draw3 = randGen_p->Uniform(0,1);

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << pop.at(0) << ", " << pop.at(1) << std::endl;
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << probRep.at(0) << ", " << probRep.at(1) << std::endl;
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << probDie.at(0) << ", " << probDie.at(1) << std::endl;

      double probRepTot = 0.;
      double probDieTot = 0.;
 
      for(unsigned int jI = 0; jI < pop.size(); ++jI){
	probDieTot += probDie.at(jI);
	if(draw1 < probDieTot){
	  pop.at(jI) -= 1;
	  break;
	}
      }

      for(unsigned int jI = 0; jI < pop.size(); ++jI){
	probRepTot += probRep.at(jI);
	if(draw2 < probRepTot){
	  pop.at(jI) += 1;
	  break;
	}
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(draw3 < mutRate && pop.at(0) != 0){
	pop.at(0) -= 1;
	pop.push_back(1);
	double sGenNew = randGen_p->Exp(s0);
	sGen_p->Fill(sGenNew);
	relFit.push_back(1.+sGenNew);
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      pops = pop.at(0);
      for(unsigned int jI = 1; jI < pop.size(); ++jI){
	if(pop.at(jI) > pops) pops = pop.at(jI);
      }
    }

    if(pop.at(0) != nPopSize){
      if(sGen_Win_.size()%10 == 0) std::cout << "Size: " << sGen_Win_.size() << std::endl;

      for(unsigned int jI = 1; jI < pop.size(); ++jI){
	if(pop.at(jI) == nPopSize){
	  sGen_Win_.push_back(relFit.at(jI)-1.);
	  break;
	}
      }
    }
  }

  for(unsigned int i = 0; i < sGen_Win_.size(); ++i){
    sGen_Win_p->Fill(sGen_Win_.at(i));
  }

  outFile_p->cd();

  sGen_p->Scale(1./sGen_p->Integral());
  sGen_Win_p->Scale(1./sGen_Win_p->Integral());

  sGen_p->Write("", TObject::kOverwrite);
  sGen_Win_p->Write("", TObject::kOverwrite);

  delete sGen_p;
  delete sGen_Win_p;

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 2 && argc != 3){
    std::cout << "Usage: ./moranModel.exe <nPopSize-optional> <mutRate-optional>" << std::endl;
  }

  int retVal = 0;
  if(argc == 1) retVal += moranModel();
  else if(argc == 2) retVal += moranModel(std::stoi(argv[1]));
  else if(argc == 3) retVal += moranModel(std::stoi(argv[1]), std::stof(argv[2]));
  return retVal;
}
