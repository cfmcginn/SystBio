#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TDatime.h"
#include "TLatex.h"

#include "include/doGlobalDebug.h"
#include "include/plotUtilities.h"
#include "include/getLinBins.h"
#include "include/kirchnerPalette.h"

int moranModel(const unsigned int nSim = 500, const int nPopSize=10, const double mutRate=0.0)
{
  TRandom3* randGen_p = new TRandom3(0);
  const double s0 = .02;

  TDatime* date = new TDatime();

  TFile* outFile_p = new TFile(("output/moranModel_" + std::to_string(date->GetDate()) + ".root").c_str(), "UPDATE");
  TH1F* sGen_p = new TH1F(("sGen_nSim" + std::to_string(nSim) + "_NPop" + std::to_string(nPopSize) + "_mutRate" + prettyString(mutRate, 6, true) + "_h").c_str(), ";s_{Gen.};Counts", 100, 0, s0*10.);
  TH1F* sGen_Win_p = new TH1F(("sGen_Win_nSim " + std::to_string(nSim) + " _NPop" + std::to_string(nPopSize) + "_mutRate" + prettyString(mutRate, 6, true) + "_h").c_str(), ";s_{Gen. Win};Counts", 100, 0, s0*10.);
  sGen_p->Sumw2();
  sGen_Win_p->Sumw2();

  std::vector<double> sGen_Win_;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  bool didGraph = false;

  while(nSim > sGen_Win_.size()){
    std::vector<int> pop;
    std::vector<double> relFit;

    pop.push_back(nPopSize-1);
    pop.push_back(1);
    relFit.push_back(1.);
    double sGenStart = randGen_p->Exp(s0);
    sGen_p->Fill(sGenStart);
    relFit.push_back(1. + sGenStart);

    std::vector<int> pop1;
    std::vector<int> pop2;
    pop1.push_back(nPopSize-1);
    pop2.push_back(1);

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    //    std::cout << "Starting fitness: " << sGenStart << std::endl;

    int pops = 1;
    while(pops != nPopSize){
      unsigned int pos = 0;
      while(pos < pop.size()){
	if(pos == 0) ++pos;
	else if(pop.at(pos) == 0){
	  pop.erase(pop.begin() + pos);
	  relFit.erase(relFit.begin() + pos);
	}
	else ++pos;
      }

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

      const Int_t nMatrix = probRep.size();
      Double_t matrix[nMatrix][nMatrix];
      Double_t newDenom = 0.;

      for(Int_t mI = 0; mI < nMatrix; ++mI){
	for(Int_t nI  = 0; nI < nMatrix; ++nI){
	  matrix[mI][nI] = probRep.at(mI)*probDie.at(nI);
	  if(mI != nI) newDenom += matrix[mI][nI];
	}
      }
      for(Int_t mI = 0; mI < nMatrix; ++mI){
	for(Int_t nI  = 0; nI < nMatrix; ++nI){
	  matrix[mI][nI] /= newDenom;
	}
      }

      double draw1 = randGen_p->Uniform(0,1);
      
      double fullProb = 0;
      int repPos = -1;
      int deathPos = -1;
      bool isFound = false;

      for(Int_t mI = 0; mI < nMatrix; ++mI){
	for(Int_t nI  = 0; nI < nMatrix; ++nI){
	  if(mI == nI) continue;

	  fullProb += matrix[mI][nI];
	  if(draw1 < fullProb){
	    repPos = mI;
	    deathPos = nI;
	    isFound = true;
	    break;
	  }
	}
	if(isFound) break;
      }
					       
      pop.at(deathPos) -= 1;
      pop.at(repPos) += 1;

      pop1.push_back(pop.at(0));
      pop2.push_back(pop.at(1));

      double draw3 = randGen_p->Uniform(0,1);

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

      if(mutRate < 0.0000000001 && nPopSize == 49 && !didGraph && pop1.size() < 200){
	didGraph = true;

	for(unsigned int time = 0; time < pop1.size(); ++time){
	  TCanvas* canv_c = new TCanvas("canv_c", "canv_c", 1000, 1000);
	  const Double_t valLow = -3.5;
	  const Double_t valHi = 3.5;

	  kirchnerPalette col;

	  TH1F* dummyHist_h = new TH1F("dummyHist_h", ";;", 10, valLow, valHi);
	  TGraph* pop1_g = new TGraph();
	  TGraph* pop2_g = new TGraph();

	  TLatex* label_p = new TLatex();
	  label_p->SetNDC();
	  label_p->SetTextFont(43);
	  label_p->SetTextSize(30);


	  const Int_t nBinsTemp = 7;
	  Double_t xVals[nPopSize];
	  Double_t yVals[nPopSize];

	  for(Int_t pI = 0; pI < nPopSize; ++pI){
	    xVals[pI] = valLow + .5 + pI%nBinsTemp;
	    yVals[pI] = valLow + .5 + pI/nBinsTemp;
	  }

	  for(Int_t pI = 0; pI < nPopSize; ++pI){
	    if(pop1.at(time) > pI) pop1_g->SetPoint(pI, xVals[pI], yVals[pI]);
	    else pop2_g->SetPoint(pI - pop1.at(time), xVals[pI], yVals[pI]);
	  }

	  pop1_g->SetMarkerSize(2.);
	  pop2_g->SetMarkerSize(2.);
	  pop1_g->SetMarkerStyle(20);
	  pop2_g->SetMarkerStyle(20);

	  pop1_g->SetMarkerColor(col.getColor(2));
	  pop2_g->SetMarkerColor(col.getColor(0));

	  canv_c->cd();
	  dummyHist_h->SetMaximum(valHi);
	  dummyHist_h->SetMinimum(valLow);
	  dummyHist_h->DrawCopy("P");

	  gStyle->SetOptStat(0);

	  if(time != pop1.size()-1) pop1_g->Draw("P");
	  pop2_g->Draw("P");

	  std::string timeStr = std::to_string(time);
	  if(time < 10) timeStr = "00" + timeStr;
	  else if(time < 100) timeStr = "0" + timeStr;

	  label_p->DrawLatex(.15, .95, ("sGen=" + prettyString(1-relFit.at(1), 3, false) + "; Time=" + std::to_string(time)).c_str());

	  canv_c->SaveAs(("gifDir/takeOver_" + timeStr + ".gif").c_str());

	  delete label_p;
	  delete pop1_g;
	  delete pop2_g;
	  delete dummyHist_h;
	  delete canv_c;
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

  delete date;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 2 && argc != 3 && argc != 4){
    std::cout << "Usage: ./moranModel.exe <nSim-optional> <nPopSize-optional> <mutRate-optional>" << std::endl;
  }

  int retVal = 0;
  if(argc == 1) retVal += moranModel();
  else if(argc == 2) retVal += moranModel(std::stoi(argv[1]));
  else if(argc == 3) retVal += moranModel(std::stoi(argv[1]), std::stoi(argv[2]));
  else if(argc == 4) retVal += moranModel(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));
  return retVal;
}
