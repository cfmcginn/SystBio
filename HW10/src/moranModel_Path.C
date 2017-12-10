#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TDatime.h"
#include "TLatex.h"
#include "TLine.h"

#include "include/doGlobalDebug.h"
#include "include/plotUtilities.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/kirchnerPalette.h"

int moranModel_Path(const unsigned int nSim = 500, const int nPopSize=10)
{
  TRandom3* randGen_p = new TRandom3(0);
  const double s0A = .05;
  const double s0B = 0.2;
  const double s0AB = 0.5;

  const Double_t mutRateLow = 0.00001;
  const Double_t mutRateHi = 1.;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  const Int_t nMuteRate = 40;

  Double_t mutRate[nMuteRate+1];
  getLogBins(mutRateLow, mutRateHi, nMuteRate, mutRate);

  TDatime* date = new TDatime();

  TFile* outFile_p = new TFile(("output/moranModel_Path_" + std::to_string(date->GetDate()) + ".root").c_str(), "UPDATE");
  TH1F* sGen_Win_p[nMuteRate];

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t muI = 0; muI < nMuteRate; ++muI){
    sGen_Win_p[muI] = new TH1F(("sGen_Win_nSim" + std::to_string(nSim) + "_NPop" + std::to_string(nPopSize) + "_mutRate" + prettyString(mutRate[muI], 6, true) + "_h").c_str(), ";s_{Gen. Win};Counts", 2, s0A - (s0B - s0A)/2., s0B + (s0B - s0A)/2.);
    sGen_Win_p[muI]->Sumw2();
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t muI = 0; muI < nMuteRate; ++muI){
    std::cout << "Running mutation: " << mutRate[muI] << std::endl;

    std::vector<double> sGen_Win_;
    
    int counter = 0;

    while(nSim > sGen_Win_.size()){
      std::vector<int> pop;
      std::vector<double> relFit;
      std::vector<std::string> popPath;
      
      pop.push_back(nPopSize-1);
      pop.push_back(0);
      pop.push_back(0);
      pop.push_back(0);
      pop.push_back(0);

      if(counter%2 == 0) pop.at(1) = 1;
      else pop.at(2) = 1;
      ++counter;
      
      relFit.push_back(1.);
      relFit.push_back(1. + s0A);
      relFit.push_back(1. + s0B);
      relFit.push_back(1. + s0AB);
      relFit.push_back(1. + s0AB);

      popPath.push_back("nom");
      popPath.push_back("a");
      popPath.push_back("b");
      popPath.push_back("ab");
      popPath.push_back("ba");
      
      Int_t pops = 1;

      while(pops != nPopSize){
	double denom = 0.;
	for(unsigned int jI = 0; jI < pop.size(); ++jI){
	  denom += pop.at(jI)*relFit.at(jI);
	}
	
	std::vector<double> probRep;
	std::vector<double> probDie;

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

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
	double draw2 = randGen_p->Uniform(0,1);
      
	if(draw2 < mutRate[muI] && repPos != 3 && repPos != 4){
	  if(repPos == 1) pop.at(3) += 1;
	  else if(repPos == 2) pop.at(4) += 1;
	  else{
	    double draw3 = randGen_p->Uniform(0, 1);	    
	    if(draw3 < .5) pop.at(2) += 1;
	    else pop.at(1) += 1;
	  }
	}
	else pop.at(repPos) += 1;
      	
	pops = pop.at(0);
	for(unsigned i = 1; i < pop.size(); ++i){
	  if(pop.at(i) > pops) pops = pop.at(i);
	}
      }

      if(pop.at(1) == nPopSize || pop.at(3) == nPopSize) sGen_Win_.push_back(relFit.at(1)-1.);
      else if(pop.at(2) == nPopSize || pop.at(4) == nPopSize) sGen_Win_.push_back(relFit.at(2)-1.);
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


    for(unsigned int i = 0; i < sGen_Win_.size(); ++i){
      sGen_Win_p[muI]->Fill(sGen_Win_.at(i));
    }
  }

  outFile_p->cd();

  TCanvas* canv_c = new TCanvas("canv_c", "canv_c", 500, 500);
  prettyCanv(canv_c);
  TH1F* dummyHist_p = new TH1F("dummyHist_h", ";Mutation Rate;Probability Path(s0B)", 10, mutRate[0]/10., mutRate[nMuteRate-1]*10);
  dummyHist_p->GetXaxis()->CenterTitle();
  dummyHist_p->GetYaxis()->CenterTitle();

  dummyHist_p->SetMaximum(1.0);
  dummyHist_p->SetMinimum(0);
  dummyHist_p->DrawCopy();
  gStyle->SetOptStat(0);
  gPad->SetLogx();

  TGraphErrors* graph = new TGraphErrors();

  for(Int_t muI = 0; muI < nMuteRate; ++muI){
    double s0ATemp = sGen_Win_p[muI]->GetBinContent(1);
    double s0BTemp = sGen_Win_p[muI]->GetBinContent(2);

    std::cout << muI << ", " << s0BTemp << ", " << s0ATemp << std::endl;

    graph->SetPoint(muI, mutRate[muI], s0BTemp/(s0ATemp+s0BTemp));
    graph->SetPointError(muI, 0, TMath::Sqrt(s0BTemp)/(s0ATemp+s0BTemp));

    sGen_Win_p[muI]->Scale(1./sGen_Win_p[muI]->Integral());
    sGen_Win_p[muI]->Write("", TObject::kOverwrite);
    delete sGen_Win_p[muI];
  }

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.);

  graph->Draw("P");

  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);
  line_p->DrawLine(mutRate[0]/10., s0B/(s0B+s0A), mutRate[nMuteRate-1]*10., s0B/(s0B+s0A));
  delete line_p;

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);

  label_p->DrawLatex(0.15, .95, ("N_{pop}=" + std::to_string(nPopSize) + "; Simulate " + std::to_string(nSim) + " times").c_str());
  label_p->DrawLatex(.2, .4, ("s0A = " + prettyString(s0A, 3, false)).c_str());
  label_p->DrawLatex(.2, .35, ("s0B = " + prettyString(s0B, 3, false)).c_str());
  label_p->DrawLatex(.2, .3, ("s0AB = " + prettyString(s0AB, 3, false)).c_str());

  canv_c->SaveAs("pdfDir/moranModel_Path.pdf");

  delete canv_c;
  delete graph;
  delete dummyHist_p;


  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  delete date;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 2 && argc != 3){
    std::cout << "Usage: ./moranModel_Path.exe <nSim-optional> <nPopSize-optional>" << std::endl;
  }

  int retVal = 0;
  if(argc == 1) retVal += moranModel_Path();
  else if(argc == 2) retVal += moranModel_Path(std::stoi(argv[1]));
  else if(argc == 3) retVal += moranModel_Path(std::stoi(argv[1]), std::stoi(argv[2]));
  return retVal;
}
