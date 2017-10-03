#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TRandom3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "include/plotUtilities.h"

int geneDuplication(const double prob)
{
  TRandom3* randGen_p = new TRandom3(0);
  
  std::cout << "Running probability p=" << prob << "..." << std::endl;

  std::map<unsigned long long, std::vector<unsigned long long> > nodeNet;
  
  std::vector<unsigned long long> temp;
  temp.push_back(1);
  
  nodeNet[0] = temp;
  temp.clear();
  temp.push_back(0);

  nodeNet[1] = temp;
  temp.clear();

  unsigned long long max = 1;

  for(unsigned long long i = 2; i < 1002; ++i){
    unsigned long long keyVal = 10000;

    //Done this way as one can generate a node that had been removed
    while(nodeNet.find(keyVal) == nodeNet.end()){
      keyVal = (unsigned long long)randGen_p->Integer(i);
    }

    temp = nodeNet[keyVal];

    std::vector<unsigned long long> newNodeSet;
    for(unsigned long long j = 0; j < temp.size(); ++j){
      if(randGen_p->Uniform(0., 1.) < prob){
	newNodeSet.push_back(temp.at(j));
      }
    }

    if(newNodeSet.size() > 0){
      if(newNodeSet.size() > max) max = newNodeSet.size();
      nodeNet[i] = newNodeSet;

      temp.push_back(i);
      nodeNet[keyVal] = temp;
      if(temp.size() > max) max = temp.size();
    }

    temp.clear();
    newNodeSet.clear();
  }

  std::cout << "Maximum degree k_{Max}=" << max << "." << std::endl;

  TFile* outFile_p = new TFile("output/nodeHists.root", "RECREATE");

  const std::string histName = "degreeHist_Prob" + prettyString(prob, 3, true) + "_h";
  TH1F* degreeHist_p = new TH1F(histName.c_str(), ";k;Events", max+1, 0.5, max+1.5);

  std::map<unsigned long long, std::vector<unsigned long long> >::iterator it = nodeNet.begin();

  while(it != nodeNet.end()){
    degreeHist_p->Fill(it->second.size());
    ++it;
  }

  TDatime* date = new TDatime();
  const std::string canvName = "degreeHist_Prob" + prettyString(prob, 3, true) + "_c";
  const std::string saveName = "pdfDir/degreeHist_Prob" + prettyString(prob, 3, true) + "_" + std::to_string(date->GetDate()) + ".pdf";
  delete date;
  
  TCanvas* canv_p = new TCanvas(canvName.c_str(), canvName.c_str(), 500, 500);
  prettyCanv(canv_p);
  
  outFile_p->cd();

  canv_p->cd();
  degreeHist_p->DrawCopy("HIST E1");
  degreeHist_p->Write("", TObject::kOverwrite);

  canv_p->SaveAs(saveName.c_str());

  delete canv_p;
  delete degreeHist_p;

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./geneDuplication.exe <probability>" << std::endl;
    return 0;
  }

  int retVal = 0;
  retVal += geneDuplication(std::stod(argv[1]));
  return retVal;
}
