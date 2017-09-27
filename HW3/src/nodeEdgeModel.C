#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

int nodeEdgeModel()
{
  const int timeSteps = 20;

  Double_t bins[timeSteps+1];
  for(int i = 0; i < timeSteps+1; ++i){
    if(i == 0) bins[i] = 1;
    else bins[i] = (TMath::Power(2.,i) + TMath::Power(2.,i+1))/2.;
  }

  TFile* nodeEdgeFile_p = new TFile("nodeEdgeFile.root", "RECREATE");
  TH1F* nodeHist_p = new TH1F("nodeHist_h", ";Time (steps);Nodes", timeSteps, -0.5, timeSteps-.5);
  TH1F* edgeHist_p = new TH1F("edgeHist_h", ";Time (steps);Edges", timeSteps, -0.5, timeSteps-.5);

  TH1F* edgesPerNode_p = new TH1F("edgesPerNode_h", ";Number of Edges;Nodes", timeSteps, bins);

  std::vector<unsigned long long> edgesPerNode;

  nodeHist_p->SetBinContent(1, 2);
  nodeHist_p->SetBinError(1, 0);

  edgeHist_p->SetBinContent(1, 1);
  edgeHist_p->SetBinError(1, 0);

  edgesPerNode.push_back(1);
  edgesPerNode.push_back(1);

  for(int i = 1; i < timeSteps; ++i){
    nodeHist_p->SetBinContent(i+1, edgeHist_p->GetBinContent(i) + nodeHist_p->GetBinContent(i));
    nodeHist_p->SetBinError(i+1, 0);

    edgeHist_p->SetBinContent(i+1, 3*edgeHist_p->GetBinContent(i));
    edgeHist_p->SetBinError(i+1, 0);

    for(unsigned int j = 0; j < edgesPerNode.size(); ++j){
      edgesPerNode.at(j) = edgesPerNode.at(j)*2;
    }

    for(int j = 0; j < edgeHist_p->GetBinContent(i); ++j){
      edgesPerNode.push_back(2);
    }
  }

  nodeEdgeFile_p->cd();

  nodeHist_p->Write("", TObject::kOverwrite);
  edgeHist_p->Write("", TObject::kOverwrite);

  for(unsigned int i = 0; i < edgesPerNode.size(); ++i){
    edgesPerNode_p->Fill(edgesPerNode.at(i));
  }

  edgesPerNode_p->Write("", TObject::kOverwrite);

  delete nodeHist_p;
  delete edgeHist_p;

  delete edgesPerNode_p;
  edgesPerNode.clear();

  nodeEdgeFile_p->Close();
  delete nodeEdgeFile_p;

  return 0;
}

int main()
{

  int retVal = 0;
  retVal += nodeEdgeModel();
  return retVal;
}
