#include <iostream>
#include <map>
#include <vector>

#include "TFile.h"
#include "TH2I.h"
#include "TRandom3.h"
#include "TMath.h"

int randomProbability(const int nSteps, const int nSimul)
{
  TRandom3* randGen_p = new TRandom3(0);

  TFile* outFile_p = new TFile("randomProb.root", "RECREATE");
  const int nGrid = 2*nSteps+1;
  TH2F* probability_h = new TH2F("probability_h", ";x;y", nGrid, -nSteps-.5, nSteps+.5, nGrid, -nSteps-.5, nSteps+.5);
  probability_h->Sumw2();

  std::map<int, std::vector<int> > pathsToPosMap; 

  for(int simIter = 0; simIter < nSimul; ++simIter){

    int posX = 0;
    int posY = 0;
    int path = 0;

    for(int i = 0; i < nSteps; ++i){
      bool isX = true;
      if(randGen_p->Uniform(0,1) > .5) isX = false;

      bool isPos = true;
      if(randGen_p->Uniform(0,1) > .5) isPos = false;

      if(isX && isPos){
	++posX;
	path += 1*TMath::Power(10, i);
      }
      else if(isX && !isPos){
	--posX;
	path += 2*TMath::Power(10, i);
      }
      else if(!isX && isPos){
	++posY;
	path += 3*TMath::Power(10, i);
      }
      else if(!isX && !isPos){
	--posY;
	path += 4*TMath::Power(10, i);
      }
    }

    probability_h->Fill(posX, posY);
    
    int mapPos = posX*100 + posY;

    std::vector<int> tempVect = pathsToPosMap[mapPos];
    
    bool isPathNew = true;
    for(unsigned int i = 0; i < tempVect.size(); ++i){
      if(path == tempVect.at(i)){
	isPathNew = false;
	break;
      }
    }

    if(isPathNew) tempVect.push_back(path);
    pathsToPosMap[mapPos] = tempVect;
  }

  int total = 0;

  for(std::map<int, std::vector<int>>::iterator it = pathsToPosMap.begin(); it != pathsToPosMap.end(); ++it){
    total += it->second.size();
  }

  int zeroVals = pathsToPosMap[0].size();
  std::vector<int> check = pathsToPosMap[0];

  for(unsigned int i =0; i < check.size(); ++i){
    std::cout << check.at(i) << std::endl;
  }


  std::cout << zeroVals << "/" << total << std::endl;

  outFile_p->cd();
  probability_h->Scale(1./probability_h->GetEntries());
  probability_h->Write("", TObject::kOverwrite);
  delete probability_h;
  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./randomProbability.exe <inNumberOfSteps> <inNumberSimul>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += randomProbability(std::stoi(argv[1]), std::stoi(argv[2]));
  return retVal;
}
