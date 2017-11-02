#include <iostream>
#include <map>
#include <vector>

#include "TFile.h"
#include "TH2I.h"
#include "TRandom3.h"
#include "TMath.h"

int boundaryWalk(const int nSize)
{
  if(nSize%2 == 0){
    std::cout << "Please use odd spacing for convenience" << std::endl;
    return 1;
  }

  const int nSimul = 100000;
  const int nStepsToCheck = 7;
  const int stepsToCheck[nStepsToCheck] = {10, 100, 1000, 5000, 10000, 50000, 100000};

  TRandom3* randGen_p = new TRandom3(0);

  TFile* outFile_p = new TFile("randomProb.root", "RECREATE");
  TH2F* probability_h[nStepsToCheck];
  TH2F* probability_Reduce_h[nStepsToCheck];

  for(Int_t i = 0; i < nStepsToCheck; ++i){
    probability_h[i] = new TH2F(("probability_" + std::to_string(stepsToCheck[i]) + "_h").c_str(), ";x;y", nSize, -(nSize/2) - .5, nSize/2 + .5, nSize, -(nSize/2) - .5, nSize/2 + .5);
    probability_h[i]->Sumw2();

    probability_Reduce_h[i] = new TH2F(("probability_" + std::to_string(stepsToCheck[i]) + "_Reduce_h").c_str(), ";x;y", 5, -(nSize/2) - .5, nSize/2 + .5, 5, -(nSize/2) - .5, nSize/2 + .5);
    probability_Reduce_h[i]->Sumw2();
  }

  for(int simIter = 0; simIter < nSimul; ++simIter){
    int posX = 0;
    int posY = 0;

    for(int i = 0; i < stepsToCheck[nStepsToCheck-1]+1; ++i){
      bool isX = true;
      if(randGen_p->Uniform(0,1) > .5) isX = false;

      bool isPos = true;
      if(randGen_p->Uniform(0,1) > .5) isPos = false;

      if(isX && isPos) ++posX;
      else if(isX && !isPos) --posX;
      else if(!isX && isPos) ++posY;
      else if(!isX && !isPos) --posY;
      
      if(posX == nSize/2 + 1) posX = -(nSize/2);
      else if(posX == -nSize/2 - 1) posX = (nSize/2);
      else if(posY == nSize/2 + 1) posY = -(nSize/2);
      else if(posY == -nSize/2 - 1) posY = (nSize/2);
      
      for(Int_t j = 0; j < nStepsToCheck; ++j){
	if(i == stepsToCheck[j]){
	  probability_h[j]->Fill(posX, posY);
	  probability_Reduce_h[j]->Fill(posX, posY);
	}
      }
    }
  }


  outFile_p->cd();

  for(Int_t i = 0; i < nStepsToCheck; ++i){
    probability_h[i]->Scale(1./probability_h[i]->GetEntries());
    probability_h[i]->Write("", TObject::kOverwrite);
    delete probability_h[i];

    probability_Reduce_h[i]->Scale(1./probability_Reduce_h[i]->GetEntries());
    probability_Reduce_h[i]->Write("", TObject::kOverwrite);
    delete probability_Reduce_h[i];
  }

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./boundaryWalk.exe <inSize>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += boundaryWalk(std::stoi(argv[1]));
  return retVal;
}
