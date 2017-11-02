#include <iostream>
#include <map>
#include <vector>

#include "TFile.h"
#include "TH2I.h"
#include "TRandom3.h"
#include "TMath.h"

int boundaryWalkWithTarget()
{
  const int nSize = 101;
  
  const int nSimul = 100000;
  TRandom3* randGen_p = new TRandom3(0);

  TFile* outFile_p = new TFile("output/boundaryWalkWithTarget.root", "RECREATE");
  TH1F* stepsToTarget_p = new TH1F("stepsToTarget_h", ";Steps (N);Counts", 1000, -0.5, 99999.5);

  const int printInterval = 1000;

  for(int simIter = 0; simIter < nSimul; ++simIter){
    if(simIter%printInterval == 0) std::cout << "Sim " << simIter << "/" << nSimul << std::endl;

    int posX = 0;
    int posY = 0;

    int targetX = 0;
    double targetXPull = randGen_p->Uniform(-20.5, 20.5);
    for(Int_t i = 0; i < 41; ++i){
      if(targetXPull >= -20.5 + i && targetXPull <= -20.5 + (i+1)){
	targetX = -20 + i;
	break;
      }
    }

    bool targetYIsPos = true;
    if(randGen_p->Uniform(0,1) < .5) targetYIsPos = false;
    int targetY = 20 - TMath::Abs(targetX);
    if(!targetYIsPos) targetY *= -1;

    if(simIter%printInterval == 0) std::cout << "targetX,Y: " << targetX << ", " << targetY << std::endl;

    int stepCount = 0;

    while(posX != targetX || posY != targetY){
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

      ++stepCount;
    }

    stepsToTarget_p->Fill(stepCount);
  }


  outFile_p->cd();

  stepsToTarget_p->Write("", TObject::kOverwrite);
  delete stepsToTarget_p;

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += boundaryWalkWithTarget();
  return retVal;
}
