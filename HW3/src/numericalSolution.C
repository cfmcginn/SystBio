#include <iostream>

#include "TFile.h"
#include "TH1F.h"

#include "include/getLinBins.h"

int numericalSolution(const double gammaX, const double x0 = 1., const double y0 = 0.)
{
  const double vx = .1;
  const double kx = 4.;
  const double ky = 2.;

  const int nTimeBins = 10000;
  const float timeLow = 0.;
  const float timeHigh = 10.0;
  const double timeBins[nTimeBins+1];
  getLinBins(timeLow, timeHigh, nTimeBins, timeBins);
  const double timeWidth = timeBins[1] - timeBins[0];

  TFile* outFile_p = new TFile("numericalSolution.root", "RECREATE");
  TH1F* xHist_p = new TH1F("xHist_h", ";Time;x(t)", nTimeBins, timeBins);
  TH1F* yHist_p = new TH1F("yHist_h", ";Time;y(t)", nTimeBins, timeBins);
  
  xHist_p->SetBinContent(1, 1.);
  yHist_p->SetBinContent(1, .0);

  xHist_p->SetBinError(1, 0.);
  yHist_p->SetBinError(1, 0.);
  
  for(int iter = 1; iter < nTimeBins; ++iter){
    double xVal = xHist_p->GetBinContent(iter);
    double yVal = yHist_p->GetBinContent(iter);

    double newXVal = xVal + timeWidth*gammaX*(vx + kx*xVal*xVal/((1+yVal)*(1+xVal*xVal)) - xVal);
    double newYVal = yVal + timeWidth*(ky*xVal- yVal);

    xHist_p->SetBinContent(iter+1, newXVal);
    xHist_p->SetBinError(iter+1, 0);

    yHist_p->SetBinContent(iter+1, newYVal);
    yHist_p->SetBinError(iter+1, 0);
  }

  outFile_p->cd();

  xHist_p->Write("", TObject::kOverwrite);
  yHist_p->Write("", TObject::kOverwrite);

  delete xHist_p;
  delete yHist_p;

  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main()
{
  if(argv != 2 && argc != 3 && argc != 4){
    std::cout << "Usage ./numericalSolution.exe <gammaX> <optX0> <optY0>" << std::endl;
  }
  

  int retVal = 0;
  if(argc == 2) retVal += numericalSolution(std::stof(argv[1]));
  else if(argc == 3) retVal += numericalSolution(std::stof(argv[1]), std::stof(argv[2]));
  else if(argc == 4) retVal += numericalSolution(std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]));
  return retVal;
}
