#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"

#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int stochBistable()
{
  kirchnerPalette col;

  double gammaVal = 1.;
  double k1k2_1 = 1;
  //  double k1k2_2 = 1.;
  double v0 = 12.5;
  double v1 = 200.;

  double intervalEmu = .0001;
  double maxX = 300;
  
  std::vector<double> xVals;
  std::vector<double> yValsProd;
  std::vector<double> yValsDeg;
  
  xVals.push_back(0);
  yValsProd.push_back(v0);
  yValsDeg.push_back(0);

  bool isProdTop = true;

  while(xVals.at(xVals.size()-1) < maxX){
    xVals.push_back(xVals.at(xVals.size()-1) + intervalEmu);
    
    double prodVal = (v0 + v1*k1k2_1*xVals.at(xVals.size()-1)*xVals.at(xVals.size()-1))/(1 + k1k2_1*xVals.at(xVals.size()-1)*xVals.at(xVals.size()-1));
    double degVal = gammaVal*xVals.at(xVals.size()-1);

    if(isProdTop){
      if(prodVal < degVal){
	std::cout << prodVal << ", " << degVal << std::endl;
	isProdTop = false;
      }
    }
    else{
      if(prodVal > degVal){
	std::cout << prodVal << ", " << degVal << std::endl;
	isProdTop = true;
      }
    }

    yValsProd.push_back(prodVal);
    yValsDeg.push_back(degVal);
  }

  TCanvas* canvBistable_c = new TCanvas("canvBistable_c", "canvBistable_c", 500, 500);
  prettyCanv(canvBistable_c);

  TGraph* graphBistableProd_p = new TGraph();
  TGraph* graphBistableDeg_p = new TGraph();
  
  graphBistableProd_p->SetMarkerSize(.4);
  graphBistableProd_p->SetMarkerStyle(20);
  graphBistableProd_p->SetMarkerColor(col.getColor(0));

  graphBistableDeg_p->SetMarkerSize(.4);
  graphBistableDeg_p->SetMarkerStyle(20);
  graphBistableDeg_p->SetMarkerColor(col.getColor(1));

  for(unsigned int i = 0; i < xVals.size(); ++i){
    if(i%100 == 0) graphBistableProd_p->SetPoint(i/100, xVals.at(i), yValsProd.at(i));
    if(i%100 == 0) graphBistableDeg_p->SetPoint(i/100, xVals.at(i), yValsDeg.at(i));
  }

  canvBistable_c->cd();
  graphBistableDeg_p->Draw("A P");
  graphBistableProd_p->Draw("P");

  canvBistable_c->SaveAs("pdfDir/canvBistable.pdf");

  delete graphBistableProd_p;
  delete graphBistableDeg_p;
  delete canvBistable_c;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += stochBistable();
  return retVal;
}
