#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"

#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int stochPotential()
{
  kirchnerPalette col;

  double gammaVal = 1.;
  double k1k2_1 = .0001;
  double k1k2_2 = 1.;
  double v0 = 12.5;
  double v1 = 200.;

  double intervalEmu = .000001;
  double maxX = 200;
  
  std::vector<double> xVals;
  std::vector<double> yVals_1;
  std::vector<double> yVals_2;
  
  xVals.push_back(0);
  yVals_1.push_back(v0*intervalEmu);
  yVals_2.push_back(v0*intervalEmu);

  bool isLess = true;
  while(xVals.at(xVals.size()-1) < maxX){
    double xVal = xVals.at(xVals.size()-1) + intervalEmu;
    double prodVal_1 = (v0 + v1*k1k2_1*xVal*xVal)/(1 + k1k2_1*xVal*xVal);
    double prodVal_2 = (v0 + v1*k1k2_2*xVal*xVal)/(1 + k1k2_2*xVal*xVal);
    double degVal = gammaVal*xVal;

    double val_1 = -2*intervalEmu*(prodVal_1 - degVal)/(prodVal_1 + degVal);
    double val_2 = -2*intervalEmu*(prodVal_2 - degVal)/(prodVal_2 + degVal);
    
    xVals.push_back(xVal);
    yVals_1.push_back(yVals_1.at(yVals_1.size()-1) + val_1);
    yVals_2.push_back(yVals_2.at(yVals_2.size()-1) + val_2);    

    if(isLess){
      if(yVals_1.at(yVals_1.size()-1) > yVals_1.at(yVals_1.size()-2)){
	std::cout << xVals.at(xVals.size()-1) << "-" << xVals.at(xVals.size()-2) << std::endl;
	isLess = false;
      }
    }
    else{
      if(yVals_1.at(yVals_1.size()-1) < yVals_1.at(yVals_1.size()-2)){
	std::cout << xVals.at(xVals.size()-1) << "-" << xVals.at(xVals.size()-2) << std::endl;
	isLess = true;
      }
    }
  }

  TCanvas* canvPotential_1_c = new TCanvas("canvPotential_1_c", "canvPotential_1_c", 500, 500);
  prettyCanv(canvPotential_1_c);

  TCanvas* canvPotential_2_c = new TCanvas("canvPotential_2_c", "canvPotential_2_c", 500, 500);
  prettyCanv(canvPotential_2_c);

  TGraph* graphPotential_1_p = new TGraph();
  TGraph* graphPotential_2_p = new TGraph();
  
  graphPotential_1_p->SetMarkerSize(.4);
  graphPotential_1_p->SetMarkerStyle(20);
  graphPotential_1_p->SetMarkerColor(col.getColor(0));

  graphPotential_2_p->SetMarkerSize(.4);
  graphPotential_2_p->SetMarkerStyle(20);
  graphPotential_2_p->SetMarkerColor(col.getColor(1));

  const int points = 100000;
  for(unsigned int i = 0; i < xVals.size(); ++i){
    if(i%points == 0) graphPotential_1_p->SetPoint(i/points, xVals.at(i), yVals_1.at(i));
    if(i%points == 0) graphPotential_2_p->SetPoint(i/points, xVals.at(i), yVals_2.at(i));
  }

  canvPotential_1_c->cd();
  graphPotential_1_p->Draw("A P");
  canvPotential_1_c->SaveAs("pdfDir/canvPotential_Approx.pdf");
  delete graphPotential_1_p;
  delete canvPotential_1_c;

  canvPotential_2_c->cd();
  graphPotential_2_p->Draw("A P");
  canvPotential_2_c->SaveAs("pdfDir/canvPotential_Approx_2.pdf");
  delete graphPotential_2_p;
  delete canvPotential_2_c;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += stochPotential();
  return retVal;
}
