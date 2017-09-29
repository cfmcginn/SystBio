#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"

#include "include/getLinBins.h"
#include "include/kirchnerPalette.h"

std::string prettyString(const double inVal, const int prec, const bool doDot)
{
  std::string retStr = std::to_string(inVal);
  while(retStr.find(".") < retStr.size()-1-prec){retStr.replace(retStr.size()-1, 1,"");}
  if(doDot) retStr.replace(retStr.find("."), 1, "p");
  return retStr;
}


void prettyCanv(TCanvas* canv_p)
{
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(1.3*canv_p->GetLeftMargin());
  canv_p->SetBottomMargin(canv_p->GetLeftMargin());

  return;
}


void prettyTH1(TH1* hist_p, const double size, const int style, const int col)
{
  hist_p->SetMarkerSize(size);
  hist_p->SetMarkerStyle(style);
  hist_p->SetMarkerColor(col);
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();

  return;
}


void reMaxMinTH1(TH1* hist_p)
{
  double max = hist_p->GetMaximum();
  double min = hist_p->GetMinimum();
  double interval = max-min;

  max += interval/10.;
  min -= interval/10.;
  if(min < 0) min = 0;

  hist_p->SetMaximum(max);
  hist_p->SetMinimum(min);

  return;
}


int numericalSolution(const double x0 = 1., const double y0 = 0.)
{
  kirchnerPalette col;
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);

  TFile* outFile_p = new TFile("output/numericalSolution.root", "RECREATE");

  const int nGammaX = 28;
  const double gammaX[nGammaX] = {3, 25./7. + .00001, 25./7.+.01, 25./7.+.1, 3.7, 3.8, 3.85, 3.9, 3.95, 4, 5, 6, 8, 10, 15, 20, 50, 100, 150, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000};
  double gammaXBins[nGammaX+1];
  gammaXBins[0] = 2.5;
  gammaXBins[nGammaX] = 2000000;

  for(int i = 0; i < nGammaX-1; ++i){
    gammaXBins[i+1] = (gammaX[i+1] + gammaX[i])/2.;
  }

  for(int i = 0; i < nGammaX+1; ++i){
    std::cout << gammaXBins[i] << ", ";
  }
  std::cout << std::endl;

  TCanvas* periodPerGamma_c = new TCanvas("periodPerGamma_c", "periodPerGamma_c", 500, 500);
  TCanvas* peakToPeakX_c = new TCanvas("peakToPeakX_c", "peakToPeakX_c", 500, 500);
  TCanvas* peakToPeakY_c = new TCanvas("peakToPeakY_c", "peakToPeakY_c", 500, 500);

  prettyCanv(periodPerGamma_c);
  prettyCanv(peakToPeakX_c);
  prettyCanv(peakToPeakY_c);

  TH1F* periodPerGamma_h = new TH1F("periodPerGamma_h", ";#gamma_{X};Period (Normalized to #gamma_{Y})", nGammaX, gammaXBins);
  TH1F* peakToPeakX_h = new TH1F("peakToPeakX_h", ";#gamma_{X};Peak-to-Peak Height (X)", nGammaX, gammaXBins);
  TH1F* peakToPeakY_h = new TH1F("peakToPeakY_h", ";#gamma_{X};Peak-to-Peak Height (Y)", nGammaX, gammaXBins);

  prettyTH1(periodPerGamma_h, .5, 20, col.getColor(0));
  prettyTH1(peakToPeakX_h, .5, 20, col.getColor(0));
  prettyTH1(peakToPeakY_h, .5, 20, col.getColor(0));
  
  std::cout << "A: " << periodPerGamma_h->GetMarkerSize() << std::endl;

  const double vx = .1;
  const double kx = 4.;
  const double ky = 2.;
  
  const unsigned long long nTimeSteps = 10000000;
  const unsigned long long nTimeBins = 10000;
  const float timeLow = 0.;
  const float timeHigh = 100.0;
  double timeBins[nTimeBins+1];
  getLinBins(timeLow, timeHigh, nTimeBins, timeBins);
  const double timeWidth = (timeHigh - timeLow)/double(nTimeSteps);

  for(int gI = 0; gI < nGammaX; ++gI){    
    const std::string canvName = "xy_Gamma" + prettyString(gammaX[gI], 4, true) + "_x" + prettyString(x0, 1, true) + "_y" + prettyString(y0, 1, true) + "_THigh" + std::to_string(int(timeHigh)) + "_c";
    

    TCanvas* canv_p = new TCanvas(canvName.c_str(), canvName.c_str(), 500, 500);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetBottomMargin(canv_p->GetLeftMargin());
    prettyCanv(canv_p);
    
    gStyle->SetOptStat(0);
    
    std::string titleString = ";Time [s/#gamma_{Y}];#color[" + std::to_string(col.getColor(0)) + "]{x(t)/A_{3},1 at t=0}, #color[" + std::to_string(col.getColor(2)) + "]{y(t)/A_{1}, 0 at t=0}";
    TH1F* dummyHist_p = new TH1F("dummyHist_h", titleString.c_str(), nTimeBins, timeBins);
    prettyTH1(dummyHist_p, 1, 20, col.getColor(0));

    TH1F* xHist_p = new TH1F("xHist_h", ";Time [s/#gamma_{Y}];x(t)/A_{3}", nTimeBins, timeBins);
    TH1F* yHist_p = new TH1F("yHist_h", ";Time [s/#gamma_{Y}];y(t)/A_{1}", nTimeBins, timeBins);

    prettyTH1(xHist_p, .2, 20, col.getColor(0));
    prettyTH1(yHist_p, .2, 20, col.getColor(2));
    
    double xVal = x0;
    double yVal = y0;
    
    xHist_p->SetBinContent(1, x0);
    yHist_p->SetBinContent(1, y0);
    
    xHist_p->SetBinError(1, 0.);
    yHist_p->SetBinError(1, 0.);
    
    bool isRisingX = false;
    bool isRisingY = false;
    
    std::vector<double> peakValTimesX;
    std::vector<double> peakValTimesY;
    
    double maximumX = -1;
    double minimumX = 100000.;
    
    double maximumY = -1;
    double minimumY = 100000.;
    
    for(unsigned long long iter = 1; iter < nTimeSteps; ++iter){
      double newXVal = xVal + timeWidth*gammaX[gI]*(vx + kx*xVal*xVal/((1+yVal)*(1+xVal*xVal)) - xVal);
      double newYVal = yVal + timeWidth*(ky*xVal- yVal);
      
      if(newXVal > xVal) isRisingX = true;
      if(newXVal < xVal && isRisingX){
	isRisingX = false;
	peakValTimesX.push_back(timeWidth*double(iter));
      }
      
      if(newXVal > xVal) isRisingY = true;
      if(newXVal < xVal && isRisingY){
	isRisingY = false;
	peakValTimesY.push_back(timeWidth*double(iter));
      }
      
      xVal = newXVal;
      yVal = newYVal;
      
      if(timeWidth*double(iter) > 5){
	if(xVal > maximumX) maximumX = xVal;
	if(yVal > maximumY) maximumY = yVal;
	if(xVal < minimumX) minimumX = xVal;
	if(yVal < minimumY) minimumY = yVal;
      }

      if(iter%(nTimeSteps/nTimeBins) != 0) continue;
      
      int bin = iter/(nTimeSteps/nTimeBins);
      
      xHist_p->SetBinContent(bin+1, newXVal);
      xHist_p->SetBinError(bin+1, 0);
      
      yHist_p->SetBinContent(bin+1, newYVal);
      yHist_p->SetBinError(bin+1, 0);
    }
    
    dummyHist_p->SetMaximum(2.8);
    dummyHist_p->SetMinimum(0.);
    
    std::cout << "DeltaPeak (period): " << peakValTimesX.at(peakValTimesX.size()-1) - peakValTimesX.at(peakValTimesX.size()-2) << std::endl;
    std::cout << "DeltaPeak (period): " << peakValTimesY.at(peakValTimesY.size()-1) - peakValTimesY.at(peakValTimesY.size()-2) << std::endl;

    if(gammaX[gI] > 25./7.){
      periodPerGamma_h->SetBinContent(gI+1,  peakValTimesX.at(peakValTimesX.size()-1) -  peakValTimesX.at(peakValTimesX.size()-2));
      periodPerGamma_h->SetBinError(gI+1, 0);
      
      peakToPeakX_h->SetBinContent(gI+1, maximumX - minimumX);
      peakToPeakX_h->SetBinError(gI+1, 0);
      
      peakToPeakY_h->SetBinContent(gI+1, maximumY - minimumY);
      peakToPeakY_h->SetBinError(gI+1, 0);
    }

    outFile_p->cd();
    
    canv_p->cd();
    dummyHist_p->DrawCopy();
    xHist_p->DrawCopy("P SAME");
    yHist_p->DrawCopy("P SAME");
    
    std::string xLabel = "X_{Min}=" + prettyString(minimumX,4, false) + ";X_{Max}=" + prettyString(maximumX,4, false);
    std::string yLabel = "Y_{Min}=" + prettyString(minimumY,4, false) + ";Y_{Max}=" + prettyString(maximumY,4, false);
    std::string period = "Period [s/#gamma_{Y}] = " + prettyString(peakValTimesX.at(peakValTimesX.size()-1) - peakValTimesX.at(peakValTimesX.size()-2), 4, false);
    std::string paramStr = "k_{X}=4. k_{Y}=2. v_{X}=.1 v_{Y}=0. #gamma_{X}=" + prettyString(gammaX[gI], 4, false);

    label_p->DrawLatex(2+timeLow, 2.7, xLabel.c_str());
    label_p->DrawLatex(2+timeLow, 2.6, yLabel.c_str());
    label_p->DrawLatex(timeHigh*5./11., 2.7, period.c_str());
    label_p->DrawLatex(timeHigh*5./11., 2.6, paramStr.c_str());

    canv_p->Write("", TObject::kOverwrite);
    std::string saveName = canv_p->GetName();
    saveName = "pdfDir/" +saveName + ".pdf";
    canv_p->SaveAs(saveName.c_str());
    delete canv_p;
    
    xHist_p->Write("", TObject::kOverwrite);
    yHist_p->Write("", TObject::kOverwrite);
    
    delete dummyHist_p;
    delete xHist_p;
    delete yHist_p;
  }

  outFile_p->cd();

  std::cout << "B: " << periodPerGamma_h->GetMarkerSize() << std::endl;

  periodPerGamma_c->cd();
  periodPerGamma_h->Write("", TObject::kOverwrite);
  reMaxMinTH1(periodPerGamma_h);
  periodPerGamma_h->DrawCopy("P");
  gPad->SetLogx();
  periodPerGamma_c->Write("", TObject::kOverwrite);
  std::string saveName = periodPerGamma_c->GetName();
  saveName = "pdfDir/" + saveName + ".pdf";
  periodPerGamma_c->SaveAs(saveName.c_str());

  peakToPeakX_c->cd();
  peakToPeakX_h->Write("", TObject::kOverwrite);
  reMaxMinTH1(peakToPeakX_h);
  peakToPeakX_h->DrawCopy("P");
  gPad->SetLogx();
  peakToPeakX_c->Write("", TObject::kOverwrite);
  saveName = peakToPeakX_c->GetName();
  saveName = "pdfDir/" + saveName + ".pdf";
  peakToPeakX_c->SaveAs(saveName.c_str());

  peakToPeakY_c->cd();
  peakToPeakY_h->Write("", TObject::kOverwrite);
  reMaxMinTH1(peakToPeakY_h);
  peakToPeakY_h->DrawCopy("P");
  gPad->SetLogx();
  peakToPeakY_c->Write("", TObject::kOverwrite);
  saveName = peakToPeakY_c->GetName();
  saveName = "pdfDir/" + saveName + ".pdf";
  peakToPeakY_c->SaveAs(saveName.c_str());

  delete periodPerGamma_c;
  delete peakToPeakX_c;
  delete peakToPeakY_c;

  delete periodPerGamma_h;
  delete peakToPeakX_h;
  delete peakToPeakY_h;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 2 && argc != 3){
    std::cout << "Usage ./numericalSolution.exe <optX0> <optY0>" << std::endl;
  }

  int retVal = 0;
  if(argc == 1) retVal += numericalSolution();
  else if(argc == 2) retVal += numericalSolution(std::stof(argv[1]));
  else if(argc == 3) retVal += numericalSolution(std::stof(argv[1]), std::stof(argv[2]));
  return retVal;
}
