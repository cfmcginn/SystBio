#include <iostream>

#include "TF1.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

#include "include/getLinBins.h"
#include "include/kirchnerPalette.h"

int campSim()
{
  kirchnerPalette col;

  const double beta = 4;
  const double kappa = 0.5;
  const double bigD = .0000001;
  const double sigma = 0.1;
  const double i0 = 0.1;
  const double gausHeight = i0/4.;

  const int spacePoints = 2000.;
  const int timeSteps = 200000.;

  const double xStart = -1.;
  const double xEnd = 1.;
  double points[spacePoints+1];
  double cPoints[spacePoints];
  double iPoints[spacePoints];
  getLinBins(xStart, xEnd, spacePoints, points);

  for(int i = 0; i < spacePoints; ++i){
    cPoints[i] = 0;
    iPoints[i] = i0;
  }

  const double deltaTimeStep = 0.001;
  const double deltaSpaceStep = points[1] - points[0];
  /*
  TF1* gaus_p = new TF1("gaus_p", "gaus", -1., 1.);
  gaus_p->SetParameter(0, gausHeight);
  gaus_p->SetParameter(1, 0);
  gaus_p->SetParameter(2, sigma);
  */
  TF1* gaus_p = new TF1("gaus_p", ".1/(TMath::Power((x/.1), 4) + 1)", -1., 1.);

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(20);
  label_p->SetNDC();

  for(int i = 0; i < timeSteps; ++i){
    if(i != 0){
      for(int j = 0; j < spacePoints; ++j){
	double didt = kappa*(cPoints[j] - iPoints[j]);
	if(j == 0 || j == spacePoints-1) didt = 0;
	iPoints[j] += didt*deltaTimeStep;

	double dcdt = cPoints[j]*cPoints[j]/(beta*iPoints[j]*iPoints[j] + cPoints[j]*cPoints[j])  - cPoints[j];
	double cLow = 0.;
	if(j != 0) cLow = cPoints[j-1];
	double cHi = 0.;
	if(j != spacePoints-1) cHi = cPoints[j+1];

	dcdt += bigD*(cLow + cHi - 2*cPoints[j])/(deltaSpaceStep*deltaSpaceStep);

	if(j == 0 || j == spacePoints-1) dcdt = 0;

	cPoints[j] += dcdt*deltaTimeStep;
      }
    }
    else{
      for(int j = 0; j < spacePoints; ++j){
	cPoints[j] = gaus_p->Eval((points[j] + points[j+1])/2.);
      }
    }

    if(i%1000 == 0){
      TCanvas* canv_p = new TCanvas("canv", "canv", 500, 500);
      TH1F* cHist_p = new TH1F("cHist_h", ";x;C or I", spacePoints, points);
      TH1F* iHist_p = new TH1F("iHist_h", ";x;C or I", spacePoints, points);
      
      cHist_p->SetMarkerStyle(20);
      cHist_p->SetMarkerSize(0.4);
      cHist_p->SetMarkerColor(col.getColor(0));

      iHist_p->SetMarkerStyle(20);
      iHist_p->SetMarkerSize(0.4);
      iHist_p->SetMarkerColor(col.getColor(2));

      for(Int_t j = 0; j < spacePoints; ++j){
	cHist_p->SetBinContent(j+1, cPoints[j]);
	cHist_p->SetBinError(j+1, 0);
	
	iHist_p->SetBinContent(j+1, iPoints[j]);
	iHist_p->SetBinError(j+1, 0);
      }
      
      canv_p->cd();
      gStyle->SetOptStat(0);
      cHist_p->SetMaximum(0.6);
      cHist_p->SetMinimum(0.0);
      cHist_p->GetXaxis()->CenterTitle();
      cHist_p->GetYaxis()->CenterTitle();
      cHist_p->DrawCopy("P");
      iHist_p->DrawCopy("P SAME");

      label_p->DrawLatex(.7, .95, ("t=" + std::to_string(int(deltaTimeStep*i)) + " (s)").c_str());

      TLegend* leg_p = new TLegend(0.3,.92, .6, .99);
      leg_p->SetFillColor(0);
      leg_p->SetBorderSize(0);
      leg_p->AddEntry(iHist_p, "I", "P");
      leg_p->AddEntry(cHist_p, "C", "P");

      leg_p->Draw("SAME");

      canv_p->SaveAs(("gifDir1/simul_" + std::to_string(i) + ".gif").c_str());
      if(i==60000)canv_p->SaveAs(("pdfDir/simul_1_" + std::to_string(i) + ".pdf").c_str());
      delete cHist_p;
      delete iHist_p;
      delete leg_p;
      
      delete canv_p;
    }
  }

  delete label_p;


  return 0;
}

int main()
{
  int retVal = 0;
  retVal += campSim();
  return retVal;
}
