#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLine.h"

#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"

int logMap()
{
  kirchnerPalette col;

  const double rParam = 4;
  const int nX0 = 2;
  const double x0s[nX0] = {0.2, 0.2000000001};
  const int nTime = 80;

  TGraph* cobWebs_p[nX0];
  TGraph* traj_p[nX0];
  TGraph* deltaTraj_p = new TGraph();
  for(int i = 0; i < nX0; ++i){
    cobWebs_p[i] = new TGraph();
    cobWebs_p[i]->SetMarkerColor(col.getColor(i));
    cobWebs_p[i]->SetMarkerStyle(20);
    cobWebs_p[i]->SetLineColor(col.getColor(i));
    cobWebs_p[i]->SetMarkerSize(.6);

    traj_p[i] = new TGraph();
    traj_p[i]->SetMarkerColor(col.getColor(i*2));
    traj_p[i]->SetMarkerStyle(20);
    traj_p[i]->SetLineColor(col.getColor(i*2));
    traj_p[i]->SetMarkerSize(.6);
  }

  double prevXVal[nX0] = {x0s[0],x0s[1]};

  double maxXVal = x0s[1];
  double maxYVal = 0;

  double minXVal = x0s[0];
  double minYVal = 0;

  for(int j = 0; j < nX0; ++j){
    cobWebs_p[j]->SetPoint(0, prevXVal[j], 0);
    traj_p[j]->SetPoint(0, 0, x0s[j]);
  }

  deltaTraj_p->SetPoint(0, 0, TMath::Abs(x0s[0] - x0s[1]));
  deltaTraj_p->SetMarkerColor(col.getColor(0));
  deltaTraj_p->SetMarkerStyle(20);
  deltaTraj_p->SetLineColor(col.getColor(0));
  deltaTraj_p->SetMarkerSize(.6);


  double maxTraj = 0;

  for(int i = 0; i < nTime; ++i){
    double delta = 0;

    for(int j = 0; j < nX0; ++j){
      double yVal = rParam*prevXVal[j]*(1-prevXVal[j]);
      cobWebs_p[j]->SetPoint(i*2+1, prevXVal[j], yVal);
      cobWebs_p[j]->SetPoint(i*2+2, yVal, yVal);
      traj_p[j]->SetPoint(i+1, i+1, yVal);

      if(j == 0) delta = yVal;
      else delta -= yVal;

      if(minXVal > prevXVal[j]) minXVal = prevXVal[j];
      if(minYVal > yVal) minYVal = yVal;

      if(maxXVal < prevXVal[j]) maxXVal = prevXVal[j];
      if(maxYVal < yVal) maxYVal = yVal;

      prevXVal[j] = yVal;
    }
    
    delta = TMath::Abs(delta);
    if(delta > maxTraj) maxTraj = delta;
    deltaTraj_p->SetPoint(i+1, i+1, delta);

    TCanvas* canv_p = new TCanvas("cobWeb_c", "cobWeb_c", 500, 500);
    prettyCanv(canv_p);
    TH1F* dummy_p = new TH1F("dummyHist_h", ";x;y", 10, TMath::Min(0., minXVal), maxXVal*1.1);
    dummy_p->SetMaximum(maxYVal*1.1);
    dummy_p->SetMinimum(minYVal);
    dummy_p->GetXaxis()->CenterTitle();
    dummy_p->GetYaxis()->CenterTitle();
    dummy_p->DrawCopy();
    gStyle->SetOptStat(0);

    delete dummy_p;

    TLine* line_p = new TLine();
    
    for(int k = 0; k < nX0; ++k){
      cobWebs_p[k]->Draw("P");
      line_p->SetLineColor(col.getColor(k));
      
      for(int j = 1; j < cobWebs_p[k]->GetN(); ++j){
	double x1,y1,x2,y2;
	
	cobWebs_p[k]->GetPoint(j-1, x1, y1);
	cobWebs_p[k]->GetPoint(j, x2, y2);
	
	line_p->DrawLine(x1,y1,x2,y2);
      }
    }
    delete line_p;
    
    canv_p->SaveAs(("pdfDir/cobWeb2_" + std::to_string(i) + ".gif").c_str());
    
    delete canv_p;


    canv_p = new TCanvas("traj_c", "traj_c", 1000, 500);
    prettyCanv(canv_p);
    dummy_p = new TH1F("dummyHist_h", ";Time;Traj", 83, -0.5, 82.5);
    dummy_p->SetMaximum(1.1);
    dummy_p->SetMinimum(0.);
    dummy_p->GetXaxis()->CenterTitle();
    dummy_p->GetYaxis()->CenterTitle();
    dummy_p->DrawCopy();
    gStyle->SetOptStat(0);

    delete dummy_p;

    line_p = new TLine();
    
    for(int k = 0; k < nX0; ++k){
      traj_p[k]->Draw("P");
      line_p->SetLineColor(col.getColor(k*2));
      
      for(int j = 1; j < traj_p[k]->GetN(); ++j){
	double x1,y1,x2,y2;
	
	traj_p[k]->GetPoint(j-1, x1, y1);
	traj_p[k]->GetPoint(j, x2, y2);
	
	line_p->DrawLine(x1,y1,x2,y2);
      }
    }
    delete line_p;
    
    canv_p->SaveAs(("pdfDir/traj_" + std::to_string(i) + ".gif").c_str());
    
    delete canv_p;
  }

  TCanvas* canv_p = new TCanvas("cobWeb_c", "cobWeb_c", 500, 500);
  prettyCanv(canv_p);
  TH1F* dummy_p = new TH1F("dummyHist_h", ";x;y", 10, TMath::Min(0., minXVal), maxXVal*1.1);
  dummy_p->SetMaximum(maxYVal*1.1);
  dummy_p->SetMinimum(minYVal);
  dummy_p->GetXaxis()->CenterTitle();
  dummy_p->GetYaxis()->CenterTitle();
  dummy_p->DrawCopy();
  gStyle->SetOptStat(0);

  delete dummy_p;

  TLine* line_p = new TLine();

  for(int i = 0; i < nX0; ++i){
    cobWebs_p[i]->Draw("P");
    line_p->SetLineColor(col.getColor(i));

    for(int j = 1; j < cobWebs_p[i]->GetN(); ++j){
      double x1,y1,x2,y2;
      
      cobWebs_p[i]->GetPoint(j-1, x1, y1);
      cobWebs_p[i]->GetPoint(j, x2, y2);

      line_p->DrawLine(x1,y1,x2,y2);
    }
  }
  delete line_p;

  canv_p->SaveAs("pdfDir/cobWeb2.pdf");

  delete canv_p;

  canv_p = new TCanvas("deltaTraj_c", "deltaTraj_c", 500, 500);  
  prettyCanv(canv_p);
  dummy_p = new TH1F("dummyHist_h", ";Time Steps;#DeltaTrajectory", 83, -0.5, 82.5);
  dummy_p->SetMaximum(maxTraj*5.);
  dummy_p->SetMinimum((x0s[1]-x0s[0])/5.); 
  dummy_p->GetXaxis()->CenterTitle();
  dummy_p->GetYaxis()->CenterTitle();
  dummy_p->DrawCopy();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  delete dummy_p;
  
  deltaTraj_p->Draw("P");
  canv_p->SaveAs("pdfDir/deltaTraj.pdf");

  delete canv_p;

  for(int i = 0; i < nX0; ++i){
    delete cobWebs_p[i];
  }

  delete deltaTraj_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += logMap();
  return retVal;
}
