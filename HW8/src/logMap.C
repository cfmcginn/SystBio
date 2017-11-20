#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"

#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"

int logMap()
{
  kirchnerPalette col;

  const int nRParam = 3;
  const double rParams[nRParam] = {2., 3.1, 3.5};
  const double x0 = 0.2;
  const int nTime = 20;

  TGraph* cobWebs_p[nRParam];
  TH1F* dummys_p[nRParam];

  for(int i = 0; i < nRParam; ++i){
    cobWebs_p[i] = new TGraph();
    cobWebs_p[i]->SetMarkerColor(col.getColor(i));
    cobWebs_p[i]->SetMarkerStyle(20);
    cobWebs_p[i]->SetLineColor(col.getColor(i));
    cobWebs_p[i]->SetMarkerSize(.6);

    dummys_p[i] = new TH1F(("dummys_" + std::to_string(int(rParams[i]*10)) + "_h").c_str(), ";;", 10, 0, 1);

    dummys_p[i]->SetMarkerColor(col.getColor(i));
    dummys_p[i]->SetMarkerStyle(20);
    dummys_p[i]->SetLineColor(col.getColor(i));
    dummys_p[i]->SetMarkerSize(.6);
  }

  double prevXVal[nRParam] = {x0,x0,x0};

  double maxXVal = x0;
  double maxYVal = 0;

  double minXVal = x0;
  double minYVal = 0;

  for(int j = 0; j < nRParam; ++j){
    cobWebs_p[j]->SetPoint(0, prevXVal[j], 0);
  }

  for(int i = 0; i < nTime; ++i){
    for(int j = 0; j < nRParam; ++j){
      double yVal = rParams[j]*prevXVal[j]*(1-prevXVal[j]);
      cobWebs_p[j]->SetPoint(i*2+1, prevXVal[j], yVal);
      cobWebs_p[j]->SetPoint(i*2+2, yVal, yVal);

      if(minXVal > prevXVal[j]) minXVal = prevXVal[j];
      if(minYVal > yVal) minYVal = yVal;

      if(maxXVal < prevXVal[j]) maxXVal = prevXVal[j];
      if(maxYVal < yVal) maxYVal = yVal;

      prevXVal[j] = yVal;
    }
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

  TLegend* leg_p = new TLegend(0.2, 0.65, 0.5, 0.95);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(18);

  TLine* line_p = new TLine();

  for(int i = 0; i < nRParam; ++i){
    cobWebs_p[i]->Draw("P");
    line_p->SetLineColor(col.getColor(i));

    for(int j = 1; j < cobWebs_p[i]->GetN(); ++j){
      double x1,y1,x2,y2;
      
      cobWebs_p[i]->GetPoint(j-1, x1, y1);
      cobWebs_p[i]->GetPoint(j, x2, y2);

      line_p->DrawLine(x1,y1,x2,y2);
    }

    leg_p->AddEntry(dummys_p[i], ("r=" + prettyString(rParams[i], 1, false)).c_str(), "P L");
  }

  leg_p->Draw("SAME");

  delete line_p;

  canv_p->SaveAs("pdfDir/cobWeb.pdf");

  delete leg_p;
  delete canv_p;

  for(int i = 0; i < nRParam; ++i){
    delete cobWebs_p[i];
  }

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += logMap();
  return retVal;
}
