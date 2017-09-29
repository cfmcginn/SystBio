#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TMath.h"
#include "TLatex.h"

void prettyCanv(TCanvas* canv_p)
{
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(canv_p->GetLeftMargin()*1.3);
  canv_p->SetBottomMargin(canv_p->GetLeftMargin());
  canv_p->SetTopMargin(canv_p->GetLeftMargin()/2.);

  return;
}


int emulateCircadian()
{
  gStyle->SetOptStat(0);

  const double vx = 0.1;
  //  const double vy = 0.0;

  const double kx = 4.;
  const double ky = 2.;

  const double gammax = 10;

  TFile* outFile_p = new TFile("output/outFile.root", "RECREATE");

  TCanvas* xDotNull1Canv_p = new TCanvas("xDotNull1Canv_p", "xDotNull1Canv_p", 500, 500);
  TCanvas* xDotNull2Canv_p = new TCanvas("xDotNull2Canv_p", "xDotNull2Canv_p", 500, 500);
  TCanvas* yDotNullCanv_p = new TCanvas("yDotNullCanv_p", "yDotNullCanv_p", 500, 500);

  TCanvas* canv1_p = new TCanvas("canv1_c", "canv1_c", 500, 500);
  TCanvas* canv2_p = new TCanvas("canv2_c", "canv2_c", 500, 500);

  prettyCanv(xDotNull1Canv_p);
  prettyCanv(xDotNull2Canv_p);
  prettyCanv(yDotNullCanv_p);
  prettyCanv(canv1_p);
  prettyCanv(canv2_p);

  const std::string xDotStr1 = std::to_string(gammax) + "*(" + std::to_string(vx) + " + " + std::to_string(kx) + "*1/(1+y) - x)";
  const std::string xDotStr2 = std::to_string(gammax) + "*(" + std::to_string(vx) + " + " + std::to_string(kx) + "*(1/(1+y))*(x*x/(1+x*x)) - x)";
  const std::string yDotStr = std::to_string(ky) + "*x - y";

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();

  TF2* xDot1_p = new TF2("xDot1_p", xDotStr1.c_str(), 0, 100000, 0, 100000);
  TF2* xDot2_p = new TF2("xDot2_p", xDotStr2.c_str(), 0, 100000, 0, 100000);
  TF2* yDot_p = new TF2("yDot_p", yDotStr.c_str(), 0, 100000, 0, 100000);

  const std::string xDotNullStr1 = std::to_string(kx) + "/(x - " + std::to_string(vx) + ") - 1";
  TF1* xDotNull1_p = new TF1("xDotNull1_p", xDotNullStr1.c_str(), 0, 10000);
  xDotNull1_p->SetMarkerColor(kRed);
  xDotNull1_p->SetLineColor(kRed);
  xDotNull1_p->SetMarkerSize(1);

  const std::string xDotNullStr2 = std::to_string(kx) + "*x*x/((1 + x*x)*(x - " + std::to_string(vx)  + ")) - 1";
  TF1* xDotNull2_p = new TF1("xDotNull2_p", xDotNullStr2.c_str(), 0, 10000);
  xDotNull2_p->SetMarkerColor(kRed);
  xDotNull2_p->SetLineColor(kRed);
  xDotNull2_p->SetMarkerSize(1);
  xDotNull2_p->SetLineStyle(1);

  const std::string yDotNullStr = std::to_string(ky) + "*x";
  TF1* yDotNull_p = new TF1("yDotNull_p", yDotNullStr.c_str(), 0, 10000);
  yDotNull_p->SetMarkerColor(kBlue);
  yDotNull_p->SetLineColor(kBlue);
  yDotNull_p->SetMarkerSize(1);

  TH1F* dummyHist_p = new TH1F("dummyHist_h", ";X;Y", 3, 0, 3.);
  dummyHist_p->SetMaximum(3.);
  dummyHist_p->SetMinimum(0);
  dummyHist_p->GetXaxis()->CenterTitle();
  dummyHist_p->GetYaxis()->CenterTitle();

  TH1F* dummyHist2_p = new TH1F("dummyHist2_h", ";X;Y", 2, 0, 2);
  dummyHist2_p->SetMaximum(2);
  dummyHist2_p->SetMinimum(0);
  dummyHist2_p->GetXaxis()->CenterTitle();
  dummyHist2_p->GetYaxis()->CenterTitle();

  xDotNull1Canv_p->cd();
  dummyHist_p->DrawCopy();
  label_p->DrawLatex(.35, .96, "dx/dt = 0; (A) Scenario");
  xDotNull1_p->DrawCopy("SAME");

  TArrow* arrow = new TArrow();

  for(int i = 0; i < 3; ++i){
    double xPos = 0.5 + i;

    for(int j = 0; j < 3; ++j){
      double yPos = 0.5 + j;

      std::string opt = ">";
      if(xDot1_p->Eval(xPos,yPos) < 0) opt = "<";

      arrow->DrawArrow(xPos - .2, yPos, xPos + .2, yPos, 0, opt.c_str());
    }
  }

  xDotNull2Canv_p->cd();
  dummyHist2_p->DrawCopy();
  label_p->DrawLatex(.35, .96, "dx/dt = 0; (B) Scenario");
  xDotNull2_p->DrawCopy("SAME");

  for(int i = 0; i < 4; ++i){
    double xPos = 0.5 + i;

    for(int j = 0; j < 4; ++j){
      double yPos = 0.5 + j;

      std::string opt = ">";
      if(xDot2_p->Eval(xPos,yPos) < 0) opt = "<";

      arrow->DrawArrow(xPos - .2, yPos, xPos + .2, yPos, 0, opt.c_str());
    }
  }

  yDotNullCanv_p->cd();
  dummyHist_p->DrawCopy();
  label_p->DrawLatex(.35, .96, "dy/dt = 0; (A) and (B) Scenario");
  yDotNull_p->DrawCopy("SAME");

  for(int i = 0; i < 3; ++i){
    double xPos = 0.5 + i;

    for(int j = 0; j < 3; ++j){
      double yPos = 0.5 + j;

      std::string opt = ">";
      if(yDot_p->Eval(xPos,yPos) < 0) opt = "<";

      arrow->DrawArrow(xPos, yPos - .2, xPos , yPos + .2, 0, opt.c_str());
    }
  }

  canv1_p->cd();
  dummyHist_p->DrawCopy();
  xDotNull1_p->DrawCopy("SAME");
  //  xDotNull2_p->DrawCopy("SAME");
  yDotNull_p->DrawCopy("SAME");
  label_p->DrawLatex(.35, .96, "Arrows show dx/dt, dy/dt vector, (A) Scenario");

  for(int i = 0; i < 29; ++i){
    double xPos = 0.1 + i*.1;

    for(int j = 0; j < 29; ++j){
      double yPos = 0.1 + j*.1;
      
      double xEval = xDot1_p->Eval(xPos,yPos);
      double yEval = yDot_p->Eval(xPos,yPos);

      double xMin = xPos - .025*TMath::Abs(xEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);
      double xMax = xPos + .025*TMath::Abs(xEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);

      double yMin = yPos - .025*TMath::Abs(yEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);
      double yMax = yPos + .025*TMath::Abs(yEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);
      
      if(yEval*xEval < 0){
	double tempY = yMin;
	yMin = yMax;
	yMax = tempY;
      }

      if(xPos < .5 && yPos > 3){
	std::cout << "xPos, yPos, xEval, yEval: " << xPos << ", " << yPos << ", " << xEval << ", " << yEval << std::endl;
      }

      std::string opt = ">";
      if(xEval < 0) opt = "<";

      arrow->DrawArrow(xMin, yMin, xMax, yMax, 0.01, opt.c_str());
    }
  }

  canv2_p->cd();
  dummyHist2_p->DrawCopy();
  //  xDotNull1_p->DrawCopy("SAME");
  xDotNull2_p->DrawCopy("SAME");
  yDotNull_p->DrawCopy("SAME");
  label_p->DrawLatex(.35, .96, "Arrows show dx/dt, dy/dt vector, (B) Scenario");

  for(int i = 0; i < 24; ++i){
    double xPos = 0.08 + i*.08;

    for(int j = 0; j < 24; ++j){
      double yPos = 0.08 + j*.08;
      
      double xEval = xDot2_p->Eval(xPos,yPos);
      double yEval = yDot_p->Eval(xPos,yPos);

      double xMin = xPos - .02*TMath::Abs(xEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);
      double xMax = xPos + .02*TMath::Abs(xEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);

      double yMin = yPos - .02*TMath::Abs(yEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);
      double yMax = yPos + .02*TMath::Abs(yEval)/TMath::Sqrt(xEval*xEval + yEval*yEval);
      
      if(yEval*xEval < 0){
	double tempY = yMin;
	yMin = yMax;
	yMax = tempY;
      }

      if(xPos < .5 && yPos > 3){
	std::cout << "xPos, yPos, xEval, yEval: " << xPos << ", " << yPos << ", " << xEval << ", " << yEval << std::endl;
      }

      std::string opt = ">";
      if(xEval < 0) opt = "<";

      arrow->DrawArrow(xMin, yMin, xMax, yMax, 0.01, opt.c_str());
    }
  }


  outFile_p->cd();

  xDotNull1Canv_p->Write("", TObject::kOverwrite);
  xDotNull2Canv_p->Write("", TObject::kOverwrite);
  yDotNullCanv_p->Write("", TObject::kOverwrite);

  canv1_p->Write("", TObject::kOverwrite);
  canv2_p->Write("", TObject::kOverwrite);


  xDotNull1Canv_p->SaveAs("pdfDir/xDotNull1Canv.pdf");
  xDotNull2Canv_p->SaveAs("pdfDir/xDotNull2Canv.pdf");
  yDotNullCanv_p->SaveAs("pdfDir/yDotNullCanv.pdf");
  canv1_p->SaveAs("pdfDir/canv1.pdf");
  canv2_p->SaveAs("pdfDir/canv2.pdf");

  delete xDotNull1Canv_p;
  delete xDotNull2Canv_p;
  delete yDotNullCanv_p;

  delete canv1_p;
  delete canv2_p;

  delete xDotNull1_p;
  delete xDotNull2_p;
  delete yDotNull_p;
  
  delete dummyHist_p;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal = emulateCircadian();
  return retVal;
}
