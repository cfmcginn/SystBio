#include <iostream>
#include <vector>

#include "TMath.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include "include/plotUtilities.h"
#include "include/kirchnerPalette.h"

int stochProteins()
{
  kirchnerPalette col;

  TLegend* leg_p = new TLegend(0.2, 0.5, 0.5, 0.9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(18);

  TLegend* leg2_p = new TLegend(0.2, 0.7, 0.5, 0.9);
  leg2_p->SetBorderSize(0);
  leg2_p->SetFillStyle(0);
  leg2_p->SetTextFont(43);
  leg2_p->SetTextSize(18);

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  label_p->SetNDC();

  double aVal10 = 10;
  double aVal100 = 100;
  double gamma = 1;

  double kVal = 100;
  double bVal = (300.*300. - 100.*100.)/400.;

  double finalTime10 = 10000;
  double finalTime100 = 10000;
  double finalTime100K = 10000;

  int numberX0 = 3;
  double time = 0;

  std::vector<double> times10;
  std::vector<int> number10;
  std::vector<double> number10True;
  std::vector<double> times100;
  std::vector<int> number100;
  std::vector<double> number100True;
  std::vector<double> times100K;
  std::vector<int> number100K;

  times10.push_back(time);
  number10.push_back(numberX0);
  times100.push_back(time);
  number100.push_back(numberX0);
  times100K.push_back(time);
  number100K.push_back(numberX0);

  TRandom3* randGen_p = new TRandom3(0);

  while(times10.at(times10.size()-1) < finalTime10){
    double prodTime10 = randGen_p->Exp(1./aVal10);
    double degTime10 = randGen_p->Exp(1./(number10.at(number10.size()-1)*gamma));
    
    if(prodTime10 < degTime10){
      times10.push_back(times10.at(times10.size()-1) + prodTime10);
      number10.push_back(number10.at(number10.size()-1) + 1);

      number10True.push_back((10 - (10 - 3)*TMath::Exp(-times10.at(times10.size()-1))));
    }
    else{
      times10.push_back(times10.at(times10.size()-1) + degTime10);
      number10.push_back(number10.at(number10.size()-1) - 1);

      number10True.push_back((10 - (10 - 3)*TMath::Exp(-times10.at(times10.size()-1))));
    }
  }

  while(times100.at(times100.size()-1) < finalTime100){
    double prodTime100 = randGen_p->Exp(1./aVal100);
    double degTime100 = randGen_p->Exp(1./(number100.at(number100.size()-1)*gamma));
    
    if(prodTime100 < degTime100){
      times100.push_back(times100.at(times100.size()-1) + prodTime100);
      number100.push_back(number100.at(number100.size()-1) + 1);

      number100True.push_back((100 - (100 - 3)*TMath::Exp(-times100.at(times100.size()-1))));
    }
    else{
      times100.push_back(times100.at(times100.size()-1) + degTime100);
      number100.push_back(number100.at(number100.size()-1) - 1);

      number100True.push_back((100 - (100 - 3)*TMath::Exp(-times100.at(times100.size()-1))));
    }

  }

  while(times100K.at(times100K.size()-1) < finalTime100K){
    double rate = bVal*kVal/(kVal + number100K.at(number100K.size()-1));
    double prodTime100K = randGen_p->Exp(1./rate);
    double degTime100K = randGen_p->Exp(1./(number100K.at(number100K.size()-1)*gamma));
    
    if(prodTime100K < degTime100K){
      times100K.push_back(times100K.at(times100K.size()-1) + prodTime100K);
      number100K.push_back(number100K.at(number100K.size()-1) + 1);
    }
    else{
      times100K.push_back(times100K.at(times100K.size()-1) + degTime100K);
      number100K.push_back(number100K.at(number100K.size()-1) - 1);
    }
  }

  gStyle->SetOptStat(0);
  TCanvas* canv10_p = new TCanvas("canv10_c", "canv10_c", 500, 500);
  prettyCanv(canv10_p);
  TCanvas* canv10_2_p = new TCanvas("canv10_2_c", "canv10_2_c", 500, 500);
  prettyCanv(canv10_2_p);

  TGraph* graph10_p = new TGraph();
  TGraph* graph10True_p = new TGraph();
  TH1F* hist10_p = new TH1F("hist10_h", ";x;Counts (Constitutive)", 35, -0.5, 34.5);
  TH1F* dummy10_p = new TH1F("dummy10_h", ";Time (s);Counts (Constitutive)", 10, 0, 100);

  prettyTH1(hist10_p, 0.01, 20, 1);
  hist10_p->GetXaxis()->SetTitleFont(43);
  hist10_p->GetYaxis()->SetTitleFont(43);
  hist10_p->GetXaxis()->SetTitleSize(18);
  hist10_p->GetYaxis()->SetTitleSize(18);

  prettyTH1(dummy10_p, 0.01, 20, 1);
  dummy10_p->GetXaxis()->SetTitleFont(43);
  dummy10_p->GetYaxis()->SetTitleFont(43);
  dummy10_p->GetXaxis()->SetTitleSize(18);
  dummy10_p->GetYaxis()->SetTitleSize(18);

  double max10 = 0;
  double min10 = 500;

  for(unsigned int i = 0; i < times10.size(); ++i){
    if(times10.at(i) < 100){
      graph10_p->SetPoint(i, times10.at(i), number10.at(i));
      graph10True_p->SetPoint(i, times10.at(i), number10True.at(i));
    }
    hist10_p->Fill(number10.at(i));

    if(number10.at(i) > max10) max10 = number10.at(i);
    if(number10.at(i) < min10) min10 = number10.at(i);
  }

  canv10_p->cd();
  graph10_p->SetMarkerColor(col.getColor(0));
  graph10_p->SetMarkerSize(.4);
  graph10_p->SetMarkerStyle(20);

  graph10True_p->SetMarkerColor(col.getColor(2));
  graph10True_p->SetMarkerSize(.4);
  graph10True_p->SetMarkerStyle(20);

  dummy10_p->SetMaximum(max10);
  dummy10_p->SetMinimum(min10);
  dummy10_p->DrawCopy();
  graph10_p->Draw("P");
  graph10True_p->Draw("P");

  TH1F* leg1 = new TH1F();
  leg1->SetMarkerColor(col.getColor(0));
  leg1->SetMarkerStyle(20);
  leg1->SetMarkerSize(2);

  TH1F* leg2 = new TH1F();
  leg2->SetMarkerColor(col.getColor(2));
  leg2->SetMarkerStyle(20);
  leg2->SetMarkerSize(2);

  leg2_p->AddEntry(leg1, "Stoch.", "P");
  leg2_p->AddEntry(leg2, "Det.", "P");

  leg2_p->Draw("SAME");

  canv10_p->SaveAs("pdfDir/canv10.pdf");

  delete canv10_p;

  canv10_2_p->cd();
  hist10_p->DrawCopy("HIST E1");
  TF1 *f1_10_p = new TF1("f1_10_p","[1]*TMath::Poisson(x,[0])",4.5, 20.5);
  f1_10_p->SetParameter(0, 10.);
  f1_10_p->SetParameter(1, hist10_p->GetMaximum());

  hist10_p->Fit("f1_10_p", "M N", "", 4.5, 20.5);
  f1_10_p->SetMarkerColor(col.getColor(0));
  f1_10_p->SetLineColor(col.getColor(0));
  f1_10_p->SetMarkerSize(1);
  f1_10_p->SetLineWidth(2);
  f1_10_p->DrawCopy("SAME");  
  label_p->DrawLatex(.65, .8, ("#lambda_{Fit}=" + prettyString(f1_10_p->GetParameter(0), 2, false) + "#pm" + prettyString(f1_10_p->GetParError(0), 2, false)).c_str());
  label_p->DrawLatex(.65, .7, ("#chi^{2}/NDF=" + prettyString(f1_10_p->GetChisquare(), 2, false) + "/" + std::to_string((int)f1_10_p->GetNDF())).c_str());
  canv10_2_p->SaveAs("pdfDir/canv10_2.pdf");

  delete hist10_p;
  delete canv10_2_p;

  TCanvas* canv100_p = new TCanvas("canv100_c", "canv100_c", 500, 500);
  prettyCanv(canv100_p);
  TCanvas* canv100_2_p = new TCanvas("canv100_2_c", "canv100_2_c", 500, 500);
  prettyCanv(canv100_2_p);

  TGraph* graph100_p = new TGraph();
  TGraph* graph100True_p = new TGraph();
  TH1F* hist100_p = new TH1F("hist100_h", ";x;Counts (Constitutive)", 100, 49.5, 149.5);
  TH1F* dummy100_p = new TH1F("dummy100_h", ";Time (s);Counts (Constitutive)", 10, 0, 100);
  prettyTH1(hist100_p, 0.01, 20, 1);
  hist100_p->GetXaxis()->SetTitleFont(43);
  hist100_p->GetYaxis()->SetTitleFont(43);
  hist100_p->GetXaxis()->SetTitleSize(18);
  hist100_p->GetYaxis()->SetTitleSize(18);

  prettyTH1(dummy100_p, 0.01, 20, 1);
  dummy100_p->GetXaxis()->SetTitleFont(43);
  dummy100_p->GetYaxis()->SetTitleFont(43);
  dummy100_p->GetXaxis()->SetTitleSize(18);
  dummy100_p->GetYaxis()->SetTitleSize(18);

  double max100 = 0;
  double min100 = 500;

  for(unsigned int i = 0; i < times100.size(); ++i){
    if(times100.at(i) < 100){
      graph100_p->SetPoint(i, times100.at(i), number100.at(i));
      graph100True_p->SetPoint(i, times100.at(i), number100True.at(i));
    }
    hist100_p->Fill(number100.at(i));

    if(number100.at(i) > max100) max100 = number100.at(i);
    if(number100.at(i) < min100) min100 = number100.at(i);
  }

  canv100_p->cd();
  graph100_p->SetMarkerColor(col.getColor(0));
  graph100_p->SetMarkerSize(.4);
  graph100_p->SetMarkerStyle(20);

  graph100True_p->SetMarkerColor(col.getColor(2));
  graph100True_p->SetMarkerSize(.4);
  graph100True_p->SetMarkerStyle(20);

  dummy100_p->SetMaximum(max100);
  dummy100_p->SetMinimum(min100);
  dummy100_p->DrawCopy();
  graph100_p->Draw("P");
  graph100True_p->Draw("P");

  leg2_p->SetX1(.5);
  leg2_p->SetY1(.2);
  leg2_p->SetX2(.9);
  leg2_p->SetY2(.5);
  leg2_p->Draw("SAME");
  

  canv100_p->SaveAs("pdfDir/canv100.pdf");


  delete graph10_p;
  delete graph10True_p;
  delete graph100_p;
  delete graph100True_p;
  delete canv100_p;

  canv100_2_p->cd();
  hist100_p->DrawCopy("HIST E1");
  TF1 *f1_100_p = new TF1("f1_100_p","[1]*TMath::Poisson(x,[0])", 79.5, 120.5);
  f1_100_p->SetParameter(0,100.);
  f1_100_p->SetParameter(1, hist100_p->GetMaximum());

  hist100_p->Fit("f1_100_p", "M N", "", 79.5, 120.5);
  f1_100_p->SetMarkerColor(col.getColor(0));
  f1_100_p->SetMarkerSize(1);
  f1_100_p->SetLineWidth(2);
  f1_100_p->SetLineColor(col.getColor(0));
  f1_100_p->DrawCopy("SAME");
  label_p->DrawLatex(.65, .8, ("#lambda_{Fit}=" + prettyString(f1_100_p->GetParameter(0), 2, false) + "#pm" + prettyString(f1_100_p->GetParError(0), 2, false)).c_str());
  label_p->DrawLatex(.65, .7, ("#chi^{2}/NDF=" + prettyString(f1_100_p->GetChisquare(), 2, false) + "/" + std::to_string((int)f1_100_p->GetNDF())).c_str());
  canv100_2_p->SaveAs("pdfDir/canv100_2.pdf");

  delete canv100_2_p;

  TCanvas* canv100K_p = new TCanvas("canv100K_c", "canv100K_c", 500, 500);
  prettyCanv(canv100K_p);
  TCanvas* canv100K_2_p = new TCanvas("canv100K_2_c", "canv100K_2_c", 500, 500);
  prettyCanv(canv100K_2_p);

  TGraph* graph100K_p = new TGraph();
  TH1F* hist100K_p = new TH1F("hist100K_h", ";x;Counts", 100, 49.5, 149.5);
  TH1F* dummy100K_p = new TH1F("dummy100K_h", ";Time (s);Counts (Repressed)", 10, 0, 100);
  prettyTH1(hist100K_p, 0.01, 20, 1);
  hist100K_p->GetXaxis()->SetTitleFont(43);
  hist100K_p->GetYaxis()->SetTitleFont(43);
  hist100K_p->GetXaxis()->SetTitleSize(18);
  hist100K_p->GetYaxis()->SetTitleSize(18);

  prettyTH1(dummy100K_p, 0.01, 20, 1);
  dummy100K_p->GetXaxis()->SetTitleFont(43);
  dummy100K_p->GetYaxis()->SetTitleFont(43);
  dummy100K_p->GetXaxis()->SetTitleSize(18);
  dummy100K_p->GetYaxis()->SetTitleSize(18);

  double max100K = 0;
  double min100K = 500;

  for(unsigned int i = 0; i < times100K.size(); ++i){
    if(times100K.at(i) < 100) graph100K_p->SetPoint(i, times100K.at(i), number100K.at(i));
    hist100K_p->Fill(number100K.at(i));

    if(number100K.at(i) > max100K) max100K = number100K.at(i);
    if(number100K.at(i) < min100K) min100K = number100K.at(i);
  }

  canv100K_p->cd();
  graph100K_p->SetMarkerColor(col.getColor(0));
  graph100K_p->SetMarkerSize(.4);
  graph100K_p->SetMarkerStyle(20);

  dummy100K_p->SetMaximum(max100K);
  dummy100K_p->SetMinimum(min100K);
  dummy100K_p->DrawCopy();
  graph100K_p->Draw("P");
  canv100K_p->SaveAs("pdfDir/canv100K.pdf");

  delete graph100K_p;
  delete canv100K_p;

  canv100K_2_p->cd();
  hist100K_p->DrawCopy("HIST E1");
  hist100_p->SetLineColor(col.getColor(0));
  hist100_p->DrawCopy("HIST E1 SAME");

  hist100_p->SetLineWidth(3);
  hist100K_p->SetLineWidth(3);

  leg_p->AddEntry(hist100_p, "Constitutive", "P L");
  leg_p->AddEntry(hist100K_p, "Repressed", "P L");

  leg_p->Draw("SAME");

  canv100K_2_p->SaveAs("pdfDir/canv100K_2.pdf");

  delete hist100_p;
  delete hist100K_p;
  delete canv100K_2_p;

  delete randGen_p;
  
  return 0;
}

int main()
{
  int retVal = 0;
  retVal += stochProteins();
  return retVal;
}
