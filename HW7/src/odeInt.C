#include <iostream>
#include <string>
#include <fstream>

#include <boost/numeric/odeint.hpp>

#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TStyle.h"

#include "include/kirchnerPalette.h"

typedef std::vector<double> state_type;
typedef boost::numeric::odeint::runge_kutta_dopri5< state_type > stepper_type;

const double Kappa = .01;
const double Beta = 1.;

void rhs(const state_type& x, state_type& dxdt, const double t)
{
  dxdt.at(0) = x.at(0)*x.at(0)/(Beta*x.at(1)*x.at(1) + x.at(0)*x.at(0)) - x.at(0);
  dxdt.at(1) = Kappa*(x.at(0) - x.at(1));
}


void write_cout(const state_type& x, const double t)
{
  std::ofstream file("odeInt.txt", std::ios::app);
  file << t << ", " << x.at(0) << ", " << x.at(1) << std::endl;
  file.close();
}


int odeInt(const std::string inFileName)
{
  std::cout << inFileName << std::endl;
  
  std::ofstream file("odeInt.txt");
  file.close();

  state_type x = {1.0, 0.0};    
  boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1E-12, 1E-12, stepper_type() ), rhs, x, 0., 1000.0, 0.1, write_cout);

  std::ifstream fileIn("odeInt.txt");

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(20);
  label_p->SetNDC();

  TGraph* graphA_p = new TGraph();
  TGraph* graphB_p = new TGraph();
  TH1F* hist_p = new TH1F("hist", ";Time (s);C or I", 10, 0, 1000);
  TH1F* hist2_p = new TH1F("hist2", ";Time (s);C or I", 10, 0, 1000);
  hist_p->SetMaximum(1.1);
  hist_p->SetMinimum(0.0);

  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();

  graphA_p->SetMarkerStyle(20);
  graphB_p->SetMarkerStyle(20);

  graphA_p->SetMarkerSize(.3);
  graphB_p->SetMarkerSize(.3);

  kirchnerPalette col;
  graphA_p->SetMarkerColor(col.getColor(0));
  graphB_p->SetMarkerColor(col.getColor(2));

  std::string getLine;
  int points = 0;
  while(std::getline(fileIn,getLine)){
    double time = std::stof(getLine.substr(0,getLine.find(",")));
    getLine.replace(0, getLine.find(",")+1, "");
    double c = std::stof(getLine.substr(0,getLine.find(",")));
    getLine.replace(0, getLine.find(",")+1, "");
    double i = std::stof(getLine);

    std::cout << time << ", " << c << ", " << i << std::endl;

    graphA_p->SetPoint(points, time, c);
    graphB_p->SetPoint(points, time, i);
    points++;
  }

  TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);

 

  hist_p->DrawCopy();
  graphA_p->Draw("P");
  graphB_p->Draw("P");

  gStyle->SetOptStat(0);

  label_p->DrawLatex(.2, .95, ("#Beta=" + std::to_string(Beta) + ", #Kappa="+std::to_string(Kappa)).c_str());

  canv_p->SaveAs("pdfDir/timeLine.pdf");

  delete graphA_p;
  delete graphB_p;
  delete canv_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./odeInt.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += odeInt(argv[1]);
  return retVal;
}
