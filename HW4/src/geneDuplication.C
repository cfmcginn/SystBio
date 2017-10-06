#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TRandom3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"

#include "include/plotUtilities.h"
#include "include/getLogBins.h"
#include "include/getLinBins.h"
#include "include/kirchnerPalette.h"

double getRandExcludeCenter(TRandom3* randGen_p, double low, double high)
{
  double retVal = 0;
  while(TMath::Abs(retVal) < .01){retVal = randGen_p->Uniform(low, high);}
  return retVal;
}

int geneDuplication(const double prob, const bool doGammaPlot = false)
{
  kirchnerPalette col;
  TDatime* date = new TDatime();
  TLine* line_p = new TLine();

  TRandom3* randGen_p = new TRandom3(0);
  
  std::cout << "Running probability p=" << prob << "..." << std::endl;

  std::map<unsigned long long, std::vector<unsigned long long> > nodeNet;
  std::map<unsigned long long, double> nodeX;
  std::map<unsigned long long, double> nodeY;
  
  std::vector<unsigned long long> temp;
  temp.push_back(1);
  
  nodeNet[0] = temp;
  nodeX[0] = getRandExcludeCenter(randGen_p, -2.8, 2.8);
  nodeY[0] = getRandExcludeCenter(randGen_p, -2.8, 2.8);
  temp.clear();
  temp.push_back(0);

  nodeNet[1] = temp;
  nodeX[1] = getRandExcludeCenter(randGen_p, -2.8, 2.8);
  nodeY[1] = getRandExcludeCenter(randGen_p, -2.8, 2.8);
  temp.clear();

  unsigned long long max = 1;

  TFile* outFile_p = new TFile("output/nodeHists.root", "RECREATE");

  for(unsigned long long i = 2; i < 1002; ++i){
    unsigned long long keyVal = 2000000;

    //Done this way as one can generate a node that had been removed
    while(nodeNet.find(keyVal) == nodeNet.end()){
      keyVal = (unsigned long long)randGen_p->Integer(i);
    }

    temp = nodeNet[keyVal];

    std::vector<unsigned long long> newNodeSet;
    for(unsigned long long j = 0; j < temp.size(); ++j){
      if(randGen_p->Uniform(0., 1.) < prob){
	newNodeSet.push_back(temp.at(j));
	nodeNet[temp.at(j)].push_back(i);

	if(nodeNet[temp.at(j)].size() > max) max = nodeNet[temp.at(j)].size();
      }
    }

    if(newNodeSet.size() > 0){
      if(newNodeSet.size() > max) max = newNodeSet.size();
      nodeNet[i] = newNodeSet;
      nodeX[i] = getRandExcludeCenter(randGen_p, -2.8, 2.8);
      nodeY[i] = getRandExcludeCenter(randGen_p, -2.8, 2.8);
    }

    if(i < 50){
      TGraph* tempGraph_p = new TGraph();
      TGraph* tempGraph2_p = new TGraph();
      TGraph* tempGraph3_p = new TGraph();
    
      tempGraph2_p->SetPoint(0, nodeX[keyVal], nodeY[keyVal]);
    
      if(nodeX.find(i) != nodeX.end()) tempGraph3_p->SetPoint(0, nodeX[i], nodeY[i]);
      else tempGraph3_p->SetPoint(0, getRandExcludeCenter(randGen_p, -2.8, 2.8), getRandExcludeCenter(randGen_p, -2.8, 2.8));

      for(unsigned long long j = 0; j <= i; ++j){
	if(nodeNet.find(j) == nodeNet.end()) continue;

	float xpos = nodeX[j];
	float ypos = nodeY[j];
	//	std::cout << "J,x,y: " << j << ", " << xpos << ", " << ypos << std::endl;
	tempGraph_p->SetPoint(j, xpos, ypos);
      }

      int pos = 0;
      while(pos < tempGraph_p->GetN()){
	double x,y;
	tempGraph_p->GetPoint(pos, x, y);
	if(TMath::Abs(x) < .001 && TMath::Abs(y) < .001) tempGraph_p->RemovePoint(pos);
	else ++pos;
      }

      pos = 0;
      while(pos < tempGraph2_p->GetN()){
	double x,y;
	tempGraph2_p->GetPoint(pos, x, y);
	if(TMath::Abs(x) < .001 && TMath::Abs(y) < .001) tempGraph2_p->RemovePoint(pos);
	else ++pos;
      }

      pos = 0;
      while(pos < tempGraph3_p->GetN()){
	double x,y;
	tempGraph3_p->GetPoint(pos, x, y);
	if(TMath::Abs(x) < .001 && TMath::Abs(y) < .001) tempGraph3_p->RemovePoint(pos);
	else ++pos;
      }

      const std::string canvGraphName = "graph_Prob" + prettyString(prob, 3, true) + "_Step" + std::to_string(i) + "_c";
      const std::string saveGraphName = "pdfDir/graph_Prob" + prettyString(prob, 3, true) + "_Step" + std::to_string(i) + "_" + std::to_string(date->GetDate()) + "_" + ".pdf";
      TCanvas* canvGraph_p = new TCanvas(canvGraphName.c_str(), canvGraphName.c_str(), 500, 500);
      prettyCanv(canvGraph_p);
      TH1F* dummyHist_p = new TH1F("dummyHist_h", ";;", 10, -3, 3);
      dummyHist_p->SetLineColor(0);
      dummyHist_p->SetMaximum(3);
      dummyHist_p->SetMinimum(-3);
      
      dummyHist_p->DrawCopy();
      tempGraph_p->SetMarkerStyle(20);
      tempGraph_p->SetMarkerSize(.5);
      tempGraph_p->SetMarkerColor(1);

      tempGraph2_p->SetMarkerStyle(20);
      tempGraph2_p->SetMarkerSize(1);
      tempGraph2_p->SetMarkerColor(col.getColor(2));

      tempGraph3_p->SetMarkerStyle(20);
      tempGraph3_p->SetMarkerSize(1);
      tempGraph3_p->SetMarkerColor(col.getColor(0));

      tempGraph_p->Draw("P");
      tempGraph2_p->Draw("P");
      tempGraph3_p->Draw("P");

      gStyle->SetOptStat(0);
      line_p->SetLineColor(1);
      line_p->SetLineWidth(1);

      for(unsigned long long j = 0; j <= i; ++j){
	if(nodeNet.find(j) == nodeNet.end()) continue;

	temp = nodeNet[j];
        double xpos = nodeX[j];
	double ypos = nodeY[j];

	for(unsigned long long k = 0; k < temp.size(); ++k){
	  double xpos2 = nodeX[temp.at(k)];
	  double ypos2 = nodeY[temp.at(k)];

	  line_p->DrawLine(xpos, ypos, xpos2, ypos2);
	}
      }
    
      double xpos = nodeX[keyVal];
      double ypos = nodeY[keyVal];
      temp = nodeNet[keyVal];
      line_p->SetLineColor(col.getColor(2));
      line_p->SetLineWidth(3);

      for(unsigned long long k = 0; k < temp.size(); ++k){
	double xpos2 = nodeX[temp.at(k)];
	double ypos2 = nodeY[temp.at(k)];
	
	line_p->DrawLine(xpos, ypos, xpos2, ypos2);
      }

      
      tempGraph3_p->GetPoint(0, xpos, ypos);
      line_p->SetLineColor(col.getColor(0));

      //      std::cout << "New node lines: " << newNodeSet.size() << std::endl;

      for(unsigned int k = 0; k < newNodeSet.size(); ++k){
	double xpos2 = nodeX[newNodeSet.at(k)];
	double ypos2 = nodeY[newNodeSet.at(k)];
	
	line_p->DrawLine(xpos, ypos, xpos2, ypos2);
      }
      

      //      canvGraph_p->SaveAs(saveGraphName.c_str());
    
      delete canvGraph_p;
      delete dummyHist_p;
      delete tempGraph_p;
      delete tempGraph2_p;
      delete tempGraph3_p;
    }

    temp.clear();
    newNodeSet.clear();
  }

  //  std::cout << "Maximum degree k_{Max}=" << max << "." << std::endl;


  const std::string histName = "degreeHist_Prob" + prettyString(prob, 3, true) + "_h";

  float binMinTemp = 0.5;
  float binMaxTemp = 200.5;

  if(prob > 0.7){
    binMinTemp = .5;
    binMaxTemp = 1000.5;
  }

  const int nBins = 100;
  const float binMin = binMinTemp;
  const float binMax = binMaxTemp;
  double bins[nBins+1];

  getLinBins(binMin, binMax, nBins, bins);

  TH1F* degreeHist_p = new TH1F(histName.c_str(), ";k (Number of edges);#frac{1}{N} #frac{dN}{dk} (N is number of nodes)", nBins, bins);
  //TH1F* degreeHist_p = new TH1F(histName.c_str(), ";k;Number of Nodes", nBins, bins);
  prettyTH1(degreeHist_p, 1.5, 20, col.getColor(0));
  degreeHist_p->SetLineColor(col.getColor(0));
  degreeHist_p->GetXaxis()->SetTitleFont(43);
  degreeHist_p->GetXaxis()->SetTitleSize(20);

  std::map<unsigned long long, std::vector<unsigned long long> >::iterator it = nodeNet.begin();

  while(it != nodeNet.end()){
    degreeHist_p->Fill(it->second.size());
    ++it;
  }

  const std::string canvName = "degreeHist_Prob" + prettyString(prob, 3, true) + "_c";
  const std::string saveName = "pdfDir/degreeHist_Prob" + prettyString(prob, 3, true) + "_" + std::to_string(date->GetDate()) + "_" + std::to_string(date->GetTime()) + ".pdf";
  
  TCanvas* canv_p = new TCanvas(canvName.c_str(), canvName.c_str(), 500, 500);
  prettyCanv(canv_p);
  
  outFile_p->cd();

  canv_p->cd();
  
  degreeHist_p->Scale(1./degreeHist_p->Integral());
  for(int i = 0; i < degreeHist_p->GetNbinsX(); ++i){
    degreeHist_p->SetBinContent(i+1, degreeHist_p->GetBinContent(i+1)/degreeHist_p->GetBinWidth(i+1));
    degreeHist_p->SetBinError(i+1, degreeHist_p->GetBinError(i+1)/degreeHist_p->GetBinWidth(i+1));
  }
  
  //  degreeHist_p->SetMaximum(1.2);
  //  degreeHist_p->SetMinimum(0.000001);

  TF1* f1_p = new TF1("f1_p", "[0]*TMath::Power(x, [1])", 0.5, 50);
  f1_p->SetParameter(0, max);
  f1_p->SetParameter(1, -3);

  TF1* f2_p = new TF1("f2_p", "[0]*TMath::Power(x, [1])", 5, 50);
  f2_p->SetParameter(0, max);
  f2_p->SetParameter(1, -3);

  degreeHist_p->Fit("f1_p", "M N", "", 0.5, 50);
  degreeHist_p->Fit("f2_p", "M N", "", 5, 50);

  TFile* outFile2_p = new TFile("output/fitFile_Compile.root", "UPDATE");
  TH1F* fitHist1_p = 0;
  TH1F* fitHist2_p = 0;

  if(outFile2_p->GetListOfKeys()->Contains("fitHist1_p")) fitHist1_p = (TH1F*)outFile2_p->Get("fitHist1_p");
  else fitHist1_p = new TH1F("fitHist1_p", ";#gamma;Counts", 41, 0.5, 4.5);

  if(outFile2_p->GetListOfKeys()->Contains("fitHist2_p")) fitHist2_p = (TH1F*)outFile2_p->Get("fitHist2_p");
  else fitHist2_p = new TH1F("fitHist2_p", ";#gamma;Counts", 41, 0.5, 4.5);
  
  fitHist1_p->Fill(-f1_p->GetParameter(1));
  fitHist2_p->Fill(-f2_p->GetParameter(1));

  fitHist1_p->Write("", TObject::kOverwrite);
  fitHist2_p->Write("", TObject::kOverwrite);

  if(doGammaPlot){
    TCanvas* gammaCanv_p = new TCanvas("gammaCanv", "gammaCanv", 500, 500);
    prettyCanv(gammaCanv_p);
    TH1F* dummyHist2_p = new TH1F("dummyHist2_p", ";Fitted #gamma_{1,2};Counts", 41, 0.5, 4.5);
    prettyTH1(dummyHist2_p, 1, 20, 1);
    dummyHist2_p->GetXaxis()->SetTitleFont(43);
    dummyHist2_p->GetYaxis()->SetTitleFont(43);

    dummyHist2_p->GetXaxis()->SetTitleSize(24);
    dummyHist2_p->GetYaxis()->SetTitleSize(20);

    dummyHist2_p->SetMaximum(TMath::Max(fitHist1_p->GetMaximum()*1.1, fitHist2_p->GetMaximum()*1.1));
    dummyHist2_p->DrawCopy();

    fitHist1_p->SetMarkerColor(col.getColor(2));
    fitHist1_p->SetFillColor(col.getColor(2));
    fitHist1_p->SetLineColor(col.getColor(2));

    fitHist2_p->SetMarkerColor(col.getColor(3));
    fitHist2_p->SetFillColor(col.getColor(3));
    fitHist2_p->SetLineColor(col.getColor(3));

    fitHist1_p->DrawCopy("SAME E1 HIST");
    fitHist2_p->DrawCopy("SAME E1 HIST");

   
    TLatex* label_p = new TLatex();
    label_p->SetTextFont(43);
    label_p->SetTextSize(18);
    label_p->SetNDC();

    label_p->DrawLatex(.6, .95, ("N_{Fills} = " + std::to_string((int)fitHist1_p->GetEntries())).c_str());
    label_p->DrawLatex(.6, .89, ("#color[" + std::to_string(col.getColor(2)) + "]{#mu_{#gamma_{1}}=" + prettyString(fitHist1_p->GetMean(), 4, false) + "#pm" + prettyString(fitHist1_p->GetMeanError(), 4, false) + "}").c_str());
    label_p->DrawLatex(.6, .83, ("#color[" + std::to_string(col.getColor(2)) + "]{#sigma_{#gamma_{1}}=" + prettyString(fitHist1_p->GetStdDev(), 4, false) + "#pm" + prettyString(fitHist1_p->GetStdDevError(), 4, false) + "}").c_str());

    label_p->DrawLatex(.6, .77, ("#color[" + std::to_string(col.getColor(3)) + "]{#mu_{#gamma_{2}}=" + prettyString(fitHist2_p->GetMean(), 4, false) + "#pm" + prettyString(fitHist2_p->GetMeanError(), 4, false) + "}").c_str());
    label_p->DrawLatex(.6, .71, ("#color[" + std::to_string(col.getColor(3)) + "]{#sigma_{#gamma_{2}}=" + prettyString(fitHist2_p->GetStdDev(), 4, false) + "#pm" + prettyString(fitHist2_p->GetStdDevError(), 4, false) + "}").c_str());

    std::string gammaSaveName = gammaCanv_p->GetName();
    gammaSaveName = "pdfDir/" + gammaSaveName + ".pdf";

    std::string gammaSaveNameLOG = gammaCanv_p->GetName();
    gammaSaveNameLOG = "pdfDir/" + gammaSaveNameLOG + "_LOG.pdf";

    gammaCanv_p->SaveAs(gammaSaveName.c_str());
    
    gammaCanv_p->cd();
    gPad->SetLogy();
    gammaCanv_p->SaveAs(gammaSaveNameLOG.c_str());
    
    delete dummyHist2_p;
    delete gammaCanv_p;
  }

  outFile2_p->Close();
  delete outFile2_p;
  
  outFile_p->cd();

  degreeHist_p->DrawCopy("P E1");  


  f1_p->SetMarkerColor(col.getColor(2));
  f1_p->SetLineColor(col.getColor(2));
  if(prob < .7) f1_p->Draw("SAME");

  f2_p->SetMarkerColor(col.getColor(3));
  f2_p->SetLineColor(col.getColor(3));
  if(prob < .7) f2_p->Draw("SAME");

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  label_p->SetNDC();

  std::string gammaFit1 = "#color[" + std::to_string(col.getColor(2)) + "]{#gamma_{1}=" + prettyString(-f1_p->GetParameter(1), 3, false) + "}";
  std::string gammaFit2 = "#color[" + std::to_string(col.getColor(3)) + "]{#gamma_{2}=" + prettyString(-f2_p->GetParameter(1), 3, false) + "}";

  label_p->DrawLatex(.48, .94, ("1000 time steps, p_{edge}=" + prettyString(prob, 2, false)).c_str());
  if(prob < .7){
    label_p->DrawLatex(.48, .86, "Form #frac{1}{N}#frac{dN(k)}{dk}=k_{0}*k^{-#gamma} (k_{0}, #gamma fit)");
    label_p->DrawLatex(.7, .78, gammaFit1.c_str());
    label_p->DrawLatex(.7, .70, gammaFit2.c_str());
  }

  gStyle->SetOptStat(0);
  gPad->SetLogx();
  gPad->SetLogy();
  degreeHist_p->Write("", TObject::kOverwrite);

  canv_p->SaveAs(saveName.c_str());

  delete canv_p;
  delete degreeHist_p;
  delete f1_p;

  outFile_p->Close();
  delete outFile_p;

  delete date;
  delete line_p;
  delete randGen_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./geneDuplication.exe <probability> <doGammaPlot>" << std::endl;
    return 0;
  }

  int retVal = 0;
  if(argc == 2) retVal += geneDuplication(std::stod(argv[1]));
  else if(argc == 3) retVal += geneDuplication(std::stod(argv[1]), std::stoi(argv[2]));
  return retVal;
}
