#include <string>

#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"


void fitData()
{
  const Int_t nDataPoints = 20;
  const Double_t dataV[nDataPoints] = {.15, .34, .35, .42, .48, .60, .52, .63, .63, .63, .60, .66, .69, .63, .73, .69, .74, .77, .72, .75};
  const Double_t dataS[nDataPoints] = {.10, .15, .19, .24, .29, .34, .38, .43, .48, .53, .57, .62, .67, .72, .76, .81, .86, .91, .95, 1.00};
  
  TFile* outFile_p = new TFile("fitDataSystBio.root", "RECREATE");
  TGraph* dataGraph_p = new TGraph(nDataPoints, dataV, dataS);

  std::string funcForm1 = "[0]*(x/([1] + x))";
  TF1* f1_p = new TF1("f1_p", funcForm.c_str(), 0, 1.1);
  dataGraph_p->Fit("f1_p", "W M", "", 0, 1.1);

  TCanvas* canv_p = new TCanvas("data1_c", "data1_c", 500, 500);
  dataGraph_p->Draw();
  

  outFile_p->cd();

  canv_p->Write("", TObject::kOverwrite);
  delete;

  dataGraph_p->Write("", TObject::kOverwrite);
  delete dataGraph_p;
  delete f1_p;

  outFile_p->Close();
  delete outFile_p;

  return;
}
