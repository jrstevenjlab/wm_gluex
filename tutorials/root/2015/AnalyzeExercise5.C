#define AnalyzeExercise5_cxx

#include "AnalyzeExercise5.h"
#include <TH2.h>
#include <TStyle.h>

//******** Definition section *********
TCanvas* ebeamCanvas = NULL;
TCanvas* chi2Canvas = NULL;
TCanvas* scatterCanvas = NULL;

TH1* ebeamHist = NULL;
TH1* chi2Hist = NULL;
TH1* scatterHist = NULL;

void AnalyzeExercise5::Begin(TTree * /*tree*/)
{
  TString option = GetOption();

  //******** Initialization section *********
  ebeamCanvas = new TCanvas("c1", "ebeam canvas",10,10,700,500);   // Exercise 3
  ebeamHist = new TH1F("ebeam","Histogram of ebeam",100,149,151);
  ebeamHist->GetXaxis()->SetTitle("ebeam [GeV]");
  ebeamHist->GetYaxis()->SetTitle("number of events");

  chi2Canvas = new TCanvas("c2", "chi2 canvas",50,50,700,500);     // Exercise 3
  chi2Hist = new TH1F("chi2","Histogram of Chi2",100,0,20);
  chi2Hist->GetXaxis()->SetTitle("chi2");
  chi2Hist->GetYaxis()->SetTitle("number of events");

  scatterCanvas = new TCanvas("c3", "scatterplot",90,90,700,500);  // Exercise 5
  scatterHist = new TH2F("scatter","Chi2 vs Ebeam",100,0,20,100,149,151);
  scatterHist->GetXaxis()->SetTitle("chi2");
  scatterHist->GetYaxis()->SetTitle("ebeam [GeV]");
}

void AnalyzeExercise5::SlaveBegin(TTree* tree) {}

Bool_t AnalyzeExercise5::Process(Long64_t entry)
{
  fReader.SetEntry(entry);

  //******** Loop section *********
  //* You will probably want to put a GetEntry here. 
  GetEntry(entry);

  ebeamHist->Fill(*ebeam);         // Exercise 3
  chi2Hist->Fill(*chi2);
  scatterHist->Fill(*chi2,*ebeam); // Exercise 5

  return kTRUE;
}

void AnalyzeExercise5::SlaveTerminate() {}

void AnalyzeExercise5::Terminate()
{
  //******** Wrap-up section *********
  ebeamHist->Fit("gaus");          // Exercise 4
  gStyle->SetOptFit();             // Exercise 4
  ebeamCanvas->cd();               // Exercise 3
  ebeamHist->Draw("e1"); 

  chi2Canvas->cd();                // Exercise 3
  chi2Hist->Draw("e1");            // Exercise 2

  scatterCanvas->cd();             // Exercise 5
  scatterHist->Draw();
}
