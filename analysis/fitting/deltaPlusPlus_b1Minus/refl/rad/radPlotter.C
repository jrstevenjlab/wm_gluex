/* This file will read in a Tree of fit parameters and create a set of
 * plots showing how the parameters of interest are affected by seed
 * and event number.
 * Original goal of this was to determine severity of m-projection
 * leakage in vector-pseudoscalar input/output test.
 */

#include <fstream>
#include <iostream>
#include "glueXstyle.C"
#include "makeTree.C"


void plotSimple(TTree *fitParameters, TString subDir);
void plotParallelCoords(TTree *fitParameters, TString subDir);
void fitGaussToPullDist(TTree *fitParameters, TString subDir);

void radPlotter(TString subDir="./") { 
  gluex_style(); // call in this TStyle first

  // if file doesn't exist, make it (argument has weird inverse logic)
  if(gSystem->AccessPathName(subDir+"seedTree.root")) {
    cout << "Tree file does not exist, making one now" << "\n";
    makeTree(subDir);
  }
  TFile *treeFile = TFile::Open(subDir+"/seedTree.root");
  TTree *fitParameters = (TTree*)(treeFile->Get("fitParameters"));
  if(!fitParameters) {
    cout << "Tree could not be obtained from treeFile" << "\n";
    return;
  }

  plotSimple(fitParameters, subDir);
  plotParallelCoords(fitParameters, subDir);
  fitGaussToPullDist(fitParameters, subDir);
  return;
}

/*************
 * FUNCTIONS *
 *************
 */

void plotSimple(TTree *fitParameters, TString subDir) {
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  // 1D hists of each s-wave (and dsRatio + total reflectivities)
  fitParameters->Draw("dsRatio>>h_dsRatio(60, 0.24, 0.30)");  

  // 2D hists for correlations of interest
  fitParameters->Draw("p1pps:p1pms>>"
		      "p1ppsVSp1pms(15, 0.0, 0.15, 15, 0.0, 0.15)",
		      "", "colz");
  // Note, to change markerStyle/Size of 2D TGraph, change it for the 
  //   tree not the graph i.e. "fitParameters->SetMarkerStyle(7)"

  // Below grab hists if fixing style (axis titles, etc)
  TH1D *h_dsRatio = (TH1D*)gDirectory->Get("h_dsRatio");
  TH1D *h_m1pms = (TH1D*)gDirectory->Get("h_m1pms");
  TH1D *h_m1pps = (TH1D*)gDirectory->Get("h_m1pps");
  TH1D *h_p1pms = (TH1D*)gDirectory->Get("h_p1pms");
  TH1D *h_p1pps = (TH1D*)gDirectory->Get("h_p1pps");
  TH1D *h_m1p = (TH1D*)gDirectory->Get("h_m1p");
  TH1D *h_p1p = (TH1D*)gDirectory->Get("h_p1p");

  delete c1;

  return;
}

/* Fits a gaussian distribution to every pull distribution whose
   expected value is > 0. 
   Ideal case is centered at 0, width of 1. Pull is value - expected, 
   so the mean to the left(right) corresponds to fit parameter being 
   over(under) estimated. Width greather than 1 corresponds to being 
   outside 1 sigma of fit uncertainty
 */
void fitGaussToPullDist(TTree *fitParameters, TString subDir) {
  gStyle->SetOptFit();
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  TObjArray *branchList = fitParameters->GetListOfBranches();
  for(int i=0; i<branchList->GetEntries(); i++) {
    TString branchName = branchList->At(i)->GetName();
      
    if(!branchName.EndsWith("pull")) {continue;}

    // draw branch with custom binning. 
    TString min = std::to_string(fitParameters->
				 GetMinimum(branchName));
    TString max = std::to_string(fitParameters->
				 GetMaximum(branchName));
    fitParameters->Draw(Form("%s>>htemp(10, %s, %s)", 
			     branchName.Data(),
			     min.Data(), max.Data()));
    TH1D *htemp = (TH1D*)gPad->GetPrimitive("htemp");
    
    // skip pullDists bound to be >0
    if(htemp->GetBinLowEdge(1) >= 0.0) {continue;}

    cout << Form("\nPerforming gaus fit to %s",branchName.Data()) 
	 << "\n";
    htemp->Fit("gaus");

    c1->Update();
    c1->Print(subDir + branchName + "_gaus.pdf");
  }

  delete c1;
  return; 
}

/* Visualize multi-dimensional data by plotting dimensions as a set 
 * of parallel strings. A single event will appear as a line 
 * connecting the coordinates in each dimension. 
 * https://root.cern/doc/master/classTParallelCoord.html
 * NOTE: The more values in a dimension that a single point in another
 *   dimension links to, the more cluttered these plots become. 
 *   Ex: a 2D plot of x vs y that is a vertical bar, will mean that in 
 *   parallel coordinates a single x point will "spray" on the y-axis
 *   Also first selection is always color blue. Is a known issue 
 *   that is too old to be fixed
 */

void plotParallelCoords(TTree *fitParameters, TString subDir) {
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);  

  // Normal Parallel Coords 
  fitParameters->Draw("m1pps:m1pms:p1pps:p1pms:dsRatio:p1p0s:m1p0s", 
		      "" , "para");
  TParallelCoord *paraTemp = (TParallelCoord*)gPad->
    GetListOfPrimitives()->FindObject("ParaCoord");
  
  // add some selections
  TParallelCoordVar* firstaxis = (TParallelCoordVar*)paraTemp->
    GetVarList()->FindObject("m1pps");
  firstaxis->AddRange(new TParallelCoordRange(firstaxis,0.08,0.12));
  paraTemp->AddSelection("violet");
  paraTemp->GetCurrentSelection()->SetLineColor(kViolet);
  
  firstaxis->AddRange(new TParallelCoordRange(firstaxis,0.0,0.013));
  paraTemp->AddSelection("voilet");
  paraTemp->GetCurrentSelection()->SetLineColor(kViolet);
  gPad->Modified();
  
  c2->Print(subDir + "paraCoord.pdf");
  
  // PULL DISTRIBUTIONS
  fitParameters->Draw("m1pps_pull:m1pms_pull:"
		      "p1pps_pull:p1pms_pull:"
		      "dsRatio_pull:"
		      "p1p0s_pull:m1p0s_pull", "" , "para");
  TParallelCoord *paraPull = (TParallelCoord*)gPad->
    GetListOfPrimitives()->FindObject("ParaCoord");
  //paraPull->SetDotsSpacing(1);
  c2->Print(subDir+"paraPull.pdf");

  // remove pulls outside some value
  TCut cutPulls = "abs(m1pps_pull) < 8 && abs(p1pps_pull) < 8 &&"
                  "abs(m1p0s_pull) < 8 && abs(p1p0s_pull) < 8 &&"
                  "abs(dsRatio_pull) < 8 &&"
                  "abs(m1pms_pull) < 8 && abs(p1pms_pull) < 8";
  fitParameters->Draw("m1pps_pull:m1pms_pull:"
		      "p1pps_pull:p1pms_pull:"
		      "dsRatio_pull:"
		      "p1p0s_pull:m1p0s_pull", cutPulls , "para");
  TParallelCoord *paraPull_cut = (TParallelCoord*)gPad->
    GetListOfPrimitives()->FindObject("ParaCoord");
  //paraPull_cut->SetDotsSpacing(1);
  paraPull_cut->SetGlobalScale(true);
  
  // add some selections and save
  TParallelCoordVar* parVar = (TParallelCoordVar*)paraPull_cut->
    GetVarList()->FindObject("m1pps_pull");
  parVar->AddRange(new TParallelCoordRange(parVar,2.0,4.81138));
  paraPull_cut->AddSelection("blue");
  paraPull_cut->GetCurrentSelection()->SetLineColor(kBlue);
  gPad->Modified();
  c2->Print(subDir+"paraPull_cut.pdf");
  
  return;
}
