/* This file will read in a Tree of fit parameters and create a set of
 * plots showing how the parameters of interest are affected by seed
 * and event number.
 */

/* FUTURE WORK
 * PlotParallelCoordPermutations:
 *   Determine if "cuts" branch names can be turned into a TCut on 
 *   the tree. Downside could be it makes unwanted cuts to the tree
 *   before plotting
 */

#include <fstream>
#include <iostream>
#include "glueXstyle.C"
#include "MakeTree.C"
#include "MyParaCoord.C"


void PlotSimple(TTree *fitParameters, TString subDir);
void FitGaussToPullDist(TTree *fitParameters, TString subDir);
void PlotParallelCoordPermutations(TTree *fitParameters,
				   TString subDir);
void PlotParallelCoords(TTree *fitParameters, TString subDir);
map<TString, int> GetPhaseMap(TTree *fitParameters, 
			      TString branchAmp);


void RadPlotter(TString subDir="./") { 
  if(!subDir.EndsWith("/")) { subDir += "/"; }

  gluex_style(); // call in this TStyle first

  // if file doesn't exist, make it (argument has weird inverse logic)
  if(gSystem->AccessPathName(subDir+"seedTree.root")) {
    cout << "Tree file does not exist, making one now" << "\n";
    MakeTree(subDir);
  }
  TFile *treeFile = TFile::Open(subDir+"seedTree.root");
  TTree *fitParameters = (TTree*)(treeFile->Get("fitParameters"));
  if(!fitParameters) {
    cout << "Tree could not be obtained from treeFile" << "\n";
    return;
  }
  
  PlotSimple(fitParameters, subDir);
  //FitGaussToPullDist(fitParameters, subDir);
  //PlotParallelCoordPermutations(fitParameters, subDir);
  PlotParallelCoords(fitParameters, subDir);

  return;
}

/*************
 * FUNCTIONS *
 *************
 */


/* Phase differences are stored in a vector, whose index corresponds
 * to the amplitude name in vector "phaseNames". This creates a 
 * mapping between the name and index, and a map on top of this to 
 * store that map for every amplitude 
 * Ex: To get index of phase difference between p1pps and p1pms wave,
 * use i = phaseMap["p1pps"]["p1pms"], then p1pps->phaseDiff[i]
 *   (All because trees can't handle maps...)
 */
std::map<TString, std::map<TString, int>> 
	    GetPhaseMap(TTree *fitParameters) {
  std::map<TString, std::map<TString, int>> phaseMap;
  TObjArray *branchList = fitParameters->GetListOfBranches();
  std::vector<string>::iterator it;
  std::string className;
  TString branchName;

  for(int i=0; i<branchList->GetEntries(); i++) {
    branchName = branchList->At(i)->GetName();     

    TBranch* b = fitParameters->GetBranch(branchName);
    className = b->GetClassName();
    if(className != "Amplitude") continue;

    Amplitude* pAmp = 0;
    std::vector<string> *phaseNames = {};
    fitParameters->SetBranchAddress(branchName, &pAmp);
    fitParameters->GetEntry(1); // needs to be called to access vector
    phaseNames = &(pAmp->phaseNames);
    for(TString ampName : *phaseNames) {      
      it = std::find(phaseNames->begin(),
		     phaseNames->end(),
		     ampName);
      int index = it - phaseNames->begin();
      phaseMap[branchName][ampName] = index;
    }
  } 
  return phaseMap;
}

void PlotSimple(TTree *fitParameters, TString subDir) {
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  // 1D hists 
  TH1D *htemp;
  fitParameters->Draw("m1pps->cs>>h_m1pps(100, 0.0, 0.17)");  
  htemp = (TH1D*)gDirectory->Get("h_m1pps");
  htemp->GetXaxis()->SetTitle("m1pps");
  c1->Update();
  c1->Print(subDir + "m1pps.pdf");

  fitParameters->Draw("m1pms->cs>>h_m1pms(100, 0.61, 0.84)");  
  htemp = (TH1D*)gDirectory->Get("h_m1pms");
  htemp->GetXaxis()->SetTitle("m1pms");
  c1->Update();
  c1->Print(subDir + "m1pms.pdf");

  // 2D hists 
  TH2D *h2temp;

  fitParameters->SetMarkerStyle(8);  
  fitParameters->Draw("m1pps->cs:m1pms->cs>>"
		      "m1ppsVSm1pms",
		      "", "SCAT");
  h2temp = (TH2D*)gDirectory->Get("m1ppsVSm1pms");
  h2temp->SetTitle(";m1pms;m1pps");
  c1->Update();
  c1->Print(subDir + "m1ppsVSm1pms.pdf");

  fitParameters->Draw("m1pms->phaseDiff[4]:"
		      "p1pms->phaseDiff[10]>>"
		      "h_pd(100, -3.14, 3.14, 100, -3.14, 3.14)",
		      "", "SCAT");
  h2temp = (TH2D*)gDirectory->Get("h_pd");
  h2temp->SetTitle(";#phi_p1pms_p1pps;#phi_m1pms_m1pps");
  c1->Update();
  c1->Print(subDir + "phase_2d.pdf");
    
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
void FitGaussToPullDist(TTree *fitParameters, TString subDir) {
  gStyle->SetOptFit();
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  TObjArray *branchList = fitParameters->GetListOfBranches();
  for(int i=0; i<branchList->GetEntries(); i++) {
    // get amplitude name and add "->pull" to extract pull dist
    TString branchName = branchList->At(i)->GetName();
    branchName += "->pull"; 

    // draw branch with custom binning
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

    float binWidth = htemp->GetBinWidth(1);

    htemp->GetXaxis()->SetTitle(Form("%s", branchName.Data()));    
    htemp->GetYaxis()->SetTitle(Form("Seeds / %.3f", binWidth));

    c1->Update();
    c1->Print(subDir + branchName + "_gaus.pdf");
  }

  delete c1;
  return; 
}

/* FUNCTION NOT COMPLETE
 * Multi-dimensional data can be plotted where every dimension is a 
 * parallel line. Single event will be a line connecting the 
 * coordinates in each dimension. 
 * https://root.cern/doc/master/classTParallelCoord.html
 * Ordering of dimensions can reduce clutter greatly, and so a first 
 * step is to try ALL unique orderings, as is done in this function
 * NOTE: The more values in 1 dimension that a single point links to 
 *       in another, a "spray" effect occurs between axes.       
 */
void PlotParallelCoordPermutations(TTree *fitParameters,
				   TString subDir) {
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  TObjArray *tempBranchList = fitParameters->GetListOfBranches();
  std::vector<TString> branchList;  

  // convert TObjArray to vector, to access next_permutation ability
  // below for loop contains different "cuts" to change
  for(int i=0; i<tempBranchList->GetEntries(); i++) {
    TString branchName = tempBranchList->At(i)->GetName();    

    // skip non s-waves (conveniently keeps dsRatio)
    if(!branchName.Contains("s")) {continue;}
    if(branchName.Contains("0")) {continue;} // skip m=0
    
    branchList.push_back(branchName);
  }
  
  // With amplitudes of interest now selected, plot paraCoord of each
  //   using their parameters

  return;
}

/* Plot parallel coordinates in ordering best determined (manually)
 * by using the PlotParallelCoordPermutations function
 */
void PlotParallelCoords(TTree *fitParameters, TString subDir) {
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);  

  // Example Plots
  fitParameters->Draw("m1pps->cs:m1pms->cs",
  		      "" , "para");
  MyParaCoord pcEx;
  pcEx.RoundRange(2);
  pcEx.CleanupAxesNames();
  c2->Print(subDir + "paraExample.pdf");

  // COHERENT SUM PARACOORD  
  fitParameters->Draw("m1pms->cs:m1pps->cs:p1pps->cs:p1pms->cs:"
  		      "dsRatio->val",
  		      "" , "para");
  MyParaCoord pcSum;
  pcSum.RoundRange(2);
  pcSum.AddSelection("m1pms->cs", 0.8, 0.84, kViolet);
  pcSum.AddSelection("m1pms->cs", 0.61, 0.7, kViolet);
  pcSum.CleanupAxesNames();
  c2->Print(subDir + "paraCoord.pdf");
  // for presentation below
  fitParameters->Draw("m1pms->cs:m1pps->cs:p1pps->cs:p1pms->cs:"
		      "m1p0s->cs:p1p0s->cs:"
  		      "dsRatio->val",
  		      "" , "para");
  MyParaCoord pcAllM;
  pcAllM.RoundRange(2);
  pcAllM.CleanupAxesNames();
  c2->Print(subDir + "paraCoord_allm.pdf");

  // PULL DISTRIBUTIONS
  fitParameters->Draw("m1pms->pull:m1pps->pull:"
		      "p1pps->pull:p1pms->pull:"
		      "dsRatio->pull",
		      "", "para");
  MyParaCoord pcPull;
  pcPull.para_->SetGlobalScale(true);
  pcPull.RoundRange(2);
  pcPull.CleanupAxesNames();
  c2->Print(subDir+"paraPull.pdf");

  // remove pulls outside some value
  TCut cutPulls = "abs(m1pps->pull) < 8 && abs(p1pps->pull) < 8 &&"
                  "abs(m1p0s->pull) < 8 && abs(p1p0s->pull) < 8 &&"
                  "abs(dsRatio->pull) < 8 &&"
                  "abs(m1pms->pull) < 8 && abs(p1pms->pull) < 8";

  fitParameters->Draw("m1pms->pull:m1pps->pull:"
		      "p1pps->pull:p1pms->pull:"
		      "dsRatio->pull",
		      cutPulls , "para");
  MyParaCoord pcPullCut;
  pcPullCut.para_->SetGlobalScale(true);
  pcPullCut.AddSelection("m1pms->pull", -1.0, -0.04, kBlue);
  pcPullCut.RoundRange(2);
  pcPullCut.CleanupAxesNames();
  c2->Print(subDir+"paraPull_cut.pdf");

  // PHASE DIFFERENCES
  std::map<TString, std::map<TString, int>> phaseMap = 
    GetPhaseMap(fitParameters);
  
  // refl- phase diffs with main amplitude
  fitParameters->Draw(Form("m1pms->cs:"
			   "m1pms->phaseDiff[%i]:"
			   "m1pms->phaseDiff[%i]",
			   phaseMap["m1pms"]["m1pps"],
			   phaseMap["m1pms"]["m1p0s"]),
		      "", "para");
  MyParaCoord pcPhase_reflm;
  pcPhase_reflm.RoundRange(2);
  pcPhase_reflm.AddSelection("m1pms->cs", 0.8, 0.84, kViolet);
  pcPhase_reflm.AddSelection("m1pms->cs", 0.61, 0.7, kViolet);
  pcPhase_reflm.CleanupAxesNames(phaseMap);

  c2->Print(subDir + "phase_reflm.pdf");

  // refl+ phase diffs with main amplitude
  fitParameters->Draw(Form("p1pms->cs:"
			   "p1pms->phaseDiff[%i]:"
			   "p1pms->phaseDiff[%i]",
			   phaseMap["p1pms"]["p1pps"],
			   phaseMap["p1pms"]["p1p0s"]),
		      "", "para");
  MyParaCoord pcPhase_reflp;
  pcPhase_reflp.RoundRange(2);
  pcPhase_reflp.AddSelection("p1pms->cs", 0.07, 0.11, kViolet);
  pcPhase_reflp.AddSelection("p1pms->cs", 0.0, 0.04, kViolet);
  pcPhase_reflp.CleanupAxesNames(phaseMap);

  c2->Print(subDir + "phase_reflp.pdf");

  return;
}

