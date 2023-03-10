/* This file will read in a Tree of fit parameters and create a set of
 * plots showing how the parameters of interest are affected by seed
 * and event number.
 * Original goal of this was to determine severity of m-projection
 * leakage in vector-pseudoscalar input/output test.
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
  
  //PlotSimple(fitParameters, subDir);
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
 * to the amplitude name in vector "phaseNames". This creates the
 * mapping between the name and index, for easier drawing later.
 *   (All because trees can't handle maps...)
 */
map<TString, int> GetPhaseMap(TTree *fitParameters, 
			      TString branchAmp) {
  map<TString, int> phaseMap;
  Amplitude* pAmp = 0;
  fitParameters->SetBranchAddress(branchAmp, &pAmp);
  fitParameters->GetEntry(1); // needs to be called to access vector
  vector<string> *phaseNames = &(pAmp->phaseNames);
  // loop through phase difference amplitudes, and create map between
  // them and associated phase difference index 
  for(TString ampName : *phaseNames) {
    std::vector<string>::iterator it = std::find(phaseNames->begin(),
					       phaseNames->end(),
					       ampName);
    int index = it - phaseNames->begin();    
    phaseMap[ampName] = index;
  }
  
  return phaseMap;
}

void PlotSimple(TTree *fitParameters, TString subDir) {
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  // 1D hists of each s-wave (and dsRatio + total reflectivities)
  fitParameters->Draw("dsRatio->val>>h_dsRatio(60, 0.24, 0.30)");  
  c1->Print("TEST_dsRatio.pdf");

  // 2D hists for correlations of interest
  fitParameters->Draw("p1pps->cs:p1pms->cs>>"
		      "p1ppsVSp1pms(15, 0.0, 0.15, 15, 0.0, 0.15)",
		      "", "colz");
  c1->Print("TEST_p1ppsVsp1pms.pdf");

  // Note, to change markerStyle/Size of 2D TGraph, change it for the 
  //   tree not the graph i.e. "fitParameters->SetMarkerStyle(7)"

  // Grab made hists if fixing style (axis titles, etc)
  TH1D *h_dsRatio = (TH1D*)gDirectory->Get("h_dsRatio");
  
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

/* Rewrites the variable names with some custom filters depending
   on the name. Avoids overlapping, and makes the phase differnces
   understandable.
 */
void CleanupParaAxes(TParallelCoord* para, 
		     map<TString,int> phaseMap = {}) {
  TList *varList = para->GetVarList();
  std::string newName;
  int start, end;
  for(TObject *var : *varList) {
    TParallelCoordVar *paraVar = (TParallelCoordVar*)var;    
    newName = paraVar->GetName();
    // remove labels (phase difference case handled below)
    if(newName.find("->") != std::string::npos &&
       newName.find("[") == std::string::npos) {
      end = newName.find("->");
      newName = newName.substr(0,end);
      paraVar->SetName(newName.c_str());
    }
    // rename phase difference to be clearer
    if(newName.find("[") != std::string::npos &&
       phaseMap.size() != 0) {      
    start = newName.find("[") + 1;
      end = newName.find("]") - 1;
      int i = std::stoi(newName.substr(start,end));
      for(auto &it : phaseMap) {
	if(it.second == i) newName = "#phi_" + it.first;
      }
      paraVar->SetName(newName.c_str());
    }
  }    
  return;
}

/* Adds a selection range to variable "varName", from rangeStart to
 * rangeEnd. A bug exists in ROOT src that first selection is always
 * blue, and any subsequent selection receives color of previous one.
 */
void AddParaSelection(TParallelCoord* para, TString varName,
		      double rangeStart, double rangeEnd,
		      Color_t myColor) {
  TParallelCoordVar* axis = (TParallelCoordVar*)para->
    GetVarList()->FindObject(varName);
  axis->AddRange(new TParallelCoordRange(axis, rangeStart, rangeEnd));
  para->AddSelection("");
  para->GetCurrentSelection()->SetLineColor(myColor);
  return;
}


/* Plot parallel coordinates in ordering best determined (manually)
 * by using the PlotParallelCoordPermutations function
 */
void PlotParallelCoords(TTree *fitParameters, TString subDir) {
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);  

  // 
  map<TString, int> map_m1pms = GetPhaseMap(fitParameters, "m1pms");

  fitParameters->Draw(Form("m1pms->cs:m1pps->cs:m1pms->phaseDiff[%i]:"
			   "p1pps->cs:m1pms->phaseDiff[%i]:"
			   "p1pms->cs:m1pms->phaseDiff[%i]",
			   map_m1pms["m1pps"], map_m1pms["p1pps"],
			   map_m1pms["p1pms"] 
			   ),
		      "", "para");
  TParallelCoord *pcMixed = (TParallelCoord*)gPad->
    GetListOfPrimitives()->FindObject("ParaCoord");  
  AddParaSelection(pcMixed, "m1pms->cs", 0.61578, 0.7, kViolet);
  AddParaSelection(pcMixed, "m1pms->cs", 0.8, 0.8345, kViolet);
  CleanupParaAxes(pcMixed, map_m1pms);
  c2->Print("Mixed.pdf");
    
  // COHERENT SUM PARACOORD
  fitParameters->Draw("m1pps->cs:m1pms->cs:p1pps->cs:p1pms->cs:"
		      "dsRatio->val:p1p0s->cs:m1p0s->cs", 
		      "" , "para");

  TParallelCoord *paraTemp = (TParallelCoord*)gPad->
    GetListOfPrimitives()->FindObject("ParaCoord");  

  // add some selections
  TParallelCoordVar* firstaxis = (TParallelCoordVar*)paraTemp->
    GetVarList()->FindObject("m1pps->cs");
  firstaxis->AddRange(new TParallelCoordRange(firstaxis,0.08,0.12));
  paraTemp->AddSelection("violet");
  paraTemp->GetCurrentSelection()->SetLineColor(kViolet);

  paraTemp->GetCurrentSelection()->SetLineColor(kViolet);
  firstaxis->AddRange(new TParallelCoordRange(firstaxis,0.0,0.013));
  paraTemp->AddSelection("voilet");
  paraTemp->GetCurrentSelection()->SetLineColor(kViolet);

  CleanupParaAxes(paraTemp);
  gPad->Modified();
  
  c2->Print(subDir + "paraCoord.pdf");
  
  // PULL DISTRIBUTIONS
  fitParameters->Draw("m1pps->pull:m1pms->pull:"
		      "p1pps->pull:p1pms->pull:"
		      "dsRatio->pull:"
		      "p1p0s->pull:m1p0s->pull", "" , "para");
  TParallelCoord *paraPull = (TParallelCoord*)gPad->
    GetListOfPrimitives()->FindObject("ParaCoord");
  c2->Print(subDir+"paraPull.pdf");

  // remove pulls outside some value
  TCut cutPulls = "abs(m1pps->pull) < 8 && abs(p1pps->pull) < 8 &&"
                  "abs(m1p0s->pull) < 8 && abs(p1p0s->pull) < 8 &&"
                  "abs(dsRatio->pull) < 8 &&"
                  "abs(m1pms->pull) < 8 && abs(p1pms->pull) < 8";
  fitParameters->Draw("m1pps->pull:m1pms->pull:"
		      "p1pps->pull:p1pms->pull:"
		      "dsRatio->pull:"
		      "p1p0s->pull:m1p0s->pull", cutPulls , "para");
  TParallelCoord *paraPull_cut = (TParallelCoord*)gPad->
    GetListOfPrimitives()->FindObject("ParaCoord");
  paraPull_cut->SetGlobalScale(true);
  
  // add some selections and save
  TParallelCoordVar* parVar = (TParallelCoordVar*)paraPull_cut->
    GetVarList()->FindObject("m1pps->pull");
  parVar->AddRange(new TParallelCoordRange(parVar,2.0,4.81138));
  paraPull_cut->AddSelection("blue");
  paraPull_cut->GetCurrentSelection()->SetLineColor(kBlue);
  gPad->Modified();

  CleanupParaAxes(paraPull_cut);

  c2->Print(subDir+"paraPull_cut.pdf");
  
  return;
}

