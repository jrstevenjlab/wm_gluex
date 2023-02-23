/* This file will read in a set of fit fractions (and other parameters)
 * produced by the gen_vec_ps fitter, and create a set of plots showing
 * how the parameters of interest are affected by seed and event number
 * Original goal of this was to determine severity of m-projection
 * leakage in vector-pseudoscalar input/output test.
 */


/* FUTURE WORK:
 * 1. Change arrays to vectors and use iterators of those vectors, will make 
 *    current naming of array iterators more clear
 * 2. Determine why "SetDotsSpacing" for ParaCoords dissapears for half of 
 *    the first axis plotted. Also why LineColors seem pre-fixed unchangeable
 */

#include <fstream>
#include <iostream>
#include "glueXstyle.C"

int getNumberOfSeeds();
void fillValues(int seedNum, double *categories[], int numberOfCategories,
		double expectedValues[]);
void plotSingleSeeds();
void plotManySeeds(TTree *fitParameters);
void plotParallelCoords(TTree *fitParameters);


void radPlotter() { 
  gluex_style(); // call in this TStyle first thing
  
  // determine # of seeds we're working with
  int numSeeds = getNumberOfSeeds();
  if(numSeeds == 0) {
    cout << "No seed directories found" << "\n"; 
    return;
  }

  cout << Form("\nFound %i seed directories", numSeeds) << "\n";
  
  // Create TTree to store each parameter with its error
  TTree *fitParameters = new TTree("fitParameters", "Fit Parameters");

  // I'd reeeally like to clean this all up below, as scaling this will quickly
  //   become a nightmare

  // Make categories for the tree (parameters from fitPars files)
  double dsRatio, dsRatio_err, dsRatio_pull,
    p1pps, p1pps_err, p1pps_pull, m1pps, m1pps_err, m1pps_pull, 
    p1p0s, p1p0s_err, p1p0s_pull, m1p0s, m1p0s_err, m1p0s_pull,
    p1pms, p1pms_err, p1pms_pull, m1pms, m1pms_err, m1pms_pull,
    p1ppd, p1ppd_err, p1ppd_pull, m1ppd, m1ppd_err, m1ppd_pull,
    p1p0d, p1p0d_err, p1p0d_pull, m1p0d, m1p0d_err, m1p0d_pull,
    p1pmd, p1pmd_err, p1pmd_pull, m1pmd, m1pmd_err, m1pmd_pull,
    p1p, p1p_err, p1p_pull, m1p, m1p_err, m1p_pull; 
  double *categories[] = {&dsRatio, &dsRatio_err, &dsRatio_pull,
			  &p1pps, &p1pps_err, &p1pps_pull, 
			  &m1pps, &m1pps_err, &m1pps_pull,
			  &p1p0s, &p1p0s_err, &p1p0s_pull,
			  &m1p0s, &m1p0s_err, &m1p0s_pull,
			  &p1pms, &p1pms_err, &p1pms_pull,
			  &m1pms, &m1pms_err, &m1pms_pull,
			  &p1ppd, &p1ppd_err, &p1ppd_pull,
			  &m1ppd, &m1ppd_err, &m1ppd_pull,
			  &p1p0d, &p1p0d_err, &p1p0d_pull,
			  &m1p0d, &m1p0d_err, &m1p0d_pull,
			  &p1pmd, &p1pmd_err, &p1pmd_pull,
			  &m1pmd, &m1pmd_err, &m1pmd_pull,
			  &p1p, &p1p_err, &p1p_pull, 
			  &m1p, &m1p_err, &m1p_pull};
  TString categoryNames[] = {"dsRatio", "dsRatio_err", "dsRatio_pull",
			     "p1pps", "p1pps_err", "p1pps_pull", 
			     "m1pps", "m1pps_err", "m1pps_pull",
			     "p1p0s", "p1p0s_err", "p1p0s_pull",
			     "m1p0s", "m1p0s_err", "m1p0s_pull",
			     "p1pms", "p1pms_err", "p1pms_pull",
			     "m1pms", "m1pms_err", "m1pms_pull",
			     "p1ppd", "p1ppd_err", "p1ppd_pull",
			     "m1ppd", "m1ppd_err", "m1ppd_pull",
			     "p1p0d", "p1p0d_err", "p1p0d_pull",
			     "m1p0d", "m1p0d_err", "m1p0d_pull",
			     "p1pmd", "p1pmd_err", "p1pmd_pull",
			     "m1pmd", "m1pmd_err", "m1pmd_pull",
			     "p1p", "p1p_err", "p1p_pull",
			     "m1p", "m1p_err", "m1p_pull"};
  // Expected values for pull parameters. Edit this for each new batch
  //   The total reflectivity contributions are easy to know, but we need to 
  //   account for fractions splitting between s/d wave
  //   i.e. solving "s+d=total, sqrt(d/s) = 0.27"
  double sWaveFactor = 1 + 0.27*0.27; 
  double dWaveFactor = 1 + 1/(0.27*0.27);
  double expectedValues[] = {0.27,
			     0., 0., 0., 0., 0.1/sWaveFactor, 0.9/sWaveFactor,
			     0., 0., 0., 0., 0.1/dWaveFactor, 0.9/dWaveFactor,
			     0.1, 0.9};
  
  
  int numberOfCategories = sizeof(categories) / sizeof(*categories);

  // Create branches for each category
  for(int i=0; i<numberOfCategories; i++) {
    fitParameters->Branch(categoryNames[i], *(categories+i), 
			  categoryNames[i]+"/D");
  }
  
  // fill TTree with values from fit of each seed
  for(int seedNum=1; seedNum <= numSeeds; seedNum++) {
    fillValues(seedNum, categories, numberOfCategories, expectedValues);    
    fitParameters->Fill();
  } 

  // Plots change style past 10 seeds, as plotting by individual seed
  // number becomes impractical 
  if(numSeeds <= 10) {
    TCanvas *cc = new TCanvas("cc", "cc", 800, 600);
    /* Go through tree line by line, and somehow access the array of
     * params. Then add each point of interest to graphs corresponding
     * to the line (or seed #). May move this to above for-loop
     * since each array already has these points built in
     */
  }

  if(numSeeds > 10) {    
    plotManySeeds(fitParameters);
    plotParallelCoords(fitParameters);
  }


  // Save filled tree to file 
  TString fileOutName = "seedTree.root";
  TFile outFile(fileOutName, "RECREATE");
  fitParameters->Write();
  cout << Form("TTree saved to %s", fileOutName.Data()) << "\n";
  outFile.Close();

  return;
}

/*************
 * FUNCTIONS *
 *************
 */

/* Loop through files in current directory to count # of directories
 * starting with "seed_". Assumes same dir as "run_rad.csh" script
 */
int getNumberOfSeeds() {
  int seedCount{0};
  
  // Find directories named "seed_*"
  TSystemDirectory dir("myDir", "./");
  TList *files = dir.GetListOfFiles();

  if(!files) {
    cout << "No files found" << "\n";
    return seedCount;
  }

  TSystemFile *file;
  TString fileName;
  TIter next(files);
  while ( (file=(TSystemFile*)next()) ) {
    fileName = file->GetName();
    if(file->IsDirectory() && fileName.BeginsWith("seed_"))
      seedCount += 1;
  }
  return seedCount;
}

/* Opens fit parameters file produced by vecps, and stores fit
 * fractions (and other params) in array that points to Tree entries
 *
 * Uses the same naming conventions as cfg files for amplitudes,
 * and assumes that amplitude fit fracs are written such that 
 * combinations are varied from left to right below:
 * 	reflectivity : J : p : m : l
 * i.e. amps go in order of p1pps, m1pps, p1p0s, m1p0s, ...
 * 
 * In future, try using some sort of dictionary to match the line
 * entries to the associated amplitude column in the TTree
 */
void fillValues(int seedNum, double *categories[], int numberOfCategories,
		double expectedValues[]) {
  TString fileName = Form("seed_%i/vecps_fitPars_reflRatio_rad.root",
			  seedNum); // generalize this later  
  fstream dataFile;
  dataFile.open(fileName, std::fstream::in);
  
  if(!dataFile.is_open()) {
    cout << Form("Skipping %s, file could not be opened\n", fileName.Data());
    return;
  }

  cout << Form("\n____Opening File: %s____\n", fileName.Data());
  string line;
  int arrayPosition = 3; // skips dsRatio value/err/pull  
  int arrayLimit = numberOfCategories - arrayPosition;

  // read opened file line by line
  while (std::getline(dataFile, line)) {
    std::istringstream iss(line);
      
    // for some reason vec_ps spits out dsratio with tab delim,
    // so handle that special case here. Could break in future!
    if(line.find("\t") < 15) {
      string ratio, err;       
      iss >> ratio >> err;
      **categories     = std::stod(ratio);
      **(categories+1) = std::stod(err);
      **(categories+2) = (std::stod(ratio) - 
			  expectedValues[arrayPosition/3-1]) / std::stod(err);
      continue;
    }

    // convert line into vector, broken up by spaces
    vector<string> lineVec;
    string entry;
    while (std::getline(iss, entry, ' ')) {lineVec.push_back(entry); }

    // check if we're on a coherent sum line
    if(lineVec[2] != "(coherent") {continue;}
    if(arrayPosition > arrayLimit) {
      cout << "Number of coherent sums has exceeded number of "
	      "given parameters" << "\n";
      continue;
    }

    cout << line << "\n";
    // below assumes val/err always in same line position
    **(categories + arrayPosition) = std::stod(lineVec[7]);
    **(categories + arrayPosition + 1) = std::stod(lineVec[9]);
    **(categories + arrayPosition + 2) = (std::stod(lineVec[7]) - 
					  expectedValues[arrayPosition/3]
					  ) / std::stod(lineVec[9]);
    
    arrayPosition+=3; // skip over error and pull categories
  }
  dataFile.close();
}

/**********************
 * PLOTTING FUNCTIONS *
 **********************
 */

void plotManySeeds(TTree *fitParameters) {
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  
  // 1D hists of each s-wave (and dsRatio + total reflectivities)
  fitParameters->Draw("dsRatio>>h_dsRatio(60, 0.24, 0.30)");  
  fitParameters->Draw("p1pps>>h_p1pps(15, 0.0, 0.15)");
  fitParameters->Draw("m1pps>>h_m1pps(15, 0.0, 0.15)");
  fitParameters->Draw("p1p0s>>h_p1p0s(15, 0.0, 0.15)");
  fitParameters->Draw("m1p0s>>h_m1p0s(15, 0.0, 0.15)");
  fitParameters->Draw("p1pms>>h_p1pms(15, 0.0, 0.15)");
  fitParameters->Draw("m1pms>>h_m1pms(40, 0.6, 1.0)");
  fitParameters->Draw("m1p>>h_m1p(20, 0.8, 1.0)");
  fitParameters->Draw("p1p>>h_p1p(20, 0.0, 0.20)");

  // 2D hists for correlations of interest
  fitParameters->Draw("p1pps:p1pms>>p1ppsVSp1pms(15, 0.0, 0.15, 15, 0.0, 0.15)",
		      "", "colz");
  // Note, to change markerStyle/Size of 2D TGraph, change it for the tree,
  //   not the graph i.e. "fitParameters->SetMarkerStyle(7)"

  // Below grab hists if fixing style (axis titles, etc)
  TH1D *h_dsRatio = (TH1D*)gDirectory->Get("h_dsRatio");
  TH1D *h_m1pms = (TH1D*)gDirectory->Get("h_m1pms");
  TH1D *h_m1pps = (TH1D*)gDirectory->Get("h_m1pps");
  TH1D *h_p1pms = (TH1D*)gDirectory->Get("h_p1pms");
  TH1D *h_p1pps = (TH1D*)gDirectory->Get("h_p1pps");
  TH1D *h_m1p = (TH1D*)gDirectory->Get("h_m1p");
  TH1D *h_p1p = (TH1D*)gDirectory->Get("h_p1p");

  return;
}

/* Visualize multi-dimensional data by plotting dimensions as a set of parallel
 * strings. A single event will appear as a line connecting the coordinates in 
 * each dimension. https://root.cern/doc/master/classTParallelCoord.html
 * NOTE: The more values in a dimension that a single point in another dimension
 *   links to, the more cluttered these plots become. Ex: a 2D plot of x vs y 
 *   that is a vertical bar, will mean that in parallel coordinates a single 
 *   x point will "spray" onto the y-axis
 */

void plotParallelCoords(TTree *fitParameters) {
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);  

  // Normal Parallel Coords
  fitParameters->Draw("p1pps:m1pps:p1p0s:m1p0s:p1pms:m1pms", "" , "para");
  c2->Print("para_m+0-.pdf");
  fitParameters->Draw("p1pps:p1p0s:p1pms:m1pps:m1p0s:m1pms", "" , "para");
  c2->Print("para_e+-.pdf");
  fitParameters->Draw("m1pps:m1pms:p1pps:p1pms:p1p0s:m1p0s", "" , "para");
  c2->Print("para_m+-0.pdf");
  fitParameters->Draw("p1pps:m1p0s:p1p0s:m1pms:m1pps:p1pms", "" , "para");
  c2->Print("para_mix.pdf");
  

  fitParameters->Draw("m1pps:m1pms:p1pps:p1pms:p1p0s:m1p0s", "" , "para");
  TParallelCoord *paraTemp = 
    (TParallelCoord*)gPad->GetListOfPrimitives()->FindObject("ParaCoord");
  TParallelCoordVar* firstaxis = 
    (TParallelCoordVar*)paraTemp->GetVarList()->FindObject("m1pps");
  firstaxis->AddRange(new TParallelCoordRange(firstaxis,0.08,0.12));
  paraTemp->AddSelection("violet");
  paraTemp->GetCurrentSelection()->SetLineColor(kViolet);
  
  firstaxis->AddRange(new TParallelCoordRange(firstaxis,0.0,0.013));
  paraTemp->AddSelection("blue");
  paraTemp->GetCurrentSelection()->SetLineColor(kBlue);


  // Pull Distributions
  fitParameters->Draw("m1pps_pull:m1pms_pull:"
		      "p1pps_pull:p1pms_pull:"
		      "p1p0s_pull:m1p0s_pull", "" , "para");
  TParallelCoord *paraPull = 
    (TParallelCoord*)gPad->GetListOfPrimitives()->FindObject("ParaCoord");
  //paraPull->SetDotsSpacing(1);
  c2->Print("paraPull.pdf");

  // remove pulls outside some value
  TCut cutPulls = "abs(m1pps_pull) < 8 && abs(p1pps_pull) < 8 &&"
                  "abs(m1p0s_pull) < 8 && abs(p1p0s_pull) < 8 &&"
                  "abs(m1pms_pull) < 8 && abs(p1pms_pull) < 8";
  fitParameters->Draw("m1pps_pull:m1pms_pull:"
		      "p1pps_pull:p1pms_pull:"
		      "p1p0s_pull:m1p0s_pull", cutPulls , "para");
  TParallelCoord *paraPull_cut = 
    (TParallelCoord*)gPad->GetListOfPrimitives()->FindObject("ParaCoord");
  //paraPull_cut->SetDotsSpacing(1);
  paraPull_cut->SetGlobalScale(true);

  TParallelCoordVar* parVar = 
    (TParallelCoordVar*)paraPull_cut->GetVarList()->FindObject("m1pps_pull");
  parVar->AddRange(new TParallelCoordRange(parVar,2.0,4.81138));
  paraPull_cut->AddSelection("violet");
  paraPull_cut->GetCurrentSelection()->SetLineColorAlpha(kViolet, 0.5);
  // MAKE ROOT FORUM POST ABOUT LINE COLOR NOT WORKING. 
  // ALSO WHY SETDOTSSPACING REMOVE DOTS IN FIRST AXIS
  gPad->Modified();
  c2->Print("paraPull_cut.pdf");
  
  return;
}
