/* This file will read in a set of fit fractions (and other parameters)
 * produced by the gen_vec_ps fitter, and create a set of plots showing
 * how the parameters of interest are affected by seed and event number
 * Original goal of this was to determine severity of m-projection
 * leakage in vector-pseudoscalar input/output test.
 */

#include <fstream>
#include <iostream>

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
 * fractions (and other params) in the array that link to the TTrees
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
void fillArrays(int seedNum, double *categories[], int numberOfCategories) {
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
  int arrayPosition = 2; // start at 2 to skip dsRatio value/error
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

    arrayPosition+=2;
  }
  dataFile.close();
}

/*| *****************
  | * MAIN FUNCTION * 
  | *****************
 */


void radPlot_fitFrac() { 
  // may have user pass generated spin proj's here
  // determine # of seeds we're working with
  int numSeeds = getNumberOfSeeds();
  if(numSeeds == 0) {
    cout << "No seed directories found" << "\n"; 
    return;
  }

  cout << Form("\nFound %i seed directories", numSeeds) << "\n";
  
  // Create NTuple of type double for each parameter and its error
  TTree *fitParameters = new TTree("fitParameters", "Fit Parameters");

  // Make categories for the tree (parameters from fitPars files)
  double dsRatio, dsRatio_err, 
    p1pps, p1pps_err, m1pps, m1pps_err, p1p0s, p1p0s_err, m1p0s, m1p0s_err,
    p1pms, p1pms_err, m1pms, m1pms_err, p1ppd, p1ppd_err, m1ppd, m1ppd_err, 
    p1p0d, p1p0d_err, m1p0d, m1p0d_err, p1pmd, p1pmd_err, m1pmd, m1pmd_err, 
    p1p, p1p_err, m1p, m1p_err; 
  double *categories[] = {&dsRatio, &dsRatio_err, 
			  &p1pps, &p1pps_err, &m1pps, &m1pps_err, 
			  &p1p0s, &p1p0s_err, &m1p0s, &m1p0s_err, 
			  &p1pms, &p1pms_err, &m1pms, &m1pms_err, 
			  &p1ppd, &p1ppd_err, &m1ppd, &m1ppd_err, 
			  &p1p0d, &p1p0d_err, &m1p0d, &m1p0d_err, 
			  &p1pmd, &p1pmd_err, &m1pmd, &m1pmd_err, 
			  &p1p, &p1p_err, &m1p, &m1p_err};
  TString categoryNames[] = {"dsRatio", "dsRatio_err", 
			  "p1pps", "p1pps_err", "m1pps", "m1pps_err", 
			  "p1p0s", "p1p0s_err", "m1p0s", "m1p0s_err", 
			  "p1pms", "p1pms_err", "m1pms", "m1pms_err", 
			  "p1ppd", "p1ppd_err", "m1ppd", "m1ppd_err", 
			  "p1p0d", "p1p0d_err", "m1p0d", "m1p0d_err", 
			  "p1pmd", "p1pmd_err", "m1pmd", "m1pmd_err", 
			  "p1p", "p1p_err", "m1p", "m1p_err"};
  
  int numberOfCategories = sizeof(categories) / sizeof(*categories);

  // Create branches for each category
  for(int i=0; i<numberOfCategories; i++) {
    fitParameters->Branch(categoryNames[i], *(categories+i), 
			  categoryNames[i]+"/D");
  }
  
  // fill TTree with values from fit of each seed
  for(int seedNum=1; seedNum <= numSeeds; seedNum++) {
    fillArrays(seedNum, categories, numberOfCategories);    
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
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    
    // below is all testing
    fitParameters->Draw("dsRatio>>h_dsRatio(60, 0.24, 0.30)");
    c1->Print("dsRatio.pdf"); c1->Clear();

    fitParameters->Draw("m1pms>>m1pms(30, 0.6, 0.9)");
    c1->Print("m1pms.pdf"); c1->Clear();
    fitParameters->Draw("m1pps>>m1pps(15, 0.0, 0.15)");
    c1->Print("m1pps.pdf"); c1->Clear();
    fitParameters->Draw("p1pms>>p1pms(15, 0.0, 0.15)");
    c1->Print("p1pms.pdf"); c1->Clear();
    fitParameters->Draw("p1pps>>p1pps(15, 0.0, 0.15)");
    c1->Print("p1pps.pdf"); c1->Clear();
    
    fitParameters->Draw("m1p>>m1p(20, 0.8, 1.0)");
    c1->Print("m1p.pdf"); c1->Clear();
    fitParameters->Draw("p1p>>p1p(20, 0.0, 0.20)");
    c1->Print("p1p.pdf"); c1->Clear();
    
    
    TH1F *h_dsRatio = (TH1F*)gDirectory->Get("h_dsRatio");
  }

  TFile outFile("seedTree.root", "RECREATE");
  fitParameters->Write();
  outFile.Close();

  return;
}
