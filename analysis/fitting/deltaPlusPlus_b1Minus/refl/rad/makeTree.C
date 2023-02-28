/* This file reads in a set of fit fractions and other fit parameters
 * produced by the gen_vec_ps fitter, and creates a TTree whose 
 * branches are each fit parameter e.g. D/S ratio, 1ppms
 * This Tree can then be read in by the radPlotter.C macro to produce
 * plots of interest
 */

#include <fstream>
#include <iostream>

int getNumberOfSeeds(TString subDir);
void fillValues(int seedNum, TString subDir, 
		std::vector<double*> branches, 
		std::vector<double> expectedValues);

void makeTree(TString subDir="./") {
  
  // determine # of seeds we're working with
  int numSeeds = getNumberOfSeeds(subDir);
  if(numSeeds == 0) {
    cout << "No seed directories found" << "\n"; 
    return;
  }

  cout << Form("\nFound %i seed directories", numSeeds) << "\n";
  
  // Create TTree to store each parameter with its error
  TTree *fitParameters = new TTree("fitParameters", "Fit Parameters");

  // Make branches for the tree (parameters from fitPars files)
  double dsRatio, dsRatio_err, dsRatio_pull,
    p1pps, p1pps_err, p1pps_pull, m1pps, m1pps_err, m1pps_pull, 
    p1p0s, p1p0s_err, p1p0s_pull, m1p0s, m1p0s_err, m1p0s_pull,
    p1pms, p1pms_err, p1pms_pull, m1pms, m1pms_err, m1pms_pull,
    p1ppd, p1ppd_err, p1ppd_pull, m1ppd, m1ppd_err, m1ppd_pull,
    p1p0d, p1p0d_err, p1p0d_pull, m1p0d, m1p0d_err, m1p0d_pull,
    p1pmd, p1pmd_err, p1pmd_pull, m1pmd, m1pmd_err, m1pmd_pull,
    p1p, p1p_err, p1p_pull, m1p, m1p_err, m1p_pull; 
  std::vector<double*> branches{
    &dsRatio, &dsRatio_err, &dsRatio_pull,
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
      &m1p, &m1p_err, &m1p_pull
      };
  std::vector<TString> branchNames{
    "dsRatio", "dsRatio_err", "dsRatio_pull",
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
      "m1p", "m1p_err", "m1p_pull"
      };

  // Expected values for pull parameters. Edit this for each new batch
  //   The total reflectivity contributions are easy to know, but we 
  //   need to account for fractions splitting between s/d wave
  //   i.e. solving "s+d=total, sqrt(d/s) = 0.27"
  double sWaveFactor = 1 + 0.27*0.27; 
  double dWaveFactor = 1 + 1/(0.27*0.27);
  std::vector<double> expectedValues{
    0.27,
      0., 0., 0., 0., 0.1/sWaveFactor, 0.9/sWaveFactor,
      0., 0., 0., 0., 0.1/dWaveFactor, 0.9/dWaveFactor,
      0.1, 0.9
      };

  if(branches.size() != branchNames.size() ||
     branches.size()/3 != expectedValues.size()) {
    cout << "Mismatch in size of branch vectors! Exiting" << "\n";
    return;
  }

  // Create branches for each category
  for(int i=0; i<branches.size(); i++) {
    fitParameters->Branch(branchNames[i], branches[i],
			  branchNames[i]+"/D");
  }
  // fill Tree with values from fit of each seed
  for(int seedNum=1; seedNum<=numSeeds; seedNum++) {
    fillValues(seedNum, subDir, branches, expectedValues);
    fitParameters->Fill();
  }

  // Save filled tree to file 
  TString fileOutName = "seedTree.root";
  TFile outFile(subDir+fileOutName, "RECREATE");
  fitParameters->Write();
  cout << Form("\nTree saved to %s", fileOutName.Data()) << "\n";
  outFile.Close();
  
  return;
}


/* Loop through files in current directory to count # of directories
 * starting with "seed_". Assumes same dir as "run_rad.csh" script
 */
int getNumberOfSeeds(TString subDir) {
  int seedCount = 0;
  
  // Find directories named "seed_*"
  TSystemDirectory dir("myDir", "./" + subDir);
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
void fillValues(int seedNum, TString subDir, 
		std::vector<double*> branches, 
		std::vector<double> expectedValues) {
  TString fileName = Form(subDir + 
			  "seed_%i/vecps_fitPars_reflRatio_rad.root",
			  seedNum); // generalize this later  
  fstream dataFile;
  dataFile.open(fileName, std::fstream::in);
  
  if(!dataFile.is_open()) {
    cout << Form("Skipping %s, file could not be opened", 
		 fileName.Data()) << "\n";
    return;
  }

  cout << Form("\n____Opening File: %s____\n", fileName.Data());
  string line;
  int i = 3; // set to 3 to skip dsRatio value/err/pull  
  int iterLimit = branches.size() - i;
  
  // read opened file line by line
  while (std::getline(dataFile, line)) {    
    std::istringstream iss(line);
      
    // for some reason vec_ps spits out dsratio with tab delim,
    // so handle that special case here. Could break in future!
    if(line.find("\t") < 15) {
      string ratio, err;       
      iss >> ratio >> err;
      *(branches[0]) = std::stod(ratio);
      *(branches[1]) = std::stod(err);
      *(branches[2]) = (std::stod(ratio) - expectedValues[0]) / 
	               std::stod(err);
      continue;
    }

    // convert line into vector, broken up by spaces
    vector<string> lineVec;
    string entry;
    while (std::getline(iss, entry, ' ')) {lineVec.push_back(entry); }

    // check if we're on a coherent sum line
    if(lineVec[2] != "(coherent") {
      continue;
    }

    if(i > iterLimit) {
      cout << "\nNumber of coherent sums has exceeded number of "
	"given parameters: "<< i << " " << iterLimit << "\n";
      return;
    }
    
    //cout << line << "\n";
    // below assumes val/err always in same line position
    *(branches[i]) = std::stod(lineVec[7]);
    *(branches[i+1]) = std::stod(lineVec[9]);
    *(branches[i+2]) = (std::stod(lineVec[7]) - expectedValues[i/3]) /
                       std::stod(lineVec[9]);
    i+=3; // skip over error and pull categories
  }
  dataFile.close();
}
