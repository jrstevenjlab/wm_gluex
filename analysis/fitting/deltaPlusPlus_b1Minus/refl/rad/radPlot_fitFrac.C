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
int getNumSeeds() {
  int countSeeds{0};
  // Find directories named "seed_*"
  TSystemDirectory dir("myDir", "./");
  TList *files = dir.GetListOfFiles();
  if(files) {
    TSystemFile *file;
    TString fileName;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fileName = file->GetName();
      if(file->IsDirectory() && fileName.BeginsWith("seed_"))
	countSeeds += 1;
    }
  }
  else {
    cout << "No files found" << "\n";
    exit(1);
  }
  return countSeeds;
}

/* Opens fit parameters file produced by vecps, and stores fit
 * fractions (and other specified params) in TTrees.
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
void fillTrees(int seedNum, double *valArray[], double *errArray[],
	       int arraySize) {
  TString fileName = Form("seed_%i/vecps_fitPars_reflRatio_rad.root",
			  seedNum); // generalize this later  
  fstream dataFile;
  dataFile.open(fileName, std::fstream::in);
  
  if(dataFile.is_open()) {
    cout << Form("\n____Opening File: %s____\n", fileName.Data());
    string line;
    int arrayPosition = 1; // start at 1 to skip dsRatio position
    int arrayLimit = arraySize - arrayPosition;

    // read opened file line by line
    while (std::getline(dataFile, line)) {
      std::istringstream iss(line);
      vector<string> lineVec;
      string entry;
      
      // for some reason vec_ps spits out dsratio with tab delim,
      // so handle that special case here. Could break in future!
      if(line.find("\t") < 15) {
	int i = 0;
	while(std::getline(iss, entry, '\t')) {
	  if (i==0) **valArray = std::stod(entry);
	  if (i==1) **errArray = std::stod(entry);
	  i++;
	}
      }
      else {	
	// convert line into vector, broken up by spaces
	while (std::getline(iss, entry, ' '))
	  lineVec.push_back(entry);
	// find "coherent sum" line and get value with error
	if(lineVec[2] == "(coherent") {
	  if(arrayPosition > arrayLimit) {
	    cout << "Number of coherent sums has exceeded number of "
	      "given parameters" << "\n";
	    continue;
	  }
	  cout << line << "\n";
	  // this assumes val and err will always be in same position
	  **(valArray + arrayPosition) = std::stod(lineVec[7]);
	  **(errArray + arrayPosition) = std::stod(lineVec[9]);
	  arrayPosition++;
	}
      }
    } // end getline
  } // end dataFile.is_open
  else {
    cout << "Could not open data file!" << "\n";
    exit(1);
  }    

  dataFile.close();
}

/*| *****************
  | * MAIN FUNCTION * 
  | *****************
 */


void radPlot_fitFrac() { 
  // first, determine # of seeds we're working with
  TString genSpinProjs; // later, expand to vec to handle multiple m's
  int numSeeds = getNumSeeds();
  if(numSeeds == 0) {
    cout << "No seed directories found" << "\n"; 
    return;
  }
  else 
    cout << Form("\nFound %i seed directories", numSeeds) << "\n";
  
  // Create NTuple of type double for each parameter and its error
  TTree *paramVals = new TTree("paramVals", "Parameter Values");
  TTree *paramErrs = new TTree("paramErrs", "Parameter Errors");

  // param list is hard coded for now
  double dsRatio, p1pps, m1pps, p1p0s, m1p0s, p1pms, m1pms,
    p1ppd, m1ppd, p1p0d, m1p0d, p1pmd, m1pmd, p1p, m1p; 
  double dsRatio_err, p1pps_err, m1pps_err, p1p0s_err, m1p0s_err, 
    p1pms_err, m1pms_err, p1ppd_err, m1ppd_err, p1p0d_err, m1p0d_err,
    p1pmd_err, m1pmd_err, p1p_err, m1p_err; 
  double *valArray[] = {&dsRatio, &p1pps, &m1pps, &p1p0s, &m1p0s, &p1pms,
			 &m1pms, &p1ppd, &m1ppd, &p1p0d, &m1p0d, &p1pmd,
			 &m1pmd, &p1p, &m1p};
  double *errArray[] = {&dsRatio_err, &p1pps_err, &m1pps_err, &p1p0s_err,
		       &m1p0s_err, &p1pms_err, &m1pms_err, &p1ppd_err, 
		       &m1ppd_err, &p1p0d_err, &m1p0d_err, &p1pmd_err,
		       &m1pmd_err, &p1p_err, &m1p_err};

  // awaiting ROOT forums to tell me how to fix this nightmare below
  paramVals->Branch("dsRatio", &dsRatio, "dsRatio/D");
  paramVals->Branch("p1pps", &p1pps, "p1pps/D");
  paramVals->Branch("m1pps", &m1pps, "m1pps/D");
  paramVals->Branch("p1p0s", &p1p0s, "p1p0s/D");
  paramVals->Branch("m1p0s", &m1p0s, "m1p0s/D");
  paramVals->Branch("p1pms", &p1pms, "p1pms/D");
  paramVals->Branch("m1pms", &m1pms, "m1pms/D");
  paramVals->Branch("p1ppd", &p1ppd, "p1ppd/D");
  paramVals->Branch("m1ppd", &m1ppd, "m1ppd/D");
  paramVals->Branch("p1p0d", &p1p0d, "p1p0d/D");
  paramVals->Branch("m1p0d", &m1p0d, "m1p0d/D");
  paramVals->Branch("p1pmd", &p1pmd, "p1pmd/D");
  paramVals->Branch("m1pmd", &m1pmd, "m1pmd/D");
  paramVals->Branch("p1p", &p1p, "p1p/D");
  paramVals->Branch("m1p", &m1p, "m1p/D");
  
  paramErrs->Branch("dsRatio_err", &dsRatio_err, "dsRatio_err/D");
  paramErrs->Branch("p1pps_err", &p1pps_err, "p1pps_err/D");
  paramErrs->Branch("m1pps_err", &m1pps_err, "m1pps_err/D");
  paramErrs->Branch("p1p0s_err", &p1p0s_err, "p1p0s_err/D");
  paramErrs->Branch("m1p0s_err", &m1p0s_err, "m1p0s_err/D");
  paramErrs->Branch("p1pms_err", &p1pms_err, "p1pms_err/D");
  paramErrs->Branch("m1pms_err", &m1pms_err, "m1pms_err/D");
  paramErrs->Branch("p1ppd_err", &p1ppd_err, "p1ppd_err/D");
  paramErrs->Branch("m1ppd_err", &m1ppd_err, "m1ppd_err/D");
  paramErrs->Branch("p1p0d_err", &p1p0d_err, "p1p0d_err/D");
  paramErrs->Branch("m1p0d_err", &m1p0d_err, "m1p0d_err/D");
  paramErrs->Branch("p1pmd_err", &p1pmd_err, "p1pmd_err/D");
  paramErrs->Branch("m1pmd_err", &m1pmd_err, "m1pmd_err/D");
  paramErrs->Branch("p1p_err", &p1p_err, "p1p_err/D");
  paramErrs->Branch("m1p_err", &m1p_err, "m1p_err/D");
  
  int arraySize = sizeof(valArray) / sizeof(*valArray);
  if(arraySize != sizeof(errArray) / sizeof(*errArray)) {
    cout << "Mismatch in # of Parameter Values and Errors!" << "\n";
    exit(1);
  }
    
  for(int seedNum=1; seedNum <= numSeeds; seedNum++) {
    fillTrees(seedNum, valArray, errArray, arraySize);

    paramVals->Fill();
    paramErrs->Fill();
    /* 
       In loop here make the plots using the array of pars & errs
       above. Seperate the plots by reflectivity, since 50/50 
       would cause a headache. 

       The TTree has a built in hist method, that won't work for 
       low stats like how I want to show it, but seems like its 
       perfect for the high # of seeds case
   */
  }  

  TFile outFile("testTree.root", "RECREATE");
  paramVals->Write();
  paramErrs->Write();
  outFile.Close();

  return;
}
