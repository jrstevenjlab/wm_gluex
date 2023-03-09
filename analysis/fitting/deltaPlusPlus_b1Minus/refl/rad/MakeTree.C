/* Idea with this is to have each parameter and amplitude be a struct,
 * that contains associated values (coherent sum, pull values, etc)
 * Some useful links below
 */

#include <fstream>
#include <iostream>

class Amplitude {
public:
  TString name;
  double cs; // coherent sum
  double cs_err;
  double pull; 
  // unfortunately trees don't support std::map correctly, so two
  // vectors are needed
  std::vector<string> phaseNames; 
  std::vector<double> phaseDiff;  

  Amplitude() : name(""), cs(), cs_err(), pull() {};

  Amplitude(std::string str,
	    double pullReference, 
	    std::vector<std::string> vecNames) {
    cs_expected_ = pullReference;
    name = str;
    InitPhaseDif(vecNames);
  }
  void SetPull() {
    pull = (cs - cs_expected_) / cs_err;
  }

private:
  double cs_expected_; 
  // initialize phase differences between an amp and current amp
  //   i.e. for an ampNmae "p1pps", the phaseDiff is between it and
  //   the current amplitude stored in "name"
  // Note: ALL amplitudes included, so non-sensical phase differences
  //   between opposite reflectivites are 0, but should be ignored
  void InitPhaseDif(std::vector<std::string> vecNames) {
    for(std::string ampName : vecNames) {
      if(ampName == name) continue;
      phaseNames.push_back(ampName);
      phaseDiff.push_back(0.0);
    }
    std::sort(phaseNames.begin(), phaseNames.end());
  }
};

class Parameter {
public:
  TString name;
  double val;
  double val_err;
  double pull;

  Parameter() : name(""), val(), val_err(), pull() {};

  Parameter(std::string str, double valReference) {
    name = str;
    val_expected_ = valReference;
  }
  void SetPull() {
    pull = (val - val_expected_) / val_err;
  }

private:
  double val_expected_;
};

int GetNumberOfSeeds(TString subDir);
bool FillValues(int seedNum, TString subDir,
		std::vector<Parameter*>, std::vector<Amplitude*>);

void MakeTree(TString subDir="./") {
  
  // determine # of seeds we're working with
  int numSeeds = GetNumberOfSeeds(subDir);
  if(numSeeds == 0) {
    cout << "No seed directories found" << "\n"; 
    return;
  }
  cout << Form("\nFound %i seed directories", numSeeds) << "\n";
  
  // This tree will hold all outputs of interest from vecFitPars file
  TTree *fitParameters = new TTree("fitParameters", "Fit Parameters");
  
  // values to account for splitting between d/s wave
  double sWaveFactor = 1 + 0.27*0.27; 
  double dWaveFactor = 1 + 1/(0.27*0.27);

  std::vector<std::string> amplitudes{
    "p1pps", "p1ppd", "p1p0s", "p1p0d", "p1pms", "p1pmd",
    "m1pps", "m1ppd", "m1p0s", "m1p0d", "m1pms", "m1pmd",
  };  

  // change 2nd paramater (expected value) depending on gen.cfg
  Parameter dsRatio("dsRatio", 0.27); 
  Amplitude p1pps("p1pps", 0.0, amplitudes);  
  Amplitude p1ppd("p1ppd", 0.0, amplitudes);
  Amplitude p1p0s("p1p0s", 0.0, amplitudes);
  Amplitude p1p0d("p1p0d", 0.0, amplitudes);
  Amplitude p1pms("p1pms", 0.1/sWaveFactor, amplitudes);
  Amplitude p1pmd("p1pmd", 0.1/dWaveFactor, amplitudes);
  Amplitude m1pps("m1pps", 0.0, amplitudes);
  Amplitude m1ppd("m1ppd", 0.0, amplitudes);
  Amplitude m1p0s("m1p0s", 0.0, amplitudes);
  Amplitude m1p0d("m1p0d", 0.0, amplitudes);
  Amplitude m1pms("m1pms", 0.9/sWaveFactor, amplitudes);
  Amplitude m1pmd("m1pmd", 0.9/dWaveFactor, amplitudes);
  
  std::vector<Parameter*> myParameters;
  std::vector<Amplitude*> myAmplitudes;
  myParameters.push_back(&dsRatio);
  myAmplitudes.push_back(&p1pps);
  myAmplitudes.push_back(&p1ppd);
  myAmplitudes.push_back(&p1p0s);
  myAmplitudes.push_back(&p1p0d);
  myAmplitudes.push_back(&p1pms);
  myAmplitudes.push_back(&p1pmd);
  myAmplitudes.push_back(&m1pps);
  myAmplitudes.push_back(&m1ppd);
  myAmplitudes.push_back(&m1p0s);
  myAmplitudes.push_back(&m1p0d);
  myAmplitudes.push_back(&m1pms);
  myAmplitudes.push_back(&m1pmd);
  
  // create branches
  for(Parameter* param : myParameters) {
    fitParameters->Branch(param->name, param);
  }
  for(Amplitude* amp : myAmplitudes) {
    fitParameters->Branch(amp->name, amp);
  }

  // fill Tree with values from fit of each seed
  for(int seedNum=1; seedNum<=numSeeds; seedNum++) {
    bool isFilled = FillValues(seedNum, subDir, myParameters, myAmplitudes);
    if(isFilled) fitParameters->Fill();
  }


  // Save filled tree to file 
  TString fileOutName = "seedTree.root";
  TFile outFile(subDir+fileOutName, "RECREATE");
  fitParameters->Write();
  cout << Form("\nTree saved to %s", 
	       (subDir+fileOutName).Data()) << "\n";
  outFile.Close();
  
  return;
}


/* Loop through files in current directory to count # of directories
 * starting with "seed_". Assumes same dir as "run_rad.csh" script
 */
int GetNumberOfSeeds(TString subDir) {
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


bool FillValues(int seedNum, TString subDir,
		std::vector<Parameter*> myParameters, 
		std::vector<Amplitude*> myAmplitudes) {
  // open fitPars file
  TString fileName = Form(subDir + 
			  "seed_%i/vecps_fitPars_reflRatio_rad.root",
			  seedNum); 
  fstream dataFile;
  dataFile.open(fileName, std::fstream::in);
  
  if(!dataFile.is_open()) {
    cout << Form("Skipping %s, file could not be opened", 
		 fileName.Data()) << "\n";
    return false;
  }
  cout << Form("\n____Opening File: %s____\n", fileName.Data());
  
  // read opened file line by line
  string line;
  while (std::getline(dataFile, line)) {    
    std::istringstream iss(line);
    
    // for some reason vec_ps spits out dsratio with tab delim,
    // so handle that special case here. Could break in future!
    if(line.find("\t") < 15) {
      string ratio, err;       
      iss >> ratio >> err;
      // only one param. In future may have to search like below
      myParameters[0]->val = std::stod(ratio);
      myParameters[0]->val_err = std::stod(err);
      myParameters[0]->SetPull();
      continue;
    }

    // convert line into vector, broken up by spaces
    vector<string> lineVec;
    string entry;
    TString refl, amplitude;
    TString amp1, amp2;
    int index;
    std::vector<string>::iterator it;
    while (std::getline(iss, entry, ' ')) {lineVec.push_back(entry); }

    // store coherent sum lines
    if(lineVec[2] == "(coherent") {
      refl = (lineVec[4] == "PosRefl") ? "p" : "m";
      amplitude = refl + lineVec[5];
      
      for(Amplitude* amp : myAmplitudes) {
	if(amp->name != amplitude) continue;
	amp->cs = std::stod(lineVec[7]);
	amp->cs_err = std::stod(lineVec[9]);
	amp->SetPull();
      }
    }
    // store phase differences
    if(lineVec[0] == "PHASE") {
      refl = (lineVec[2].find("ImagPosSign") != std::string::npos)
	? "p" : "m";
      amp1 = refl + lineVec[2].substr(lineVec[2].length() - 4);
      amp2 = refl + lineVec[3].substr(lineVec[3].length() - 4);

      // since map not supported, get index of amplitude in 
      //   "phaseNames" and then store value in phaseDiff at the 
      //   associated index
      for(Amplitude* amp : myAmplitudes) {
	if(amp->name == amp1) {
	  it = std::find(amp->phaseNames.begin(), 
			      amp->phaseNames.end(), amp2);
	  index = it - amp->phaseNames.begin();
	  amp->phaseDiff[index] = std::stod(lineVec[4]);
	}
	if(amp->name == amp2) {
	  it = std::find(amp->phaseNames.begin(), 
			      amp->phaseNames.end(), amp1);
	  index = it - amp->phaseNames.begin();
	  amp->phaseDiff[index] = std::stod(lineVec[4]);
	}      
      }
    }
  }
  return true;
}

