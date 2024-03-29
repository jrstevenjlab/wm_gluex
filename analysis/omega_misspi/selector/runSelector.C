// macro to process analysis TTree with TSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
   

void runSelector(TString runNumber = "30496", TString myPath = "/sciclone/gluex10/gluex_simulations/REQUESTED_MC/2017_ver03_28_bggen_batch01/tree_pi0pipmisspim__B1_T1_U1_M7_Effic/", TString myTree = "pi0pipmisspim__B1_T1_U1_M7_Effic", TString myOption = "")
{
  bool mc;
  if(myPath.Contains("REQUESTED_MC")) {
    mc = true;
  }
  else {
    mc = false;  
  }
  // Load DSelector library
  gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
  int Proof_Nthreads = 8;

  // process signal 
  TString sampleDir = myPath; // + "/tree_" + myTree + "/";
  //sampleDir += Form("0%s/", runNumber.Data());
  cout<<"running selector on files in: "<<sampleDir.Data()<<endl;
  
  // name for output histograms
  TString sampleName = "misspip";
  if(sampleDir.Contains("misspim")) sampleName = "misspim"; 

  TChain *chain = new TChain(myTree + "_Tree");
  TSystemDirectory dir(sampleDir, sampleDir);
  TList *files = dir.GetListOfFiles();
  int ifile = 0;
  if(files) {
	  TSystemFile *file;
	  TString fileName;
	  TIter next(files);
	  
	  // loop over files
	  while ((file=(TSystemFile*)next())) {
		  fileName = file->GetName();
		  if(fileName.Contains(runNumber)) {
			  cout<<fileName.Data()<<endl;
			  
			  // check if file corrupted
			  TFile f(sampleDir+fileName);
			  if(f.TestBit(TFile::kRecovered)) {
				  cout<<"file corrupted -> skipping"<<endl;
				  continue;
			  }
			  if(f.IsZombie()) {
				  cout<<"file is a Zombie -> skipping"<<endl;
				  continue;
			  }
			  
			  // add file to chain
			  chain->Add(sampleDir+fileName);
			  ifile++;
		  }
	  }

	  cout<<"total entries in TChain = "<<chain->GetEntries()<<" from "<<ifile<<" files"<<endl;
	  if(mc == true){
	    DPROOFLiteManager::Process_Chain(chain, "DSelector_omega_misspi.C+", Proof_Nthreads, Form("hist_omega_%s_gen_%s.acc.root", sampleName.Data(), runNumber.Data()), "", myOption.Data());
	  }
	  else {
	    DPROOFLiteManager::Process_Chain(chain, "DSelector_omega_misspi.C+", Proof_Nthreads, Form("hist_omega_%s_data_%s.acc.root", sampleName.Data(), runNumber.Data()), "", myOption.Data());
	  }
  }

  return;
}
