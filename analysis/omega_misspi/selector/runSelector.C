// macro to process analysis TTree with TSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
   

void runSelector(TString runNumber = "30730", TString myPath = "/sciclone/gluex10/RunPeriod-2017-01/analysis/ver33/tree_pi0pipmisspim__B1_T1_U1_M7_Effic/merged/", TString myTreeName = "pi0pipmisspim__B1_T1_U1_M7_Effic_Tree") //use for mc
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
  int Proof_Nthreads = 1;

  // process signal 
  TString sampleDir = myPath;
  //sampleDir += Form("0%s/", runNumber.Data());
  cout<<"running selector on files in: "<<sampleDir.Data()<<endl;
  
  TChain *chain = new TChain(myTreeName);
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
	    DPROOFLiteManager::Process_Chain(chain, "DSelector_omega_misspi.C+", Proof_Nthreads, Form("hist_omega_misspi_gen_%s.acc.root", runNumber.Data()));
	  }
	  else {
	    DPROOFLiteManager::Process_Chain(chain, "DSelector_omega_misspi.C+", Proof_Nthreads, Form("hist_omega_misspi_data_%s.acc.root", runNumber.Data()));
	  }
  }

  return;
}
