// macro to process analysis TTree with TSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
   
void runSelector(TString runNumber = "31057", TString myPath = "/sciclone/data10/jrstevens01/RunPeriod-2017-01/analysis/ver08/tree_jpsi__M83/merged/") 
//void runSelector(TString runNumber = "11366", TString myPath = "/sciclone/data10/jrstevens01/RunPeriod-2016-02/analysis/ver06/tree_lowmass_jpsi/")
{
  // Load DSelector library
  gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
  int Proof_Nthreads = 8;

  // process signal 
  TString sampleDir = myPath;
  //sampleDir += Form("0%s/", runNumber.Data());
  cout<<"running selector on files in: "<<sampleDir.Data()<<endl;
  
  TChain *chain = new TChain("jpsi__M83_Tree"); //"lowmass_electrons_robison_Tree"); //
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
	  DPROOFLiteManager::Process_Chain(chain, "DSelector_tcs.C+", Proof_Nthreads, Form("hist_tcs_%s.acc.root", runNumber.Data()));
  }

  return;
}
