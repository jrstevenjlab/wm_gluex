// macro to process analysis TTree with TSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
   
void runSelector(TString runNumber = "51384", TString myPath = "/sciclone/gluex10/RunPeriod-2018-08/analysis/ver05/tree_gg__B4/merged/")
{
  // Load DSelector library
  gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
  int Proof_Nthreads = 8;

  // process signal 
  TString sampleDir = myPath;
  cout<<"running selector on files in: "<<sampleDir.Data()<<endl;
  
  TChain *chain = new TChain("gg__B4_Tree");
  TSystemDirectory dir(sampleDir, sampleDir);
  TList *files = dir.GetListOfFiles();
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
		  }
	  }

	  cout<<"total entries in TChain = "<<chain->GetEntries()<<endl;
	  DPROOFLiteManager::Process_Chain(chain, "DSelector_p2gamma.C+", Proof_Nthreads, Form("hist_p2gamma_%s.acc.root", runNumber.Data()), "", "DEFAULTFLATOFF"); //Form("tree_gg__B4_0%s.acc.root", runNumber.Data()), Form("tree_flat_p2gamma_%s.acc.root", runNumber.Data())
  }

  return;
}
