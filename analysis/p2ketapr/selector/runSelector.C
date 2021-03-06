// macro to process analysis TTree with TSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
   
void runSelector(TString runNumber = "307", TString myPath = "/sciclone/gluex10/RunPeriod-2017-01/analysis/ver39/tree_kpkmetapr__B4_M35_M17/merged/")
{
  // Load DSelector library
  gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
  int Proof_Nthreads = 8;

  // process signal 
  TString sampleDir = myPath;
  //sampleDir += Form("0%s/", runNumber.Data());
  cout<<"running selector on files in: "<<sampleDir.Data()<<endl;
  
  TChain *chain = new TChain("kpkmetapr__B4_M35_M17_Tree");
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
	  DPROOFLiteManager::Process_Chain(chain, "DSelector_p2ketapr.C+", Proof_Nthreads, Form("hist_p2ketapr_%s.acc.root", runNumber.Data())); 
  }

  return;
}
