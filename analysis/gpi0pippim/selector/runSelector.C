// macro to process analysis TTree with TSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

// replace "*" with number 1-4, also files without _M7 flag exist for no pi0 mass constraint
// OLD MC path: 		/sciclone/gluex10/gluex_simulations/REQUESTED_MC/F2018_ver02_18_bggen_batch01/tree_gpi0pippim__B4_M7/ 
// New MC (pi0 M unconstrained):/sciclone/gluex10/gluex_simulations/REQUESTED_MC/F2018_ver02_21_bggen_batch01/tree_gpi0pippim__B4_M7/ 51117, 51314
// New MC (pi0 M csontrained):	/sciclone/gluex10/gluex_simulations/REQUESTED_MC/F2018_ver02_21_bggen_batch01/tree_gpi0pippim__B4/ 
// Thrown Paths: 		/sciclone/gluex10/gluex_simulations/REQUESTED_MC/F2018_ver02_21_bggen_batch03/tree_thrown/merged/ 51030

// 2017 data path: 		/sciclone/gluex10/RunPeriod-2017-01/analysis/ver46/tree_gpi0pippim__B4_M7/merged/
// 2018 datasets: 		/sciclone/gluex10/RunPeriod-2018-01/analysis/ver15/tree_gpi0pippim__B4_M7/merged/ 51287, 51582, 51593, 51628
// 				/sciclone/gluex10/RunPeriod-2018-08/analysis/ver15/tree_gpi0pippim__B4_M7/merged/ 50704
// test bkg data subset:	/sciclone/gluex10/RunPeriod-2017-01/analysis/ver46/tree_gpi0pippim__B4_M7_original/merged/
void runSelector(TString runNumber = "51314", TString myPath = "/sciclone/gluex10/gluex_simulations/REQUESTED_MC/F2018_ver02_21_bggen_batch01/tree_gpi0pippim__B4_M7/") 
{
  // Load DSelector library
  gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
  int Proof_Nthreads = 8;

  // process signal 
  TString sampleDir = myPath;
  cout<<"running selector on files in: "<<sampleDir.Data()<<endl;
  
  TChain *chain = new TChain("gpi0pippim__B4_M7_Tree");	// pi0 mass unconstrained
  //TChain *chain = new TChain("gpi0pippim__B4_Tree"); 	// pi0 mass constrained
  //TChain *chain = new TChain("Thrown_Tree"); 		// Thrown Trees
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
	  DPROOFLiteManager::Process_Chain(chain, "DSelector_gpi0pippim.C+", Proof_Nthreads, Form("hist_gpi0pippim_%s.acc.root", runNumber.Data()));
  }

  return;
}
