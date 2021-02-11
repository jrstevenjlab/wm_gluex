// macro to process analysis TTree with TSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
//GlueX 2017 data located at: /sciclone/gluex10/RunPeriod-2017-01/analysis/ver36/tree_pi0pippippimpim__B4_M7/merged/
//GlueX 2018-01 data located at: /sciclone/gluex10/RunPeriod-2018-01/analysis/ver08/tree_pi0pippippimpim__B4_M7/merged/
//GlueX 2018-08 data located at: /sciclone/gluex10/RunPeriod-2018-08/analysis/ver06/tree_pi0pippippimpim__B4_M7/merged/
//Simulations located at: /sciclone/gluex10/amschertz/pomega2pi_omega3pi/simulations/test##/root/trees/
//bggen located at /sciclone/gluex10/gluex_simulations/REQUESTED_MC/2017_bggen_p1/tree_pi0pippippimpim__B4_M7/
   

void runSelector(TString runNumber = "3", TString myPath = "/sciclone/gluex10/amschertz/deltaPlusPlus_b1Minus/refl/test03/phasespace/root/trees/", TString myOption = "") 

{
  // Load DSelector library
  gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
  int Proof_Nthreads = 8;

  // process signal 
  TString sampleDir = myPath;
  //sampleDir += Form("0%s/", runNumber.Data());
  cout<<"running selector on files in: "<<sampleDir.Data()<<endl;
  
  TChain *chain = new TChain("pi0pippippimpim__B4_M7_Tree");
  if(myPath.Contains("tree_thrown")) treeName = "Thrown_Tree";
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
	  DPROOFLiteManager::Process_Chain(chain, "DSelector_pomega2pi_omega3pi.C+", Proof_Nthreads, Form("hist_pomega2pi_omega3pi_%s.acc.root", runNumber.Data()), Form("tree_pomega2pi_omega3pi_%s.acc.root", runNumber.Data()), Form("%s DEFAULTFLATOFF", myOption.Data()));
  }

  return;
}
