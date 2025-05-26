#define Analyze_cxx

#include "Analyze.h"
#include <TH2.h>
#include <TStyle.h>

//******** Definition section *********

void Analyze::Begin(TTree * /*tree*/)
{
  TString option = GetOption();

  //******** Initialization section *********
}

void Analyze::SlaveBegin(TTree *) {}

Bool_t Analyze::Process(Long64_t entry)
{
  // Don’t delete this line! Without it the program will crash. 
  fReader.SetEntry(entry);

  //******** Loop section *********
  // You probably want GetEntry(entry) here.

  return kTRUE;
}

void Analyze::SlaveTerminate() {}

void Analyze::Terminate()
{
  //******** Wrap-up section *********
}
