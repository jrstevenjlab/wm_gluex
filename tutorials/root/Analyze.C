#define Analyze_cxx

#include "Analyze.h"
#include <TH2.h>
#include <TStyle.h>

//******** Definition section *********

void Analyze::Begin(TTree * /*tree*/)
{
  TString option = GetOption();

  //******** Set-up section *********
}

void Analyze::SlaveBegin(TTree* tree) {}

Bool_t Analyze::Process(Long64_t entry)
{
  //******** Loop section *********
  //* You will probably want to put a GetEntry here. 

  return kTRUE;
}

void Analyze::SlaveTerminate() {}

void Analyze::Terminate()
{
  //******** Wrap-up section *********
}
