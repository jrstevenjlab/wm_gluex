//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  7 12:48:43 2017 by ROOT version 5.34/34
// from TTree Thrown_Tree/Thrown_Tree
// found on file: tree_thrown_30496.root
//////////////////////////////////////////////////////////

#ifndef Thrown_534_h
#define Thrown_534_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class Thrown_534 : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   
   TFile          *outfile;

   // Declaration of leaf types
   UInt_t          RunNumber;
   ULong64_t       EventNumber;
   Int_t           ThrownBeam__PID;
   TLorentzVector  *ThrownBeam__X4;
   TLorentzVector  *ThrownBeam__P4;
   ULong64_t       NumPIDThrown_FinalState;
   ULong64_t       PIDThrown_Decaying;
   Float_t         MCWeight;
   UInt_t          NumThrown;
   Int_t           Thrown__ParentIndex[7];   //[NumThrown]
   Int_t           Thrown__PID[7];   //[NumThrown]
   TClonesArray    *Thrown__X4;
   TClonesArray    *Thrown__P4;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_ThrownBeam__PID;   //!
   TBranch        *b_ThrownBeam__X4;   //!
   TBranch        *b_ThrownBeam__P4;   //!
   TBranch        *b_NumPIDThrown_FinalState;   //!
   TBranch        *b_PIDThrown_Decaying;   //!
   TBranch        *b_MCWeight;   //!
   TBranch        *b_NumThrown;   //!
   TBranch        *b_Thrown__ParentIndex;   //!
   TBranch        *b_Thrown__PID;   //!
   TBranch        *b_Thrown__X4;   //!
   TBranch        *b_Thrown__P4;   //!

   Thrown_534(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Thrown_534() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Thrown_534,0);
};

#endif

#ifdef Thrown_534_cxx
void Thrown_534::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ThrownBeam__X4 = 0;
   ThrownBeam__P4 = 0;
   Thrown__X4 = 0;
   Thrown__P4 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("ThrownBeam__PID", &ThrownBeam__PID, &b_ThrownBeam__PID);
   fChain->SetBranchAddress("ThrownBeam__X4", &ThrownBeam__X4, &b_ThrownBeam__X4);
   fChain->SetBranchAddress("ThrownBeam__P4", &ThrownBeam__P4, &b_ThrownBeam__P4);
   fChain->SetBranchAddress("NumPIDThrown_FinalState", &NumPIDThrown_FinalState, &b_NumPIDThrown_FinalState);
   fChain->SetBranchAddress("PIDThrown_Decaying", &PIDThrown_Decaying, &b_PIDThrown_Decaying);
   fChain->SetBranchAddress("MCWeight", &MCWeight, &b_MCWeight);
   fChain->SetBranchAddress("NumThrown", &NumThrown, &b_NumThrown);
   fChain->SetBranchAddress("Thrown__ParentIndex", Thrown__ParentIndex, &b_Thrown__ParentIndex);
   fChain->SetBranchAddress("Thrown__PID", Thrown__PID, &b_Thrown__PID);
   fChain->SetBranchAddress("Thrown__X4", &Thrown__X4, &b_Thrown__X4);
   fChain->SetBranchAddress("Thrown__P4", &Thrown__P4, &b_Thrown__P4);

   outfile = new TFile("hist_thrown.root", "recreate");
}

Bool_t Thrown_534::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Thrown_534_cxx
