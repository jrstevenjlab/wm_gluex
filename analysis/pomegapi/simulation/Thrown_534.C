#define Thrown_534_cxx
// The class definition in Thrown_534.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("Thrown_534.C")
// Root > T->Process("Thrown_534.C","some options")
// Root > T->Process("Thrown_534.C+")
//

#include "Thrown_534.h"
#include <TH2.h>
#include <TStyle.h>

TH1F* dHist_MissingMassSquared;
TH1F* dHist_BeamEnergy;
TH1F* dHist_3PiMass_Measured;
TH1F* dHist_MM2_Weighted;
TH1F* dHist_lambda_peak;
TH1F* dHist_lambda_wings;
TH1F* dHist_lambda_uncut;
TH1F* dHist_4PiMass;
TH1F* dHist_OmegaPiMass;
TH2F* dHist_3vs4;
TH1F* dHist_Man_t;
TH1F* dHist_costheta;
TH1F* dHist_phi;
TH1F* dHist_costhetaH;
TH1F* dHist_phiH;
TH2F* dHist_CosThetaVsMass;
TH2F* dHist_PhiVsMass;
TH2F* dHist_CosThetaHVsMass;
TH2F* dHist_PhiHVsMass;
TH1F* dHist_CosTheta_t1;
TH1F* dHist_Phi_t1;
TH1F* dHist_CosThetaH_t1;
TH1F* dHist_PhiH_t1;
TH1F* dHist_CosTheta_t2;
TH1F* dHist_Phi_t2;
TH1F* dHist_CosThetaH_t2;
TH1F* dHist_PhiH_t2;
TH2F* dHist_CosThetaVsMass_t1;
TH2F* dHist_PhiVsMass_t1;
TH2F* dHist_CosThetaHVsMass_t1;
TH2F* dHist_PhiHVsMass_t1;
TH2F* dHist_CosThetaVsMass_t2;
TH2F* dHist_PhiVsMass_t2;
TH2F* dHist_CosThetaHVsMass_t2;
TH2F* dHist_PhiHVsMass_t2;


void Thrown_534::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   dHist_MissingMassSquared = new TH1F("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
   dHist_BeamEnergy = new TH1F("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
   dHist_3PiMass_Measured = new TH1F("3PiMass_thrown", ";3 Pion Mass (GeV/c^{2})", 600, 0.6, 1.1);
   dHist_MM2_Weighted = new TH1F("MM2_Weighted", ";Weighted Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
   dHist_lambda_peak = new TH1F("lambda_peak", ";{lambda} (peak)", 600, 0, 1);
   dHist_lambda_wings = new TH1F("lambda_wings", ";{lambda} (wings)", 600, 0, 1);
   dHist_lambda_uncut = new TH1F("lambda_uncut", ";{lambda} (uncut)", 600, 0, 1);
   dHist_4PiMass = new TH1F("4PiMass", ";Pi+Pi-Pi0Pi0 Mass (GeV/c^{2})", 600, 0, 3);
   dHist_OmegaPiMass = new TH1F("OmegaPiMass", ";Omega Pi0 Mass (GeV)", 600, 0, 3);
   dHist_3vs4 = new TH2F("3vs4", ";3 Pion Mass vs 4 Pion Mass", 600, 0, 3, 600, 0, 3);
   dHist_Man_t = new TH1F("Man_t", ";Four-Momentum Transfer Squared (GeV)^{2}", 600, 0.0, 1.0);
   //Decay Angles
   dHist_costheta = new TH1F("costheta", ";Cos(theta)", 600, -1.0, 1.0);
   dHist_phi = new TH1F("phi", ";phi (radians)", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
   dHist_costhetaH = new TH1F("costhetaH", ";Cos(theta_H)", 600, -1.0, 1.0);
   dHist_phiH = new TH1F("phiH", ";phi_H (rad)", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
   //Angles vs mass
   dHist_CosThetaVsMass = new TH2F("CosThetaVsMass", ";Cos(theta) vs Omega Pi Mass", 600, -1.0, 1.0, 20, 1.0, 3.0);
   dHist_PhiVsMass = new TH2F("PhiVsMass", ";Phi vs Omega Pi Mass", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
   dHist_CosThetaHVsMass = new TH2F("CosThetaHVsMass", ";Cos(theta_H) vs Omega Pi Mass", 600, -1.0, 1.0, 20, 1.0, 3.0);
   dHist_PhiHVsMass = new TH2F("PhiHVsMass", ";Phi vs Omega Pi Mass", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
   //Angles in t range [0.1, 0.3]
   dHist_CosTheta_t1 = new TH1F("CosTheta_t1", ";Cos(theta) with t[0.1, 0.3]", 600, -1.0, 1.0);
   dHist_Phi_t1 = new TH1F("Phi_t1", ";Phi with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
   dHist_CosThetaH_t1 = new TH1F("CosThetaH_t1", ";Cos(theta_H) with t[0.1, 0.3]", 600, -1.0, 1.0);
   dHist_PhiH_t1 = new TH1F("PhiH_t1", ";Phi_H with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
   //Angles in t range [0.3, 1.0]
   dHist_CosTheta_t2 = new TH1F("CosTheta_t2", ";Cos(theta) with t[0.3, 1.0]", 600, -1.0, 1.0);
   dHist_Phi_t2 = new TH1F("Phi_t2", ";Phi with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
   dHist_CosThetaH_t2 = new TH1F("CosThetaH_t2", ";Cos(theta_H) with t[0.3, 1.0]", 600, -1.0, 1.0);
   dHist_PhiH_t2 = new TH1F("PhiH_t2", ";Phi_H with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
   //Angles vs mass in t range [0.1, 0.3]
   dHist_CosThetaVsMass_t1 = new TH2F("CosThetaVsMass_t1", ";Cos(theta) vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0, 1.0, 20, 1.0, 3.0);
   dHist_PhiVsMass_t1 = new TH2F("PhiVsMass_t1", ";Phi vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
   dHist_CosThetaHVsMass_t1 = new TH2F("CosThetaHVsMass_t1", ";Cos(theta_H) vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0, 1.0, 20, 1.0, 3.0);
   dHist_PhiHVsMass_t1 = new TH2F("PhiHVsMass_t1", ";Phi vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
   //Angles vs mass in t range [0.3, 1.0]
   dHist_CosThetaVsMass_t2 = new TH2F("CosThetaVsMass_t2", ";Cos(theta) vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0, 1.0, 20, 1.0, 3.0);
   dHist_PhiVsMass_t2 = new TH2F("PhiVsMass_t2", ";Phi vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
   dHist_CosThetaHVsMass_t2 = new TH2F("CosThetaHVsMass_t2", ";Cos(theta_H) vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0, 1.0, 20, 1.0, 3.0);
   dHist_PhiHVsMass_t2 = new TH2F("PhiHVsMass_t2", ";Phi_H vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);


}

void Thrown_534::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Thrown_534::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Thrown_534::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
  
  fChain->GetEntry(entry);
  
  //Get P4's:
  //Step 0
  TLorentzVector dTargetP4(0, 0, 0, 0.9382720813);
  TLorentzVector *locBeamP4_Measured = (TLorentzVector*)ThrownBeam__P4;
  TLorentzVector *locPiPlusP4_Measured = static_cast<TLorentzVector*>(Thrown__P4->At(0));
  TLorentzVector *locPiMinusP4_Measured = static_cast<TLorentzVector*>(Thrown__P4->At(1));
  TLorentzVector *locProtonP4_Measured = static_cast<TLorentzVector*>(Thrown__P4->At(6));
  //Step 1
  TLorentzVector *locPhoton1P4_Measured = static_cast<TLorentzVector*>(Thrown__P4->At(2));
  TLorentzVector *locPhoton2P4_Measured = static_cast<TLorentzVector*>(Thrown__P4->At(3));
  //Step 2
  TLorentzVector *locPhoton3P4_Measured = static_cast<TLorentzVector*>(Thrown__P4->At(4));
  TLorentzVector *locPhoton4P4_Measured = static_cast<TLorentzVector*>(Thrown__P4->At(5));

  /****************************************** COMBINE FOUR-MOMENTA *******************************************************/

  //Copy/pasted from DSelector_pomegapi.C
  TLorentzVector locPi01P4_Measured = *locPhoton1P4_Measured + *locPhoton2P4_Measured;
  TLorentzVector loc3Pi1P4_Measured = *locPiPlusP4_Measured + *locPiMinusP4_Measured + locPi01P4_Measured;
  TLorentzVector locPi02P4_Measured = *locPhoton3P4_Measured + *locPhoton4P4_Measured;
  TLorentzVector loc3Pi2P4_Measured = *locPiPlusP4_Measured + *locPiMinusP4_Measured + locPi02P4_Measured;
  TLorentzVector loc4PiP4_Measured = *locPiPlusP4_Measured + *locPiMinusP4_Measured + locPi01P4_Measured + locPi02P4_Measured;
  TLorentzVector locSqrt_t = loc4PiP4_Measured - *locBeamP4_Measured;
  
  //Define omega-pi0 decay angles
  //Boost to gamma-p rest frame
  TLorentzVector locGammapP4 = *locBeamP4_Measured + dTargetP4; 
  TVector3 locGammapBoost = locGammapP4.BoostVector();     //create boost vector from gamma-p rest frame to lab frame - use negative of this to boost to gamma-p rest frame
  TLorentzVector locBeamP4_gpRest = *locBeamP4_Measured;     //boost the beam to the gamma-p rest frame
  locBeamP4_gpRest.Boost(-1.0*locGammapBoost);
  TLorentzVector locOmegaPi0P4_gpRest = loc4PiP4_Measured;    //boost omega-pi0 to the gamma-p rest frame
  locOmegaPi0P4_gpRest.Boost(-1.0*locGammapBoost);
  TLorentzVector locOmega1P4_gpRest = loc3Pi1P4_Measured;   //boost omega to gamma-p rest frame - first omega
  locOmega1P4_gpRest.Boost(-1.0*locGammapBoost);
  TLorentzVector locOmega2P4_gpRest = loc3Pi2P4_Measured;   //boost omega to gamma-p rest frame -second omega
  locOmega2P4_gpRest.Boost(-1.0*locGammapBoost);
  
  //Define unit 3-vectors - turns out we could have combined a lot of these steps
  TVector3 locBeamP3_gpRest = locBeamP4_gpRest.Vect();
  TVector3 locOmegaPi0P3_gpRest = locOmegaPi0P4_gpRest.Vect();
  TVector3 locOmega1P3_gpRest = locOmega1P4_gpRest.Vect();
  TVector3 locOmega2P3_gpRest = locOmega2P4_gpRest.Vect();
  TVector3 lock = locBeamP3_gpRest.Unit();
  TVector3 locz = locOmegaPi0P3_gpRest.Unit();
  TVector3 lockcrossz = lock.Cross(locz);
  TVector3 locy = lockcrossz.Unit();
  TVector3 locx = locy.Cross(locz);
  
  //Boost omega from gamma-p to omegapi0 rest frame
  TVector3 locOmegapiBoost = locOmegaPi0P4_gpRest.BoostVector();  //create boost vector from omega-pi rest frame to gamma-p rest frame
  TLorentzVector locOmega1P4_opiRest = locOmega1P4_gpRest;
  locOmega1P4_opiRest.Boost(-1.0*locOmegapiBoost);
  TVector3 locOmega1P3 = locOmega1P4_opiRest.Vect();
  TLorentzVector locOmega2P4_opiRest = locOmega2P4_gpRest;
  locOmega2P4_opiRest.Boost(-1.0*locOmegapiBoost);
  TVector3 locOmega2P3 = locOmega2P4_opiRest.Vect();
  
  //Define helicity unit vectors
  TVector3 loczH1 = locOmega1P3.Unit();
  TVector3 locyH1 = locz.Cross(loczH1).Unit();
  TVector3 locxH1 = locyH1.Cross(loczH1).Unit();
  TVector3 loczH2 = locOmega2P3.Unit();
  TVector3 locyH2 = locz.Cross(loczH2).Unit();
  TVector3 locxH2 = locyH2.Cross(loczH2).Unit();

  //Boost charged pions from lab to Omega rest frame
  //First Omega combo
  TVector3 locOmega1Boost = locOmega1P4_opiRest.BoostVector();   //create boost vector from omega rest frame to omega-pi rest frame
  TLorentzVector locPiPlusP4_o1Rest = *locPiPlusP4_Measured;
  locPiPlusP4_o1Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
  locPiPlusP4_o1Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
  locPiPlusP4_o1Rest.Boost(-1.0*locOmega1Boost);     //Boost to omega rf
  TVector3 locPiPlusP3_o1Rest = locPiPlusP4_o1Rest.Vect();
  TLorentzVector locPiMinusP4_o1Rest = *locPiMinusP4_Measured;
  locPiMinusP4_o1Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
  locPiMinusP4_o1Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
  locPiMinusP4_o1Rest.Boost(-1.0*locOmega1Boost);     //Boost to omega rf
  TVector3 locPiMinusP3_o1Rest = locPiMinusP4_o1Rest.Vect();
  //Second Omega combo
  TVector3 locOmega2Boost = locOmega2P4_opiRest.BoostVector();   //create boost vector from omega rest frame to omega-pi rest frame
  TLorentzVector locPiPlusP4_o2Rest = *locPiPlusP4_Measured;
  locPiPlusP4_o2Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
  locPiPlusP4_o2Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
  locPiPlusP4_o2Rest.Boost(-1.0*locOmega2Boost);     //Boost to omega rf
  TVector3 locPiPlusP3_o2Rest = locPiPlusP4_o2Rest.Vect();
  TLorentzVector locPiMinusP4_o2Rest = *locPiMinusP4_Measured;
  locPiMinusP4_o2Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
  locPiMinusP4_o2Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
  locPiMinusP4_o2Rest.Boost(-1.0*locOmega2Boost);     //Boost to omega rf
  TVector3 locPiMinusP3_o2Rest = locPiMinusP4_o2Rest.Vect();
  
  //Normal vector to decay plane
  TVector3 locnormal1 = locPiPlusP3_o1Rest.Cross(locPiMinusP3_o1Rest).Unit();
  TVector3 locnormal2 = locPiPlusP3_o2Rest.Cross(locPiMinusP3_o2Rest).Unit();
  

  // Combine 4-vectors
  TLorentzVector locMissingP4_Measured = *locBeamP4_Measured + dTargetP4;
  locMissingP4_Measured -= *locPiPlusP4_Measured + *locPiMinusP4_Measured + *locProtonP4_Measured + *locPhoton1P4_Measured + *locPhoton2P4_Measured + *locPhoton3P4_Measured + *locPhoton4P4_Measured;
  
  /************************************************ TIME FOR SOME HISTOGRAMS **************************************************/
  
  /************************************************ HISTOGRAM BEAM ENERGY ****************************************************/
  dHist_BeamEnergy->Fill((*locBeamP4_Measured).E());

  /*********************************************** HISTOGRAM 3PI MASS *******************************************************/
  double loc3PiMass1_Measured = loc3Pi1P4_Measured.M();
  double loc3PiMass2_Measured = loc3Pi2P4_Measured.M();

  dHist_3PiMass_Measured->Fill(loc3PiMass1_Measured);
  dHist_3PiMass_Measured->Fill(loc3PiMass2_Measured);

  /********************************************** HISTOGRAM MISSING MASS SQUARED ******************************************/
  double locMissingMassSquared = locMissingP4_Measured.M2();

  dHist_MissingMassSquared->Fill(locMissingMassSquared);

  /********************************************** HISTOGRAM WEIGHTED MISSING MASS SQUARED *********************************/
  double weight1;
  if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
  else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
  else weight1 = 0;
  dHist_MM2_Weighted->Fill(locMissingMassSquared, weight1);

  double weight2;
  if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
  else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
  else weight2 = 0;
  dHist_MM2_Weighted->Fill(locMissingMassSquared, weight2);

  /********************************************** HISTOGRAM DECAY MATRIX ELEMENT SQUARED (LAMBDA) **************************/
  //Boost Pi+ and Pi- to the 3Pi rest frame
  //Sum the 3Pi 4-Momenta
  TLorentzVector loc3PiP4_1 = *locPiPlusP4_Measured + *locPiMinusP4_Measured + locPi01P4_Measured;
  TVector3 b1 = loc3PiP4_1.BoostVector();
  TLorentzVector locPiPlusP4_Boosted1 = *locPiPlusP4_Measured;
  locPiPlusP4_Boosted1.Boost(-1*b1);
  TLorentzVector locPiMinusP4_Boosted1 = *locPiMinusP4_Measured;
  locPiMinusP4_Boosted1.Boost(-1*b1);
  
  //Extract 3-Momenta from 4-Momenta
  TVector3 locPiPlusP3_1 = locPiPlusP4_Boosted1.Vect();
  TVector3 locPiMinusP3_1 = locPiMinusP4_Boosted1.Vect();
  TVector3 locPlusCrossMinus1 = locPiPlusP3_1.Cross(locPiMinusP3_1);		

  //Decay matrix element squared
  double loclambda1 = 4/3. * fabs(locPlusCrossMinus1.Dot(locPlusCrossMinus1))/TMath::Power((1/9. * loc3Pi1P4_Measured.M2() - locPi01P4_Measured.M2()), 2.);

  //Boost Pi+ and Pi- to the 3Pi rest frame
  //Sum the 3Pi 4-Momenta
  TLorentzVector loc3PiP4_2 = *locPiPlusP4_Measured + *locPiMinusP4_Measured + locPi02P4_Measured;
  TVector3 b2 = loc3PiP4_2.BoostVector();
  TLorentzVector locPiPlusP4_Boosted2 = *locPiPlusP4_Measured;
  locPiPlusP4_Boosted2.Boost(-1*b2);
  TLorentzVector locPiMinusP4_Boosted2 = *locPiMinusP4_Measured;
  locPiMinusP4_Boosted2.Boost(-1*b2);

  //Extract 3-Momenta from 4-Momenta
  TVector3 locPiPlusP3_2 = locPiPlusP4_Boosted2.Vect();
  TVector3 locPiMinusP3_2 = locPiMinusP4_Boosted2.Vect();
  TVector3 locPlusCrossMinus2 = locPiPlusP3_2.Cross(locPiMinusP3_2);		

  //Decay matrix element squared
  double loclambda2 = 4/3. * fabs(locPlusCrossMinus2.Dot(locPlusCrossMinus2))/TMath::Power((1/9. * loc3Pi2P4_Measured.M2() - locPi02P4_Measured.M2()), 2.);
  
  //Histogram lambda at the mass peak
  double weightpeak1;
  if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weightpeak1 = +1;
  else weightpeak1 = 0;
  dHist_lambda_peak->Fill(loclambda1, weightpeak1);
  double weightpeak2;
  if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weightpeak2 = +1;
  else weightpeak2 = 0;
  dHist_lambda_peak->Fill(loclambda2, weightpeak2);

  //Histogram lambda at the mass wings
  double weightwings1;
  if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weightwings1 = +1;
  else weightwings1 = 0;
  dHist_lambda_wings->Fill(loclambda1, weightwings1);
  double weightwings2;
  if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weightwings2 = +1;
  else weightwings2 = 0;
  dHist_lambda_wings->Fill(loclambda2, weightwings2);

  //Histogram lambda without mass cuts
  dHist_lambda_uncut->Fill(loclambda1);
  dHist_lambda_uncut->Fill(loclambda2);

  /******************************************** HISTOGRAM 4PI MASS SPECTRUM ***************************************************/
  double loc4PiMass = loc4PiP4_Measured.M();
  dHist_4PiMass->Fill(loc4PiMass);
  

  /*********************************** HISTOGRAM OMEGAPI0 MASS SPECTRUM **************************************************/
  //OmegaPi mass is 4Pi mass, with mass cuts, BUT: The code below was written earlier and is sloppy
  double locOmegaPiMass = loc4PiP4_Measured.M();
  dHist_OmegaPiMass->Fill(locOmegaPiMass, weight1);
  dHist_OmegaPiMass->Fill(locOmegaPiMass, weight2);

  /************************************* HISTOGRAM 3PION MASS VS 4PION MASS ***************************/
  //3 Pion and 4 Pion masses have been declared previously
  dHist_3vs4->Fill(loc3PiMass1_Measured, loc4PiMass);
  dHist_3vs4->Fill(loc3PiMass2_Measured, loc4PiMass);

  /************************************* HISTOGRAM 4-MOMENTUM TRANSFER SQUARED (t) ********************/
  double locMan_t = fabs(locSqrt_t.Dot(locSqrt_t));
  dHist_Man_t->Fill(locMan_t, weight1);
  dHist_Man_t->Fill(locMan_t, weight2);


	/********************************************* HISTOGRAM COS(THETA)  *****************************************************/
  //need proton momentum in omega rest frame
  TLorentzVector locProtonP4_OmegaPiCM = *locProtonP4_Measured;
  locProtonP4_OmegaPiCM.Boost(-1.0*locGammapBoost);
  locProtonP4_OmegaPiCM.Boost(-1.0*locOmegapiBoost);
  TVector3 locProtonP3_OmegaPiCM = locProtonP4_OmegaPiCM.Vect();

  //Boost helicity unit vectors to omega rf
  TVector3 locHelicityZAxis_OmegaPiCM = -1.0*locProtonP3_OmegaPiCM.Unit();
  TVector3 locHelicityYAxis_OmegaPiCM = -1.0*locBeamP4_gpRest.Vect().Cross(locProtonP3_OmegaPiCM).Unit();
  TVector3 locHelicityXAxis_OmegaPiCM = locHelicityYAxis_OmegaPiCM.Cross(locHelicityZAxis_OmegaPiCM).Unit();
		
  //Project the omega momentum onto these axes and read off the angles
  TVector3 locOmega1P3_Angles(locOmega1P3.Dot(locHelicityXAxis_OmegaPiCM),locOmega1P3.Dot(locHelicityYAxis_OmegaPiCM),locOmega1P3.Dot(locHelicityZAxis_OmegaPiCM));
  TVector3 locOmega2P3_Angles(locOmega2P3.Dot(locHelicityXAxis_OmegaPiCM),locOmega2P3.Dot(locHelicityYAxis_OmegaPiCM),locOmega2P3.Dot(locHelicityZAxis_OmegaPiCM));
  double loccostheta1 = locOmega1P3_Angles.CosTheta();
  double loccostheta2 = locOmega2P3_Angles.CosTheta();

  dHist_costheta->Fill(loccostheta1, weight1);
  dHist_costheta->Fill(loccostheta2, weight2);

  /*********************************************** HISTOGRAM PHI *************************************************************/
		
  double locphi1 = TMath::ATan2(locOmega1P3.Dot(locy), locOmega1P3.Dot(locx));
  double locphi2 = TMath::ATan2(locOmega2P3.Dot(locy), locOmega2P3.Dot(locx));
  dHist_phi->Fill(locphi1, weight1);
  dHist_phi->Fill(locphi2, weight2);
		
  /****************************************** HISTOGRAM COS(THETA_H) *****************************************************/
  double loccosthetaH1 = locnormal1.Dot(loczH1);
  double loccosthetaH2 = locnormal2.Dot(loczH2);
  dHist_costhetaH->Fill(loccosthetaH1, weight1);
  dHist_costhetaH->Fill(loccosthetaH2, weight2);



  /*********************************************** HISTOGRAM PHI_H *************************************************************/
  double locphiH1 = TMath::ATan2(locnormal1.Dot(locyH1), locnormal1.Dot(locxH1));
  double locphiH2 = TMath::ATan2(locnormal2.Dot(locyH2), locnormal2.Dot(locxH2));
  dHist_phiH->Fill(locphiH1, weight1);
  dHist_phiH->Fill(locphiH2, weight2);

  /************************************************* HISTOGRAM COS(THETA_H) IN MASS BINS ****************************************************************/
  //loccosthetaH and locOmegaPiMass have already been declared
  dHist_CosThetaHVsMass->Fill(loccosthetaH1,locOmegaPiMass, weight1);
  dHist_CosThetaHVsMass->Fill(loccosthetaH2,locOmegaPiMass, weight2);
	
  /************************************************* HISTOGRAM PHI_H IN MASS BINS ********************************************************/
  //locphiH and locOmegaPiMass have already been declared
  dHist_PhiHVsMass->Fill(locphiH1,locOmegaPiMass, weight1);
  dHist_PhiHVsMass->Fill(locphiH2,locOmegaPiMass, weight2);


  /************************************************* HISTOGRAM COS(THETA) IN MASS BINS ****************************************************************/
  //loccostheta and locOmegaPiMass have already been declared
  dHist_CosThetaVsMass->Fill(loccostheta1,locOmegaPiMass, weight1);
  dHist_CosThetaVsMass->Fill(loccostheta2,locOmegaPiMass, weight2);

  /************************************************* HISTOGRAM PHI IN MASS BINS ********************************************************/
  //locphi and locOmegaPiMass have already been declared
  dHist_PhiVsMass->Fill(locphi1,locOmegaPiMass, weight1);
  dHist_PhiVsMass->Fill(locphi2,locOmegaPiMass, weight2);		 

  /***************************************** HISTOGRAM COS(THETA) IN t [0.1, 0.3] BIN ***********************************************/
  //loccostheta, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      //Apply omega mass cuts and histogram results
      dHist_CosTheta_t1->Fill(loccostheta1, weight1);
      dHist_CosTheta_t1->Fill(loccostheta2, weight2);
    }
		  
  /***************************************** HISTOGRAM COS(THETA) IN t [0.3, 1.0] BIN ***********************************************/
  //loccostheta, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      //Apply omega mass cuts and histogram results
      dHist_CosTheta_t2->Fill(loccostheta1, weight1);
      dHist_CosTheta_t2->Fill(loccostheta2, weight2);
    }		 

  /***************************************** HISTOGRAM PHI IN t [0.1, 0.3] BIN ***********************************************/
  //locphi, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      //Apply omega mass cuts and histogram results
      dHist_Phi_t1->Fill(locphi1, weight1);
      dHist_Phi_t1->Fill(locphi2, weight2);
    }		  

  /***************************************** HISTOGRAM PHI IN t [0.3, 1.0] BIN ***********************************************/
  //locphi, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      //Apply omega mass cuts and histogram results
      dHist_Phi_t2->Fill(locphi1, weight1);
      dHist_Phi_t2->Fill(locphi2, weight2);
    }
		  
  /***************************************** HISTOGRAM COS(THETA_H) IN t [0.1, 0.3] BIN ***********************************************/
  //loccosthetaH, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      //Apply omega mass cuts and histogram results
      dHist_CosThetaH_t1->Fill(loccosthetaH1, weight1);
      dHist_CosThetaH_t1->Fill(loccosthetaH2, weight2);
    }

  /***************************************** HISTOGRAM COS(THETA_H) IN t [0.3, 1.0] BIN ***********************************************/
  //loccosthetaH, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      //Apply omega mass cuts and histogram results
      dHist_CosThetaH_t2->Fill(loccosthetaH1, weight1);
      dHist_CosThetaH_t2->Fill(loccosthetaH2, weight2);
    }
		  

  /***************************************** HISTOGRAM PHI_H IN t [0.1, 0.3] BIN ***********************************************/
  //locphiH, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      //Apply omega mass cuts and histogram results
      dHist_PhiH_t1->Fill(locphiH1, weight1);
      dHist_PhiH_t1->Fill(locphiH2, weight2);
    }
		  

  /***************************************** HISTOGRAM PHI_H IN t [0.3, 1.0] BIN ***********************************************/
  //locphiH, locMan_t, etc have been declared earlier
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      //Apply omega mass cuts and histogram results
      dHist_PhiH_t2->Fill(locphiH1, weight1);
      dHist_PhiH_t2->Fill(locphiH2, weight2);
    }
  
  /**************************************** HISTOGRAM ANGLES VS MASS IN t[0.1, 0.3] ******************************************/
  /************************************************* HISTOGRAM COS(THETA_H) IN MASS BINS ****************************************************************/
  //loccosthetaH and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      //apply omega mass cuts
      dHist_CosThetaHVsMass_t1->Fill(loccosthetaH1,locOmegaPiMass, weight1);
      dHist_CosThetaHVsMass_t1->Fill(loccosthetaH2,locOmegaPiMass, weight2);
    }

  /************************************************* HISTOGRAM PHI_H IN MASS BINS ********************************************************/
  //locphiH and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      dHist_PhiHVsMass_t1->Fill(locphiH1,locOmegaPiMass, weight1);
      dHist_PhiHVsMass_t1->Fill(locphiH2,locOmegaPiMass, weight2);
    }

  /************************************************* HISTOGRAM COS(THETA) IN MASS BINS ****************************************************************/
  //loccostheta and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      dHist_CosThetaVsMass_t1->Fill(loccostheta1,locOmegaPiMass, weight1);
      dHist_CosThetaVsMass_t1->Fill(loccostheta2,locOmegaPiMass, weight2);
    }

  /************************************************* HISTOGRAM PHI IN MASS BINS ********************************************************/
  //locphi and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.1 && locMan_t < 0.3)
    {
      dHist_PhiVsMass_t1->Fill(locphi1,locOmegaPiMass, weight1);
      dHist_PhiVsMass_t1->Fill(locphi2,locOmegaPiMass, weight2);
    }

  /**************************************** HISTOGRAM ANGLES VS MASS IN t[0.3, 1.0] ******************************************/
  /************************************************* HISTOGRAM COS(THETA_H) IN MASS BINS ****************************************************************/
  //loccosthetaH and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      dHist_CosThetaHVsMass_t2->Fill(loccosthetaH1,locOmegaPiMass, weight1);
      dHist_CosThetaHVsMass_t2->Fill(loccosthetaH2,locOmegaPiMass, weight2);
    }

  /************************************************* HISTOGRAM PHI_H IN MASS BINS ********************************************************/
  //locphiH and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      dHist_PhiHVsMass_t2->Fill(locphiH1,locOmegaPiMass, weight1);
      dHist_PhiHVsMass_t2->Fill(locphiH2,locOmegaPiMass, weight2);
    }

  /************************************************* HISTOGRAM COS(THETA) IN MASS BINS ****************************************************************/
  //loccostheta and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      dHist_CosThetaVsMass_t2->Fill(loccostheta1,locOmegaPiMass, weight1);
      dHist_CosThetaVsMass_t2->Fill(loccostheta2,locOmegaPiMass, weight2);
    }

  /************************************************* HISTOGRAM PHI IN MASS BINS ********************************************************/
  //locphi and locOmegaPiMass have already been declared
  //Make Mandelstam t cuts
  if(locMan_t > 0.3 && locMan_t < 1.0)
    {
      dHist_PhiVsMass_t2->Fill(locphi1,locOmegaPiMass, weight1);
      dHist_PhiVsMass_t2->Fill(locphi2,locOmegaPiMass, weight2);
    }





  return kTRUE;
}

void Thrown_534::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Thrown_534::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  dHist_BeamEnergy->Write();
  dHist_MissingMassSquared->Write();
  dHist_3PiMass_Measured->Write();
  dHist_MM2_Weighted->Write();
  dHist_lambda_peak->Write();
  dHist_lambda_wings->Write();
  dHist_lambda_uncut->Write();
  dHist_4PiMass->Write();
  dHist_OmegaPiMass->Write();
  dHist_3vs4->Write();
  dHist_Man_t->Write();
  dHist_costheta->Write();
  dHist_phi->Write();
  dHist_costhetaH->Write();
  dHist_phiH->Write();
  dHist_CosThetaVsMass->Write();
  dHist_PhiVsMass->Write();
  dHist_CosThetaHVsMass->Write();
  dHist_PhiHVsMass->Write();
  dHist_CosTheta_t1->Write();
  dHist_Phi_t1->Write();
  dHist_CosThetaH_t1->Write();
  dHist_PhiH_t1->Write();
  dHist_CosTheta_t2->Write();
  dHist_Phi_t2->Write();
  dHist_CosThetaH_t2->Write();
  dHist_PhiH_t2->Write();
  dHist_CosThetaVsMass_t1->Write();
  dHist_PhiVsMass_t1->Write();
  dHist_CosThetaHVsMass_t1->Write();
  dHist_PhiHVsMass_t1->Write();
  dHist_CosThetaVsMass_t2->Write();
  dHist_PhiVsMass_t2->Write();
  dHist_CosThetaHVsMass_t2->Write();
  dHist_PhiHVsMass_t2->Write();

  outfile->Write();
  outfile->Close();

}
