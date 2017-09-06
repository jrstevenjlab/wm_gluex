#include "DSelector_pomegapi.h"

void DSelector_pomegapi::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "pomegapi.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none

	//DO THIS NEXT
	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	//THEN THIS
	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_dEdxProton(dComboWrapper, false, Proton, SYS_CDC));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, Proton, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, Proton, SYS_BCAL));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiPlus, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, PiPlus, SYS_BCAL));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiMinus, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, PiMinus, SYS_BCAL));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, Gamma, SYS_BCAL));

        //MASSES
        deque<Particle_t> locPi0PIDs;  locPi0PIDs.push_back(Gamma); locPi0PIDs.push_back(Gamma);
        deque<Particle_t> locOmegaPIDs;  locOmegaPIDs.push_back(Pi0); locOmegaPIDs.push_back(PiPlus); locOmegaPIDs.push_back(PiMinus);
        deque<Particle_t> locb1PIDs;  locb1PIDs.push_back(PiPlus); locb1PIDs.push_back(PiMinus); locb1PIDs.push_back(Pi0); locb1PIDs.push_back(Pi0);
        double minPi0 = 0.09; double maxPi0 = 0.18;
        double minOmega = 0.5; double maxOmega = 1.1;
        double minb1 = 0.5; double maxb1 = 2.0;
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Pi0, 100, minPi0, maxPi0, "Pi0_NoCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 100, minOmega, maxOmega, "omega_NoCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PIDs, 150, minb1, maxb1, "b1_NoCut"));
        //dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, false, 1, locPi0PIDs, 0, locOmegaPIDs, 100, minPi0, maxPi0, 100, minOmega, maxOmega, "Omega_Pi0_NoCut"));
        //dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 0, locb1PIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1_NoCut"));
        dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

        //KINFIT RESULTS
        dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

        //CUT MISSING MASS
        dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));
        dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.00001));

        //MASSES (after KINFIT CL cut)
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Pi0, 100, minPi0, maxPi0, "Pi0_KinCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 100, minOmega, maxOmega, "omega_KinCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PIDs, 150, minb1, maxb1, "b1_KinCut"));
        //dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, false, 1, locPi0PIDs, 0, locOmegaPIDs, 100, minPi0, maxPi0, 100, minOmega, maxOmega, "Omega_Pi0_KinCut"));
        //dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 0, locb1PIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1_KinCut"));

        //CUT PI0 MASS
        dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, false, Pi0, 0.114, 0.156));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 100, minOmega, maxOmega, "omega_Pi0Cut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PIDs, 150, minb1, maxb1, "b1_Pi0Cut"));
        //dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 0, locb1PIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1_Pi0Cut"));

	//CUT OMEGA MASS
        dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 0.76, 0.81));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PIDs, 150, minb1, maxb1, "b1_OmegaCut"));

        //CUT BEAM ENERGY
        dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, true));
        dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, true, 7.0, 12.0));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PIDs, 150, minb1, maxb1, "b1_BeamCut"));

        //KINEMATICS
        dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, true));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
	dHist_3PiMass_Measured = new TH1I("3PiMass_Measured", ";3 Pion Mass (GeV/c^{2})", 600, 0.6, 1.1);
	dHist_MM2_Weighted = new TH1I("MM2_Weighted", ";Weighted Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_lambda_peak = new TH1I("lambda_peak", ";{lambda} (peak)", 600, 0, 1);
       	dHist_lambda_wings = new TH1I("lambda_wings", ";{lambda} (wings)", 600, 0, 1);
       	dHist_lambda_uncut = new TH1I("lambda_uncut", ";{lambda} (uncut)", 600, 0, 1);
	dHist_4PiMass = new TH1I("4PiMass", ";Pi+Pi-Pi0Pi0 Mass (GeV/c^{2})", 600, 0, 3);
	dHist_OmegaPiMass = new TH1I("OmegaPiMass", ";Omega Pi0 Mass (GeV)", 600, 0, 3);
	dHist_3vs4 = new TH2I("3vs4", ";3 Pion Mass vs 4 Pion Mass", 600, 0, 3, 600, 0, 3);
	dHist_Man_t = new TH1I("Man_t", ";Four-Momentum Transfer Squared (GeV)^{2}", 600, 0.0, 1.0);
	//Decay Angles
	dHist_costheta = new TH1I("costheta", ";Cos(theta)", 600, -1.0, 1.0);
	dHist_phi = new TH1I("phi", ";phi (radians)", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
	dHist_costhetaH = new TH1I("costhetaH", ";Cos(theta_H)", 600, -1.0, 1.0);
	dHist_phiH = new TH1I("phiH", ";phi_H (rad)", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
	//Angles vs mass
	dHist_CosThetaVsMass = new TH2I("CosThetaVsMass", ";Cos(theta) vs Omega Pi Mass", 600, -1.0, 1.0, 20, 1.0, 3.0);
	dHist_PhiVsMass = new TH2I("PhiVsMass", ";Phi vs Omega Pi Mass", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
	dHist_CosThetaHVsMass = new TH2I("CosThetaHVsMass", ";Cos(theta_H) vs Omega Pi Mass", 600, -1.0, 1.0, 20, 1.0, 3.0);
	dHist_PhiHVsMass = new TH2I("PhiHVsMass", ";Phi vs Omega Pi Mass", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
	//Angles in t range [0.1, 0.3]
	dHist_CosTheta_t1 = new TH1I("CosTheta_t1", ";Cos(theta) with t[0.1, 0.3]", 600, -1.0, 1.0);
	dHist_Phi_t1 = new TH1I("Phi_t1", ";Phi with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
	dHist_CosThetaH_t1 = new TH1I("CosThetaH_t1", ";Cos(theta_H) with t[0.1, 0.3]", 600, -1.0, 1.0);
	dHist_PhiH_t1 = new TH1I("PhiH_t1", ";Phi_H with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
	//Angles in t range [0.3, 1.0]
	dHist_CosTheta_t2 = new TH1I("CosTheta_t2", ";Cos(theta) with t[0.3, 1.0]", 600, -1.0, 1.0);
	dHist_Phi_t2 = new TH1I("Phi_t2", ";Phi with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
	dHist_CosThetaH_t2 = new TH1I("CosThetaH_t2", ";Cos(theta_H) with t[0.3, 1.0]", 600, -1.0, 1.0);
	dHist_PhiH_t2 = new TH1I("PhiH_t2", ";Phi_H with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi());
	//Angles vs mass in t range [0.1, 0.3]
	dHist_CosThetaVsMass_t1 = new TH2I("CosThetaVsMass_t1", ";Cos(theta) vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0, 1.0, 20, 1.0, 3.0);
	dHist_PhiVsMass_t1 = new TH2I("PhiVsMass_t1", ";Phi vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
	dHist_CosThetaHVsMass_t1 = new TH2I("CosThetaHVsMass_t1", ";Cos(theta_H) vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0, 1.0, 20, 1.0, 3.0);
	dHist_PhiHVsMass_t1 = new TH2I("PhiHVsMass_t1", ";Phi vs Omega Pi Mass with t[0.1, 0.3]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
	//Angles vs mass in t range [0.3, 1.0]
	dHist_CosThetaVsMass_t2 = new TH2I("CosThetaVsMass_t2", ";Cos(theta) vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0, 1.0, 20, 1.0, 3.0);
	dHist_PhiVsMass_t2 = new TH2I("PhiVsMass_t2", ";Phi vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
	dHist_CosThetaHVsMass_t2 = new TH2I("CosThetaHVsMass_t2", ";Cos(theta_H) vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0, 1.0, 20, 1.0, 3.0);
	dHist_PhiHVsMass_t2 = new TH2I("PhiHVsMass_t2", ";Phi_H vs Omega Pi Mass with t[0.3, 1.0]", 600, -1.0*TMath::Pi(), 1.0*TMath::Pi(), 20, 1.0, 3.0);
	

	/***************************************** ADVANCED: CHOOSE BRANCHES TO READ ****************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_pomegapi::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_3PiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MM2_Weighted;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_lambda_peak;
       	set<map<Particle_t, set<Int_t> > > locUsedSoFar_lambda_wings;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_lambda_uncut;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_4PiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaPiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_3vs4;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Man_t;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_costheta;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_phi;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_costhetaH;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_phiH;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaVsMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiVsMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaHVsMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiHVsMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosTheta_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Phi_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaH_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiH_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosTheta_t2;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Phi_t2;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaH_t2;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiH_t2;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaVsMass_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiVsMass_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaHVsMass_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiHVsMass_t1;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaVsMass_t2;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiVsMass_t2;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CosThetaHVsMass_t2;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PhiHVsMass_t2;


	

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		//Step 2
		Int_t locPhoton3NeutralID = dPhoton3Wrapper->Get_NeutralID();
		Int_t locPhoton4NeutralID = dPhoton4Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locDecayingPi01P4 = dDecayingPi01Wrapper->Get_P4();
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
		//Step 2
		TLorentzVector locDecayingPi02P4 = dDecayingPi02Wrapper->Get_P4();
		TLorentzVector locPhoton3P4 = dPhoton3Wrapper->Get_P4();
		TLorentzVector locPhoton4P4 = dPhoton4Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();
		//Step 2
		TLorentzVector locPhoton3P4_Measured = dPhoton3Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton4P4_Measured = dPhoton4Wrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE
		TLorentzVector locPi01P4_Measured = locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector loc3Pi1P4_Measured = locPiPlusP4_Measured + locPiMinusP4_Measured + locPi01P4_Measured;
		TLorentzVector locPi02P4_Measured = locPhoton3P4_Measured + locPhoton4P4_Measured;
		TLorentzVector loc3Pi2P4_Measured = locPiPlusP4_Measured + locPiMinusP4_Measured + locPi02P4_Measured;
		TLorentzVector loc4PiP4_Measured = locPiPlusP4_Measured + locPiMinusP4_Measured + locPi01P4_Measured + locPi02P4_Measured;
		TLorentzVector locSqrt_t = loc4PiP4_Measured - locBeamP4_Measured;

		//Define omega-pi0 decay angles
		//Boost to gamma-p rest frame
		TLorentzVector locGammapP4 = locBeamP4_Measured + dTargetP4; 
		TVector3 locGammapBoost = locGammapP4.BoostVector();     //create boost vector from gamma-p rest frame to lab frame - use negative of this to boost to gamma-p rest frame
		TLorentzVector locBeamP4_gpRest = locBeamP4_Measured;     //boost the beam to the gamma-p rest frame
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
		TLorentzVector locPiPlusP4_o1Rest = locPiPlusP4_Measured;
		locPiPlusP4_o1Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
		locPiPlusP4_o1Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
		locPiPlusP4_o1Rest.Boost(-1.0*locOmega1Boost);     //Boost to omega rf
		TVector3 locPiPlusP3_o1Rest = locPiPlusP4_o1Rest.Vect();
		TLorentzVector locPiMinusP4_o1Rest = locPiMinusP4_Measured;
		locPiMinusP4_o1Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
		locPiMinusP4_o1Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
		locPiMinusP4_o1Rest.Boost(-1.0*locOmega1Boost);     //Boost to omega rf
		TVector3 locPiMinusP3_o1Rest = locPiMinusP4_o1Rest.Vect();
		//Second Omega combo
		TVector3 locOmega2Boost = locOmega2P4_opiRest.BoostVector();   //create boost vector from omega rest frame to omega-pi rest frame
		TLorentzVector locPiPlusP4_o2Rest = locPiPlusP4_Measured;
		locPiPlusP4_o2Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
		locPiPlusP4_o2Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
		locPiPlusP4_o2Rest.Boost(-1.0*locOmega2Boost);     //Boost to omega rf
		TVector3 locPiPlusP3_o2Rest = locPiPlusP4_o2Rest.Vect();
		TLorentzVector locPiMinusP4_o2Rest = locPiMinusP4_Measured;
		locPiMinusP4_o2Rest.Boost(-1.0*locGammapBoost);    //Boost to gamma-p rf		
		locPiMinusP4_o2Rest.Boost(-1.0*locOmegapiBoost);   //Boost to Omega-pi rf
		locPiMinusP4_o2Rest.Boost(-1.0*locOmega2Boost);     //Boost to omega rf
		TVector3 locPiMinusP3_o2Rest = locPiMinusP4_o2Rest.Vect();

		//Normal vector to decay plane
		TVector3 locnormal1 = locPiPlusP3_o1Rest.Cross(locPiMinusP3_o1Rest).Unit();
		TVector3 locnormal2 = locPiPlusP3_o2Rest.Cross(locPiMinusP3_o2Rest).Unit();
		

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured + locPhoton4P4_Measured;

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E());
			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}
		/******************************************* HISTOGRAM 3PI MASS *****************************************************/
		double loc3PiMass1_Measured = loc3Pi1P4_Measured.M();
		double loc3PiMass2_Measured = loc3Pi2P4_Measured.M();
	
		//Uniqueness Tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_3PiMass;
		locUsedThisCombo_3PiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_3PiMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_3PiMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_3PiMass[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_3PiMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_3PiMass[PiMinus].insert(locPiMinusTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_3PiMass.find(locUsedThisCombo_3PiMass) == locUsedSoFar_3PiMass.end())
		  {
		    //unique missing mass combo: histogram it, and register this combo of particles
		    dHist_3PiMass_Measured->Fill(loc3PiMass1_Measured);
		    dHist_3PiMass_Measured->Fill(loc3PiMass2_Measured);
		    locUsedSoFar_3PiMass.insert(locUsedThisCombo_3PiMass);
		  }

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton4NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared);
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}

		/*************************************** HISTOGRAM WEIGHTED MISSING MASS SQUARED *********************************/
	       	//Missing Mass Squared
		double locMM2_Weighted = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
		map<Particle_t, set<Int_t> > locUsedThisCombo_MM2_Weighted;
		locUsedThisCombo_MM2_Weighted[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MM2_Weighted[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MM2_Weighted[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_MM2_Weighted[Proton].insert(locProtonTrackID);
		locUsedThisCombo_MM2_Weighted[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_MM2_Weighted[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_MM2_Weighted[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_MM2_Weighted[Gamma].insert(locPhoton4NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_MM2_Weighted.find(locUsedThisCombo_MM2_Weighted) == locUsedSoFar_MM2_Weighted.end())
		{
		  //unique missing mass combo: histogram it, and register this combo of particles
		  //Assign weight for first 3Pi combo
		  double weight1;
		  if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
		  else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
		  else weight1 = 0;
		   dHist_MM2_Weighted->Fill(locMM2_Weighted, weight1);
		  
		  //Assign weight for second 3Pi combo
		  double weight2;
		  if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
		  else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
		  else weight2 = 0;
		   dHist_MM2_Weighted->Fill(locMM2_Weighted, weight2);
		   locUsedSoFar_MM2_Weighted.insert(locUsedThisCombo_MM2_Weighted);
		}


		/**************************************** HISTOGRAM DECAY MATRIX ELEMENT SQUARED (LAMBDA)(PEAK) ************************/
		//Boost Pi+ and Pi- to the 3Pi rest frame
		//Sum the 3Pi 4-Momenta
		TLorentzVector loc3PiP4_1 = locPiPlusP4_Measured + locPiMinusP4_Measured + locPi01P4_Measured;
		TVector3 b1 = loc3PiP4_1.BoostVector();
		TLorentzVector locPiPlusP4_Boosted1 = locPiPlusP4_Measured;
		locPiPlusP4_Boosted1.Boost(-1*b1);
		TLorentzVector locPiMinusP4_Boosted1 = locPiMinusP4_Measured;
		locPiMinusP4_Boosted1.Boost(-1*b1);

		//Extract 3-Momenta from 4-Momenta
	       	TVector3 locPiPlusP3_1 = locPiPlusP4_Boosted1.Vect();
		TVector3 locPiMinusP3_1 = locPiMinusP4_Boosted1.Vect();
		TVector3 locPlusCrossMinus1 = locPiPlusP3_1.Cross(locPiMinusP3_1);		

		//Decay matrix element squared
	       	double loclambda_peak1 = 4/3. * fabs(locPlusCrossMinus1.Dot(locPlusCrossMinus1))/powf((1/9. * loc3Pi1P4_Measured.M2() - locPi01P4_Measured.M2()), 2.);

		//Boost Pi+ and Pi- to the 3Pi rest frame
		//Sum the 3Pi 4-Momenta
		TLorentzVector loc3PiP4_2 = locPiPlusP4_Measured + locPiMinusP4_Measured + locPi02P4_Measured;
		TVector3 b2 = loc3PiP4_2.BoostVector();
		TLorentzVector locPiPlusP4_Boosted2 = locPiPlusP4_Measured;
		locPiPlusP4_Boosted2.Boost(-1*b2);
		TLorentzVector locPiMinusP4_Boosted2 = locPiMinusP4_Measured;
		locPiMinusP4_Boosted2.Boost(-1*b2);

		//Extract 3-Momenta from 4-Momenta
	       	TVector3 locPiPlusP3_2 = locPiPlusP4_Boosted2.Vect();
		TVector3 locPiMinusP3_2 = locPiMinusP4_Boosted2.Vect();
		TVector3 locPlusCrossMinus2 = locPiPlusP3_2.Cross(locPiMinusP3_2);		

		//Decay matrix element squared
	       	double loclambda_peak2 = 4/3. * fabs(locPlusCrossMinus2.Dot(locPlusCrossMinus2))/powf((1/9. * loc3Pi2P4_Measured.M2() - locPi02P4_Measured.M2()), 2.);
		
		//Uniqueness Tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_lambda_peak;
		locUsedThisCombo_lambda_peak[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_lambda_peak[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_lambda_peak[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_lambda_peak[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_lambda_peak[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_lambda_peak[Gamma].insert(locPhoton4NeutralID);

		//Compare to what's been used so far
		if(locUsedSoFar_lambda_peak.find(locUsedThisCombo_lambda_peak) == locUsedSoFar_lambda_peak.end())
		  {
		    //unique lambda combo: histogram it, and register this combo of particles
		    double weightpeak1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weightpeak1 = +1;
		    else weightpeak1 = 0;
		    dHist_lambda_peak->Fill(loclambda_peak1, weightpeak1);
		    double weightpeak2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weightpeak2 = +1;
		    else weightpeak2 = 0;
		    dHist_lambda_peak1->Fill(loclambda_peak2, weightpeak2);

		    locUsedSoFar_lambda_peak.insert(locUsedThisCombo_lambda_peak);
		  }
		


	/**************************************** HISTOGRAM DECAY MATRIX ELEMENT SQUARED (LAMBDA)(WINGS, FIRST PI0) ************************/
		//Decay matrix element squared
	       	double loclambda_wings1 = 4/3. * fabs(locPlusCrossMinus1.Dot(locPlusCrossMinus1))/powf((1/9. * loc3Pi1P4_Measured.M2() - locPi01P4_Measured.M2()), 2.);
	       	double loclambda_wings2 = 4/3. * fabs(locPlusCrossMinus2.Dot(locPlusCrossMinus2))/powf((1/9. * loc3Pi2P4_Measured.M2() - locPi02P4_Measured.M2()), 2.);
		
		//Uniqueness Tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_lambda_wings;
		locUsedThisCombo_lambda_wings[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_lambda_wings[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_lambda_wings[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_lambda_wings[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_lambda_wings[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_lambda_wings[Gamma].insert(locPhoton4NeutralID);

		//Compare to what's been used so far
		if(locUsedSoFar_lambda_wings.find(locUsedThisCombo_lambda_wings) == locUsedSoFar_lambda_wings.end())
		  {
		    //unique lambda combo: histogram it, and register this combo of particles
		    double weightwings1;
		    if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weightwings1 = +1;
		    else weightwings1 = 0;
		    dHist_lambda_wings->Fill(loclambda_wings1, weightwings1);
		    double weightwings2;
		    if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weightwings2 = +1;
		    else weightwings2 = 0;
		    dHist_lambda_wings->Fill(loclambda_wings2, weightwings2);

		    locUsedSoFar_lambda_wings.insert(locUsedThisCombo_lambda_wings);
		  }
		
		/*********************************** HISTOGRAM AN UNCUT LAMBDA **************************************************/
		//Decay matrix element squared
	       	double loclambda_uncut = 4/3. * fabs(locPlusCrossMinus2.Dot(locPlusCrossMinus2))/powf((1/9. * loc3Pi2P4_Measured.M2() - locPi02P4_Measured.M2()), 2.);
		
		//Uniqueness Tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_lambda_uncut;
		locUsedThisCombo_lambda_uncut[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_lambda_uncut[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_lambda_uncut[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_lambda_uncut[Gamma].insert(locPhoton4NeutralID);

		//Compare to what's been used so far
		if(locUsedSoFar_lambda_uncut.find(locUsedThisCombo_lambda_uncut) == locUsedSoFar_lambda_uncut.end())
		  {
		    //unique lambda combo: histogram it, and register this combo of particles
		    dHist_lambda_uncut->Fill(loclambda_uncut);
		    locUsedSoFar_lambda_uncut.insert(locUsedThisCombo_lambda_uncut);
		  }

		/*********************************** HISTOGRAM 4PI MASS SPECTRUM **************************************************/
		double loc4PiMass = loc4PiP4_Measured.M();
		
		//Uniqueness tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_4PiMass;
		locUsedThisCombo_4PiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_4PiMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_4PiMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_4PiMass[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_4PiMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_4PiMass[PiMinus].insert(locPiMinusTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_4PiMass.find(locUsedThisCombo_4PiMass) == locUsedSoFar_4PiMass.end())
		  {
		    //unique missing mass combo: histogram it, and register this combo of particles
		    dHist_4PiMass->Fill(loc4PiMass);
		    locUsedSoFar_4PiMass.insert(locUsedThisCombo_4PiMass);
		  }

		/*********************************** HISTOGRAM OMEGAPI0 MASS SPECTRUM **************************************************/
		double locOmegaPiMass = loc4PiP4_Measured.M();
		
		//Uniqueness tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_OmegaPiMass;
		locUsedThisCombo_OmegaPiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_OmegaPiMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_OmegaPiMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_OmegaPiMass[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_OmegaPiMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_OmegaPiMass[PiMinus].insert(locPiMinusTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_OmegaPiMass.find(locUsedThisCombo_OmegaPiMass) == locUsedSoFar_OmegaPiMass.end())
		  {
		    //unique missing mass combo: histogram it, and register this combo of particles
		    //Assign weight for first omega combo
		    double omega1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) omega1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) omega1 = -1;
		    else omega1 = 0;
		    dHist_OmegaPiMass->Fill(locOmegaPiMass, omega1);

		    //Assign weight for second omega combo
		    double omega2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) omega2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) omega2 = -1;
		    else omega2 = 0;
		    dHist_OmegaPiMass->Fill(locOmegaPiMass, omega2);
		    locUsedSoFar_OmegaPiMass.insert(locUsedThisCombo_OmegaPiMass);
		  }
		
		/************************************* HISTOGRAM 3PION MASS VS 4PION MASS ***************************/
		//I think we've already declared loc4PiMass etc. Let's see if it works
		//Uniqueness tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_3vs4;
		locUsedThisCombo_3vs4[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_3vs4[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_3vs4[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_3vs4[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_3vs4[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_3vs4[PiMinus].insert(locPiMinusTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_3vs4.find(locUsedThisCombo_3vs4) == locUsedSoFar_3vs4.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    dHist_3vs4->Fill(loc3PiMass1_Measured, loc4PiMass);
		    locUsedSoFar_3vs4.insert(locUsedThisCombo_3vs4);
		  }
		
		/************************************* HISTOGRAM 4-MOMENTUM TRANSFER SQUARED (t) ********************/
		double locMan_t = fabs(locSqrt_t.Dot(locSqrt_t));

		//Uniqueness tracking:
		map<Particle_t, set<Int_t> > locUsedThisCombo_Man_t;
		locUsedThisCombo_Man_t[Unknown].insert(locBeamID);
		locUsedThisCombo_Man_t[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_Man_t[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_Man_t[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_Man_t[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_Man_t[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_Man_t[Gamma].insert(locPhoton4NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_Man_t.find(locUsedThisCombo_Man_t) == locUsedSoFar_Man_t.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //Assign weight for first 3Pi combo
		    double weightt1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weightt1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weightt1 = -1;
		    else weightt1 = 0;
		    dHist_Man_t->Fill(locMan_t, weightt1);
		  
		    //Assign weight for second 3Pi combo
		    double weightt2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weightt2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weightt2 = -1;
		    else weightt2 = 0;
		    dHist_Man_t->Fill(locMan_t, weightt2);
		    locUsedSoFar_Man_t.insert(locUsedThisCombo_Man_t);
		  }

		/********************************************* HISTOGRAM COS(THETA)  *****************************************************/
		//double loccostheta1 = (locOmega1P3.Dot(locz)) / (locOmega1P3.Mag());
		//double loccostheta2 = (locOmega2P3.Dot(locz)) / (locOmega2P3.Mag());

		//need proton momentum in omega rest frame
		TLorentzVector locProtonP4_OmegaPiCM = locProtonP4_Measured;
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
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_costheta;
		locUsedThisCombo_costheta[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_costheta[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_costheta[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_costheta[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_costheta[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_costheta[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_costheta[Proton].insert(locProtonTrackID);
		locUsedThisCombo_costheta[Unknown].insert(locBeamID);	

		//compare to what's been used so far
		if(locUsedSoFar_costheta.find(locUsedThisCombo_costheta) == locUsedSoFar_costheta.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //Apply omega mass cuts
		    double thetaweight1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) thetaweight1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) thetaweight1 = -1;
		    else thetaweight1 = 0;
		    dHist_costheta->Fill(loccostheta1, thetaweight1);
		    double thetaweight2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) thetaweight2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) thetaweight2 = -1;
		    else thetaweight2 = 0;
		    dHist_costheta->Fill(loccostheta2, thetaweight2);
		    locUsedSoFar_costheta.insert(locUsedThisCombo_costheta);
		  }

		/*********************************************** HISTOGRAM PHI *************************************************************/
		
		double locphi1 = TMath::ATan2(locOmega1P3.Dot(locy), locOmega1P3.Dot(locx));
		double locphi2 = TMath::ATan2(locOmega2P3.Dot(locy), locOmega2P3.Dot(locx));
		

		//Uniqueness tracking
	       	map<Particle_t, set<Int_t> > locUsedThisCombo_phi;
		locUsedThisCombo_phi[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_phi[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_phi[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_phi[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_phi[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_phi[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_phi[Proton].insert(locProtonTrackID);
		locUsedThisCombo_phi[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_phi.find(locUsedThisCombo_phi) == locUsedSoFar_phi.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //Apply omega mass cuts
		    double phiweight1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) phiweight1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) phiweight1 = -1;
		    else phiweight1 = 0;
		       dHist_phi->Fill(locphi1, phiweight1);
		    //Apply omega mass cuts
		    double phiweight2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) phiweight2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) phiweight2 = -1;
		    else phiweight2 = 0;
		       dHist_phi->Fill(locphi2, phiweight2);
		    locUsedSoFar_phi.insert(locUsedThisCombo_phi);
		  }
		
		
		/****************************************** HISTOGRAM COS(THETA_H) *****************************************************/
		double loccosthetaH1 = locnormal1.Dot(loczH1);
		double loccosthetaH2 = locnormal2.Dot(loczH2);
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_costhetaH;
		locUsedThisCombo_costhetaH[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_costhetaH[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_costhetaH[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_costhetaH[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_costhetaH[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_costhetaH[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_costhetaH[Proton].insert(locProtonTrackID);
		locUsedThisCombo_costhetaH[Unknown].insert(locBeamID);	

		//compare to what's been used so far
		if(locUsedSoFar_costhetaH.find(locUsedThisCombo_costhetaH) == locUsedSoFar_costhetaH.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //Apply omega mass cuts
		    double thetaHweight1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) thetaHweight1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) thetaHweight1 = -1;
		    else thetaHweight1 = 0;
		       dHist_costhetaH->Fill(loccosthetaH1, thetaHweight1);
		    //Apply omega mass cuts
		    double thetaHweight2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) thetaHweight2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) thetaHweight2 = -1;
		    else thetaHweight2 = 0;
		       dHist_costhetaH->Fill(loccosthetaH2, thetaHweight2);
		    locUsedSoFar_costhetaH.insert(locUsedThisCombo_costhetaH);
		  }		
		
		/*********************************************** HISTOGRAM PHI_H *************************************************************/
		double locphiH1 = TMath::ATan2(locnormal1.Dot(locyH1), locnormal1.Dot(locxH1));
		double locphiH2 = TMath::ATan2(locnormal2.Dot(locyH2), locnormal2.Dot(locxH2));

		//Uniqueness tracking
	       	map<Particle_t, set<Int_t> > locUsedThisCombo_phiH;
		locUsedThisCombo_phiH[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_phiH[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_phiH[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_phiH[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_phiH[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_phiH[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_phiH[Proton].insert(locProtonTrackID);
		locUsedThisCombo_phiH[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_phiH.find(locUsedThisCombo_phiH) == locUsedSoFar_phiH.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //Apply omega mass cuts
		    double phiHweight1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) phiHweight1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) phiHweight1 = -1;
		    else phiHweight1 = 0;
		       dHist_phiH->Fill(locphiH1, phiHweight1);
		    double phiHweight2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) phiHweight2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) phiHweight2 = -1;
		    else phiHweight2 = 0;
		       dHist_phiH->Fill(locphiH2, phiHweight2);
		    locUsedSoFar_phiH.insert(locUsedThisCombo_phiH);
		  }
		

		/************************************************* HISTOGRAM COS(THETA_H) IN MASS BINS ****************************************************************/
		//loccosthetaH and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaHVsMass;
		locUsedThisCombo_CosThetaHVsMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaHVsMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaHVsMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaHVsMass[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaHVsMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaHVsMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaHVsMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaHVsMass[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaHVsMass.find(locUsedThisCombo_CosThetaHVsMass) == locUsedSoFar_CosThetaHVsMass.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //apply omega mass cuts
		    double wH1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) wH1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) wH1 = -1;
		    else wH1 = 0;
		    dHist_CosThetaHVsMass->Fill(loccosthetaH1,locOmegaPiMass, wH1);
		    double wH2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) wH2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) wH2 = -1;
		    else wH2 = 0;
		    dHist_CosThetaHVsMass->Fill(loccosthetaH2,locOmegaPiMass, wH2);
		      locUsedSoFar_CosThetaHVsMass.insert(locUsedThisCombo_CosThetaHVsMass);
		  }

		/************************************************* HISTOGRAM PHI_H IN MASS BINS ********************************************************/
		//locphiH and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiHVsMass;
		locUsedThisCombo_PhiHVsMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiHVsMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiHVsMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiHVsMass[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiHVsMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiHVsMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiHVsMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiHVsMass[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiHVsMass.find(locUsedThisCombo_PhiHVsMass) == locUsedSoFar_PhiHVsMass.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //apply omega mass cuts
		    double pH1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) pH1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) pH1 = -1;
		    else pH1 = 0;
		    dHist_PhiHVsMass->Fill(locphiH1,locOmegaPiMass, pH1);
		    double pH2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) pH2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) pH2 = -1;
		    else pH2 = 0;
		    dHist_PhiHVsMass->Fill(locphiH2,locOmegaPiMass, pH2);
		      locUsedSoFar_PhiHVsMass.insert(locUsedThisCombo_PhiHVsMass);
		  }

		/************************************************* HISTOGRAM COS(THETA) IN MASS BINS ****************************************************************/
		//loccostheta and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaVsMass;
		locUsedThisCombo_CosThetaVsMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaVsMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaVsMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaVsMass[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaVsMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaVsMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaVsMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaVsMass[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaVsMass.find(locUsedThisCombo_CosThetaVsMass) == locUsedSoFar_CosThetaVsMass.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //apply omega mass cuts
		    double w1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) w1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) w1 = -1;
		    else w1 = 0;
		    dHist_CosThetaVsMass->Fill(loccostheta1,locOmegaPiMass, w1);
		    double w2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) w2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) w2 = -1;
		    else w2 = 0;
		    dHist_CosThetaVsMass->Fill(loccostheta2,locOmegaPiMass, w2);
		      locUsedSoFar_CosThetaVsMass.insert(locUsedThisCombo_CosThetaVsMass);
		  }

		/************************************************* HISTOGRAM PHI IN MASS BINS ********************************************************/
		//locphi and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiVsMass;
		locUsedThisCombo_PhiVsMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiVsMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiVsMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiVsMass[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiVsMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiVsMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiVsMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiVsMass[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiVsMass.find(locUsedThisCombo_PhiVsMass) == locUsedSoFar_PhiVsMass.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    //apply omega mass cuts
		    double p1;
		    if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) p1 = +1;
		    else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) p1 = -1;
		    else p1 = 0;
		    dHist_PhiVsMass->Fill(locphi1,locOmegaPiMass, p1);
		    double p2;
		    if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) p2 = +1;
		    else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) p2 = -1;
		    else p2 = 0;
		    dHist_PhiVsMass->Fill(locphi2,locOmegaPiMass, p2);
		      locUsedSoFar_PhiVsMass.insert(locUsedThisCombo_PhiVsMass);
		  }

		/***************************************** HISTOGRAM COS(THETA) IN t [0.1, 0.3] BIN ***********************************************/
		//loccostheta, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosTheta_t1;
		locUsedThisCombo_CosTheta_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosTheta_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosTheta_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosTheta_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosTheta_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosTheta_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosTheta_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosTheta_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosTheta_t1.find(locUsedThisCombo_CosTheta_t1) == locUsedSoFar_CosTheta_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_CosTheta_t1->Fill(loccostheta1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_CosTheta_t1->Fill(loccostheta2, weight2);
			locUsedSoFar_CosTheta_t1.insert(locUsedThisCombo_CosTheta_t1);
		      }
		  }

		/***************************************** HISTOGRAM COS(THETA) IN t [0.3, 1.0] BIN ***********************************************/
		//loccostheta, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosTheta_t2;
		locUsedThisCombo_CosTheta_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosTheta_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosTheta_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosTheta_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosTheta_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosTheta_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosTheta_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosTheta_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosTheta_t2.find(locUsedThisCombo_CosTheta_t2) == locUsedSoFar_CosTheta_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_CosTheta_t2->Fill(loccostheta1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_CosTheta_t2->Fill(loccostheta2, weight2);
			locUsedSoFar_CosTheta_t2.insert(locUsedThisCombo_CosTheta_t2);
		      }
		  }

		/***************************************** HISTOGRAM PHI IN t [0.1, 0.3] BIN ***********************************************/
		//locphi, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_Phi_t1;
		locUsedThisCombo_Phi_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_Phi_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_Phi_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_Phi_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_Phi_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_Phi_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_Phi_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_Phi_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_Phi_t1.find(locUsedThisCombo_Phi_t1) == locUsedSoFar_Phi_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_Phi_t1->Fill(locphi1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_Phi_t1->Fill(locphi2, weight2);
			locUsedSoFar_Phi_t1.insert(locUsedThisCombo_Phi_t1);
		      }
		  }

		/***************************************** HISTOGRAM PHI IN t [0.3, 1.0] BIN ***********************************************/
		//locphi, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_Phi_t2;
		locUsedThisCombo_Phi_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_Phi_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_Phi_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_Phi_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_Phi_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_Phi_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_Phi_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_Phi_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_Phi_t2.find(locUsedThisCombo_Phi_t2) == locUsedSoFar_Phi_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_Phi_t2->Fill(locphi1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_Phi_t2->Fill(locphi2, weight2);
			locUsedSoFar_Phi_t2.insert(locUsedThisCombo_Phi_t2);
		      }
		  }

		/***************************************** HISTOGRAM COS(THETA_H) IN t [0.1, 0.3] BIN ***********************************************/
		//loccosthetaH, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaH_t1;
		locUsedThisCombo_CosThetaH_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaH_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaH_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaH_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaH_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaH_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaH_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaH_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaH_t1.find(locUsedThisCombo_CosThetaH_t1) == locUsedSoFar_CosThetaH_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_CosThetaH_t1->Fill(loccosthetaH1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_CosThetaH_t1->Fill(loccosthetaH2, weight2);
			locUsedSoFar_CosThetaH_t1.insert(locUsedThisCombo_CosThetaH_t1);
		      }
		  }

		/***************************************** HISTOGRAM COS(THETA_H) IN t [0.3, 1.0] BIN ***********************************************/
		//loccosthetaH, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaH_t2;
		locUsedThisCombo_CosThetaH_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaH_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaH_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaH_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaH_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaH_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaH_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaH_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaH_t2.find(locUsedThisCombo_CosThetaH_t2) == locUsedSoFar_CosThetaH_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_CosThetaH_t2->Fill(loccosthetaH1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_CosThetaH_t2->Fill(loccosthetaH2, weight2);
			locUsedSoFar_CosThetaH_t2.insert(locUsedThisCombo_CosThetaH_t2);
		      }
		  }

		/***************************************** HISTOGRAM PHI_H IN t [0.1, 0.3] BIN ***********************************************/
		//locphiH, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiH_t1;
		locUsedThisCombo_PhiH_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiH_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiH_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiH_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiH_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiH_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiH_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiH_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiH_t1.find(locUsedThisCombo_PhiH_t1) == locUsedSoFar_PhiH_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_PhiH_t1->Fill(locphiH1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_PhiH_t1->Fill(locphiH2, weight2);
			locUsedSoFar_PhiH_t1.insert(locUsedThisCombo_PhiH_t1);
		      }
		  }

		/***************************************** HISTOGRAM PHI_H IN t [0.3, 1.0] BIN ***********************************************/
		//locphiH, locMan_t, etc have been declared earlier
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiH_t2;
		locUsedThisCombo_PhiH_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiH_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiH_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiH_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiH_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiH_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiH_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiH_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiH_t2.find(locUsedThisCombo_PhiH_t2) == locUsedSoFar_PhiH_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//Apply omega mass cuts and histogram results
			double weight1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) weight1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) weight1 = -1;
			else weight1 = 0;
			dHist_PhiH_t2->Fill(locphiH1, weight1);
			double weight2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) weight2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) weight2 = -1;
			else weight2 = 0;
			dHist_PhiH_t2->Fill(locphiH2, weight2);
			locUsedSoFar_PhiH_t2.insert(locUsedThisCombo_PhiH_t2);
		      }
		  }


		/**************************************** HISTOGRAM ANGLES VS MASS IN t[0.1, 0.3] ******************************************/
		/************************************************* HISTOGRAM COS(THETA_H) IN MASS BINS ****************************************************************/
		//loccosthetaH and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaHVsMass_t1;
		locUsedThisCombo_CosThetaHVsMass_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaHVsMass_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaHVsMass_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaHVsMass_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaHVsMass_t1.find(locUsedThisCombo_CosThetaHVsMass_t1) == locUsedSoFar_CosThetaHVsMass_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double wH1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) wH1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) wH1 = -1;
			else wH1 = 0;
			dHist_CosThetaHVsMass_t1->Fill(loccosthetaH1,locOmegaPiMass, wH1);
			double wH2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) wH2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) wH2 = -1;
			else wH2 = 0;
			dHist_CosThetaHVsMass_t1->Fill(loccosthetaH2,locOmegaPiMass, wH2);
			locUsedSoFar_CosThetaHVsMass_t1.insert(locUsedThisCombo_CosThetaHVsMass_t1);
		      }
		  }

		/************************************************* HISTOGRAM PHI_H IN MASS BINS ********************************************************/
		//locphiH and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiHVsMass_t1;
		locUsedThisCombo_PhiHVsMass_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiHVsMass_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiHVsMass_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiHVsMass_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiHVsMass_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiHVsMass_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiHVsMass_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiHVsMass_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiHVsMass_t1.find(locUsedThisCombo_PhiHVsMass_t1) == locUsedSoFar_PhiHVsMass_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double pH1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) pH1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) pH1 = -1;
			else pH1 = 0;
			dHist_PhiHVsMass_t1->Fill(locphiH1,locOmegaPiMass, pH1);
			double pH2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) pH2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) pH2 = -1;
			else pH2 = 0;
			dHist_PhiHVsMass_t1->Fill(locphiH2,locOmegaPiMass, pH2);
			locUsedSoFar_PhiHVsMass_t1.insert(locUsedThisCombo_PhiHVsMass_t1);
		      }
		  }

		/************************************************* HISTOGRAM COS(THETA) IN MASS BINS ****************************************************************/
		//loccostheta and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaVsMass_t1;
		locUsedThisCombo_CosThetaVsMass_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaVsMass_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaVsMass_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaVsMass_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaVsMass_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaVsMass_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaVsMass_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaVsMass_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaVsMass_t1.find(locUsedThisCombo_CosThetaVsMass_t1) == locUsedSoFar_CosThetaVsMass_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double w1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) w1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) w1 = -1;
			else w1 = 0;
			dHist_CosThetaVsMass_t1->Fill(loccostheta1,locOmegaPiMass, w1);
			double w2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) w2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) w2 = -1;
			else w2 = 0;
			dHist_CosThetaVsMass_t1->Fill(loccostheta2,locOmegaPiMass, w2);
			locUsedSoFar_CosThetaVsMass_t1.insert(locUsedThisCombo_CosThetaVsMass_t1);
		      }
		  }

		/************************************************* HISTOGRAM PHI IN MASS BINS ********************************************************/
		//locphi and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiVsMass_t1;
		locUsedThisCombo_PhiVsMass_t1[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiVsMass_t1[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiVsMass_t1[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiVsMass_t1[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiVsMass_t1[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiVsMass_t1[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiVsMass_t1[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiVsMass_t1[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiVsMass_t1.find(locUsedThisCombo_PhiVsMass_t1) == locUsedSoFar_PhiVsMass_t1.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.1 && locMan_t < 0.3)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double p1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) p1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) p1 = -1;
			else p1 = 0;
			dHist_PhiVsMass_t1->Fill(locphi1,locOmegaPiMass, p1);
			double p2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) p2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) p2 = -1;
			else p2 = 0;
			dHist_PhiVsMass_t1->Fill(locphi2,locOmegaPiMass, p2);
			locUsedSoFar_PhiVsMass_t1.insert(locUsedThisCombo_PhiVsMass_t1);
		      }
		  }

		/**************************************** HISTOGRAM ANGLES VS MASS IN t[0.3, 1.0] ******************************************/
		/************************************************* HISTOGRAM COS(THETA_H) IN MASS BINS ****************************************************************/
		//loccosthetaH and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaHVsMass_t2;
		locUsedThisCombo_CosThetaHVsMass_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaHVsMass_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaHVsMass_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaHVsMass_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaHVsMass_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaHVsMass_t2.find(locUsedThisCombo_CosThetaHVsMass_t2) == locUsedSoFar_CosThetaHVsMass_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double wH1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) wH1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) wH1 = -1;
			else wH1 = 0;
			dHist_CosThetaHVsMass_t2->Fill(loccosthetaH1,locOmegaPiMass, wH1);
			double wH2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) wH2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) wH2 = -1;
			else wH2 = 0;
			dHist_CosThetaHVsMass_t2->Fill(loccosthetaH2,locOmegaPiMass, wH2);
			locUsedSoFar_CosThetaHVsMass_t2.insert(locUsedThisCombo_CosThetaHVsMass_t2);
		      }
		  }

		/************************************************* HISTOGRAM PHI_H IN MASS BINS ********************************************************/
		//locphiH and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiHVsMass_t2;
		locUsedThisCombo_PhiHVsMass_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiHVsMass_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiHVsMass_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiHVsMass_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiHVsMass_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiHVsMass_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiHVsMass_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiHVsMass_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiHVsMass_t2.find(locUsedThisCombo_PhiHVsMass_t2) == locUsedSoFar_PhiHVsMass_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double pH1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) pH1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) pH1 = -1;
			else pH1 = 0;
			dHist_PhiHVsMass_t2->Fill(locphiH1,locOmegaPiMass, pH1);
			double pH2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) pH2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) pH2 = -1;
			else pH2 = 0;
			dHist_PhiHVsMass_t2->Fill(locphiH2,locOmegaPiMass, pH2);
			locUsedSoFar_PhiHVsMass_t2.insert(locUsedThisCombo_PhiHVsMass_t2);
		      }
		  }

		/************************************************* HISTOGRAM COS(THETA) IN MASS BINS ****************************************************************/
		//loccostheta and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_CosThetaVsMass_t2;
		locUsedThisCombo_CosThetaVsMass_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_CosThetaVsMass_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_CosThetaVsMass_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_CosThetaVsMass_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_CosThetaVsMass_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CosThetaVsMass_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CosThetaVsMass_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_CosThetaVsMass_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_CosThetaVsMass_t2.find(locUsedThisCombo_CosThetaVsMass_t2) == locUsedSoFar_CosThetaVsMass_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double w1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) w1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) w1 = -1;
			else w1 = 0;
			dHist_CosThetaVsMass_t2->Fill(loccostheta1,locOmegaPiMass, w1);
			double w2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) w2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) w2 = -1;
			else w2 = 0;
			dHist_CosThetaVsMass_t2->Fill(loccostheta2,locOmegaPiMass, w2);
			locUsedSoFar_CosThetaVsMass_t2.insert(locUsedThisCombo_CosThetaVsMass_t2);
		      }
		  }

		/************************************************* HISTOGRAM PHI IN MASS BINS ********************************************************/
		//locphi and locOmegaPiMass have already been declared

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_PhiVsMass_t2;
		locUsedThisCombo_PhiVsMass_t2[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_PhiVsMass_t2[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_PhiVsMass_t2[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_PhiVsMass_t2[Gamma].insert(locPhoton4NeutralID);
		locUsedThisCombo_PhiVsMass_t2[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_PhiVsMass_t2[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_PhiVsMass_t2[Proton].insert(locProtonTrackID);
		locUsedThisCombo_PhiVsMass_t2[Unknown].insert(locBeamID);

		//compare to what's been used so far
		if(locUsedSoFar_PhiVsMass_t2.find(locUsedThisCombo_PhiVsMass_t2) == locUsedSoFar_PhiVsMass_t2.end())
		  {
		    //Make Mandelstam t cuts
		    if(locMan_t > 0.3 && locMan_t < 1.0)
		      {
			//unique combo: histogram it, and register this combo of particles
			//apply omega mass cuts
			double p1;
			if(loc3PiMass1_Measured > 0.733 && loc3PiMass1_Measured < 0.833) p1 = +1;
			else if((loc3PiMass1_Measured > 0.683 && loc3PiMass1_Measured < 0.733) || (loc3PiMass1_Measured > 0.833 && loc3PiMass1_Measured < 0.883)) p1 = -1;
			else p1 = 0;
			dHist_PhiVsMass_t2->Fill(locphi1,locOmegaPiMass, p1);
			double p2;
			if(loc3PiMass2_Measured > 0.733 && loc3PiMass2_Measured < 0.833) p2 = +1;
			else if((loc3PiMass2_Measured > 0.683 && loc3PiMass2_Measured < 0.733) || (loc3PiMass2_Measured > 0.833 && loc3PiMass2_Measured < 0.883)) p2 = -1;
			else p2 = 0;
			dHist_PhiVsMass_t2->Fill(locphi2,locOmegaPiMass, p2);
			locUsedSoFar_PhiVsMass_t2.insert(locUsedThisCombo_PhiVsMass_t2);
		      }
		  }


	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();


	
	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/
	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		FillOutputTree();
*/

	return kTRUE;
}

void DSelector_pomegapi::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
