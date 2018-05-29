#include "DSelector_pomega2pi_omega3pi.h"

void DSelector_pomega2pi_omega3pi::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "pomega2pi_omega3pi.root"; //"" for none
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
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, Proton, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, Proton, SYS_BCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiPlus, SYS_TOF));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, PiPlus, SYS_BCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiMinus, SYS_TOF));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, PiMinus, SYS_BCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, Gamma, SYS_BCAL));

        //MASSES
	deque<Particle_t> locPi0PIDs;  locPi0PIDs.push_back(Gamma); locPi0PIDs.push_back(Gamma);
	deque<Particle_t> locRhoPIDs;  locRhoPIDs.push_back(PiPlus); locRhoPIDs.push_back(PiMinus);
	deque<Particle_t> locOmegaPIDs;  locOmegaPIDs.push_back(Pi0); locOmegaPIDs.push_back(PiPlus); locOmegaPIDs.push_back(PiMinus);
        deque<Particle_t> locb1PlusPIDs;  locb1PlusPIDs.push_back(Pi0); locb1PlusPIDs.push_back(PiPlus); locb1PlusPIDs.push_back(PiMinus); locb1PlusPIDs.push_back(PiPlus);
	deque<Particle_t> locb1MinusPIDs;  locb1MinusPIDs.push_back(Pi0); locb1MinusPIDs.push_back(PiPlus); locb1MinusPIDs.push_back(PiMinus); locb1MinusPIDs.push_back(PiMinus);
	deque<Particle_t> locDeltaPlusPlusPIDs;  locDeltaPlusPlusPIDs.push_back(PiPlus); locDeltaPlusPlusPIDs.push_back(Proton);
	deque<Particle_t> locDelta0PIDs;  locDelta0PIDs.push_back(PiMinus); locDelta0PIDs.push_back(Proton);
	deque<Particle_t> locXPIDs;  locXPIDs.push_back(omega); locXPIDs.push_back(PiPlus); locXPIDs.push_back(PiMinus);
	double minPi0 = 0.09; double maxPi0 = 0.18;
	double minOmega = 0.5; double maxOmega = 1.1;
	double minb1 = 0.5; double maxb1 = 2.0;
	double minDelta = 1.0; double maxDelta = 3.0;
	double minX = 1.0; double maxX = 3.0;
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, Pi0, 100, minPi0, maxPi0, "Pi0_NoCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 100, minOmega, maxOmega, "omega_NoCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PlusPIDs, 150, minb1, maxb1, "b1Plus_NoCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1MinusPIDs, 150, minb1, maxb1, "b1Minus_NoCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, 150, minDelta, maxDelta, "DeltaPlusPlus_NoCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDelta0PIDs, 150, minDelta, maxDelta, "Delta0_NoCut"));
	//dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 1, locPi0PIDs, 0, locOmegaPIDs, 100, minPi0, maxPi0, 100, minOmega, maxOmega, "Omega_Pi0_NoCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, locb1PlusPIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1Plus_NoCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, locb1MinusPIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1Minus_NoCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, locb1MinusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "DeltaPlusPlus_b1Minus_NoCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDelta0PIDs, locb1PlusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "Delta0_b1Plus_NoCut"));
        dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

        //KINFIT RESULTS
        dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

        //CUT MISSING MASS
        dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));
        dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.01));

        //MASSES (after KINFIT CL cut)
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, Pi0, 100, minPi0, maxPi0, "Pi0_KinCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 100, minOmega, maxOmega, "omega_KinCut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PlusPIDs, 150, minb1, maxb1, "b1Plus_KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1MinusPIDs, 150, minb1, maxb1, "b1Minus_KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, 150, minDelta, maxDelta, "DeltaPlusPlus_KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDelta0PIDs, 150, minDelta, maxDelta, "Delta0_KinCut"));
	//dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 1, locPi0PIDs, 0, locOmegaPIDs, 100, minPi0, maxPi0, 100, minOmega, maxOmega, "Omega_Pi0_KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, locb1PlusPIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1Plus_KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, locb1MinusPIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1Minus_KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, locb1MinusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "DeltaPlusPlus_b1Minus_KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDelta0PIDs, locb1PlusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "Delta0_b1Plus_KinCut"));

	//CUT PI0 MASS
        dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, Pi0, 0.114, 0.156));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 100, minOmega, maxOmega, "omega_Pi0Cut"));
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PlusPIDs, 150, minb1, maxb1, "b1Plus_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1MinusPIDs, 150, minb1, maxb1, "b1Minus_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, 150, minDelta, maxDelta, "DeltaPlusPlus_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDelta0PIDs, 150, minDelta, maxDelta, "Delta0_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locRhoPIDs, 150, 0.0, 1.5, "Rho_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, locb1PlusPIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1Plus_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, locb1MinusPIDs, 100, minOmega, maxOmega, 150, minb1, maxb1, "Omega_b1Minus_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locOmegaPIDs, locb1MinusPIDs, 100, minOmega, maxOmega, 150, 0., 1.5, "Omega_Rho_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, locb1MinusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "DeltaPlusPlus_b1Minus_Pi0Cut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDelta0PIDs, locb1PlusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "Delta0_b1Plus_Pi0Cut"));

        //CUT OMEGA MASS
        dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 0.71, 0.85));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PlusPIDs, 150, minb1, maxb1, "b1Plus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1MinusPIDs, 150, minb1, maxb1, "b1Minus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, 150, minDelta, maxDelta, "DeltaPlusPlus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDelta0PIDs, 150, minDelta, maxDelta, "Delta0_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locRhoPIDs, 150, 0.0, 1.5, "Rho_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, locb1MinusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "DeltaPlusPlus_b1Minus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDelta0PIDs, locb1PlusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "Delta0_b1Plus_OmegaCut"));

	//CUT AWAY DELTA CONTRIBUTION
	dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, 1.4, 1000.));
	dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, 0, locDelta0PIDs, 1.4, 1000.));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locRhoPIDs, 150, 0.0, 1.5, "Rho_DeltaCut"));
	dAnalysisActions.push_back(new DHistogramAction_Dalitz(dComboWrapper, true, 0, locb1PlusPIDs, locb1MinusPIDs, 150, minb1*minb1, maxb1*maxb1, 150, minb1*minb1, maxb1*maxb1, "b1Plus_b1Minus_Dalitz_DeltaCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locb1PlusPIDs, locXPIDs, 150, minb1, maxb1, 200, minX, maxX, "X_b1Plus_DeltaCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locb1MinusPIDs, locXPIDs, 150, minb1, maxb1, 200, minX, maxX, "X_b1Minus_DeltaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locXPIDs, 200, minX, maxX, "X_DeltaCut"));

	//CUT BEAM ENERGY
        dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, true));
	dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, true, 7.0, 12.0));
	dAnalysisActions.push_back(new DHistogramAction_Dalitz(dComboWrapper, true, 0, locb1PlusPIDs, locb1MinusPIDs, 150, minb1*minb1, maxb1*maxb1, 150, minb1*minb1, maxb1*maxb1, "b1Plus_b1Minus_Dalitz_BeamCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locb1PlusPIDs, locXPIDs, 150, minb1, maxb1, 200, minX, maxX, "X_b1Plus_BeamCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locb1MinusPIDs, locXPIDs, 150, minb1, maxb1, 200, minX, maxX, "X_b1Minus_BeamCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locXPIDs, 200, minX, maxX, "X_BeamCut"));

        //KINEMATICS
        dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, true));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

	/***************************************** ADVANCED: CHOOSE BRANCHES TO READ ****************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_pomega2pi_omega3pi::Process(Long64_t locEntry)
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
		Int_t locPiPlus1TrackID = dPiPlus1Wrapper->Get_TrackID();
		Int_t locPiMinus1TrackID = dPiMinus1Wrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();
		Int_t locPiPlus2TrackID = dPiPlus2Wrapper->Get_TrackID();
		Int_t locPiMinus2TrackID = dPiMinus2Wrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiPlus1P4 = dPiPlus1Wrapper->Get_P4();
		TLorentzVector locPiMinus1P4 = dPiMinus1Wrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		TLorentzVector locPiPlus2P4 = dPiPlus2Wrapper->Get_P4();
		TLorentzVector locPiMinus2P4 = dPiMinus2Wrapper->Get_P4();
		//Step 1
		TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPiPlus1P4_Measured = dPiPlus1Wrapper->Get_P4_Measured();
		TLorentzVector locPiMinus1P4_Measured = dPiMinus1Wrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		TLorentzVector locPiPlus2P4_Measured = dPiPlus2Wrapper->Get_P4_Measured();
		TLorentzVector locPiMinus2P4_Measured = dPiMinus2Wrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPiPlus1P4_Measured + locPiMinus1P4_Measured + locProtonP4_Measured + locPiPlus2P4_Measured + locPiMinus2P4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

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

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton2NeutralID);

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

void DSelector_pomega2pi_omega3pi::Finalize(void)
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
