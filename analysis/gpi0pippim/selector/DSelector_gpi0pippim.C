#include "DSelector_gpi0pippim.h"

void DSelector_gpi0pippim::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "gpi0pippim.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.3, Proton, SYS_TOF));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.6, Gamma, SYS_FCAL));
	dAnalysisActions.push_back(new DCutAction_dEdx(dComboWrapper, false, Proton, SYS_CDC));

	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));

	//CUT ON SHOWER QUALITY
	dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.0, 12.0));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));
	
	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_Pi0Mass = new TH1F("Pi0Mass", ";#gamma#gamma Mass (GeV)", 600, 0.07, 0.20);
	dHist_Pi02Mass = new TH1F("Pi02Mass", ";#gamma#gamma Mass (GeV)", 600, 0.0, 1.0);
	dHist_Pi03Mass = new TH1F("Pi03Mass", ";#gamma#gamma Mass (GeV)", 600, 0.0, 1.0);

	dHist_3PiMass = new TH1F("3PiMass", ";#pi^{+}#pi^{-}#gamma#gamma Mass (GeV)", 600, 0.1, 1.3);
	dHist_3Pi2Mass = new TH1F("3Pi2Mass", ";#pi^{+}#pi^{-}#gamma#gamma Mass (GeV)", 600, 0.1, 1.3);
	dHist_3Pi3Mass = new TH1F("3Pi3Mass", ";#pi^{+}#pi^{-}#gamma#gamma Mass (GeV)", 600, 0.1, 1.3);

	dHist_ProtonPiPlusMass = new TH1F("ProtonPiPlusMass", ";p#pi^{+} Mass (GeV)", 600, 1.0, 2.5);

	dHist_Pi0GammaVsPi0Mass = new TH2F("Pi0GammaVsPi0", "; #gamma#gamma Mass (GeV); #gamma#gamma#gamma Mass (GeV)", 100, 0.07, 0.20, 100, 0.1, 1.3);
	dHist_Pi0GammaMass = new TH1F("Pi0GammaMass", ";#pi^{0}#gamma Mass (GeV)", 600, 0.1, 1.3);

	dHist_Pi0PiMinusMass = new TH1F("Pi0PiMinusMass", ";#pi^{0}#pi^{-} Mass (GeV)", 600, 0.1, 1.3);
	dHist_Pi0PiMinusVsPi0Gamma = new TH2F("Pi0PiMinusVsPi0Gamma", "; #pi^{0}#gamma Mass (GeV); #pi^{0}#pi^{-} Mass (GeV)", 100, 0.1, 1.3, 100, 0.1, 1.3);
	dHist_KinFitChiSqVsPi0GammaMass = new TH2F("KinFitChiSqVsPi0GammaMass", "; #pi^{0}#gamma Mass (GeV); KinFit #chi^{2}", 100, 0.1, 1.3, 200, 0, 20);

	dHist_Pi0PiMinusGammaVsPi0Gamma = new TH2F("Pi0PiMinusGammaVsPi0Gamma", "; #pi^{0}#gamma Mass (GeV); #pi^{0}#pi^{-}#gamma Mass (GeV)", 100, 0.1, 1.3, 100, 0.1, 2.1);
	dHist_Pi0PiMinusGammaVsPi0PiMinus = new TH2F("Pi0PiMinusGammaVsPi0PiMinus", "; #pi^{0}#pi^{-} Mass (GeV); #pi^{0}#pi^{-}#gamma Mass (GeV)", 100, 0.1, 1.3, 100, 0.1, 2.1);
	
	dHist_OmegaPiMinusMass = new TH1F("OmegaPiMinusMass", ";#omega#pi^{-} Mass (GeV)", 600, 0.9, 2.1);

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

}

Bool_t DSelector_gpi0pippim::Process(Long64_t locEntry)
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
	//TLorentzVector locProductionX4 = Get_X4_Production();

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
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();
		Int_t locPhoton3NeutralID = dPhoton3Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
		TLorentzVector locPhoton3P4 = dPhoton3Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton3P4_Measured = dPhoton3Wrapper->Get_P4_Measured();

		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t locNumOutOfTimeBunchesInTree = 4; //YOU need to specify this number
			//Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 

		Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
		if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		 	dComboWrapper->Set_IsComboCut(true); 
		 	continue; 
		} 

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPhoton1P4_Measured + locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured;

		TLorentzVector locPi0P4 = locPhoton2P4 + locPhoton3P4;
		TLorentzVector locPi0GammaP4 = locPi0P4 + locPhoton1P4;
		TLorentzVector loc3PiP4 = locPi0P4 + locPiPlusP4 + locPiMinusP4;
		TLorentzVector loc3Pi2P4 = locPhoton1P4 + locPhoton2P4 + locPiPlusP4 + locPiMinusP4;
		TLorentzVector loc3Pi3P4 = locPhoton1P4 + locPhoton3P4 + locPiPlusP4 + locPiMinusP4;

		TLorentzVector locOmegaPiMinusP4 = locPi0GammaP4 + locPiMinusP4;
		TLorentzVector locDeltaPlusPlus = locProtonP4 + locPiPlusP4;

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		//if the active combo fails a cut, IsComboCutFlag automatically set
		if(!Execute_Actions()) 
			continue;

		// t cut
		double t = ( (locProtonP4+locPiPlusP4) - dTargetP4 ).M2();
		if( fabs(t) > 0.5 ) continue;

		// kinematic fit cut
		double kinFitChiSq = dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit();
		if(kinFitChiSq > 5.0) continue;

		// possible cut to remove events with additional unused showers
		double energyUnusedShowers = dComboWrapper->Get_Energy_UnusedShowers();
		double numUnusedShowers = dComboWrapper->Get_NumUnusedShowers();
		//if(numUnusedShowers > 0) continue;

		// possible cuts to select photons with good quality (reduce split-offs, hadronic deposit)
		//if( dPhoton1Wrapper->Get_Shower_Quality() < 0.75 ) continue;
		//if( dPhoton1Wrapper->Get_Energy_BCAL() <= 0.0 && dPhoton1Wrapper->Get_Energy_FCAL() < 0.5 ) continue;
		//if( dPhoton2Wrapper->Get_Energy_BCAL() <= 0.0 && dPhoton2Wrapper->Get_Energy_FCAL() < 0.25 ) continue;
		//if( dPhoton3Wrapper->Get_Energy_BCAL() <= 0.0 && dPhoton3Wrapper->Get_Energy_FCAL() < 0.25 ) continue;

		dHist_Pi0Mass->Fill(locPi0P4.M(), locHistAccidWeightFactor);
		dHist_Pi02Mass->Fill((locPhoton1P4+locPhoton2P4).M(), locHistAccidWeightFactor);
		dHist_Pi03Mass->Fill((locPhoton1P4+locPhoton3P4).M(), locHistAccidWeightFactor);

		// veto wrong combination pi0s
		if( fabs((locPhoton1P4+locPhoton2P4).M() - 0.135) < 0.035 ) continue;
		if( fabs((locPhoton1P4+locPhoton3P4).M() - 0.135) < 0.035 ) continue;

		dHist_3PiMass->Fill(loc3PiP4.M(), locHistAccidWeightFactor);
		dHist_3Pi2Mass->Fill(loc3Pi2P4.M(), locHistAccidWeightFactor);
		dHist_3Pi3Mass->Fill(loc3Pi3P4.M(), locHistAccidWeightFactor);
		
		// veto omega->3pi decay
		if( loc3PiP4.M() < 0.9 ) continue; 

		// plot and select Delta++ 
		dHist_ProtonPiPlusMass->Fill(locDeltaPlusPlus.M(), locHistAccidWeightFactor);
		if( locDeltaPlusPlus.M() > 1.3 ) continue; 

		dHist_Pi0Mass->Fill(locPi0P4.M(), locHistAccidWeightFactor);
		dHist_Pi0GammaVsPi0Mass->Fill(locPi0P4.M(), locPi0GammaP4.M(), locHistAccidWeightFactor);

		// select pi0
		if( fabs(locPi0P4.M() - 0.135) > 0.015 ) continue;

		// plot and select events by Kinematic Fit cut
		dHist_KinFitChiSqVsPi0GammaMass->Fill(locPi0GammaP4.M(), kinFitChiSq, locHistAccidWeightFactor);
		if(kinFitChiSq > 2.0) continue;

		// Lots of mass distributions...
		dHist_Pi0GammaMass->Fill(locPi0GammaP4.M(), locHistAccidWeightFactor); // omega radiative
		dHist_Pi0PiMinusMass->Fill((locPi0P4+locPiMinusP4).M(), locHistAccidWeightFactor); // rho -> pi0 pi-
		dHist_Pi0PiMinusVsPi0Gamma->Fill(locPi0GammaP4.M(), (locPi0P4+locPiMinusP4).M(), locHistAccidWeightFactor);
		
		dHist_Pi0PiMinusGammaVsPi0Gamma->Fill(locPi0GammaP4.M(), (locPi0GammaP4+locPiMinusP4).M(), locHistAccidWeightFactor);
		dHist_Pi0PiMinusGammaVsPi0PiMinus->Fill((locPi0P4+locPiMinusP4).M(), (locPi0GammaP4+locPiMinusP4).M(), locHistAccidWeightFactor);

		// fill sideband weighted omega-pi spectrum
		if( fabs(locPi0GammaP4.M() - 0.782) < 0.04 )
			dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusP4.M(), locHistAccidWeightFactor);
		else if( fabs(locPi0GammaP4.M() - 0.670) < 0.02 )
			dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusP4.M(), -1.0*locHistAccidWeightFactor);
		else if( fabs(locPi0GammaP4.M() - 0.894) < 0.02 )
			dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusP4.M(), -1.0*locHistAccidWeightFactor);

	
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	return kTRUE;
}

void DSelector_gpi0pippim::Finalize(void)
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
