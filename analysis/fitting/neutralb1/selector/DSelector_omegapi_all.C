#include "DSelector_omegapi_all.h"

void DSelector_omegapi_all::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "omegapi.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = "AmpToolsInputTree.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = "kin"; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	if(!(dTreeInterface->Get_Branch("NumCombos") == NULL)) Get_ComboWrappers();
	dPreviousRunNumber = 0;

	if(!(dTreeInterface->Get_Branch("NumCombos") == NULL)) {
		Get_ComboWrappers();

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

		//KINFIT RESULTS
        	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

        	//CUT MISSING MASS
	        dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));
	        dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.01));

		//INITIALIZE ACTIONS
	        //If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	        Initialize_Actions();
	}

	dHist_3piPiPlus1Mass = new TH1F("3piPiPlus1Mass", ";Invariant Mass (GeV)", 200, 0.5, 1.5);
        dHist_3piPiPlus2Mass = new TH1F("3piPiPlus2Mass", ";Invariant Mass (GeV)", 200, 0.5, 1.5);
        dHist_3piPiPlus1MassCorr = new TH2F("3piPiPlus1MassCorr", ";Invariant Mass (GeV); Invariant Mass (GeV)", 100, 0.5, 1.5, 100, 0.5, 1.5);
        dHist_4piMassSum = new TH1F("4piMassSum", ";Invariant Mass (GeV)", 200, 1.0, 2.0);
        dHist_3piMassSum = new TH1F("3piMassSum", ";Invariant Mass (GeV)", 200, 0.5, 1.5);
        dHist_2gammaMassSum = new TH1F("2gammaMassSum", ";Invariant Mass (GeV)", 500, 0.07, 0.20);
        dHist_KinFitChiSq = new TH1F("KinFitChiSq", ";#chi^{2}/NDF", 100, 0, 25);

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

	SetupAmpTools_FlatTree();

        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("t");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M3Pi");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M4Pi");
        dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Phi_Prod");
}

Bool_t DSelector_omegapi_all::Process(Long64_t locEntry)
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

	vector<TLorentzVector> locFinalStateThrownP4;
	vector<TLorentzVector> locPi0ThrownP4;
        if(Get_NumThrown()>0) {
		for (UInt_t i=0;i<Get_NumThrown();i++){
			dThrownWrapper->Set_ArrayIndex(i);
                        Int_t locThrownPID = dThrownWrapper->Get_PID();
                        if(locThrownPID==7) locPi0ThrownP4.push_back(dThrownWrapper->Get_P4()); // get Pi0s first to set proper order for AmpTools	
		}
		
		for (UInt_t i=0;i<Get_NumThrown();i++){
                        //Set branch array indices corresponding to this particle
                        dThrownWrapper->Set_ArrayIndex(i);
                        Int_t locThrownPID = dThrownWrapper->Get_PID();
                        if(locThrownPID==8 || locThrownPID==9 || locThrownPID==14)
				locFinalStateThrownP4.push_back(dThrownWrapper->Get_P4());
			if(locFinalStateThrownP4.size() == 1) {
				locFinalStateThrownP4.push_back(locPi0ThrownP4[1]);
				locFinalStateThrownP4.push_back(locPi0ThrownP4[0]);
			}
                        if(locFinalStateThrownP4.size() == 5) break;
                }
        }
	
	// fill generated from Thrown_Tree
	if(dOption.Contains("thrown") && dTreeInterface->Get_Branch("NumCombos") == NULL) {
                dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", 1.0);
                TLorentzVector locBeamP4 = dThrownBeam->Get_P4();
		//for(uint i=0; i<locFinalStateThrownP4.size(); i++)
		//	cout<<i<<" "<<locFinalStateThrownP4[i].M()<<endl;
                FillAmpTools_FlatTree(locBeamP4, locFinalStateThrownP4);
                Fill_FlatTree();
		
                return kTRUE;
        }

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

		if(dIsMC) {
			// Keep only generator beam photons for phasespace MC
			Bool_t locIsGeneratorFlag = (dThrownBeam->Get_P4().E() == dComboBeamWrapper->Get_P4().E() && fabs(dThrownBeam->Get_X4().T() - dComboBeamWrapper->Get_X4().T()) < 2.004) ? kTRUE : kFALSE;
	                if(dOption.Contains("noaccidental") && !(locIsGeneratorFlag || dComboBeamWrapper->Get_IsGenerator())) {
        	                dComboWrapper->Set_IsComboCut(true);
                	        continue;
			}
                }

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

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured + locPhoton4P4_Measured;
		double locMissingMassSquared = locMissingP4_Measured.M2();
		
		//3Pi Momenta
		TLorentzVector loc3PiP4_11 = locDecayingPi01P4 + locPiPlusP4 + locPiMinusP4;
                TLorentzVector loc3PiP4_12 = locDecayingPi02P4 + locPiPlusP4 + locPiMinusP4;

		/******************************************** ACCIDENTAL SUBTRACTION INFO *******************************************/

                // measured tagger time for combo
                TLorentzVector locBeam_X4_Measured = dComboBeamWrapper->Get_X4_Measured();

                // measured RF time for combo
                double locRFTime = dComboWrapper->Get_RFTime_Measured();

                // time difference between tagger and RF (corrected for production vertex position relative to target center)
                double locBeamDeltaT = locBeam_X4_Measured.T() - (locRFTime + (locBeam_X4_Measured.Z() - dTargetCenter.Z())/29.9792458);

                // calculate accidental subtraction weight based on time difference 
                double locAccWeight = 0.; // weight to accidentally subtracted histgorams
		bool locAccid = false; // flag to fill separate prompt and accidental histograms for later subtraction
		                if(fabs(locBeamDeltaT) < 0.5*4.008) { // prompt signal receives a weight of 1
                        locAccWeight = 1.;
                        locAccid = false;
                }
                else { // accidentals recieve a weight of 1/# RF bunches included in TTree (8 in this case)
                        locAccWeight = -1./8.;
                        locAccid = true;
		}

		// Loop through the analysis actions, executing them in order for the active particle combo
                if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
                        continue;

		// Apply cuts and if combo fails cut, then use dComboWrapper->Set_IsComboCut(true)
		double locKinFit_CL = dComboWrapper->Get_ConfidenceLevel_KinFit("");
                double locKinFit_RedChiSq = dComboWrapper->Get_ChiSq_KinFit("")/dComboWrapper->Get_NDF_KinFit("");
                if(locKinFit_CL < 1e-6) {
                        dComboWrapper->Set_IsComboCut(true);
                        continue;
                }

                if(locBeamP4.E()<8.2 || locBeamP4.E()>8.8) {
                        dComboWrapper->Set_IsComboCut(true);
                        continue;
                }

		//FILL FLAT TREE
		vector<Int_t> locFinalStatePID {2212, 111, 111, 211, -211};
                vector<TLorentzVector> locPi0P4;
                locPi0P4.push_back(locDecayingPi01P4); locPi0P4.push_back(locDecayingPi02P4);
                for(int i=0; i<(int)locPi0P4.size(); i++) {
			// omega mass cut and sideband weight
			TLorentzVector loc3PiP4 = locPi0P4[i] + locPiPlusP4 + locPiMinusP4;
                        double loc3PiMass = loc3PiP4.M();
			TLorentzVector loc3PiP4_alt = locPi0P4[abs(i-1)] + locPiPlusP4 + locPiMinusP4;
                        double loc3PiMass_alt = loc3PiP4_alt.M();

			dHist_KinFitChiSq->Fill(locKinFit_RedChiSq, locAccWeight);
                        if(locKinFit_RedChiSq > 5)
                                continue;

                        dHist_3piMassSum->Fill(loc3PiMass);
		
			// t cut
			double loct = -1. * (locProtonP4 - dTargetP4).M2();
                        if(loct > 1.0)
                                continue;

                        double loc2Dweight;
                        double locLmin = 0.690;
                        double locLmax = 0.735;
                        double locomegamin = 0.760;
                        double locomegamax = 0.805;
                        double locHmin = 0.830;
                        double locHmax = 0.875;

			if((loc3PiMass > locomegamin && loc3PiMass < locomegamax) && (loc3PiMass_alt > locomegamin && loc3PiMass_alt < locomegamax)) {
				// "double-omega" events only used once, closest to true omega mass
				if(fabs(loc3PiMass-0.783) < fabs(loc3PiMass_alt-0.783)) loc2Dweight = 1.0;
				else continue;
			}
			else if((loc3PiMass > locomegamin && loc3PiMass < locomegamax)) loc2Dweight = 1.0;
			else if((loc3PiMass_alt > locomegamin && loc3PiMass_alt < locomegamax)) continue;
			// contol region overlap (cross hatched) are counted twice with a weight of -5/8
			else if((loc3PiMass > locLmin && loc3PiMass < locLmax) && (loc3PiMass_alt > locLmin && loc3PiMass_alt < locLmax)) loc2Dweight = -0.625;
			else if((loc3PiMass > locLmin && loc3PiMass < locLmax) && (loc3PiMass_alt > locHmin && loc3PiMass_alt < locHmax)) loc2Dweight = -0.625;
			else if((loc3PiMass > locHmin && loc3PiMass < locHmax) && (loc3PiMass_alt > locLmin && loc3PiMass_alt < locLmax)) loc2Dweight = -0.625;
			else if((loc3PiMass > locHmin && loc3PiMass < locHmax) && (loc3PiMass_alt > locHmin && loc3PiMass_alt < locHmax)) loc2Dweight = -0.625;
			// control regions included once with weight of -0.5
			else if((loc3PiMass > locLmin && loc3PiMass < locLmax) || (loc3PiMass > locHmin && loc3PiMass < locHmax)) loc2Dweight = -0.5;
			else continue;
				
			double loc4PiMass = (loc3PiP4 + locPi0P4[abs(i-1)]).M();
                        dHist_4piMassSum->Fill(loc4PiMass, loc2Dweight);
			
			// set ordering of pions for amplitude analysis
			vector<TLorentzVector> locFinalStateP4;
                        locFinalStateP4.push_back(locProtonP4);
                        locFinalStateP4.push_back(locPi0P4[abs(i-1)]); // Bachelor Pi0
                        locFinalStateP4.push_back(locPi0P4[i]);
                        locFinalStateP4.push_back(locPiPlusP4);
                        locFinalStateP4.push_back(locPiMinusP4);

			dFlatTreeInterface->Fill_Fundamental<Float_t>("M3Pi", loc3PiMass);
                        dFlatTreeInterface->Fill_Fundamental<Float_t>("M4Pi", loc4PiMass);
                        dFlatTreeInterface->Fill_Fundamental<Float_t>("Phi_Prod", locProtonP4.Phi());
                        dFlatTreeInterface->Fill_Fundamental<Float_t>("t", loct);

			// set weight according to omega mass
			dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", loc2Dweight*locAccWeight);
			// set ordered final state P4 for filling flat tree
			if(!dOption.Contains("accept")) 
                		FillAmpTools_FlatTree(locBeamP4, locFinalStateP4);
                        else {
				dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", 1.0);
                        	FillAmpTools_FlatTree(locBeamP4, locFinalStateThrownP4);
			}

			Fill_FlatTree(); //for the active combo
		}

	} // end of combo loop

	return kTRUE;
}

void DSelector_omegapi_all::Finalize(void)
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
