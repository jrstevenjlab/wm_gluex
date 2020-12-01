#include "DSelector_omegapi.h"

void DSelector_omegapi::Init(TTree *locTree)
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

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

	SetupAmpTools_FlatTree();

	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("t");
	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M3Pi");
	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M4Pi");
	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("MRecoil");
	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Phi_Prod");
}

Bool_t DSelector_omegapi::Process(Long64_t locEntry)
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

	// get thrown particles in case we want to test fit with truth information
	vector<TLorentzVector> locFinalStateThrownP4;
	TLorentzVector locPi0ThrownP4;
	if(Get_NumThrown()>0) {
		for (UInt_t i=0;i<Get_NumThrown();i++){
			dThrownWrapper->Set_ArrayIndex(i);
                        Int_t locThrownPID = dThrownWrapper->Get_PID();
                        if(locThrownPID==7) locPi0ThrownP4 = dThrownWrapper->Get_P4(); // get Pi0 first to set proper order for AmpTools
		}
         	for (UInt_t i=0;i<Get_NumThrown();i++){
                	//Set branch array indices corresponding to this particle
                	dThrownWrapper->Set_ArrayIndex(i);
			Int_t locThrownPID = dThrownWrapper->Get_PID();
			if(locThrownPID==14 || locThrownPID==8 || locThrownPID==9) {
				locFinalStateThrownP4.push_back(dThrownWrapper->Get_P4());
			}
			if(locFinalStateThrownP4.size() == 2) locFinalStateThrownP4.push_back(locPi0ThrownP4);
			if(locFinalStateThrownP4.size() == 6) break;
	    	}
	}

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		dComboWrapper->Set_IsComboCut(false);
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locPiPlus1TrackID = dPiPlus1Wrapper->Get_TrackID();
		Int_t locPiPlus2TrackID = dPiPlus2Wrapper->Get_TrackID();
		Int_t locPiMinus1TrackID = dPiMinus1Wrapper->Get_TrackID();
		Int_t locPiMinus2TrackID = dPiMinus2Wrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiPlus1P4 = dPiPlus1Wrapper->Get_P4();
		TLorentzVector locPiPlus2P4 = dPiPlus2Wrapper->Get_P4();
		TLorentzVector locPiMinus1P4 = dPiMinus1Wrapper->Get_P4();
		TLorentzVector locPiMinus2P4 = dPiMinus2Wrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPiPlus1P4_Measured = dPiPlus1Wrapper->Get_P4_Measured();
		TLorentzVector locPiPlus2P4_Measured = dPiPlus2Wrapper->Get_P4_Measured();
		TLorentzVector locPiMinus1P4_Measured = dPiMinus1Wrapper->Get_P4_Measured();
		TLorentzVector locPiMinus2P4_Measured = dPiMinus2Wrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// Combine 4-vectors
                TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
                locMissingP4_Measured -= locPiPlus1P4_Measured + locPiPlus2P4_Measured + locPiMinus1P4_Measured + locPiMinus2P4_Measured + locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;
		double locMissingMassSquared = locMissingP4_Measured.M2();

		// pi0 sum
		TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;

		// ProtonPiPlus Momenta:
		TLorentzVector locProtonPiPlus1P4 = locProtonP4 + locPiPlus1P4;
		TLorentzVector locProtonPiPlus2P4 = locProtonP4 + locPiPlus2P4;

		 //3Pi Momenta:
                TLorentzVector loc3PiP4_11 = locPiPlus1P4 + locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
                TLorentzVector loc3PiP4_12 = locPiPlus1P4 + locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;
                TLorentzVector loc3PiP4_21 = locPiPlus2P4 + locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
                TLorentzVector loc3PiP4_22 = locPiPlus2P4 + locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;

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

		/****************************************** Set up omega sideband weights ************************************/

		double loc3PiMass11 = loc3PiP4_11.M();
                double loc3PiMass12 = loc3PiP4_12.M();
                double loc3PiMass21 = loc3PiP4_21.M();
                double loc3PiMass22 = loc3PiP4_22.M();

                //Sidebands range from 3-4 sigma on either side
                double weight11;
                if (loc3PiMass11 > 0.7479 && loc3PiMass11 < 0.8169) weight11 = +1;
                else if ((loc3PiMass11 > 0.7134 && loc3PiMass11 < 0.7249) || (loc3PiMass11 > 0.8399 && loc3PiMass11 < 0.8514)) weight11 = -3.;
                else weight11 = 0;
                double weight12;
                if (loc3PiMass12 > 0.7479 && loc3PiMass12 < 0.8169) weight12 = +1;
                else if ((loc3PiMass12 > 0.7134 && loc3PiMass12 < 0.7249) || (loc3PiMass12 > 0.8399 && loc3PiMass12 < 0.8514)) weight12 = -3.;
                else weight12 = 0;
                double weight21;
                if (loc3PiMass21 > 0.7479 && loc3PiMass21 < 0.8169) weight21 = +1;
                else if ((loc3PiMass21 > 0.7134 && loc3PiMass21 < 0.7249) || (loc3PiMass21 > 0.8399 && loc3PiMass21 < 0.8514)) weight21 = -3.;
                else weight21 = 0;
		double weight22;
                if (loc3PiMass22 > 0.7479 && loc3PiMass22 < 0.8169) weight22 = +1;
                else if ((loc3PiMass22 > 0.7134 && loc3PiMass22 < 0.7249) || (loc3PiMass22 > 0.8399 && loc3PiMass22 < 0.8514)) weight22 = -3.;
                else weight22 = 0;
		
		// Apply cuts and if combo fails cut, then use dComboWrapper->Set_IsComboCut(true)
		double locKinFit_CL = dComboWrapper->Get_ConfidenceLevel_KinFit("");
                if(locKinFit_CL < 1e-4) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
			
		if(locBeamP4.E()<8.2 || locBeamP4.E()>8.8) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		//FILL FLAT TREE
		vector<Int_t> locFinalStatePID {2212, -211, 111, 211, -211, 211};
		vector<TLorentzVector> locPiMinusP4; 
		locPiMinusP4.push_back(locPiMinus1P4); locPiMinusP4.push_back(locPiMinus2P4);
		vector<TLorentzVector> locPiPlusP4; 
		locPiPlusP4.push_back(locPiPlus1P4); locPiPlusP4.push_back(locPiPlus2P4);
		for(int iminus=0; iminus<(int)locPiMinusP4.size(); iminus++) {
			for(int iplus=0; iplus<(int)locPiPlusP4.size(); iplus++) {

				// omega mass cut and sideband weight
				TLorentzVector loc3PiP4 = locPi0P4 + locPiPlusP4[iplus] + locPiMinusP4[iminus];
				double loc3PiMass = loc3PiP4.M();
				
				// Delta++ cut
				TLorentzVector locRecoilP4 = locProtonP4+locPiPlusP4[abs(iplus-1)];
				if(locRecoilP4.M() > 1.6) 
					continue; 

				// loose pi0 mass
				if(fabs(locPi0P4.M() - 0.135) > 0.015) continue; 

				// loose -t cut for now
				double loct = -1. * (locRecoilP4 - dTargetP4).M2();
				if(loct > 1.0) continue;

				// set weight from 2D 3pi mass correlation (see Chung et. al. 1975)...
				double weight = 1.0;
				// PLACEHOLDER NEEDS TO BE FIXED!
                                if (loc3PiMass > 0.7479 && loc3PiMass < 0.8169) weight *= 1.0;
                                else if ((loc3PiMass > 0.7134 && loc3PiMass < 0.7249) || (loc3PiMass > 0.8399 && loc3PiMass < 0.8514)) weight = -3.;
                                else continue;

				double loc4PiMass = (loc3PiP4 + locPiMinusP4[abs(iminus-1)]).M();

				// set ordering of pions for amplitude analysis
				vector<TLorentzVector> locFinalStateP4;
				locFinalStateP4.push_back(locProtonP4);
				locFinalStateP4.push_back(locPiMinusP4[abs(iminus-1)]); // Bachelor Pi-
				locFinalStateP4.push_back(locPi0P4); 
				locFinalStateP4.push_back(locPiPlusP4[iplus]); 
				locFinalStateP4.push_back(locPiMinusP4[iminus]);
				locFinalStateP4.push_back(locPiPlusP4[abs(iplus-1)]); // Delta++ recoil Pi+

				dFlatTreeInterface->Fill_Fundamental<Float_t>("M3Pi", loc3PiMass);
				dFlatTreeInterface->Fill_Fundamental<Float_t>("M4Pi", loc4PiMass);
				dFlatTreeInterface->Fill_Fundamental<Float_t>("MRecoil", locRecoilP4.M());
				dFlatTreeInterface->Fill_Fundamental<Float_t>("Phi_Prod", locRecoilP4.Phi());
				dFlatTreeInterface->Fill_Fundamental<Float_t>("t", loct);

				// set weight according to omega mass
				dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", weight*locAccWeight);
		
				// set ordered final state P4 for filling flat tree
				FillAmpTools_FlatTree(locBeamP4, locFinalStateP4); 
				//FillAmpTools_FlatTree(locBeamP4, locFinalStateThrownP4);

				Fill_FlatTree(); //for the active combo
			}
		}

		if(dOption.Contains("signal") && weight11 <= 0. && weight12 <= 0.) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		if(dOption.Contains("bkgnd") && weight11 >= 0. && weight12 >= 0.) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

	} // end of combo loop

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
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
	if(!locIsEventCut){ // && dOutputTreeFileName != "") {
		cout<<"filled tree good entry "<<locEntry<<endl;
		cout<<dOutputTreeFileName<<endl;
		//eventCounter++;
		Fill_OutputTree();
	}

	return kTRUE;
}

void DSelector_omegapi::Finalize(void)
{
	//cout<<"Total events written = "<<eventCounter<<endl;

	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
