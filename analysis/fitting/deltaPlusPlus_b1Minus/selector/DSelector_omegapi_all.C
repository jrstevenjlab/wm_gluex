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

	dHist_ProtonPiPlus1Mass = new TH1F("ProtonPiPlus1Mass", ";Invariant Mass (GeV)", 1000, 1.0, 2.5);
	dHist_ProtonPiPlus2Mass = new TH1F("ProtonPiPlus2Mass", ";Invariant Mass (GeV)", 1000, 1.0, 2.5);

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

	SetupAmpTools_FlatTree();

	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("t");
	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M3Pi");
	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M4Pi");
	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("MRecoil");
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
	TLorentzVector locPi0ThrownP4;
	TLorentzVector locOmegaP4;
	TLorentzVector locRecoilP4;
	if(Get_NumThrown()>0) {
		for (UInt_t i=0;i<Get_NumThrown();i++){
			dThrownWrapper->Set_ArrayIndex(i);
                        Int_t locThrownPID = dThrownWrapper->Get_PID();
                        if(locThrownPID==7) { // get Pi0 for later use
				locPi0ThrownP4 = dThrownWrapper->Get_P4(); 
				locOmegaP4 = locPi0ThrownP4; // add omega decay neutral pion
			}
		}
         	for (UInt_t i=0;i<Get_NumThrown();i++){
                	//Set branch array indices corresponding to this particle
                	dThrownWrapper->Set_ArrayIndex(i);
			Int_t locThrownPID = dThrownWrapper->Get_PID();
			if(locThrownPID==14 || locThrownPID==8 || locThrownPID==9) {
				locFinalStateThrownP4.push_back(dThrownWrapper->Get_P4());
				if(locThrownPID==14) locRecoilP4 += dThrownWrapper->Get_P4();
				else {
					// add omega decay charged pions
					if(locFinalStateThrownP4.size() == 4) {
						locOmegaP4 += dThrownWrapper->Get_P4(); // pi+
					}
					if(locFinalStateThrownP4.size() == 5) {
						locOmegaP4 += dThrownWrapper->Get_P4(); // pi-
					}
				}
			}

			// set pi0 in the proper order for AmpTools
			if(locFinalStateThrownP4.size() == 2) locFinalStateThrownP4.push_back(locPi0ThrownP4);
			if(locFinalStateThrownP4.size() == 6) {
				locRecoilP4 += dThrownWrapper->Get_P4();
				break;
			}
	    	}
	}

	// fill generated from Thrown_Tree
	if(dOption.Contains("thrown") && dTreeInterface->Get_Branch("NumCombos") == NULL) {
		dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", 1.0);
		dFlatTreeInterface->Fill_Fundamental<Float_t>("MRecoil", locRecoilP4.M());
		TLorentzVector locBeamP4 = dThrownBeam->Get_P4();	
		FillAmpTools_FlatTree(locBeamP4, locFinalStateThrownP4);
		Fill_FlatTree();
		
		return kTRUE;
	}

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
        //For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
        Reset_Actions_NewEvent();

	/************************************************* LOOP OVER COMBOS *************************************************/

	// keep track of combo with best ChiSquare
	double locBestChiSq = 9e9;
	vector< vector<TLorentzVector> > locFinalStateBestChiSqP4;
	vector< TLorentzVector > locRecoilBestChiSqP4, locBeamBestChiSqP4;
	vector< double > locWeightBestChiSq, locM3PiBestChiSq, locM4PiBestChiSq, loctBestChiSq; 

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		dComboWrapper->Set_IsComboCut(false);
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

		// Loop through the analysis actions, executing them in order for the active particle combo
                if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
                        continue;

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
		double locKinFit_RedChiSq = dComboWrapper->Get_ChiSq_KinFit("")/dComboWrapper->Get_NDF_KinFit("");
                if(locKinFit_CL < 1e-6) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
			
		if(locBeamP4.E()<8.2 || locBeamP4.E()>8.8) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		dHist_KinFitChiSq->Fill(locKinFit_RedChiSq, locAccWeight);
		if(locKinFit_RedChiSq > 5) continue;

		// reset best ChiSquare vectors 
		if(dOption.Contains("bestChiSq")) {
			if(locAccid) continue;
			if(locKinFit_RedChiSq < locBestChiSq) {
				locFinalStateBestChiSqP4.clear();
				locRecoilBestChiSqP4.clear(); locBeamBestChiSqP4.clear();
				locWeightBestChiSq.clear(); locM3PiBestChiSq.clear();
				locM4PiBestChiSq.clear(); loctBestChiSq.clear();
			}
		}

		dHist_3piPiPlus1Mass->Fill(loc3PiMass11);
		dHist_3piPiPlus1Mass->Fill(loc3PiMass12);
		dHist_3piPiPlus1MassCorr->Fill(loc3PiMass12, loc3PiMass11);
		dHist_3piPiPlus2Mass->Fill(loc3PiMass21);
		dHist_3piPiPlus2Mass->Fill(loc3PiMass22);

		dHist_ProtonPiPlus2Mass->Fill(locProtonPiPlus2P4.M());
                dHist_ProtonPiPlus1Mass->Fill(locProtonPiPlus1P4.M());

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
				TLorentzVector loc3PiP4_alt = locPi0P4 + locPiPlusP4[iplus] + locPiMinusP4[abs(iminus-1)];
                                double loc3PiMass_alt = loc3PiP4_alt.M();	
	
				// Delta++ cut
				TLorentzVector locRecoilP4 = locProtonP4+locPiPlusP4[abs(iplus-1)];
				if(locRecoilP4.M() > 1.6) 
					continue; 

				double weight = 0;
				dHist_2gammaMassSum->Fill(locPi0P4.M());
				if(fabs(locPi0P4.M() - 0.135) < 0.015) 
					weight = 1.0;								
				else if((locPi0P4.M() > 0.09 && locPi0P4.M() < 0.105) || (locPi0P4.M() > 0.15 && locPi0P4.M() < 0.165))
					weight = -1.0;

				dHist_3piMassSum->Fill(loc3PiMass, weight);

				// t cut
				double loct = -1. * (locRecoilP4 - dTargetP4).M2();
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
	
				if(weight < 1.0) continue;

				double loc4PiMass = (loc3PiP4 + locPiMinusP4[abs(iminus-1)]).M();
				dHist_4piMassSum->Fill(loc4PiMass, loc2Dweight);

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

				cout<<locOmegaP4.M()<<" "<<loc3PiMass<<" "<<loc3PiMass_alt<<" "<<loc2Dweight<<endl;

				// set weight according to omega mass
				dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", loc2Dweight*locAccWeight);
				// set ordered final state P4 for filling flat tree
				if(!dOption.Contains("accept")) 
					FillAmpTools_FlatTree(locBeamP4, locFinalStateP4); 
				else {
					dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", 1.0);
					FillAmpTools_FlatTree(locBeamP4, locFinalStateThrownP4);
				}
					

				if(!dOption.Contains("bestChiSq"))
					Fill_FlatTree(); //for the active combo
				else {
					locFinalStateBestChiSqP4.push_back(locFinalStateP4);
					locRecoilBestChiSqP4.push_back(locRecoilP4);
					locBeamBestChiSqP4.push_back(locBeamP4);
					locWeightBestChiSq.push_back(loc2Dweight);
					loctBestChiSq.push_back(loct);
					locM3PiBestChiSq.push_back(loc3PiMass);
					locM4PiBestChiSq.push_back(loc4PiMass);
				}
			}
		}

	} // end of combo loop

	// fill with a single combo per event with the best ChiSquare (reduce accidentals)
	if(dOption.Contains("bestChiSq")) {
		
		// loop over omega candidates within combo
		for(uint i=0; i<locFinalStateBestChiSqP4.size(); i++) {
			dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", locWeightBestChiSq[i]);
			dFlatTreeInterface->Fill_Fundamental<Float_t>("MRecoil", locRecoilBestChiSqP4[i].M());
			dFlatTreeInterface->Fill_Fundamental<Float_t>("M3Pi", locM3PiBestChiSq[i]);
			dFlatTreeInterface->Fill_Fundamental<Float_t>("M4Pi", locM4PiBestChiSq[i]);
			dFlatTreeInterface->Fill_Fundamental<Float_t>("Phi_Prod", locRecoilBestChiSqP4[i].Phi());
			dFlatTreeInterface->Fill_Fundamental<Float_t>("t", loctBestChiSq[i]);
		
			FillAmpTools_FlatTree(locBeamBestChiSqP4[i], locFinalStateBestChiSqP4[i]);
			Fill_FlatTree();
		}
		
		return kTRUE;
	}

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
	if(!locIsEventCut && dOutputTreeFileName != "") {
		cout<<"filled tree good entry "<<locEntry<<endl;
		cout<<dOutputTreeFileName<<endl;
		//eventCounter++;
		Fill_OutputTree();
	}

	return kTRUE;
}

void DSelector_omegapi_all::Finalize(void)
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
