#include "DSelector_ppipkskm.h"

void DSelector_ppipkskm::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "ppipkskm.root"; //"" for none
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
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, KPlus, SYS_BCAL));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KMinus, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, KMinus, SYS_BCAL));

	deque<Particle_t> locKstarMPIDs;  locKstarMPIDs.push_back(KMinus); locKstarMPIDs.push_back(PiPlus);
	deque<Particle_t> locXPIDs;  locXPIDs.push_back(KMinus); locXPIDs.push_back(PiPlus); locXPIDs.push_back(KShort);

	//MISSING MASS
        dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, KShort, 100, 0.3, 0.7, "Kshort_NoCut"));
	
        //KINFIT RESULTS
        dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

        //CUT MISSING MASS
        dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));
        dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.001));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, KShort, 100, 0.3, 0.7, "Kshort_Kin0.001"));

	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.01));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, KShort, 100, 0.3, 0.7, "Kshort_Kin0.01"));

	//CUT KSHORT MASS
	dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, false, KShort, 0.48, 0.52));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locKstarMPIDs, locXPIDs, 100, 0.5, 2.5, 100, 1.0, 3.0, "2D_KShortCut"));
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false, "Kshort_Cut"));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	int n2pi = 200; double min2pi = 0.3; double max2pi = 0.7;
        int nkpi = 150; double minkpi = 0.5; double maxkpi = 2.0;
        int n2kpi = 100; double min2kpi = 1.0; double max2kpi = 3.0;

	//EXAMPLE MANUAL HISTOGRAMS:
        dHist_2piMass_deltaR = new TH2I("h2piMass_deltaR", "2piMass_deltaR", 200, 0., 5., n2pi, min2pi, max2pi);
        dHist_2piMass_deltaZ = new TH2I("h2piMass_deltaZ", "2piMass_deltaZ", 200, -20., 20., n2pi, min2pi, max2pi);
        dHist_deltaR_deltaZ = new TH2I("hdeltaR_deltaZ", "deltaR_deltaZ", 200, -20., 20., 200, 0., 5.);

	dHist_2piMass_kppimMass = new TH2I("h2piMass_kppimMass", "h2piMass_kppimMass", nkpi, minkpi, maxkpi, n2pi, min2pi, max2pi);
        dHist_2kpiMass_kppimMass = new TH2I("h2kpiMass_kppimMass", "h2kpiMass_kppimMass", nkpi, minkpi, maxkpi, n2kpi, min2kpi, max2kpi);
	dHist_2kpiMass_kppimMass_SB1 = new TH2I("h2kpiMass_kppimMass_SB1", "h2kpiMass_kppimMass_SB1", nkpi, minkpi, maxkpi, n2kpi, min2kpi, max2kpi);
	dHist_2kpiMass_kppimMass_SB2 = new TH2I("h2kpiMass_kppimMass_SB2", "h2kpiMass_kppimMass_SB2", nkpi, minkpi, maxkpi, n2kpi, min2kpi, max2kpi);
	
	dHist_2piMass_kppimMass_vertex = new TH2I("h2piMass_kppimMass_vertex", "h2piMass_kppimMass_vertex", nkpi, minkpi, maxkpi, n2pi, min2pi, max2pi);
        dHist_2kpiMass_kppimMass_vertex = new TH2I("h2kpiMass_kppimMass_vertex", "h2kpiMass_kppimMass_vertex", nkpi, minkpi, maxkpi, n2kpi, min2kpi, max2kpi);
	dHist_2kpiMass_kppimMass_vertex_SB1 = new TH2I("h2kpiMass_kppimMass_vertex_SB1", "h2kpiMass_kppimMass_vertex_SB1", nkpi, minkpi, maxkpi, n2kpi, min2kpi, max2kpi);
	dHist_2kpiMass_kppimMass_vertex_SB2 = new TH2I("h2kpiMass_kppimMass_vertex_SB2", "h2kpiMass_kppimMass_vertex_SB2", nkpi, minkpi, maxkpi, n2kpi, min2kpi, max2kpi);

	/***************************************** ADVANCED: CHOOSE BRANCHES TO READ ****************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_ppipkskm::Process(Long64_t locEntry)
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
		Int_t locKMinusTrackID = dKMinusWrapper->Get_TrackID();
		Int_t locPiPlus1TrackID = dPiPlus1Wrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPiPlus2TrackID = dPiPlus2Wrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locKMinusP4 = dKMinusWrapper->Get_P4();
		TLorentzVector locPiPlus1P4 = dPiPlus1Wrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locDecayingKShortP4 = dDecayingKShortWrapper->Get_P4();
		TLorentzVector locPiPlus2P4 = dPiPlus2Wrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locPiPlus1P4_Measured = dPiPlus1Wrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPiPlus2P4_Measured = dPiPlus2Wrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

			/************************************ MY ANALYSIS STUFF ************************************/

                TLorentzVector loc2piP4_Measured = locPiMinusP4_Measured + locPiPlus2P4_Measured;

		TLorentzVector lockppimP4 = locKMinusP4 + locPiPlus1P4;
                TLorentzVector loc2kpiP4 = locKMinusP4 + locPiPlus1P4 + locDecayingKShortP4;

		TLorentzVector loc2piX4_Measured = dPiPlus1Wrapper->Get_X4_Measured() + dPiMinusWrapper->Get_X4_Measured();
                loc2piX4_Measured *= 0.5;
                TLorentzVector locProtonX4_Measured = dProtonWrapper->Get_X4_Measured();
                TLorentzVector locDeltaX4_Measured = loc2piX4_Measured - locProtonX4_Measured;

		// AFTER KINFIT CL CUT
                dHist_2piMass_deltaR->Fill(locDeltaX4_Measured.Perp(), loc2piP4_Measured.M());
                dHist_2piMass_deltaZ->Fill(locDeltaX4_Measured.Z(), loc2piP4_Measured.M());
                dHist_deltaR_deltaZ->Fill(locDeltaX4_Measured.Z(), locDeltaX4_Measured.Perp());
		
                dHist_2piMass_kppimMass->Fill(lockppimP4.M(), loc2piP4_Measured.M());
		if(loc2piP4_Measured.M() > 0.48 && loc2piP4_Measured.M() < 0.52) {
			dHist_2kpiMass_kppimMass->Fill(lockppimP4.M(), loc2kpiP4.M());
		}
		else if(loc2piP4_Measured.M() > 0.56 && loc2piP4_Measured.M() < 0.6) {
			dHist_2kpiMass_kppimMass_SB1->Fill(lockppimP4.M(), loc2kpiP4.M());
		}
		else if(loc2piP4_Measured.M() > 0.40 && loc2piP4_Measured.M() < 0.44) {
			dHist_2kpiMass_kppimMass_SB2->Fill(lockppimP4.M(), loc2kpiP4.M());
		}

		// REJECT KSHORT WITH VERTEX CUT
                if(locDeltaX4_Measured.Perp() < 0.3 && locDeltaX4_Measured.Z() < 1.5)
                        continue;

		dHist_2piMass_kppimMass_vertex->Fill(lockppimP4.M(), loc2piP4_Measured.M());
		if(loc2piP4_Measured.M() > 0.48 && loc2piP4_Measured.M() < 0.52) {
			dHist_2kpiMass_kppimMass_vertex->Fill(lockppimP4.M(), loc2kpiP4.M());
		}
		else if(loc2piP4_Measured.M() > 0.56 && loc2piP4_Measured.M() < 0.6) {
			dHist_2kpiMass_kppimMass_vertex_SB1->Fill(lockppimP4.M(), loc2kpiP4.M());
		}
		else if(loc2piP4_Measured.M() > 0.40 && loc2piP4_Measured.M() < 0.44) {
			dHist_2kpiMass_kppimMass_vertex_SB2->Fill(lockppimP4.M(), loc2kpiP4.M());
		}


	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	return kTRUE;
}

void DSelector_ppipkskm::Finalize(void)
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
