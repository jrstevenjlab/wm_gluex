#include "DSelector_p2kpi0.h"

void DSelector_p2kpi0::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "2kpi0.root"; //"" for none
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
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.3, KPlus, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, KPlus, SYS_BCAL));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.3, KMinus, SYS_TOF));
        dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, KMinus, SYS_BCAL));

        //MISSING MASS
        dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));
	dAnalysisActions.push_back(new DHistogramAction_Energy_UnusedShowers(dComboWrapper));
	dAnalysisActions.push_back(new DCutAction_Energy_UnusedShowers(dComboWrapper, 0.5));
	dAnalysisActions.push_back(new DHistogramAction_Energy_UnusedShowers(dComboWrapper, "AfterCut"));

	//KINFIT RESULTS
        dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

        //CUT MISSING MASS
        dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));
        dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.01));

        //CUT PI0 MASS
        deque<Particle_t> locPi0PIDs; locPi0PIDs.push_back(Gamma); locPi0PIDs.push_back(Gamma);
        dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Pi0, 100, 0.09, 0.18, "Pi0_KinCut"));
        dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, false, Pi0, 0.114, 0.156));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	int n2k = 500; double min2k = 0.9; double max2k = 1.9;
        int nkpi = 450; double minkpi = 0.5; double maxkpi = 2.0;
        int nt = 100; double mint = 0; double maxt = 5.0;
        int nDelta = 100; double minDelta = 1.0; double maxDelta = 3.0;
        int n2kpi = 400; double min2kpi = 1.0; double max2kpi = 3.0;

        //EXAMPLE MANUAL HISTOGRAMS:
        dHist_kppiMass_kmpiMass = new TH2I("hkppiMass_kmpiMass", "kppiMass_kmpiMass", nkpi, minkpi, maxkpi, nkpi, minkpi, maxkpi);
        dHist_2kMass_t = new TH2I("h2kMass_t", "2kMass_t", nt, mint, maxt, n2k, min2k, max2k);

        dHist_2kMass_ProtonPiMass_noKstar = new TH2I("h2kMass_ProtonPiMass_noKstar", "2kMass_ProtonPiMass_noKstar", nDelta, minDelta, maxDelta, n2k, min2k, max2k);

        dHist_2kMass_2kpiMass_noKstar = new TH2I("h2kMass_2kpimMass_noKstar", "2kMass_2kpimMass_noKstar", 200, 1., 3., n2k, min2k, max2k);

	dHist_2kpiMass_t = new TH2I("h2kpiMass_t", "2kpiMass_t", nt, mint, maxt, n2kpi, min2kpi, max2kpi);
        dHist_2kpiMass_t_SB = new TH2I("h2kpiMass_t_SB", "2kpiMass_t_SB", nt, mint, maxt, n2kpi, min2kpi, max2kpi);

        dHist_2kpiMass_Egamma = new TH2I("h2kpiMass_Egamma", "2kpiMass_Egamma", 200, 0., 12., n2kpi, min2kpi, max2kpi);
        dHist_2kpiMass_Egamma_SB = new TH2I("h2kpiMass_Egamma_SB", "2kpiMass_Egamma_SB", 200, 0., 12., n2kpi, min2kpi, max2kpi);

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

Bool_t DSelector_p2kpi0::Process(Long64_t locEntry)
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
		Int_t locKPlusTrackID = dKPlusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locKMinusP4 = dKMinusWrapper->Get_P4();
		TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locKMinusP4_Measured + locKPlusP4_Measured + locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

		/*********************************************** GET Space-Time **********************************************/
                TLorentzVector locBeam_X4_Measured = dComboBeamWrapper->Get_X4_Measured();
                TLorentzVector locProton_X4_Measured = dProtonWrapper->Get_X4_Measured();
		
		double locKinFitCL = dComboWrapper->Get_ConfidenceLevel_KinFit();	

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/************************************ MY ANALYSIS STUFF ************************************/

                TLorentzVector lockppi0P4 = locKPlusP4 + locDecayingPi0P4;
                TLorentzVector lockmpi0P4 = locKMinusP4 + locDecayingPi0P4;
                TLorentzVector loc2kP4 = locKPlusP4 + locKMinusP4;
                TLorentzVector loc2kpi0P4 = locKPlusP4 + locKMinusP4 + locDecayingPi0P4;
                TLorentzVector locProtonPi0P4 = locProtonP4 + locDecayingPi0P4;
		TLorentzVector locProtonKMinusP4 = locProtonP4 + locKMinusP4;

                double t = (dTargetP4 - locProtonP4).M2();

		// AFTER KINFIT and PI0 CUTS
                dHist_kppiMass_kmpiMass->Fill(lockmpi0P4.M(), lockppi0P4.M());
                dHist_2kMass_t->Fill(-t, loc2kP4.M());

                // REJECT KSTAR WITH MASS CUT
                //if(fabs(lockmpi0P4.M() - 0.892) < 0.1 || fabs(lockppi0P4.M() - 0.892) < 0.1)
                //        continue;

                dHist_2kMass_ProtonPiMass_noKstar->Fill(locProtonPi0P4.M(), loc2kP4.M());
		dHist_2kMass_2kpiMass_noKstar->Fill(loc2kpi0P4.M(), loc2kP4.M());

                // CUT PHI MASS (sideband subtraction)
                if(loc2kP4.M() < 1.04 && loc2kP4.M() > 1.00) {
			dHist_2kpiMass_t->Fill(-t, loc2kpi0P4.M());
			dHist_2kpiMass_Egamma->Fill(locBeamP4.E(), loc2kpi0P4.M());
                }
                else if(loc2kP4.M() > 1.08 && loc2kP4.M() < 1.12) {
			dHist_2kpiMass_t_SB->Fill(-t, loc2kpi0P4.M());
			dHist_2kpiMass_Egamma_SB->Fill(locBeamP4.E(), loc2kpi0P4.M());
                }

		/********************************************* RF time and Beam time **********************************************/
		double tRF = 0.0;
                tRF = dComboWrapper->Get_RFTime_Measured();

                TLorentzVector P4_Beam = locBeamP4_Measured;
                TLorentzVector X4_Beam = locBeam_X4_Measured;

                double locBeamDeltaT = X4_Beam.T() - (tRF + (X4_Beam.Z() - dTargetCenter.Z())/29.9792458);

                /****************************************** FILL FLAT TREE ******************************************/
		//FILL ANY CUSTOM BRANCHES FIRST!!
		double weight = 1.;

                // weights for 2 deltaT sideband bins
                Bool_t AcciWeight  = (fabs(locBeamDeltaT)<(2.5*4.008)&&fabs(locBeamDeltaT)>(0.5*4.008)) ? true : false;
                if(AcciWeight) weight = -1/2.;

                // remove any further deltaT sidebands if they exist
                 if(fabs(locBeamDeltaT) > 2.5*4.008)
                        continue;

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E());
			locUsedSoFar_BeamEnergy.insert(locBeamID);
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

void DSelector_p2kpi0::Finalize(void)
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
