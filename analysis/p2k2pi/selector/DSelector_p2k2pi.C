#include "DSelector_p2k2pi.h"

void DSelector_p2k2pi::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "p2k2pi.root"; //"" for none
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
	dAnalysisActions.push_back(new DCutAction_dEdx(dComboWrapper, false, Proton, SYS_CDC));
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

	//MISSING MASS
	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

        //KINFIT RESULTS
        dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//MASSES
        deque<Particle_t> locPhiPIDs;  locPhiPIDs.push_back(KPlus); locPhiPIDs.push_back(KMinus);
	double minPhi = 0.95; double maxPhi = 1.2;
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locPhiPIDs, 100, minPhi, maxPhi, "Phi_NoCut"));

        //CUT MISSING MASS
        dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));
        dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.001));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locPhiPIDs, 100, minPhi, maxPhi, "Phi_Kin0.001"));

	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.01));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locPhiPIDs, 100, minPhi, maxPhi, "Phi_Kin0.01"));

	//dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.1));
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locPhiPIDs, 100, minPhi, maxPhi, "Phi_Kin0.1"));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	int n2k = 90; double min2k = 0.9; double max2k = 1.2;
	int n2pi = 200; double min2pi = 0.0; double max2pi = 2.0;
	int nkpi = 150; double minkpi = 0.5; double maxkpi = 2.0;
	int nt = 100; double mint = 0; double maxt = 5.0;
	int nDelta = 100; double minDelta = 1.0; double maxDelta = 3.0;
	int n2kpi = 100; double min2kpi = 1.0; double max2kpi = 3.0;
	int nphi2pi = 100; double minphi2pi = 1.0; double maxphi2pi = 3.0;

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_2piMass_deltaR = new TH2I("h2piMass_deltaR", "2piMass_deltaR", 200, 0., 5., n2pi, min2pi, max2pi);
	dHist_2piMass_deltaZ = new TH2I("h2piMass_deltaZ", "2piMass_deltaZ", 200, -20., 20., n2pi, min2pi, max2pi);
	dHist_deltaR_deltaZ = new TH2I("hdeltaR_deltaZ", "deltaR_deltaZ", 200, -20., 20., 200, 0., 5.);

	dHist_2kMass_2piMass = new TH2I("h2kMass_2piMass", "2kMass_2piMass", n2pi, min2pi, max2pi, n2k, min2k, max2k);
	dHist_kppimMass_kmpipMass = new TH2I("hkppimMass_kmpipMass", "kppimMass_kmpipMass", nkpi, minkpi, maxkpi, nkpi, minkpi, maxkpi);
	dHist_2kMass_t = new TH2I("h2kMass_t", "2kMass_t", nt, mint, maxt, n2k, min2k, max2k);

	dHist_2kMass_2piMass_vertex = new TH2I("h2kMass_2piMass_vertex", "2kMass_2piMass_vertex", n2pi, min2pi, max2pi, n2k, min2k, max2k);
	dHist_kppimMass_kmpipMass_vertex = new TH2I("hkppimMass_kmpipMass_vertex", "kppimMass_kmpipMass_vertex", nkpi, minkpi, maxkpi, nkpi, minkpi, maxkpi);
	dHist_2kMass_t_vertex = new TH2I("h2kMass_t_vertex", "2kMass_t_vertex", nt, mint, maxt, n2k, min2k, max2k);

	dHist_2kMass_2piMass_noKstar = new TH2I("h2kMass_2piMass_noKstar", "2kMass_2piMass_noKstar", n2pi, min2pi, max2pi, n2k, min2k, max2k);
	dHist_2kMass_ProtonPimMass_noKstar = new TH2I("h2kMass_ProtonPimMass_noKstar", "2kMass_ProtonPimMass_noKstar", nDelta, minDelta, maxDelta, n2k, min2k, max2k);
	dHist_2kMass_ProtonPipMass_noKstar = new TH2I("h2kMass_ProtonPipMass_noKstar", "2kMass_ProtonPipMass_noKstar", nDelta, minDelta, maxDelta, n2k, min2k, max2k);

	dHist_2kMass_2kpimMass_noKstar = new TH2I("h2kMass_2kpimMass_noKstar", "2kMass_2kpimMass_noKstar", n2kpi, min2kpi, max2kpi, n2k, min2k, max2k);

	dHist_2kMass_2piMass_noDelta = new TH2I("h2kMass_2piMass_noDelta", "2kMass_2piMass_noDelta", n2pi, min2pi, max2pi, n2k, min2k, max2k);

	dHist_phi2piMass_2piMass_noDelta = new TH2I("hphi2piMass_2piMass_noDelta", "phi2piMass_2piMass_noDelta", n2pi, min2pi, max2pi, nphi2pi, minphi2pi, maxphi2pi);
	dHist_phi2piMass_2piMass_noDelta_SB = new TH2I("hphi2piMass_2piMass_noDelta_SB", "phi2piMass_2piMass_noDelta_SB", n2pi, min2pi, max2pi, nphi2pi, minphi2pi, maxphi2pi);
		
	dHist_2k2piMass_t = new TH2I("h2k2piMass_t", "2k2piMass_t", nt, mint, maxt, nphi2pi, minphi2pi, maxphi2pi);
	dHist_2k2piMass_t_SB = new TH2I("h2k2piMass_t_SB", "2k2piMass_t_SB", nt, mint, maxt, nphi2pi, minphi2pi, maxphi2pi);

	dHist_2k2piMass_Egamma = new TH2I("h2k2piMass_Egamma", "2k2piMass_Egamma", 200, 0., 12., nphi2pi, minphi2pi, maxphi2pi);
	dHist_2k2piMass_Egamma_SB = new TH2I("h2k2piMass_Egamma_SB", "2k2piMass_Egamma_SB", 200, 0., 12., nphi2pi, minphi2pi, maxphi2pi);

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

Bool_t DSelector_p2k2pi::Process(Long64_t locEntry)
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
		Int_t locKPlusTrackID = dKPlusWrapper->Get_TrackID();
		Int_t locKMinusTrackID = dKMinusWrapper->Get_TrackID();
		Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
		TLorentzVector locKMinusP4 = dKMinusWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locKPlusP4_Measured + locKMinusP4_Measured + locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured;

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/************************************ MY ANALYSIS STUFF ************************************/

		TLorentzVector lockppimP4_Measured = locKPlusP4_Measured + locPiMinusP4_Measured;
		TLorentzVector lockmpipP4_Measured = locKMinusP4_Measured + locPiPlusP4_Measured;
		TLorentzVector loc2kP4_Measured = locKPlusP4_Measured + locKMinusP4_Measured;
		TLorentzVector loc2piP4_Measured = locPiPlusP4_Measured + locPiMinusP4_Measured;

		TLorentzVector lockppimP4 = locKPlusP4 + locPiMinusP4;
		TLorentzVector lockmpipP4 = locKMinusP4 + locPiPlusP4;
		TLorentzVector loc2kP4 = locKPlusP4 + locKMinusP4;
		TLorentzVector loc2piP4 = locPiPlusP4 + locPiMinusP4;
		TLorentzVector loc2k2piP4 = locKPlusP4 + locKMinusP4 + locPiPlusP4 + locPiMinusP4;
		TLorentzVector locProtonPipP4 = locProtonP4 + locPiPlusP4;
		TLorentzVector locProtonPimP4 = locProtonP4 + locPiMinusP4;
		TLorentzVector loc2kpimP4 = locKPlusP4 + locKMinusP4 + locPiMinusP4;

		TLorentzVector loc2piX4_Measured = dPiPlusWrapper->Get_X4_Measured() + dPiMinusWrapper->Get_X4_Measured();
		loc2piX4_Measured *= 0.5;
		TLorentzVector locProtonX4_Measured = dProtonWrapper->Get_X4_Measured();
		TLorentzVector locDeltaX4_Measured = loc2piX4_Measured - locProtonX4_Measured;
		
		double t = (dTargetP4 - locProtonP4).M2();

		bool tofGOOD = false;
                if((locKPlusP4.Theta()*180./TMath::Pi() < 11. && locKPlusP4.P() < 1.5) || (locKMinusP4.Theta()*180/TMath::Pi() < 11. && locKMinusP4.P() < 1.5))
                        tofGOOD = true;
                //if(!tofGOOD)
                //      continue;
		
		// AFTER KINFIT CL CUT
		dHist_2piMass_deltaR->Fill(locDeltaX4_Measured.Perp(), loc2piP4_Measured.M());
		dHist_2piMass_deltaZ->Fill(locDeltaX4_Measured.Z(), loc2piP4_Measured.M());
		dHist_deltaR_deltaZ->Fill(locDeltaX4_Measured.Z(), locDeltaX4_Measured.Perp());
		dHist_2kMass_2piMass->Fill(loc2piP4.M(), loc2kP4.M());
		dHist_kppimMass_kmpipMass->Fill(lockmpipP4.M(), lockppimP4.M());
		dHist_2kMass_t->Fill(-t, loc2kP4_Measured.M());

		// REJECT KSHORT WITH VERTEX CUT
		if(locDeltaX4_Measured.Perp() > 1.0 || fabs(locDeltaX4_Measured.Z()) > 5)
			continue;

		dHist_2kMass_2piMass_vertex->Fill(loc2piP4.M(), loc2kP4.M());
		dHist_kppimMass_kmpipMass_vertex->Fill(lockmpipP4.M(), lockppimP4.M());
		dHist_2kMass_t_vertex->Fill(-t, loc2kP4.M());

		// REJECT KSTAR WITH MASS CUT
		if(fabs(lockmpipP4.M() - 0.892) < 0.1 || fabs(lockppimP4.M() - 0.892) < 0.1)
			continue;

		dHist_2kMass_2piMass_noKstar->Fill(loc2piP4.M(), loc2kP4.M());
		dHist_2kMass_ProtonPimMass_noKstar->Fill(locProtonPimP4.M(), loc2kP4.M());
		dHist_2kMass_ProtonPipMass_noKstar->Fill(locProtonPipP4.M(), loc2kP4.M());

		// REQUIRE SMALL -t OFF PROTON (meson production)
		//if(-t > 0.5)
		//	continue;
		
		// REQUIRE NO DELTA++
		if(locProtonPipP4.M() < 1.4) {
			dHist_2kMass_2kpimMass_noKstar->Fill(loc2kpimP4.M(), loc2kP4.M());
			continue;
		}

		dHist_2kMass_2piMass_noDelta->Fill(loc2piP4.M(), loc2kP4.M());

		// CUT PHI MASS (sideband subtraction)
		if(loc2kP4.M() < 1.03 && loc2kP4.M() > 1.01) {
		   dHist_phi2piMass_2piMass_noDelta->Fill(loc2piP4.M(), loc2k2piP4.M());
		   if(loc2piP4.M() > 0.9 && loc2piP4.M() < 1.05) {
			   dHist_2k2piMass_t->Fill(-t, loc2k2piP4.M());
			   dHist_2k2piMass_Egamma->Fill(locBeamP4.E(), loc2k2piP4.M());
		   }
		}
		else if(loc2kP4.M() < 1.12 && loc2kP4.M() > 1.10) {
			dHist_phi2piMass_2piMass_noDelta_SB->Fill(loc2piP4.M(), loc2k2piP4.M());
			if(loc2piP4.M() > 0.9 && loc2piP4.M() < 1.05) {
				dHist_2k2piMass_t_SB->Fill(-t, loc2k2piP4.M());
				dHist_2k2piMass_Egamma_SB->Fill(locBeamP4.E(), loc2k2piP4.M());
			}
		}

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

	return kTRUE;
}

void DSelector_p2k2pi::Finalize(void)
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
