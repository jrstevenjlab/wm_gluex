#include "DSelector_omega_misspim.h"

void DSelector_omega_misspim::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "omega_misspim.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen

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
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	//dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1F("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.16, 0.16);
	dHist_MissingMass = new TH1F("MissingMass", ";Missing Mass (GeV)", 600, 0., 0.5);
	dHist_BeamEnergy = new TH1F("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
	dHist_MissingEnergy = new TH1F("MissingEnergy", ";Missing Energy (GeV)", 600, -4., 4.);

	dHist_BeamDeltaT = new TH1F("BeamDeltaT", "; t_{Tagger} - t_{RF} (ns)", 200, -20., 20.);
	dHist_UnusedShowerEnergy = new TH1F("UnusedShowerEnergy", ";E_{Unused} (GeV)", 600, 0., 10.);
	dHist_NumTracks = new TH1F("NumTracks", ";Number Unique Tracks", 10, 0, 10);
	dHist_KinFitCL = new TH1F("KinFit_CL", ";Kinematic Fit Confidence Level", 600, 0., 1.);

	dHist_2GammaMass = new TH1F("2GammaMass", ";M_{#gamma#gamma} (GeV)", 600, 0.1, 0.17);
	dHist_OmegaMass_missing = new TH1F("OmegaMass_missing", ";Missing Mass off the Proton (GeV)", 600, 0.3, 1.3);

	for(int i = 0; i < 10; i++) {
	  dHist_3PiMass[i] = new TH1F(Form("3PiMass_0%d", i+1), ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV)", 600, 0.3, 1.3);
	  dHist_OmegaMass_1track[i] = new TH1F(Form("OmegaMass_1track_0%d", i+1), ";Missing Mass off the Proton (GeV)", 600, 0.3, 1.3);
	  dHist_MassCorrelation[i] = new TH2F(Form("MassCorrelation_0%d", i+1), ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);Missing Mass off the Proton (GeV)", 600, 0.3, 1.3, 600, 0.3, 1.3);
	}

	dHist_DeltaPhiMissing = new TH1F("DeltaPhiMissing", ";#phi_{reco} - #phi_{miss} (deg)", 600, -180., 180.);
	dHist_DeltaThetaMissing = new TH1F("DeltaThetaMissing", ";#theta_{reco} - #theta_{miss} (deg)", 600, -180., 180.);
	dHist_DeltaPMissing = new TH1F("DeltaPMissing", ";p_{reco} - p_{miss} (GeV)", 600, -10.0, 10.0);
	dHist_DeltaPhiVsDeltaTheta_missing = new TH2F("DeltaPhiVsDeltaTheta_missing", "#Delta#phi vs #Delta#theta;#phi_{reco} - #phi_{miss} (deg);#theta_{reco} - #theta_{miss} (deg)", 600, -60., 60., 600, -60., 60.);
	dHist_DeltaPVsDeltaPhi_missing = new TH2F("DeltaPVsDeltaPhi_missing", "#Delta p vs #Delta#phi;p_{reco} - p_{miss} (GeV);#phi_{reco} - #phi_{miss} (deg)", 600, -2.0, 2.0, 600, -60., 60.);
	dHist_DeltaPVsDeltaTheta_missing = new TH2F("DeltaPVsDeltaTheta_missing", "#Delta p vs #Delta#theta;p_{reco} - p_{miss} (GeV);#theta_{reco} - #theta_{miss} (deg)", 600, -2.0, 2.0, 600, -60., 60.);

	dHist_DeltaPhiTruth = new TH1F("DeltaPhiTruth", ";#phi_{reco} - #phi_{truth} (deg)", 600, -180., 180.);
	dHist_DeltaThetaTruth = new TH1F("DeltaThetaTruth", ";#theta_{reco} - #theta_{truth} (deg)", 600, -180., 180.);
	dHist_DeltaPTruth = new TH1F("DeltaPTruth", ";p_{reco} - p_{truth} (GeV)", 600, -10., 10.);

	dHist_OmegaMassVsPhi_0track = new TH2F("OmegaMassVsPhi_0tracks", ";Missing Mass off the Proton (GeV);#phi_{missing} (deg)", 600, 0.3, 1.3, 600, -180., 180.);
	dHist_OmegaMassVsTheta_0track = new TH2F("OmegaMassVsTheta_0tracks", ";Missing Mass off the Proton (GeV);#theta_{missing} (deg)", 600, 0.3, 1.3, 600, 0., 30.);
	dHist_OmegaMassVsP_0track = new TH2F("OmegaMassVsP_0tracks", ";Missing Mass off the Proton (GeV);p_{missing} (GeV)", 600, 0.3, 1.3, 600, 0., 10.);

	dHist_OmegaMassVsPhi_1track = new TH2F("OmegaMassVsPhi_1track", ";Missing Mass off the Proton (GeV);#phi_{mmop} (deg)", 600, 0.3, 1.3, 600, -180., 180.);
	dHist_OmegaMassVsTheta_1track = new TH2F("OmegaMassVsTheta_1track", ";Missing Mass off the Proton (GeV);#theta_{mmop} (deg)", 600, 0.3, 1.3, 600, 0., 30.);
	dHist_OmegaMassVsP_1track = new TH2F("OmegaMassVsP_1track", ";Missing Mass off the Proton (GeV);p_{mmop} (GeV)", 600, 0.3, 1.3, 600, 0., 10.);
	dHist_3PiMassVsPhi = new TH2F("3PiMassVsPhi", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);#phi_{mmop} (deg)", 600, 0.3, 1.3, 600, -180., 180.);
	dHist_3PiMassVsTheta = new TH2F("3PiMassVsTheta", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);#theta_{mmop} (deg)", 600, 0.3, 1.3, 600, 0., 30.);
	dHist_3PiMassVsP = new TH2F("3PiMassVsP", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);p_{mmop} (GeV)", 600, 0.3, 1.3, 600, 0., 10.);

	for(int i = 0; i < 9; i++){
	  dHist_OmegaMassVsTheta_0track_p[i] = new TH2F(Form("OmegaMassVsTheta_0track_p%d", i+1), ";Missing Mass off the Proton (GeV);#theta_{mmop} (deg)", 600, 0.3, 1.3, 600, 0., 30.);
	  dHist_OmegaMassVsTheta_1track_p[i] = new TH2F(Form("OmegaMassVsTheta_1track_p%d", i+1), ";Missing Mass off the Proton (GeV);#theta_{mmop} (deg)", 600, 0.3, 1.3, 600, 0., 30.);
	  dHist_3PiMassVsTheta_p[i] = new TH2F(Form("3PiMassVsTheta_p%d", i+1), ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);#theta_{mmop} (deg)", 600, 0.3, 1.3, 600, 0., 30.);
	}

	dHist_ThetaVsP_reco = new TH2F("ThetaVsP_reco", ";#theta_{reco} (deg);p_{reco} (GeV)", 600, 0., 30., 600, 0., 10.);
	dHist_ThetaVsP_missing = new TH2F("ThetaVsP_missing", ";#theta_{missing} (deg);p_{missing} (GeV)", 600, 0., 30., 600, 0., 10.);

	dHist_DeltaPhi3D_truth = new TH3F("DeltaPhi3D_truth", "#Delta#phi vs #theta vs p;#phi_{reco} - #phi_{truth} (deg);#theta_{truth} (deg);p_{truth} (GeV/c)", 60, -90., 90., 60, 0., 30., 60, 0., 10.);
	dHist_DeltaTheta3D_truth = new TH3F("DeltaTheta3D_truth", "#Delta#theta vs #theta vs p;#theta_{reco} - #theta_{truth} (deg);#theta_{truth} (deg);p_{truth} (GeV/c)", 60, -30., 30., 60, 0., 30., 60, 0., 10.);
	dHist_DeltaP3D_truth = new TH3F("DeltaP3D_truth", "#Delta p vs #theta vs p;p_{reco} - p_{truth} (GeV/c);#theta_{truth} (deg);p_{truth} (GeV/c)", 60, -3., 3., 60, 0., 30., 60, 0., 10.);

	dHist_DeltaPhi3D_miss = new TH3F("DeltaPhi3D_miss", "#Delta#phi vs #theta vs p;#phi_{reco} - #phi_{miss} (deg);#theta_{miss} (deg);p_{miss} (GeV)", 60, -90., 90., 60, 0., 30., 60, 0., 10.);
	dHist_DeltaTheta3D_miss = new TH3F("DeltaTheta3D_miss", "#Delta#theta vs #theta vs p;#theta_{reco} - #theta_{miss} (deg);#theta_{miss} (deg);p_{miss} (GeV)", 60, -30., 30., 60, 0., 30., 60, 0., 10.);
	dHist_DeltaP3D_miss = new TH3F("DeltaP3D_miss", "#Delta p vs #theta vs p;p_{reco} - p_{miss} (GeV);#theta_{miss} (deg);p_{miss} (GeV)", 60, -3., 3., 60, 0., 30., 60, 0., 10.);

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_omega_misspim::Process(Long64_t locEntry)
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
  // THIS CAUSES A CRASH -- IMPORTANT?
  Reset_Actions_NewEvent();
  //dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()
  
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
  
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_2GammaMass;
	
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_DeltaMissing;

  //Should've done these in an array:
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation01;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation02;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation03;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation04;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation05;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation06;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation07;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation08;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation09;
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MassCorrelation10;
  
  set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaMassVsPhi_1track;

  set<map<Particle_t, set<Int_t> > > locUsedSoFar_DeltaTruth;

  set<map<Particle_t, set<Int_t> > > locUsedSoFar_Delta3D_truth;

  set<map<Particle_t, set<Int_t> > > locUsedSoFar_Delta3D_miss;

  set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaMass;

  set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaMassVsPhi_0track;

  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMassSquared;

  set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingEnergy;

  
  //INSERT USER ANALYSIS UNIQUENESS TRACKING HERE
  
  /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

  /*
    Int_t locMyInt = 7;
    dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

    TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
    dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

    for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
    dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
  */

  /************************************************** LOOP OVER CHARGED TRACK HYPOTHESES- REJECT ANY WITH MORE THAN 3 TRACKS ********************/
  
  vector<Int_t> locUsedSoFar_TrackID;
  int itrack = 0;

  //Loop over charged track hypotheses
  for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
    {
	    
      //Set branch array indices corresponding to this particle
      dChargedHypoWrapper->Set_ArrayIndex(loc_i);
      if (std::find(locUsedSoFar_TrackID.begin(), locUsedSoFar_TrackID.end(), dChargedHypoWrapper->Get_TrackID()) == locUsedSoFar_TrackID.end()) {
	itrack++;
	locUsedSoFar_TrackID.push_back(dChargedHypoWrapper->Get_TrackID());
      }
    }

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
      Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();
      
      //Step 1
      Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
      Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();
      
      /*********************************************** GET FOUR-MOMENTUM **********************************************/
      
      // Get P4's: //is kinfit if kinfit performed, else is measured
      //dTargetP4 is target p4
      //Step 0
      TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
      TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
      TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
      TLorentzVector locMissingPiMinusP4 = dMissingPiMinusWrapper->Get_P4();
      //Step 1
      //TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
      TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
      TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
      
      // Get Measured P4's:
      //Step 0
      TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
      TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
      TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
      //Step 1
      TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
      TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();

      /******************************************* KINFIT CONFIDENCE LEVEL ******************************************************/

      double locKinFit_CL = dComboWrapper->Get_ConfidenceLevel_KinFit("");

      /******************************************** ACCIDENTAL SUBTRACTION INFO *******************************************/
      		
      // measured tagger time for combo
      TLorentzVector locBeam_X4_Measured = dComboBeamWrapper->Get_X4_Measured(); 

      // measured RF time for combo
      double locRFTime = dComboWrapper->Get_RFTime_Measured(); 

      // time difference between tagger and RF (corrected for production vertex position relative to target center)
      double locBeamDeltaT = locBeam_X4_Measured.T() - (locRFTime + (locBeam_X4_Measured.Z() - dTargetCenter.Z())/29.9792458); 
      dHist_BeamDeltaT->Fill(locBeamDeltaT);

      // calculate accidental subtraction weight based on time difference 
      double locAccWeight = 0.; // weight to accidentally subtracted histgorams
      bool locAccid = false; // flag to fill separate prompt and accidental histograms for later subtraction
      
      if(fabs(locBeamDeltaT) < 0.5*4.008) { // prompt signal recieves a weight of 1
	locAccWeight = 1.;
	locAccid = false;
      }
      else { // accidentals recieve a weight of 1/# RF bunches included in TTree (2 in this case)
	locAccWeight = -1./2.;
	locAccid = true;
      }	
      
      /********************************************* COMBINE FOUR-MOMENTUM ********************************************/
      
      // DO YOUR STUFF HERE

      // Combine 4-vectors
      TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
      locMissingP4_Measured -= locPiPlusP4_Measured + locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

      TLorentzVector locOmegaP4_mmop = locBeamP4_Measured + dTargetP4;
      locOmegaP4_mmop -= locProtonP4_Measured;
      double locOmegaMass_mmop = locOmegaP4_mmop.M();		    

      TLorentzVector loc2GammaP4 = locPhoton1P4_Measured + locPhoton2P4_Measured;
      double loc2GammaMass = loc2GammaP4.M();

      double locPiMinusPhi_mmop = locMissingP4_Measured.Phi() * 180./TMath::Pi();
      double locPiMinusTheta_mmop = locMissingP4_Measured.Theta() * 180./TMath::Pi();
      TVector3 locPiMinusP3_mmop = locMissingP4_Measured.Vect();
      double locPiMinusP_mmop = locPiMinusP3_mmop.Mag();
      
      //Loop over charged track hypotheses
      for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
	  //Set branch array indices corresponding to this particle
	  dChargedHypoWrapper->Set_ArrayIndex(loc_i);
		    
	  //Do stuff with the wrapper here ...
	  if (dChargedHypoWrapper->Get_PID() != PiMinus)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locPiPlusTrackID)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locProtonTrackID)
	    continue;
		   

	  TLorentzVector locPiMinusP4_reco = dChargedHypoWrapper->Get_P4_Measured();
	  TLorentzVector loc3PiP4_reco = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locPiMinusP4_reco;
	  double loc3PiMass_reco = loc3PiP4_reco.M();
		   
	  //Mass Correlation with no cuts
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation01;
	  locUsedThisCombo_MassCorrelation01[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation01[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation01[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation01[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation01[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation01[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation01.find(locUsedThisCombo_MassCorrelation01) == locUsedSoFar_MassCorrelation01.end())
	    {
	      dHist_3PiMass[0]->Fill(loc3PiMass_reco, locAccWeight);
	      dHist_OmegaMass_1track[0]->Fill(locOmegaMass_mmop, locAccWeight);
	      dHist_MassCorrelation[0]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight);
	      locUsedSoFar_MassCorrelation01.insert(locUsedThisCombo_MassCorrelation01);
	    }
	}

      
      /********************************************* HISTOGRAM UNUSED SHOWER ENERGY *****************************************/
      
      double locUnusedShowerEnergy = dComboWrapper->Get_Energy_UnusedShowers();
      dHist_UnusedShowerEnergy->Fill(locUnusedShowerEnergy, locAccWeight);
      if (locUnusedShowerEnergy > 1.)
	continue;
      
      //Loop over charged track hypotheses
      for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
	  //Set branch array indices corresponding to this particle
	  dChargedHypoWrapper->Set_ArrayIndex(loc_i);
		    
	  //Do stuff with the wrapper here ...
	  if (dChargedHypoWrapper->Get_PID() != PiMinus)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locPiPlusTrackID)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locProtonTrackID)
	    continue;
		  	   
	  TLorentzVector locPiMinusP4_reco = dChargedHypoWrapper->Get_P4_Measured();
	  TLorentzVector loc3PiP4_reco = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locPiMinusP4_reco;
	  double loc3PiMass_reco = loc3PiP4_reco.M();

	  //Mass Correlation with unused energy cuts
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation02;
	  locUsedThisCombo_MassCorrelation02[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation02[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation02[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation02[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation02[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation02[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation02.find(locUsedThisCombo_MassCorrelation02) == locUsedSoFar_MassCorrelation02.end())
	    {
	      dHist_3PiMass[1]->Fill(loc3PiMass_reco, locAccWeight);
	      dHist_OmegaMass_1track[1]->Fill(locOmegaMass_mmop, locAccWeight);
	      dHist_MassCorrelation[1]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight);
	      locUsedSoFar_MassCorrelation02.insert(locUsedThisCombo_MassCorrelation02);
	    }
	}

      
      //Cut on number of unused tracks
      dHist_NumTracks->Fill(itrack, locAccWeight);
      if(itrack > 3)
	continue;

      //Loop over charged track hypotheses
      for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
	  //Set branch array indices corresponding to this particle
	  dChargedHypoWrapper->Set_ArrayIndex(loc_i);
		    
	  //Do stuff with the wrapper here ...
	  if (dChargedHypoWrapper->Get_PID() != PiMinus)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locPiPlusTrackID)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locProtonTrackID)
	    continue;
		  	  
	  TLorentzVector locPiMinusP4_reco = dChargedHypoWrapper->Get_P4_Measured();
	  TLorentzVector loc3PiP4_reco = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locPiMinusP4_reco;
	  double loc3PiMass_reco = loc3PiP4_reco.M();

	  //Mass Correlation with NumTracks and unused energy cuts
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation03;
	  locUsedThisCombo_MassCorrelation03[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation03[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation03[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation03[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation03[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation03[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation03.find(locUsedThisCombo_MassCorrelation03) == locUsedSoFar_MassCorrelation03.end())
	    {
	      dHist_3PiMass[2]->Fill(loc3PiMass_reco, locAccWeight);
	      dHist_OmegaMass_1track[2]->Fill(locOmegaMass_mmop, locAccWeight);
	      dHist_MassCorrelation[2]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight);
	      locUsedSoFar_MassCorrelation03.insert(locUsedThisCombo_MassCorrelation03);
	    }
	}
     

      /************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/
      
      //Missing Mass
     double locMissingMass = locMissingP4_Measured.M();

      //Uniqueness tracking: Build the map of particles used for the missing mass
      //For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
      map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
      locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
      locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
      locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
      locUsedThisCombo_MissingMass[Gamma].insert(locPhoton1NeutralID);
      locUsedThisCombo_MissingMass[Gamma].insert(locPhoton2NeutralID);
      
      //compare to what's been used so far
      if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
	{
	  //unique missing mass combo: histogram it, and register this combo of particles
	  dHist_MissingMass->Fill(locMissingMass, locAccWeight);
	  locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
	}

      //Cut on missing mass
      if(locMissingMass > 0.25)
	continue;

      //Loop over charged track hypotheses
      for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
	  //Set branch array indices corresponding to this particle
	  dChargedHypoWrapper->Set_ArrayIndex(loc_i);
	  
	  //Do stuff with the wrapper here ...
	  if (dChargedHypoWrapper->Get_PID() != PiMinus)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locPiPlusTrackID)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locProtonTrackID)
	    continue;

	  TLorentzVector locPiMinusP4_reco = dChargedHypoWrapper->Get_P4_Measured();
	  TLorentzVector loc3PiP4_reco = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locPiMinusP4_reco;
	  double loc3PiMass_reco = loc3PiP4_reco.M();

	  //Mass Correlation with NumTracks, unused energy, and missing mass cuts
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation04;
	  locUsedThisCombo_MassCorrelation04[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation04[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation04[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation04[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation04[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation04[Gamma].insert(locPhoton2NeutralID);
	  
	  if(locUsedSoFar_MassCorrelation04.find(locUsedThisCombo_MassCorrelation04) == locUsedSoFar_MassCorrelation04.end())
	    {
	      dHist_3PiMass[3]->Fill(loc3PiMass_reco, locAccWeight);
	      dHist_OmegaMass_1track[3]->Fill(locOmegaMass_mmop, locAccWeight);
	      dHist_MassCorrelation[3]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight);
	      locUsedSoFar_MassCorrelation04.insert(locUsedThisCombo_MassCorrelation04);
	    }
	}
      
      /*************************************** HISTOGRAM KINFIT CL *************************************/
      
      dHist_KinFitCL->Fill(locKinFit_CL, locAccWeight);
      //Cut on KinFit CL
      if(locKinFit_CL < 0.1)
	continue;
      
      //Loop over charged track hypotheses
      for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
	  //Set branch array indices corresponding to this particle
	  dChargedHypoWrapper->Set_ArrayIndex(loc_i);
	  
	  //Do stuff with the wrapper here ...
	  if (dChargedHypoWrapper->Get_PID() != PiMinus)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locPiPlusTrackID)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locProtonTrackID)
	    continue;

	  TLorentzVector locPiMinusP4_reco = dChargedHypoWrapper->Get_P4_Measured();
	  TLorentzVector loc3PiP4_reco = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locPiMinusP4_reco;
	  double loc3PiMass_reco = loc3PiP4_reco.M();
	  
	  //Mass Correlation with NumTracks, unused energy, missing mass, and KinFit CL cuts
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation05;
	  locUsedThisCombo_MassCorrelation05[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation05[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation05[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation05[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation05[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation05[Gamma].insert(locPhoton2NeutralID);
	  
	  if(locUsedSoFar_MassCorrelation05.find(locUsedThisCombo_MassCorrelation05) == locUsedSoFar_MassCorrelation05.end())
	    {
	      dHist_3PiMass[4]->Fill(loc3PiMass_reco, locAccWeight);
	      dHist_OmegaMass_1track[4]->Fill(locOmegaMass_mmop, locAccWeight);
	      dHist_MassCorrelation[4]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight);
	      locUsedSoFar_MassCorrelation05.insert(locUsedThisCombo_MassCorrelation05);
	    }
	}


      /******************************************* HISTOGRAM 2 PHOTON MASS *********************************************************/
      
      map<Particle_t, set<Int_t> > locUsedThisCombo_2GammaMass;
      locUsedThisCombo_2GammaMass[Gamma].insert(locPhoton1NeutralID);
      locUsedThisCombo_2GammaMass[Gamma].insert(locPhoton2NeutralID);
      
      if(locUsedSoFar_2GammaMass.find(locUsedThisCombo_2GammaMass) == locUsedSoFar_2GammaMass.end())
	{
	  dHist_2GammaMass->Fill(loc2GammaMass, locAccWeight);
	  locUsedSoFar_2GammaMass.insert(locUsedThisCombo_2GammaMass);
	}

      //Sideband subtraction on pi0 mass
      double locPi0Weight;
      if(loc2GammaMass > 0.12 && loc2GammaMass < 0.15)
	locPi0Weight = 1.;
      else if((loc2GammaMass > 0.095 && loc2GammaMass < 0.11) || (loc2GammaMass > 0.16 && loc2GammaMass < 0.175))
	locPi0Weight = -1.;
      else
	locPi0Weight = 0.;
      
      TLorentzVector locTruthP4;
      //Loop over throwns to get true values
      for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
	  if(Get_NumThrown() == 0)
	    continue;
	  //Set branch array indices corresponding to this particle
	  dThrownWrapper->Set_ArrayIndex(loc_i);
	  
	  //Do stuff with the wrapper here ...
	  
	  if (dThrownWrapper->Get_PID() != PiMinus)
	    continue;
	  
	  locTruthP4 = dThrownWrapper->Get_P4_Measured();
	}
      
      double locTruthPhi = locTruthP4.Phi() * 180./TMath::Pi();
      double locTruthTheta = locTruthP4.Theta() * 180./TMath::Pi();
      TVector3 locTruthP3 = locTruthP4.Vect();
      double locTruthP = locTruthP3.Mag();
      
      //Loop over charged track hypotheses
      for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
	  //Set branch array indices corresponding to this particle
	  dChargedHypoWrapper->Set_ArrayIndex(loc_i);
	  
	  //Do stuff with the wrapper here ...
	  if (dChargedHypoWrapper->Get_PID() != PiMinus)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locPiPlusTrackID)
	    continue;
	  if (dChargedHypoWrapper->Get_TrackID() == locProtonTrackID)
	    continue;

	  TLorentzVector locPiMinusP4_reco = dChargedHypoWrapper->Get_P4_Measured();
	  double locPiMinusPhi_reco = locPiMinusP4_reco.Phi() * 180./TMath::Pi();
	  double locPiMinusTheta_reco = locPiMinusP4_reco.Theta() * 180./TMath::Pi();
	  TVector3 locPiMinusP3_reco = locPiMinusP4_reco.Vect();
	  double locPiMinusP_reco = locPiMinusP3_reco.Mag();

	  TLorentzVector loc3PiP4_reco = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locPiMinusP4_reco;
	  double loc3PiMass_reco = loc3PiP4_reco.M();

	  //Mass Correlation with KinFit CL, NumTracks, unused energy, missing mass, and pi0 mass cuts
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation06;
	  locUsedThisCombo_MassCorrelation06[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation06[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation06[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation06[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation06[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation06[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation06.find(locUsedThisCombo_MassCorrelation06) == locUsedSoFar_MassCorrelation06.end())
	    {
	      dHist_3PiMass[5]->Fill(loc3PiMass_reco, locAccWeight * locPi0Weight);
	      dHist_OmegaMass_1track[5]->Fill(locOmegaMass_mmop, locAccWeight * locPi0Weight);
	      dHist_MassCorrelation[5]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight * locPi0Weight);
	      locUsedSoFar_MassCorrelation06.insert(locUsedThisCombo_MassCorrelation06);
	    }
	  
	  double locDeltaPhiTruth;
	  if(locPiMinusPhi_reco - locTruthPhi > 180.)
	    locDeltaPhiTruth = locPiMinusPhi_reco - locTruthPhi - 2.0*180.;
	  else if(locPiMinusPhi_reco - locTruthPhi < -180.)
	    locDeltaPhiTruth = locPiMinusPhi_reco - locTruthPhi + 2.0*180.;
	  else 
	    locDeltaPhiTruth = locPiMinusPhi_reco - locTruthPhi;
	  double locDeltaThetaTruth = locPiMinusTheta_reco - locTruthTheta;
	  double locDeltaPTruth = locPiMinusP_reco - locTruthP;
	  
	  double locDeltaPhiMissing;
	  if(locPiMinusPhi_reco - locPiMinusPhi_mmop > 180.)
	    locDeltaPhiMissing = locPiMinusPhi_reco - locPiMinusPhi_mmop - 2.0*180.;
	  else if(locPiMinusPhi_reco - locPiMinusPhi_mmop < -180.)
	    locDeltaPhiMissing = locPiMinusPhi_reco - locPiMinusPhi_mmop + 2.0*180.;
	  else
	    locDeltaPhiMissing = locPiMinusPhi_reco - locPiMinusPhi_mmop;
	  double locDeltaThetaMissing = locPiMinusTheta_reco - locPiMinusTheta_mmop;
	  double locDeltaPMissing = locPiMinusP_reco - locPiMinusP_mmop;

	  //Histogram DeltaPhi, DeltaTheta, DeltaP (Truth)
	  map<Particle_t, set<Int_t> > locUsedThisCombo_DeltaTruth;
	  locUsedThisCombo_DeltaTruth[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_DeltaTruth[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_DeltaTruth[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_DeltaTruth[Gamma].insert(locPhoton2NeutralID);
	  locUsedThisCombo_DeltaTruth[Unknown].insert(locBeamID);
		    
	  if(locUsedSoFar_DeltaTruth.find(locUsedThisCombo_DeltaTruth) == locUsedSoFar_DeltaTruth.end())
	    {
	      dHist_DeltaPhiTruth->Fill(locDeltaPhiTruth, locAccWeight * locPi0Weight);
	      dHist_DeltaThetaTruth->Fill(locDeltaThetaTruth, locAccWeight * locPi0Weight);
	      dHist_DeltaPTruth->Fill(locDeltaPTruth, locAccWeight * locPi0Weight);
	      locUsedSoFar_DeltaTruth.insert(locUsedThisCombo_DeltaTruth);
	    }
	  
	  //Histogram DeltaPhi, DeltaTheta, DeltaP (Missing)
	  map<Particle_t, set<Int_t> > locUsedThisCombo_DeltaMissing;
	  locUsedThisCombo_DeltaMissing[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_DeltaMissing[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_DeltaMissing[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_DeltaMissing[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_DeltaMissing[Gamma].insert(locPhoton2NeutralID);
	  locUsedThisCombo_DeltaMissing[Unknown].insert(locBeamID);
	  
	  if(locUsedSoFar_DeltaMissing.find(locUsedThisCombo_DeltaMissing) == locUsedSoFar_DeltaMissing.end())
	    {
	      dHist_DeltaPhiMissing->Fill(locDeltaPhiMissing, locAccWeight * locPi0Weight);
	      dHist_DeltaThetaMissing->Fill(locDeltaThetaMissing, locAccWeight * locPi0Weight);
	      dHist_DeltaPMissing->Fill(locDeltaPMissing, locAccWeight * locPi0Weight);
	      dHist_DeltaPhiVsDeltaTheta_missing->Fill(locDeltaPhiMissing, locDeltaThetaMissing, locAccWeight * locPi0Weight);
	      dHist_DeltaPVsDeltaPhi_missing->Fill(locDeltaPMissing, locDeltaPhiMissing, locAccWeight * locPi0Weight);
	      dHist_DeltaPVsDeltaTheta_missing->Fill(locDeltaPMissing, locDeltaThetaMissing, locAccWeight * locPi0Weight );
	      locUsedSoFar_DeltaMissing.insert(locUsedThisCombo_DeltaMissing);
	    }

	  //Create Weights for Delta(reco - truth)
	  double locTruthWeight;
	  if(Get_NumThrown() == 0)
	    locTruthWeight = 1.;
	  else if(locDeltaPhiTruth > -90. && locDeltaPhiTruth < 90. && locDeltaThetaTruth > -30. && locDeltaThetaTruth < 30. && locDeltaPTruth > -3. && locDeltaPTruth < 3.)
	    locTruthWeight = 1.;
	  else
	    locTruthWeight = 0.;

	  //Mass correlation with Delta_truth cuts
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation07;
	  locUsedThisCombo_MassCorrelation07[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation07[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation07[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation07[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation07[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation07[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation07.find(locUsedThisCombo_MassCorrelation07) == locUsedSoFar_MassCorrelation07.end())
	    {
	      dHist_3PiMass[6]->Fill(loc3PiMass_reco, locAccWeight * locPi0Weight * locTruthWeight);
	      dHist_OmegaMass_1track[6]->Fill(locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight);
	      dHist_MassCorrelation[6]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight);
	      locUsedSoFar_MassCorrelation07.insert(locUsedThisCombo_MassCorrelation07);
	    }


	  //Create Weights for Delta(reco - missing)
	  double locCutWeight;
	  if (locDeltaPhiMissing > -90. && locDeltaPhiMissing < 90. && locDeltaThetaMissing > -30. && locDeltaThetaMissing < 30. && locDeltaPMissing > -3. && locDeltaPMissing < 3.)
	    locCutWeight = 1.0;
	  else
	    locCutWeight = 0.0;

	  //Mass Correlation with all cuts (except missing energy and MM^2)
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation08;
	  locUsedThisCombo_MassCorrelation08[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation08[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation08[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation08[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation08[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation08[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation08.find(locUsedThisCombo_MassCorrelation08) == locUsedSoFar_MassCorrelation08.end())
	    {
	      dHist_3PiMass[7]->Fill(loc3PiMass_reco, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_OmegaMass_1track[7]->Fill(locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_MassCorrelation[7]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      locUsedSoFar_MassCorrelation08.insert(locUsedThisCombo_MassCorrelation08);
	    }

	  //Missing Mass, Energy etc.
	  TLorentzVector locMissingP4_reco = locBeamP4_Measured + dTargetP4;
	  locMissingP4_reco -= locPiPlusP4_Measured + locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured + locPiMinusP4_reco;

	  double locMissingMassSquared = locMissingP4_reco.M2();
 
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMassSquared;
	  locUsedThisCombo_MissingMassSquared[Unknown].insert(locBeamID); //beam
	  locUsedThisCombo_MissingMassSquared[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MissingMassSquared[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MissingMassSquared[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MissingMassSquared[Gamma].insert(locPhoton2NeutralID);
      
	  //compare to what's been used so far
	  if(locUsedSoFar_MissingMassSquared.find(locUsedThisCombo_MissingMassSquared) == locUsedSoFar_MissingMassSquared.end())
	    {
	      //unique missing mass combo: histogram it, and register this combo of particles
	      dHist_MissingMassSquared->Fill(locMissingMassSquared, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      locUsedSoFar_MissingMassSquared.insert(locUsedThisCombo_MissingMassSquared);
	    }

	  //Cut on missing mass squared
	  if(locMissingMassSquared > 0.1 || locMissingMassSquared < -0.1)
	    continue;


	  //Mass Correlation with all cuts (except missing energy)
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation09;
	  locUsedThisCombo_MassCorrelation09[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation09[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation09[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation09[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation09[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation09[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation09.find(locUsedThisCombo_MassCorrelation09) == locUsedSoFar_MassCorrelation09.end())
	    {
	      dHist_3PiMass[8]->Fill(loc3PiMass_reco, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_OmegaMass_1track[8]->Fill(locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_MassCorrelation[8]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      locUsedSoFar_MassCorrelation09.insert(locUsedThisCombo_MassCorrelation09);
	    }


	  double locMissingEnergy = locMissingP4_reco.E();
 
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MissingEnergy;
	  locUsedThisCombo_MissingEnergy[Unknown].insert(locBeamID); //beam
	  locUsedThisCombo_MissingEnergy[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MissingEnergy[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MissingEnergy[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MissingEnergy[Gamma].insert(locPhoton2NeutralID);
      
	  //compare to what's been used so far
	  if(locUsedSoFar_MissingEnergy.find(locUsedThisCombo_MissingEnergy) == locUsedSoFar_MissingEnergy.end())
	    {
	      //unique missing energy combo: histogram it, and register this combo of particles
	      dHist_MissingEnergy->Fill(locMissingEnergy, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      locUsedSoFar_MissingEnergy.insert(locUsedThisCombo_MissingEnergy);
	    }

	  //Cut on missing energy
	  if(locMissingEnergy > 3. || locMissingEnergy < -3.)
	    continue;


	  //Mass Correlation with all cuts (for real)
	  map<Particle_t, set<Int_t> > locUsedThisCombo_MassCorrelation10;
	  locUsedThisCombo_MassCorrelation10[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_MassCorrelation10[Unknown].insert(locBeamID);
	  locUsedThisCombo_MassCorrelation10[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_MassCorrelation10[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_MassCorrelation10[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_MassCorrelation10[Gamma].insert(locPhoton2NeutralID);

	  if(locUsedSoFar_MassCorrelation10.find(locUsedThisCombo_MassCorrelation10) == locUsedSoFar_MassCorrelation10.end())
	    {
	      dHist_3PiMass[9]->Fill(loc3PiMass_reco, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_OmegaMass_1track[9]->Fill(locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_MassCorrelation[9]->Fill(loc3PiMass_reco, locOmegaMass_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_ThetaVsP_reco->Fill(locPiMinusTheta_reco, locPiMinusP_reco, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      locUsedSoFar_MassCorrelation10.insert(locUsedThisCombo_MassCorrelation10);
	    }


	  //Histogram Omega Mass (mmop) and 3Pi Mass vs Phi, Theta, P (mmop)
	  map<Particle_t, set<Int_t> > locUsedThisCombo_OmegaMassVsPhi_1track;
	  locUsedThisCombo_OmegaMassVsPhi_1track[PiMinus].insert(dChargedHypoWrapper->Get_TrackID());
	  locUsedThisCombo_OmegaMassVsPhi_1track[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_OmegaMassVsPhi_1track[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_OmegaMassVsPhi_1track[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_OmegaMassVsPhi_1track[Gamma].insert(locPhoton2NeutralID);
	  locUsedThisCombo_OmegaMassVsPhi_1track[Unknown].insert(locBeamID);

	  double pmin[9] = {0.5, 0.75, 1., 1.25, 1.5, 2., 3., 4., 5.};
	  double pmax[9] = {0.75, 1., 1.25, 1.5, 2., 3., 4., 5., 6.};

	  if(locUsedSoFar_OmegaMassVsPhi_1track.find(locUsedThisCombo_OmegaMassVsPhi_1track) == locUsedSoFar_OmegaMassVsPhi_1track.end())
	    {
	      dHist_OmegaMassVsPhi_1track->Fill(locOmegaMass_mmop, locPiMinusPhi_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      dHist_3PiMassVsPhi->Fill(loc3PiMass_reco, locPiMinusPhi_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      if(locPiMinusP_mmop > 0.5){
		dHist_OmegaMassVsTheta_1track->Fill(locOmegaMass_mmop, locPiMinusTheta_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
		dHist_3PiMassVsTheta->Fill(loc3PiMass_reco, locPiMinusTheta_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      }
	      if(locPiMinusTheta_mmop > 1.15){
		dHist_OmegaMassVsP_1track->Fill(locOmegaMass_mmop, locPiMinusP_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
		dHist_3PiMassVsP->Fill(loc3PiMass_reco, locPiMinusP_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
	      }
	      for(int i = 0; i < 9; i++) {
		if(locPiMinusP_mmop > pmin[i] && locPiMinusP_mmop < pmax[i]){
		  dHist_OmegaMassVsTheta_1track_p[i]->Fill(locOmegaMass_mmop, locPiMinusTheta_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
		  dHist_3PiMassVsTheta_p[i]->Fill(loc3PiMass_reco, locPiMinusTheta_mmop, locAccWeight * locPi0Weight * locTruthWeight * locCutWeight);
		}
	      }
	      locUsedSoFar_OmegaMassVsPhi_1track.insert(locUsedThisCombo_OmegaMassVsPhi_1track);
	    }
      
	  /*********************************** PLOT DELTAPHI, THETA, P VS THETA VS P -TRUTH **************************************************/
      
	  map<Particle_t, set<Int_t> > locUsedThisCombo_Delta3D_truth;
	  locUsedThisCombo_Delta3D_truth[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_Delta3D_truth[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_Delta3D_truth[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_Delta3D_truth[Gamma].insert(locPhoton2NeutralID);
	  locUsedThisCombo_Delta3D_truth[Unknown].insert(locBeamID);
		    
	  if(locUsedSoFar_Delta3D_truth.find(locUsedThisCombo_Delta3D_truth) == locUsedSoFar_Delta3D_truth.end())
	    {
	      dHist_DeltaPhi3D_truth->Fill(locDeltaPhiTruth, locTruthTheta, locTruthP, locAccWeight * locPi0Weight);
	      dHist_DeltaTheta3D_truth->Fill(locDeltaThetaTruth, locTruthTheta, locTruthP, locAccWeight * locPi0Weight);
	      dHist_DeltaP3D_truth->Fill(locDeltaPTruth, locTruthTheta, locTruthP, locAccWeight * locPi0Weight);
	      locUsedSoFar_Delta3D_truth.insert(locUsedThisCombo_Delta3D_truth);
	    }
      

	  /*********************************** MISSING - PLOT DELTAPHI, THETA, P VS THETA VS P **************************************************/
	  
	  map<Particle_t, set<Int_t> > locUsedThisCombo_Delta3D_miss;
	  locUsedThisCombo_Delta3D_miss[PiPlus].insert(locPiPlusTrackID);
	  locUsedThisCombo_Delta3D_miss[Proton].insert(locProtonTrackID);
	  locUsedThisCombo_Delta3D_miss[Gamma].insert(locPhoton1NeutralID);
	  locUsedThisCombo_Delta3D_miss[Gamma].insert(locPhoton2NeutralID);
	  locUsedThisCombo_Delta3D_miss[Unknown].insert(locBeamID);
		    
	  if(locUsedSoFar_Delta3D_miss.find(locUsedThisCombo_Delta3D_miss) == locUsedSoFar_Delta3D_miss.end())
	    {
	      dHist_DeltaPhi3D_miss->Fill(locDeltaPhiMissing, locPiMinusTheta_mmop, locPiMinusP_mmop, locAccWeight * locPi0Weight);
	      dHist_DeltaTheta3D_miss->Fill(locDeltaThetaMissing, locPiMinusTheta_mmop, locPiMinusP_mmop, locAccWeight * locPi0Weight);
	      dHist_DeltaP3D_miss->Fill(locDeltaPMissing, locPiMinusTheta_mmop, locPiMinusP_mmop, locAccWeight * locPi0Weight);
	      locUsedSoFar_Delta3D_miss.insert(locUsedThisCombo_Delta3D_miss);
	    }  
	}



      /******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/
      
      // Loop through the analysis actions, executing them in order for the active particle combo
  
      //dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
      if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
	continue;
 
      //if you manually execute any actions, and it fails a cut, be sure to call:
      //dComboWrapper->Set_IsComboCut(true);

      /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

      /*
	TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
	//for arrays below: 2nd argument is value, 3rd is array index
	//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
	//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
	dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
      */

      /**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/
      
      //Histogram beam energy (if haven't already)
      if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
	{
	  dHist_BeamEnergy->Fill(locBeamP4.E());
	  locUsedSoFar_BeamEnergy.insert(locBeamID);
	}
      

      /******************************************* HISTOGRAM OMEGA USING MISSING MASS OFF THE PROTON ******************************/
      
      map<Particle_t, set<Int_t> > locUsedThisCombo_OmegaMass;
      locUsedThisCombo_OmegaMass[Unknown].insert(locBeamID);
      locUsedThisCombo_OmegaMass[Proton].insert(locProtonTrackID);
      
      if(locUsedSoFar_OmegaMass.find(locUsedThisCombo_OmegaMass) == locUsedSoFar_OmegaMass.end())
	{
	  dHist_OmegaMass_missing->Fill(locOmegaMass_mmop, locAccWeight * locPi0Weight);
	  locUsedSoFar_OmegaMass.insert(locUsedThisCombo_OmegaMass);
	}

      /************************************** HISTOGRAM OMEGA MASS VS PHI, THETA, 3-MOMENTUM *******************/
      
      map<Particle_t, set<Int_t> > locUsedThisCombo_OmegaMassVsPhi_0track;
      locUsedThisCombo_OmegaMassVsPhi_0track[PiPlus].insert(locPiPlusTrackID);
      locUsedThisCombo_OmegaMassVsPhi_0track[Proton].insert(locProtonTrackID);
      locUsedThisCombo_OmegaMassVsPhi_0track[Gamma].insert(locPhoton1NeutralID);
      locUsedThisCombo_OmegaMassVsPhi_0track[Gamma].insert(locPhoton2NeutralID);
      locUsedThisCombo_OmegaMassVsPhi_0track[Unknown].insert(locBeamID);

      double pmin[9] = {0.5, 0.75, 1., 1.25, 1.5, 2., 3., 4., 5.};
      double pmax[9] = {0.75, 1., 1.25, 1.5, 2., 3., 4., 5., 6.};
      
      if(locUsedSoFar_OmegaMassVsPhi_0track.find(locUsedThisCombo_OmegaMassVsPhi_0track) == locUsedSoFar_OmegaMassVsPhi_0track.end())
	{
	  //loop over charged track hypothesis
	  bool found = false;
	  for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	    {
	      //Set branch array indices corresponding to this particle
	      dChargedHypoWrapper->Set_ArrayIndex(loc_i);
	      
	      //Do stuff with the wrapper here ...
	      if (dChargedHypoWrapper->Get_PID() != PiMinus)
		continue;
	      if (dChargedHypoWrapper->Get_TrackID() == locPiPlusTrackID)
		continue;
	      if (dChargedHypoWrapper->Get_TrackID() == locProtonTrackID)
		continue;
	      found = true;	
	    }
	  if(found == true)
	    continue;

	  dHist_OmegaMassVsPhi_0track->Fill(locOmegaMass_mmop, locPiMinusPhi_mmop, locAccWeight * locPi0Weight);
	  if(locPiMinusP_mmop > 0.5) {
	    dHist_OmegaMassVsTheta_0track->Fill(locOmegaMass_mmop, locPiMinusTheta_mmop, locAccWeight * locPi0Weight);
	  }
	  if(locPiMinusTheta_mmop > 1.15) {
	    dHist_OmegaMassVsP_0track->Fill(locOmegaMass_mmop, locPiMinusP_mmop, locAccWeight * locPi0Weight);
	  }
	  dHist_ThetaVsP_missing->Fill(locPiMinusTheta_mmop, locPiMinusP_mmop, locAccWeight * locPi0Weight);
	  for(int i = 0; i < 9; i++) {
	    if(locPiMinusP_mmop > pmin[i] && locPiMinusP_mmop < pmax[i]){
	      dHist_OmegaMassVsTheta_0track_p[i]->Fill(locOmegaMass_mmop, locPiMinusTheta_mmop, locAccWeight * locPi0Weight);
	    }
	  }
	  locUsedSoFar_OmegaMassVsPhi_0track.insert(locUsedThisCombo_OmegaMassVsPhi_0track);
	}	


      /****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

      /*
      //FILL ANY CUSTOM BRANCHES FIRST!!
      Int_t locMyInt_Flat = 7;
      dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);
      
      TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
      dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);
      
      for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
      {
      dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
      TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
      dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
      }
      */
      
      //FILL FLAT TREE
      //Fill_FlatTree(); //for the active combo
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
		Fill_OutputTree();
*/

  return kTRUE;
}

void DSelector_omega_misspim::Finalize(void)
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
