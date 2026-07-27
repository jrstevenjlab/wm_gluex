#include "DSelector_pi0kpkm.h"

void DSelector_pi0kpkm::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "pi0kpkm.root"; //"" for none
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
	// dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, K, 1000, 0.9, 1.9, "phiMeson"));
	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	// dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));

	deque<Particle_t> locPi0PIDs;  locPi0PIDs.push_back(Gamma); locPi0PIDs.push_back(Gamma);
	double minPi0 = 0.09; double maxPi0 = 0.18;

	deque<Particle_t> locPhiPIDs;  locPhiPIDs.push_back(KPlus); locPhiPIDs.push_back(KMinus);
	deque<Particle_t> locKPlusPIDs; locKPlusPIDs.push_back(KPlus);
	deque<Particle_t> locKMinusPIDs; locKMinusPIDs.push_back(KMinus);
	deque<Particle_t> locProtonPIDs; locProtonPIDs.push_back(Proton);

//	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Pi0, 100, minPi0, maxPi0, "Pi0_NoCut"));

	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 1, locPi0PIDs, 100, minPi0, maxPi0, "Pi0_NoCut"));

//    	dAnalysisActions.push_back(new DHistogramAction_vanHoveFour(dComboWrapper, false, locKPlusPIDs, locKMinusPIDs,  locPi0PIDs, locProtonPIDs, "Four_particle_vanHove_no_cut"));

	/*--------------------------------------------------*/
	/// Kinematics fit on 1e-5 
//  dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 1e-5));
//	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.01));
	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.001));
//	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.1));
//    dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.5));


	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));
    	dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.05, 0.05));
//	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.01));
//    	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 1e-5));
	
//	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Pi0, 100, minPi0, maxPi0, "Pi0_KinCut"));

	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 1, locPi0PIDs, 100, minPi0, maxPi0, "Pi0_KinCut"));

	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locPhiPIDs, 100, 0.9, 1.4, "Phi_KinCut"));


	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));
//
	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));


	//*************************************
	// Date: Nov 27 2019
	// WL learning comments
	//
	// True: kinematics variable, false: measured variable
	// Number is the level of the decay chain
	//

	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locPhiPIDs, 100, 0.9, 1.4, "Phi_NoCut"));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

//    	dAnalysisActions.push_back(new DHistogramAction_vanHove(dComboWrapper, false, locKPlusPIDs, locKMinusPIDs, locPi0PIDs, ""));

//    	dAnalysisActions.push_back(new DHistogramAction_vanHove(dComboWrapper, false, KPlus, KMinus, Proton, ""));

//    	dAnalysisActions.push_back(new DHistogramAction_vanHove(dComboWrapper, false, KPlus, KMinus,  Pi0, ""));

//    	dAnalysisActions.push_back(new DHistogramAction_vanHoveFour(dComboWrapper, false, locKPlusPIDs, locKMinusPIDs,  locPi0PIDs, locProtonPIDs, "Four_particle_vanHove"));



	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

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

	int n2k = 500; double min2k = 0.9; double max2k = 1.9;
        int nkpi = 750; double minkpi = 0.5; double maxkpi = 2.0;
        int nt = 100; double mint = 0; double maxt = 5.0;
        int nDelta = 100; double minDelta = 1.0; double maxDelta = 3.0;
        int n2kpi = 400; double min2kpi = 1.0; double max2kpi = 4.0;

        int n2kpi_dalitz = 40; 
	int n2k_dalitz = 50;
	int nt_dalitz = 10;


	///*--------------------------------------------------*/
	/// Flag for pi0pippiminus

	if (dKPlusWrapper->Get_P4().M() < 0.4) {	
		min2k = 0.3;
		min2kpi = 0.5;
	}


	// pi0_mass 
	hpi0_mass = new TH1I("hpi0_mass", "pi0_mass", 100, 0, 0.3);
	hpi0_mass_measured = new TH1I("hpi0_mass_measured", "pi0_mass_measured", 100, 0, 0.3);

	h2gamma_2kpi0_mass = new TH2I("2gamma_2kpi0",  "2gamma_2kpi0", 100, 0, 0.3, n2kpi, min2kpi, max2kpi);

	h2kpi0_invariantmass_below_phi = new TH1I("h2kpi0_invariantmass_below_phi", "2kpi0_invariantmass_blow_phi", n2kpi, min2kpi, max2kpi);

	///--------------------------------------
	// Invariance mass plot

	gDirectory->mkdir("Invariantmass");
	gDirectory->cd("Invariantmass");

	h2k_invariantmass = new TH1I("h2k_invariantmass", "2k_invariantmass", n2k, min2k, max2k);

	h2kpi0_invariantmass = new TH1I("h2kpi0_invariantmass", "2kpi0_invariantmass", n2kpi, min2kpi, max2kpi);

	///--------------------------------------
	//
	// Hist a: M(2K) <= 1.010
	// Hist b: 1.010 < M(2K) <= 1.025
	// Hist c: 1.025 < M(2K) <= 1.040

	h2kpi0_invariantmass_a= new TH1I("h2kpi0_invariantmass_below_phi", "2kpi0_invariantmass_below_phi", n2kpi, min2kpi, max2kpi);

	h2kpi0_invariantmass_b = new TH1I("h2kpi0_invariantmass_phi", "2kpi0_invariantmass_phi", n2kpi, min2kpi, max2kpi);

	h2kpi0_invariantmass_c = new TH1I("h2kpi0_invariantmass_above_phi", "2kpi0_invariantmass_above_phi", n2kpi, min2kpi, max2kpi);

	hkppi0_invariantmass = new TH1I("hkppi0_invariantmass", "kppi0_invariantmass", nkpi, minkpi, maxkpi);
	hkmpi0_invariantmass = new TH1I("hkmpi0_invariantmass", "kmpi0_invariantmass", nkpi, minkpi, maxkpi);

	h2kpi0_invariantmass_below_phi_pi0_cut = new TH1I("h2kpi0_invariantmass_below_phi_pi0_cut", "2kpi0_invariantmass_blow_phi_pi0_cut", n2kpi, min2kpi, max2kpi);

	h2kpi0_invariantmass_below_phi_pi0_sideband = new TH1I("h2kpi0_invariantmass_below_phi_pi0_sideband", "2kpi0_invariantmass_blow_phi_pi0_sideband", n2kpi, min2kpi, max2kpi);

	h2kpi0_invariantmass_below_phi_pi0_lo_sideband = new TH1I("h2kpi0_invariantmass_below_phi_pi0_lo_sideband", "2kpi0_invariantmass_blow_phi_pi0_sideband", n2kpi, min2kpi, max2kpi);

	h2kpi0_invariantmass_below_phi_pi0_hi_sideband = new TH1I("h2kpi0_invariantmass_below_phi_pi0_hi_sideband", "2kpi0_invariantmass_blow_phi_pi0_sideband", n2kpi, min2kpi, max2kpi);


	gDirectory->cd("/");

	//-----------------------------------
	// t disribution

	gDirectory->mkdir("t_distribution");
	gDirectory->cd("t_distribution");
	
	dHist_2kpi_t = new TH1I("dHist_2kpi_t", "dHist_2kpi_t", 200, 0, 20);
	dHist_2k_t   = new TH1I("dHist_2k_t", "dHist_2k_t", 200, 0, 20);
	dHist_2kpiMass_t = new TH2I("dHist_2kpiMass_t", "dHist_2kpiMass_t", 200, 0, 4, 200, 1, 2);
	dHist_t2k_2kpiMass = new TH2I("dHist_t2k_2kpiMass", "dHist_t2k_2kpiMass", 200, 0, 4, 200, 1, 2);
	dHist_t_t2k = new TH2I("dHist_t_t2k", "dHist_t_2k", 300, 0, 3, 300, 0,3);

	gDirectory->cd("/");

	//-----------------------------------
	// Dalitz Plot

	gDirectory->mkdir("Dalitz");
	gDirectory->cd("Dalitz");

//	dalitz_pikp_pikm_low  = new TH2I("hpikp_pikm_low",  "pikp_pikm_low", 150, 0, 0.9, 150, 0, 0.9);
//	dalitz_pikp_pikm_mid  = new TH2I("hpikp_pikm_mid",  "pikp_pikm_mid", 150, 0, 1.4, 150, 0, 1.4);
//	dalitz_pikp_pikm_high = new TH2I("hpikp_pikm_high", "pikp_pikm_high", 300, 0, 2.0, 200, 0, 2.0);
//	dalitz_pikp_pikm_ultra = new TH2I("hpikp_pikm_ultra", "pikp_pikm_ultra", 300, 0, 4.0, 200, 0, 4.0);
	


//	dalitz_kpkm_pi = new TH2I("dalitz_2k_pi",  "kpkm_kpkmpi0", 200, 0, 4, 200, 0, 4);
//	dalitz_kp_pi = new TH2I("dalitz_kp_pi",  "kpkm_kpkmpi0", 200, 0, 4, 200, 0, 4);
//	dalitz_km_pi = new TH2I("dalitz_km_pi",  "kpkm_kpkmpi0", 200, 0, 4, 200, 0, 4);


	dalitz_kppi0_kmpi0 = new TH2I("dalitz_kpi0_kmpi0",  "kppi0_kmpi0", 20, 0, 4, 20, 0, 4);
	dalitz_kppi0_kmpi0_sq = new TH2I("dalitz_kpi0_kmpi0_sq",  "kppi0_kmpi0_sq", 20, 0, 4, 20, 0, 4);

	dalitz_kpkmpi0_kpkm = new TH2I("dalitz_kpkmpi0_kpkm",  "kpkmpi0_kpkm", 20, 0, 4, 20, 0, 4);
	dalitz_kpkmpi0_kpkm_sq  = new TH2I("dalitz_kpkmpi0_kpkm_sq",  "kpkmpi0_kpkm_sq", 20, 0, 4, 20, 0, 4);

	dalitz_kpkm_kppi0 = new TH2I("dalitz_kpkm_kppi0",  "kpkm_kppi0", 20, 0, 4, 20, 0, 4);
	dalitz_kpkm_kppi0_sq = new TH2I("dalitz_kpkm_kppi0_sq",  "kpkm_kppi0_sq", 20, 0, 4, 20, 0, 4);

	dalitz_kpkm_kmpi0 = new TH2I("dalitz_kpkm_kmpi0",  "kpkm_kmpi0", 20, 0, 4, 20, 0, 4);
	dalitz_kpkm_kmpi0_sq = new TH2I("dalitz_kpkm_kmpi0_sq",  "kpkm_kppi0_sq", 20, 0, 4, 20, 0, 4);

	dalitz_kppi0_kpkmpi0 = new TH2I("dalitz_kppi0_kpkmpi0",  "kppi0_kpkmmpi0", 20, 0, 4, 20, 0, 4);
	dalitz_kppi0_kpkmpi0_sq = new TH2I("dalitz_kppi0_kpkmpi0_sq",  "kppi0_kpkmmpi0_sq", 20, 0, 4, 20, 0, 4);

	dalitz_kmpi0_kpkmpi0 = new TH2I("dalitz_kmpi0_kpkmpi0",  "kppi0_kmkmmpi0", 20, 0, 4, 20, 0, 4);
	dalitz_kmpi0_kpkmpi0_sq = new TH2I("dalitz_kmpi0_kpkmpi0_sq",  "kmpi0_kpkmmpi0_sq", 20, 0, 4, 20, 0, 4);

	dalitz_t2kpi_t2k_kpkmpi0 = new TH3I("dalitz_t2kpi_t2k_kpkmpi0",  "t2kpi_t2k_kpkmpi0", 20, 0, 20, 20, 0, 20, 20, 0, 4);

	dalitz_t2kpi_t2k_kpkm = new TH3I("dalitz_t2kpi_t2k_kpkm",  "t2kpi_t2k_kpkm", 200, 0, 20, 20, 0, 20, 20, 0, 4);

	dalitz_kpkmpi0_kpkm_t2k = new TH3I("dalitz_kpkmpi0_kpkm_t2k",  "dalitz_kpkmpi0_kpkm_t2k", n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2kpi = new TH3I("dalitz_kpkmpi0_kpkm_t2kpi", "dalitz_kpkmpi0_kpkm_t2kpi",  n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	boundary_1 = 1.6;  
        boundary_2 = 2.0;
        boundary_3 = 2.5; 

	//*****************************************************
	// Breaking things into four different regions in terms of kpkmpi0 masses
	
	dalitz_kpkmpi0_kpkm_t2k_1 = new TH3I("dalitz_kpkmpi0_kpkm_t2k_1",  "dalitz_kpkmpi0_kpkm_t2k_1", n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2k_2 = new TH3I("dalitz_kpkmpi0_kpkm_t2k_2",  "dalitz_kpkmpi0_kpkm_t2k_2", n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2k_3 = new TH3I("dalitz_kpkmpi0_kpkm_t2k_3",  "dalitz_kpkmpi0_kpkm_t2k_3", n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2k_4 = new TH3I("dalitz_kpkmpi0_kpkm_t2k_4",  "dalitz_kpkmpi0_kpkm_t2k_4", n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2kpi_1 = new TH3I("dalitz_kpkmpi0_kpkm_t2kpi_1", "dalitz_kpkmpi0_kpkm_t2kpi_1",  n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2kpi_2 = new TH3I("dalitz_kpkmpi0_kpkm_t2kpi_2", "dalitz_kpkmpi0_kpkm_t2kpi_2",  n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2kpi_3 = new TH3I("dalitz_kpkmpi0_kpkm_t2kpi_3", "dalitz_kpkmpi0_kpkm_t2kpi_3",  n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);

	dalitz_kpkmpi0_kpkm_t2kpi_4 = new TH3I("dalitz_kpkmpi0_kpkm_t2kpi_4", "dalitz_kpkmpi0_kpkm_t2kpi_4",  n2kpi_dalitz, min2kpi, max2kpi, n2k_dalitz, min2k, max2k, nt_dalitz, mint, maxt);


	//********************************************

	gDirectory->cd("/");

	//-----------------------------------
	// DIRC Plot

	gDirectory->mkdir("DIRC_PID");
	gDirectory->cd("DIRC_PID");

	hDIRC_l_photon = new TH1I("dirc_l_photon",  "", 100, 0, 1000);
	hDIRC_l_pi     = new TH1I("dirc_l_pi",  "", 300, -600, 0);
	hDIRC_l_K      = new TH1I("dirc_l_K",  "", 300, -600, 0);
	hPID_l_pi_momentum  = new TH2I("dirc_l_pi_momentum",  "", 500, 0, 10, 300, -600, 0);
	hPID_l_K_momentum  = new TH2I("dirc_l_K_momentum",  "", 500, 0, 10, 300, -600, 0);

	hPID_extrpl_track_X_Y  = new TH2I("dirc_extrapolated_X_Y", "", 100, -200, 200, 100, -200, 200);

	/// Note that rho omega containation limi is set to be 1.17 <= 2k_mass <=1.27
	hPID_l_diff_momentum_rho  = new TH2I("dirc_l_diff_momentum_rho",  "dirc_l_diff_momentum_rho", 500, 0, 10, 200, -500, 500);
	hPID_l_diff_momentum_omega  = new TH2I("dirc_l_diff_momentum_omega",  "dirc_l_diff_momentum_omega", 500, 0, 10, 200, -500, 500);

	hPID_l_diff_momentum  = new TH2I("dirc_l_diff_momentum",  "", 500, 0, 10, 200, -500, 500);

	hPID_l_diff_momentum_phi  = new TH2I("dirc_l_diff_momentum_phi",  "", 500, 0, 10, 200, -500, 500);
	hPID_l_diff_momentum_below_phi  = new TH2I("dirc_l_diff_momentum_below_phi",  "", 500, 0, 10, 200, -500, 500);
	hPID_l_diff_momentum_above_phi  = new TH2I("dirc_l_diff_momentum_above_phi",  "", 500, 0, 10, 200, -500, 500);
	hPID_momentum  = new TH1I("dirc_momentum",  "", 500, 0, 10);

	h2kpi0_invariantmass_below_phi_dirc = new TH1I("h2kpi0_invariantmass_below_phi_dirc", "2kpi0_invariantmass_blow_phi_dirc", n2kpi, min2kpi, max2kpi);

	h2k_invariantmass_dirc = new TH1I("h2k_invariantmass_dirc", "2k_invariantmass_dirc", n2k, min2k, max2k);

	h2kpi0_invariantmass_dirc = new TH1I("h2kpi0_invariantmass_dirc", "2kpi0_invariantmass_dirc", n2kpi, min2kpi, max2kpi);


	gDirectory->mkdir("km");
	gDirectory->cd("km");

	hPID_l_diff_momentum_kpcut_km  = new TH2I("dirc_l_diff_momentum_kpcut_km",  "", 500, 0, 10, 200, -500, 500);

	hPID_momentum_km  = new TH1I("dirc_momentum_km",  "", 500, 0, 10);

	hPID_l_diff_momentum_km  = new TH2I("dirc_l_diff_momentum_km",  "", 500, 0, 10, 200, -500, 500);

	hPID_l_diff_momentum_km_phi  = new TH2I("dirc_l_diff_momentum_km_phi",  "", 500, 0, 10, 200, -500, 500);
	hPID_l_diff_momentum_km_below_phi  = new TH2I("dirc_l_diff_momentum_km_below_phi",  "", 500, 0, 10, 200, -500, 500);
	hPID_l_diff_momentum_km_above_phi  = new TH2I("dirc_l_diff_momentum_km_above_phi",  "", 500, 0, 10, 200, -500, 500);

	hPID_l_momentum_Kp_Km = new TH2I("momentum_Kp_Km",  "", 500, 0, 10, 500, 0, 10);




//	exit(0);

	gDirectory->cd("/");

	//-----------------------------------
	// Dalitz Plot

	gDirectory->mkdir("Angle");
	gDirectory->cd("Angle");

	// Angle implementation

	hAngle_GJ_Phi_pi0 = new TH2D("cos_GJ_angle_phi_pi0", "cos_GJ_angle_phi_pi0", 150, 1, 4, 100, -1, 1);

	hAngle_GJ_a0_pi0  = new TH2D("cos_GJ_angle_a0_pi0", "cos_GJ_angle_a0_pi0", 150, 1, 4, 100, -1, 1);

	hAngle_GJ_Kstar_Kplus  = new TH2D("cos_GJ_angle_kstar_kplus", "cos_GJ_angle_kstar_kplus", 100, 1, 4, 100, -1, 1);

	hAngle_GJ_Kstar_Kminus = new TH2D("cos_GJ_angle_kstar_kminus", "cos_GJ_angle_kstar_kminus", 100, 1, 4, 100, -1, 1);

	hAngle_TY_Phi_pi0 = new TH2D("TY_angle_phi_pi0", "TY_angle_phi_pi0", 150, 1, 4, 361, -180, 180);


	gDirectory->cd("/");
	gDirectory->mkdir("bar_DeltaLL");
	gDirectory->cd("bar_DeltaLL");

//	UInt_t len = sizeof(hist_deltaLL_momemtum_bar)/sizeof(hist_deltaLL_momemtum_bar[0]);

	for (UInt_t i=0; i < num_bars; i++) {
//		cout << i << endl;
		TString plot_name;
		plot_name.Form("hist_deltaLL_momemtum_bar_%i", i);
		hist_deltaLL_momemtum_bar[i]  = new TH2I(plot_name,  "", 500, 0, 10, 200, -500, 500);

	}




	/*--------------------------------------------------*/
	// Test histograms

	gDirectory->cd("/");

	hist_extrapolatedY_check = new TH1I("hist_extrapolatedY_check",  "",  100, -200, 200);

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	pi_const = 3.1415926;
	





	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

/**********************************************************/

Bool_t DSelector_pi0kpkm::Process(Long64_t locEntry) {


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

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)	{

//		continue;

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
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();


//		cout << "KPlus: " << dKPlusWrapper->Get_TrackID()    << "  " <<  dKPlusWrapper->Get_P4().M()  << endl;
//		cout << "KMinus: " << dKMinusWrapper->Get_TrackID()  << "  " <<  dKMinusWrapper->Get_P4().M() << endl;
//		cout << "KProton: " << dProtonWrapper->Get_TrackID() << "  " <<  dProtonWrapper->Get_P4().M() << endl;


//		exit(0);

		
		// continue;

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
		TLorentzVector locKMinusP4 = dKMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();


		//Step 1
		// ***********************************
		/// Problem right here
//		TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
//		continue;

		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();

		TLorentzVector locDecayingPi0P4 = locPhoton1P4 + locPhoton2P4;

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();

		TLorentzVector locDecayingPi0P4_Measured = locPhoton1P4_Measured + locPhoton2P4_Measured;

		// ***********************************
		/// Intermidiate decay state and CM transformation

		TVector3 boost_vec = (locBeamP4 + dTargetP4).BoostVector();

		TLorentzVector locDecayingPhiP4 = locKPlusP4 + locKMinusP4;
		TLorentzVector locDecayinga0P4  = locKPlusP4 + locKMinusP4;

		TLorentzVector locDecayingKstarpP4 =  locKPlusP4 + locDecayingPi0P4_Measured;
		TLorentzVector locDecayingKstarmP4 =  locKMinusP4 + locDecayingPi0P4_Measured;

		TLorentzVector locDecayingPhiP4_CM = locDecayingPhiP4;
		TLorentzVector locDecayinga0P4_CM  = locDecayinga0P4;

		TLorentzVector locBeamP4_CM  = locBeamP4;

		locDecayingPhiP4_CM.Boost(-boost_vec);   
		locDecayinga0P4_CM.Boost(-boost_vec);   

		locBeamP4_CM.Boost(-boost_vec);   

		//continue;

		///**********************************************
		// X -> phi + pi0

		TLorentzVector locDecaying_X_Phi_Pi0_P4 = locDecayingPhiP4 + locDecayingPi0P4_Measured;  
		
		TVector3 gj_phi_pi0_boost_vec = locDecaying_X_Phi_Pi0_P4.BoostVector(); 

		TLorentzVector locDecayingPhiP4_GJ = locDecayingPhiP4;
		locDecayingPhiP4_GJ.Boost(-gj_phi_pi0_boost_vec);   

		TLorentzVector locDecayingPhiP4_TY_GJ = locDecayingPhiP4_GJ;  
		locDecayingPhiP4_TY_GJ.SetPz(0.0);
			
		
		TLorentzVector locDecayingPhiP4_TY = locDecayingPhiP4_TY_GJ;  
		

		locDecayingPhiP4_TY.Boost(gj_phi_pi0_boost_vec);   


		///**********************************************
		// X -> a0 + pi0
		TLorentzVector locDecaying_X_a0_Pi0_P4 = locDecayinga0P4 + locDecayingPi0P4_Measured;  

		TVector3 gj_a0_pi0_boost_vec = locDecaying_X_a0_Pi0_P4.BoostVector(); 

		TLorentzVector locDecayinga0P4_GJ = locDecayinga0P4;
		locDecayinga0P4_GJ.Boost(-gj_a0_pi0_boost_vec);   


		///**********************************************
		// X -> Kstar + K
		TLorentzVector locDecaying_X_Kstarp_Km_P4 =  locDecayingKstarpP4 + locKMinusP4;  
		TVector3 gj_X_Kstarp_Km_boost_vec = locDecaying_X_Kstarp_Km_P4.BoostVector(); 

		TLorentzVector locDecaying_X_Kstarm_Kp_P4 =  locDecayingKstarmP4 + locKPlusP4;  
		TVector3 gj_X_Kstarm_Kp_boost_vec = locDecaying_X_Kstarm_Kp_P4.BoostVector(); 

		TLorentzVector locDecayingKstarpP4_GJ = locDecayingKstarpP4;
		locDecayingKstarpP4_GJ.Boost(-gj_X_Kstarp_Km_boost_vec);   

		TLorentzVector locDecayingKstarmP4_GJ = locDecayingKstarmP4;
		locDecayingKstarmP4_GJ.Boost(-gj_X_Kstarm_Kp_boost_vec);   


//		TVector3 gj_a0_pi0_boost_vec = locDecaying_X_a0_Pi0_P4.BoostVector(); 

//		TLorentzVector locDecayinga0P4_GJ = locDecayinga0P4;
//		locDecayinga0P4_GJ.boost(-gj_a0_pi0_boost_vec);   










		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;



//		if(dKPlusWrapper->Get_Track_NumPhotons_DIRC() <= 5.0) //if the active combo fails a cut, IsComboCutFlag automatically set
//		if(dKPlusWrapper->Get_Track_NumPhotons_DIRC() <=20.0) //if the active combo fails a cut, IsComboCutFlag automatically set
// 		if(dKPlusWrapper->Get_Track_NumPhotons_DIRC() <= 1.0 || dKPlusWrapper->Get_P4_Measured().P() > 3.5 || dKPlusWrapper->Get_P4_Measured().P() < 1.5) //if the active combo fails a cut, IsComboCutFlag automatically set
// 			continue;
// 
// 
// 		if(dKMinusWrapper->Get_Track_NumPhotons_DIRC() <= 1.0 || dKMinusWrapper->Get_P4_Measured().P() > 3.5 || dKMinusWrapper->Get_P4_Measured().P() < 1.5) //if the active combo fails a cut, IsComboCutFlag automatically set
// 			continue;

// 		if(dKPlusWrapper->Get_Track_ExtrapolatedX_DIRC() > 300 || dKMinusWrapper->Get_ExtrapolatedX_DIRC() > 300) //if the active combo fails a cut, IsComboCutFlag automatically set

//		cout << dKPlusWrapper->Get_Track_ExtrapolatedX_DIRC()  << "   " << dKPlusWrapper->Get_Track_ExtrapolatedY_DIRC() << endl;

//		cout << fabs(dKPlusWrapper->Get_Track_ExtrapolatedX_DIRC()) << "   " << fabs(dKPlusWrapper->Get_Track_ExtrapolatedY_DIRC()) << endl;

//		exit(0);


		
		hist_extrapolatedY_check->Fill(dKPlusWrapper->Get_Track_ExtrapolatedY_DIRC());


 		if( fabs(dKPlusWrapper->Get_Track_ExtrapolatedX_DIRC()) > 200 && fabs(dKPlusWrapper->Get_Track_ExtrapolatedY_DIRC()) > 200) //if the active combo fails a cut, IsComboCutFlag automatically set
 			continue;

 		if( dKPlusWrapper->Get_Track_NumPhotons_DIRC() < 5 && dKPlusWrapper->Get_Track_NumPhotons_DIRC() < 5 ) //if the active combo fails a cut, IsComboCutFlag automatically set
 			continue;



//		cout << fabs(dKPlusWrapper->Get_Track_ExtrapolatedX_DIRC()) << "   " << fabs(dKPlusWrapper->Get_Track_ExtrapolatedY_DIRC()) << endl;

//		exit(0);



 
//		if(dKPlusWrapper->Get_Track_NumPhotons_DIRC() <= 1.0 || dKPlusWrapper->Get_P4_Measured().P() > 3.5 || dKPlusWrapper->Get_P4_Measured().P() < 1.5) //if the active combo fails a cut, IsComboCutFlag automatically set
// 			continue;


		//**************************************
		//**************************************
		//
		/// Next is Andrew's clean cut to reject 2pi containitation		
		/// 		

//		if ( (dKPlusWrapper->Get_Detector_System_Timing() == SYS_TOF && (locKPlusP4.Vect()).Mag() < 2.0) || (dKMinusWrapper->Get_Detector_System_Timing() == SYS_TOF && (locKMinusP4.Vect()).Mag() < 2.0) ) { 


		//if you manually execute any actions, and it fails a cut, be sure to call:
		//dComboWrapper->Set_IsComboCut(true);
		//
	        TLorentzVector lockppi0P4 = locKPlusP4 + locDecayingPi0P4;
	        TLorentzVector lockmpi0P4 = locKMinusP4 + locDecayingPi0P4;
	        TLorentzVector loc2kP4 = locKPlusP4 + locKMinusP4;
	        TLorentzVector loc2kpi0P4 = locKPlusP4 + locKMinusP4 + locDecayingPi0P4;
	        TLorentzVector locProtonPi0P4 = locProtonP4 + locDecayingPi0P4;
		TLorentzVector locProtonKMinusP4 = locProtonP4 + locKMinusP4;
	
	        double t    = (dTargetP4 - locProtonP4).M2();
	        double t_kk = (locBeamP4 - loc2kP4).M2();




		h2k_invariantmass->Fill(loc2kP4.M());
		h2kpi0_invariantmass->Fill(loc2kpi0P4.M());
	
		hkppi0_invariantmass->Fill(lockppi0P4.M());
		hkmpi0_invariantmass->Fill(lockmpi0P4.M()); 
	
		hpi0_mass->Fill(locDecayingPi0P4.M());
	
		hpi0_mass_measured->Fill(locDecayingPi0P4_Measured.M());
	
//	        if(loc2kpi0P4.M() > 1.23 && loc2kpi0P4.M() < 1.33) {
//	
//			dalitz_pikp_pikm_low->Fill( pow(lockmpi0P4.M(),2), pow(lockppi0P4.M(),2) );
//	
//		}  else if(loc2kpi0P4.M() > 1.36 && loc2kpi0P4.M() < 1.46) {
//	
//			dalitz_pikp_pikm_mid->Fill( pow(lockmpi0P4.M(),2), pow(lockppi0P4.M(),2) ); 
//			
//		}  else if(loc2kpi0P4.M() > 1.55 && loc2kpi0P4.M() < 1.75) {
//	
//			dalitz_pikp_pikm_high->Fill( pow(lockmpi0P4.M(),2), pow(lockppi0P4.M(),2) );
//	
//		}  else if(loc2kpi0P4.M() > 1.9 && loc2kpi0P4.M() < 2.8) {
//	
//			dalitz_pikp_pikm_ultra->Fill( pow(lockmpi0P4.M(),2), pow(lockppi0P4.M(),2) );
//	
//		}

		dalitz_kppi0_kmpi0->Fill(lockmpi0P4.M(), lockppi0P4.M());
		dalitz_kppi0_kmpi0_sq->Fill(pow(lockmpi0P4.M(),2), pow(lockppi0P4.M(),2));

		dalitz_kpkmpi0_kpkm->Fill(loc2kP4.M(), loc2kpi0P4.M());
		dalitz_kpkmpi0_kpkm_sq->Fill(pow(loc2kP4.M(),2), pow(loc2kpi0P4.M(),2));

		dalitz_kpkm_kmpi0->Fill(loc2kP4.M(), lockmpi0P4.M());
		dalitz_kpkm_kmpi0_sq->Fill(pow(loc2kP4.M(),2), pow(lockmpi0P4.M(),2));

		dalitz_kpkm_kppi0->Fill(loc2kP4.M(), lockppi0P4.M());
		dalitz_kpkm_kppi0_sq->Fill(pow(loc2kP4.M(),2), pow(lockppi0P4.M(),2));

		dalitz_kppi0_kpkmpi0->Fill(lockppi0P4.M(), loc2kpi0P4.M());
		dalitz_kppi0_kpkmpi0_sq->Fill(pow(lockppi0P4.M(),2), pow(loc2kpi0P4.M(),2));

		dalitz_kmpi0_kpkmpi0->Fill(lockmpi0P4.M(), loc2kpi0P4.M());
		dalitz_kmpi0_kpkmpi0_sq->Fill(pow(lockmpi0P4.M(),2), pow(loc2kpi0P4.M(),2));

		dalitz_t2kpi_t2k_kpkmpi0->Fill(-t, -t_kk, loc2kpi0P4.M());
		dalitz_t2kpi_t2k_kpkm->Fill(-t, -t_kk, loc2kP4.M());

		dalitz_kpkmpi0_kpkm_t2k->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t_kk);
		dalitz_kpkmpi0_kpkm_t2kpi->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t);


//        	if( loc2kpi0P4.M() <= boundary_1 ) { 
//	            dalitz_kpkmpi0_kpkm_t2k_1->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t_kk);
//		    dalitz_kpkmpi0_kpkm_t2kpi_1->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t);
//		} else if( loc2kpi0P4.M() > boundary_1 && loc2kpi0P4.M() <= boundary_2) {
//	            dalitz_kpkmpi0_kpkm_t2k_2->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t_kk);
//		    dalitz_kpkmpi0_kpkm_t2kpi_2->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t);
//
//		} else if(loc2kpi0P4.M() > boundary_2 && loc2kpi0P4.M() <= boundary_3) {
//		    dalitz_kpkmpi0_kpkm_t2k_3->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t_kk);
//		    dalitz_kpkmpi0_kpkm_t2kpi_3->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t);
//		} else {
//	            dalitz_kpkmpi0_kpkm_t2k_4->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t_kk);
//		    dalitz_kpkmpi0_kpkm_t2kpi_4->Fill(loc2kpi0P4.M(), loc2kP4.M(), -t);
//		}
		
        	if( loc2kP4.M() <= 1.010 ) {

			h2kpi0_invariantmass_a->Fill(loc2kpi0P4.M());
			
			h2gamma_2kpi0_mass->Fill(locDecayingPi0P4_Measured.M(), loc2kpi0P4.M());

			dHist_2kpi_t->Fill(-t);
			dHist_2k_t->Fill(-t_kk);
			
			if (locDecayingPi0P4_Measured.M() >= 0.11 && locDecayingPi0P4_Measured.M() <= 0.16){
				h2kpi0_invariantmass_below_phi_pi0_cut->Fill(loc2kpi0P4.M());

			} else if(locDecayingPi0P4_Measured.M() >= 0.08 && locDecayingPi0P4_Measured.M() <= 0.1) {
				h2kpi0_invariantmass_below_phi_pi0_sideband->Fill(loc2kpi0P4.M());
				h2kpi0_invariantmass_below_phi_pi0_hi_sideband->Fill(loc2kpi0P4.M());
			} else if(locDecayingPi0P4_Measured.M() >= 0.17 && locDecayingPi0P4_Measured.M() <= 0.19) {

				h2kpi0_invariantmass_below_phi_pi0_sideband->Fill(loc2kpi0P4.M());
				h2kpi0_invariantmass_below_phi_pi0_hi_sideband->Fill(loc2kpi0P4.M());


				Double_t theta= locDecayinga0P4_GJ.Theta();
				hAngle_GJ_a0_pi0->Fill(loc2kpi0P4.M(), cos(theta));

			}

		} else if(loc2kP4.M() > 1.010 && loc2kP4.M() <= 1.025) {

			// Selecting phi

			h2kpi0_invariantmass_b->Fill(loc2kpi0P4.M());
 			dHist_2kpiMass_t->Fill(-t, loc2kpi0P4.M());
 			dHist_t2k_2kpiMass->Fill(-t_kk, loc2kpi0P4.M());
 			dHist_t_t2k->Fill(-t, -t_kk);

//			hAngle_GJ_Phi_pi0->Fill(loc2kpi0P4.M(), cos(locDecayingPhiP4_CM.Theta()));
//
//

//			hAngle_GJ_Phi_pi0->Fill(loc2kpi0P4.M(), cos(theta));

////			cout << loc2kpi0P4.M() << "     " << cos(locDecayingPhiP4_CM.Theta()) << "   " << locDecayingPhiP4_CM.Theta() << endl;
//
//
//			cout << locDecayingPhiP4.E()  << "  ";
//			cout << locDecayingPhiP4.Px() << "  ";
//			cout << locDecayingPhiP4.Py() << "  ";
//			cout << locDecayingPhiP4.Pz() << endl;
//
//			cout << locBeamP4.E()  << "  ";
//			cout << locBeamP4.Px() << "  ";
//			cout << locBeamP4.Py() << "  ";
//			cout << locBeamP4.Pz() << endl;
//
//			cout << boost_vec.Px()  << "  ";
//			cout << boost_vec.Py()  << "  ";
//			cout << boost_vec.Pz()  << endl;
//			
//			locDecayingPhiP4_CM.Boost(-boost_vec);
//
//			cout << locDecayingPhiP4_CM.E()  << "  ";
//			cout << locDecayingPhiP4_CM.Px() << "  ";
//			cout << locDecayingPhiP4_CM.Py() << "  ";
//			cout << locDecayingPhiP4_CM.Pz() << endl;
//
//			exit(0);

//			Double_t theta= locDecayingPhiP4_CM.Angle(locBeamP4_CM.Vect());
//			Double_t theta= locDecayingPhiP4_CM.Angle(locBeamP4_CM.Vect());
//
			Double_t theta= locDecayingPhiP4_GJ.Theta();
			hAngle_GJ_Phi_pi0->Fill(loc2kpi0P4.M(), cos(theta));

			Double_t phi= locDecayingPhiP4_TY.Phi()*180.0 / pi_const;
			hAngle_TY_Phi_pi0->Fill(loc2kpi0P4.M(), phi);
			

		} else if(loc2kP4.M() > 1.025 && loc2kP4.M() <= 1.040) {

			h2kpi0_invariantmass_c->Fill(loc2kpi0P4.M());
		
		}


	
		if(lockppi0P4.M() > 0.75 && lockppi0P4.M() < 0.90) {

			Double_t theta_p= locDecayingKstarpP4_GJ.Theta();
			hAngle_GJ_Kstar_Kminus->Fill(loc2kpi0P4.M(), cos(theta_p));
		
		}


		if(lockmpi0P4.M() > 0.75 && lockmpi0P4.M() < 0.90) {

			Double_t theta_m= locDecayingKstarmP4_GJ.Theta();
			hAngle_GJ_Kstar_Kplus ->Fill(loc2kpi0P4.M(), cos(theta_m));

		}


//         if(loc2kP4.M() < 1.04 && loc2kP4.M() > 1.00) {
// 			dHist_2kpiMass_t->Fill(-t, loc2kpi0P4.M());
// 			dHist_2kpiMass_Egamma->Fill(locBeamP4.E(), loc2kpi0P4.M());
//         } else if(loc2kP4.M() > 1.08 && loc2kP4.M() < 1.12) {
// 			dHist_2kpiMass_t_SB->Fill(-t, loc2kpi0P4.M());
// 			dHist_2kpiMass_Egamma_SB->Fill(locBeamP4.E(), loc2kpi0P4.M());
//         }


//          if(loc2kP4.M() < 1.04 && loc2kP4.M() > 1.00) {
// 			dHist_2kpiMass_t->Fill(-t, loc2kpi0P4.M());
//          }
         

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

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[KPlus].insert(locKPlusTrackID);
		locUsedThisCombo_MissingMass[KMinus].insert(locKMinusTrackID);
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
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

//	   }

//		cout << "DIRC photon #: " << dKPlusWrapper->Get_Track_NumPhotons_DIRC() << endl;
//		cout << "Lpi: " << dKPlusWrapper->Get_Track_Lpi_DIRC() << endl;
//		cout << "Lk: " << dKPlusWrapper->Get_Track_Lk_DIRC() << endl;




//	} // end of combo loop




// 		if(dKPlusWrapper->Get_Track_NumPhotons_DIRC() > 0.0){
// 
// 			cout << "Lpi: " << dKPlusWrapper->Get_Track_Lpi_DIRC() << endl;
// 			cout << "Lk: " << dKPlusWrapper->Get_Track_Lk_DIRC() << endl;
// 
// 		}








///*--------------------------------------------------*/
//	//---------------------------------
//	// DIRC related
//
//
	hDIRC_l_photon->Fill(dKPlusWrapper->Get_Track_NumPhotons_DIRC());
//
//	//------------------------------------
//	// kp 
//
		if (dKPlusWrapper->Get_Track_NumPhotons_DIRC() > 5.0) {



			hPID_extrpl_track_X_Y->Fill(dKPlusWrapper->Get_Track_ExtrapolatedX_DIRC(), dKPlusWrapper->Get_Track_ExtrapolatedY_DIRC());

			hDIRC_l_pi->Fill(dKPlusWrapper->Get_Track_Lpi_DIRC());
			hDIRC_l_K->Fill(dKPlusWrapper->Get_Track_Lk_DIRC());
	
			hPID_l_pi_momentum->Fill(dKPlusWrapper->Get_P4_Measured().P() , dKPlusWrapper->Get_Track_Lpi_DIRC());
			hPID_l_K_momentum->Fill(dKPlusWrapper->Get_P4_Measured().P(), dKPlusWrapper->Get_Track_Lk_DIRC());
	
			hPID_momentum->Fill(dKPlusWrapper->Get_P4_Measured().P());

			Float_t lk_diff = dKPlusWrapper->Get_Track_Lpi_DIRC() - dKPlusWrapper->Get_Track_Lk_DIRC();
			
			hPID_l_diff_momentum->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);

			if (loc2kP4.M() > 1.006 && loc2kP4.M() < 1.034) {
				hPID_l_diff_momentum_phi->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
			} else if (loc2kP4.M() <= 1.006) {
				hPID_l_diff_momentum_below_phi->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
			} else {
				hPID_l_diff_momentum_above_phi->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
			}

			if (loc2kP4.M() <= 1.27 && loc2kP4.M() >= 1.17) {
				hPID_l_diff_momentum_rho->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
			}

			if (loc2kpi0P4.M() <= 1.5 && loc2kpi0P4.M() >= 1.3) {
			
				hPID_l_diff_momentum_omega->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
			}

			/*--------------------------------------------------*/
			// Flag for pi0pippiminus


			if (dKPlusWrapper->Get_P4().M() < 0.4) {	

 				if (loc2kP4.M() > 1.006 && loc2kP4.M() < 1.034) {
 					hPID_l_diff_momentum_phi->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
 				} else if (loc2kP4.M() <= 1.006) {
 					hPID_l_diff_momentum_below_phi->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
 				} else {

					hPID_l_diff_momentum_above_phi->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
				}

				if (loc2kP4.M() <= 0.810 && loc2kP4.M() >= 0.675) {
					hPID_l_diff_momentum_rho->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
				}

				if (loc2kpi0P4.M() <= 0.85 && loc2kpi0P4.M() >= 0.7) {
				
					hPID_l_diff_momentum_omega->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);
				}

			}




// 			for ((Get_Track_ExtrapolatedY_DIRC() + 200))
 
			/*--------------------------------------------------*/
			/// which bar?
			Float_t box_boundary = 100;	
		
			
			Int_t det_bar_num = (dKPlusWrapper->Get_Track_ExtrapolatedY_DIRC()+ box_boundary) / (2*box_boundary/num_bars) ;

///			cout << det_bar_num << endl;
//
			if (det_bar_num < 0) {
				det_bar_num = 0;
 			} else if (det_bar_num >= num_bars) {
				det_bar_num = num_bars - 1;
			}
 			hist_deltaLL_momemtum_bar[det_bar_num]->Fill(dKPlusWrapper->Get_P4_Measured().P(), lk_diff);

//			exit(0);

 	//------------------------------------
 	// km 
 
 			if (dKMinusWrapper->Get_Track_NumPhotons_DIRC() > 5.0) {

				Float_t lk_diff_km = dKMinusWrapper->Get_Track_Lpi_DIRC() - dKMinusWrapper->Get_Track_Lk_DIRC();
 
 				hPID_momentum_km->Fill(dKMinusWrapper->Get_P4_Measured().P());
 
  				Float_t lk_diff = dKMinusWrapper->Get_Track_Lpi_DIRC() - dKMinusWrapper->Get_Track_Lk_DIRC();

				hPID_l_diff_momentum_km->Fill(dKMinusWrapper->Get_P4_Measured().P(), lk_diff_km);

//  				if (loc2kP4.M() > 1.006 && loc2kP4.M() < 1.034) {
//  					hPID_l_diff_momentum_km_phi->Fill(dKMinusWrapper->Get_P4_Measured().P(), lk_diff);
//  				} else if (loc2kP4.M() <= 1.006) {
//  					hPID_l_diff_momentum_km_below_phi->Fill(dKMinusWrapper->Get_P4_Measured().P(), lk_diff);
//  				} else {
//  					hPID_l_diff_momentum_km_above_phi->Fill(dKMinusWrapper->Get_P4_Measured().P(), lk_diff);
//  				}
 
 			}

//			/*--------------------------------------------------*/
//
//			if(lk_diff < -20) {
//	 		       if(loc2kP4.M() < 1.006 ) {
//					
//					if (dKPlusWrapper->Get_P4_Measured().P() < 3.2) {
//						h2kpi0_invariantmass_below_phi_dirc->Fill(loc2kpi0P4.M());
//					}
//				}
//			
//				Float_t lk_diff_km = dKMinusWrapper->Get_Track_Lpi_DIRC() - dKMinusWrapper->Get_Track_Lk_DIRC();
//				
//				
//				h2k_invariantmass_dirc->Fill(loc2kP4.M());
//				h2kpi0_invariantmass_dirc->Fill(loc2kpi0P4.M());
//
//			} else {
//			
//			
//	 		       	if(loc2kP4.M() < 1.006) {
//					if (dKPlusWrapper->Get_P4_Measured().P() < 3.2) {
////			       			h2kpi0_invariantmass_below_phi_pion->Fill(loc2kpi0P4.M());
////			       			h2kpi0_invariantmass_below_phi_pion->Fill(loc2kpi0P4.M());
//					}
//				}
//			}
//

		} else {

		/// Outside of the DIRC


		}

// 	//------------------------------------
// 	// km 
// 
// 		if (dKMinusWrapper->Get_Track_NumPhotons_DIRC() > 0.0) {
// 
// 			hPID_momentum_km->Fill(dKMinusWrapper->Get_P4_Measured().P());
// 
// 			Float_t lk_diff = dKMinusWrapper->Get_Track_Lpi_DIRC() - dKMinusWrapper->Get_Track_Lk_DIRC();
// 
// 			if (loc2kP4.M() > 1.006 && loc2kP4.M() < 1.034) {
// 
// 				hPID_l_diff_momentum_km_phi->Fill(dKMinusWrapper->Get_P4_Measured().P(), lk_diff);
// 
// 			} else if (loc2kP4.M() <= 1.006) {
// 
// 				hPID_l_diff_momentum_km_below_phi->Fill(dKMinusWrapper->Get_P4_Measured().P(), lk_diff);
// 
// 			} else {
// 
// 				hPID_l_diff_momentum_km_above_phi->Fill(dKMinusWrapper->Get_P4_Measured().P(), lk_diff);
// 
// 			}
// 
// 
// 		}
// 
// 
// 		if ( dKPlusWrapper->Get_Track_NumPhotons_DIRC() > 0.0 && dKMinusWrapper->Get_Track_NumPhotons_DIRC() > 0.0) {
// 
// 			hPID_l_momentum_Kp_Km->Fill(dKPlusWrapper->Get_P4_Measured().P(), dKMinusWrapper->Get_P4_Measured().P());		
// 
// 		}



//
//	//FILL HISTOGRAMS: Num combos / events surviving actions
//	Fill_NumCombosSurvivedHists();
//
//	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
///*
//	//Thrown beam: just use directly
//	if(dThrownBeam != NULL)
//		double locEnergy = dThrownBeam->Get_P4().E();
//
//	//Loop over throwns
//	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
//	{
//		//Set branch array indices corresponding to this particle
//		dThrownWrapper->Set_ArrayIndex(loc_i);
//
//		//Do stuff with the wrapper here ...
//	}
//*/
//	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
///*
//	//Loop over beam particles (note, only those appearing in combos are present)
//	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
//	{
//		//Set branch array indices corresponding to this particle
//		dBeamWrapper->Set_ArrayIndex(loc_i);
//
//		//Do stuff with the wrapper here ...
//	}
//
//	//Loop over charged track hypotheses
//	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
//	{
//		//Set branch array indices corresponding to this particle
//		dChargedHypoWrapper->Set_ArrayIndex(loc_i);
//
//		//Do stuff with the wrapper here ...
//	}
//
//	//Loop over neutral particle hypotheses
//	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
//	{
//		//Set branch array indices corresponding to this particle
//		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);
//
//		//Do stuff with the wrapper here ...
//	}
//*/
//
//	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
///*
//	Bool_t locIsEventCut = true;
//	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
//		//Set branch array indices for combo and all combo particles
//		dComboWrapper->Set_ComboIndex(loc_i);
//		// Is used to indicate when combos have been cut
//		if(dComboWrapper->Get_IsComboCut())
//			continue;
//		locIsEventCut = false; // At least one combo succeeded
//		break;
//	}
//	if(!locIsEventCut && dOutputTreeFileName != "")
//		Fill_OutputTree();
//*/
//
//
//
////	dKPlusWrapper->Get_Track_Lpi_DIRC();
//
//
//
//
//
//
//
//


	} // end of combo loop

	return kTRUE;
}

void DSelector_pi0kpkm::Finalize(void)
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
