#include "DSelector_pomega2pi_omega3pi.h"

void DSelector_pomega2pi_omega3pi::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "pomega2pi_omega3pi.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = "AmpToolsInputTree.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = "kin"; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.

	//DO THIS NEXT
	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	//THEN THIS
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

	SetupAmpTools_FlatTree();

  	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("t");
 	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M3Pi");
 	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("M4Pi");
 	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("MRecoil");
 	dFlatTreeInterface->Create_Branch_Fundamental<Float_t>("Phi_Prod");

        if((dTreeInterface->Get_Branch("NumCombos") == NULL)) return;

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
        //dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, Pi0, 0.114, 0.156));
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
        //dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, 0, locOmegaPIDs, 0.71, 0.85));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1PlusPIDs, 150, minb1, maxb1, "b1Plus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locb1MinusPIDs, 150, minb1, maxb1, "b1Minus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, 150, minDelta, maxDelta, "DeltaPlusPlus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locDelta0PIDs, 150, minDelta, maxDelta, "Delta0_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locRhoPIDs, 150, 0.0, 1.5, "Rho_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, locb1MinusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "DeltaPlusPlus_b1Minus_OmegaCut"));
	dAnalysisActions.push_back(new DHistogramAction_2DInvariantMass(dComboWrapper, true, 0, locDelta0PIDs, locb1PlusPIDs, 100, minDelta, maxDelta, 150, minb1, maxb1, "Delta0_b1Plus_OmegaCut"));

	//CUT AWAY DELTA CONTRIBUTION
	//dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, 0, locDeltaPlusPlusPIDs, 1.4, 1000.));
	//dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, 0, locDelta0PIDs, 1.4, 1000.));
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
	
	dHist_ThrownTopologies = new TH1F("hThrownTopologies","hThrownTopologies", 10, -0.5, 9.5);

	vector<TString> locThrownTopologies;
	locThrownTopologies.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]");
	locThrownTopologies.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#omega]");		
	locThrownTopologies.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#eta]");
	locThrownTopologies.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0},#omega]");
	locThrownTopologies.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0}]");
	locThrownTopologies.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#phi]");
	locThrownTopologies.push_back("2#gamma#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}p[#pi^{0}]");		
	locThrownTopologies.push_back("#gamma2#pi^{#plus}2#pi^{#minus}p[#eta']");
	locThrownTopologies.push_back("3#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#omega,#eta']");		
	locThrownTopologies.push_back("6#gamma2#pi^{#plus}2#pi^{#minus}p[3#pi^{0},#eta]");
	locThrownTopologies.push_back("2#gamma#pi^{#plus}2#pi^{#minus}K^{#plus}p[#pi^{0},#Lambda]");
	for(uint i=0; i<locThrownTopologies.size(); i++) {
		dHist_InvariantMass_ThrownTopology[locThrownTopologies[i]] = new TH1F(Form("hInvariantMass_ThrownTopology_%d", i),Form("Invariant Mass Topology: %s", locThrownTopologies[i].Data()), 1000, 1.0, 3.0);
	}

	//after omega cut
	dHist_ThrownTopologies_omegacut = new TH1F("hThrownTopologies_omegacut","hThrownTopologies_omegacut", 10, -0.5, 9.5);

	vector<TString> locThrownTopologies_omegacut;
	locThrownTopologies_omegacut.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#omega]");
	locThrownTopologies_omegacut.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0},#omega]");		
	locThrownTopologies_omegacut.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#eta]");
	locThrownTopologies_omegacut.push_back("3#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#eta,#omega]");
	locThrownTopologies_omegacut.push_back("3#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#omega,#eta']");
	locThrownTopologies_omegacut.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0}]");
	locThrownTopologies_omegacut.push_back("4#gamma#pi^{#plus}#pi^{#minus}p[2#pi^{0},#omega]");		
	locThrownTopologies_omegacut.push_back("2#gamma#pi^{#plus}2#pi^{#minus}K^{#plus}p[#pi^{0},#Lambda,#omega]");
	locThrownTopologies_omegacut.push_back("3#gammae^{#plus}e^{#minus}#pi^{#plus}#pi^{#minus}p[2#pi^{0}]");		
	locThrownTopologies_omegacut.push_back("3#gammae^{#plus}e^{#minus}#pi^{#plus}#pi^{#minus}p[2#pi^{0},#omega]");
	locThrownTopologies_omegacut.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0},#eta]");
	for(uint i=0; i<locThrownTopologies_omegacut.size(); i++) {
		dHist_InvariantMass_ThrownTopology_omegacut[locThrownTopologies_omegacut[i]] = new TH1F(Form("hInvariantMass_ThrownTopology_omegacut_%d", i),Form("Invariant Mass Topology: %s", locThrownTopologies_omegacut[i].Data()), 1000, 1.0, 3.0);
	}

	//after omega cut w/out sideband subtraction
	dHist_ThrownTopologies_omegacut_nosideband = new TH1F("hThrownTopologies_omegacut_nosideband","hThrownTopologies_omegacut_nosideband", 10, -0.5, 9.5);

	vector<TString> locThrownTopologies_omegacut_nosideband;
	locThrownTopologies_omegacut_nosideband.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#omega]"); //update these
	locThrownTopologies_omegacut_nosideband.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0},#omega]");		
	locThrownTopologies_omegacut_nosideband.push_back("2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#eta]");
	locThrownTopologies_omegacut_nosideband.push_back("3#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#eta,#omega]");
	locThrownTopologies_omegacut_nosideband.push_back("3#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0},#omega,#eta']");
	locThrownTopologies_omegacut_nosideband.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0}]");
	locThrownTopologies_omegacut_nosideband.push_back("4#gamma#pi^{#plus}#pi^{#minus}p[2#pi^{0},#omega]");		
	locThrownTopologies_omegacut_nosideband.push_back("2#gamma#pi^{#plus}2#pi^{#minus}K^{#plus}p[#pi^{0},#Lambda,#omega]");
	locThrownTopologies_omegacut_nosideband.push_back("3#gammae^{#plus}e^{#minus}#pi^{#plus}#pi^{#minus}p[2#pi^{0}]");		
	locThrownTopologies_omegacut_nosideband.push_back("3#gammae^{#plus}e^{#minus}#pi^{#plus}#pi^{#minus}p[2#pi^{0},#omega]");
	locThrownTopologies_omegacut_nosideband.push_back("4#gamma2#pi^{#plus}2#pi^{#minus}p[2#pi^{0},#eta]");
	for(uint i=0; i<locThrownTopologies_omegacut_nosideband.size(); i++) {
		dHist_InvariantMass_ThrownTopology_omegacut_nosideband[locThrownTopologies_omegacut_nosideband[i]] = new TH1F(Form("hInvariantMass_ThrownTopology_omegacut_nosideband_%d", i),Form("Invariant Mass Topology: %s", locThrownTopologies_omegacut_nosideband[i].Data()), 1000, 1.0, 3.0);
	}

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1F("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1F("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
	dHist_BeamDeltaT = new TH1F("BeamDeltaT", "; t_{Tagger} - t_{RF} (ns)", 200, -20., 20.);
	dHist_Phi_lab_ProtonPiPlus_vs_Ebeam = new TH2F("hPhi_lab_ProtonPiPlus_vs_Ebeam", ";E_{#gamma} (GeV);#Phi_{lab,p#pi^{+}} (deg)", 600, 0., 12., 600, -180., 180.);

	
	//MY HISTOGRAMS
	dHist_2PhotonMass = new TH1F("h2PhotonMass", ";M_{#gamma#gamma} (GeV)", 600, 0.1, 0.17);
	dHist_2PhotonVs3PiMass = new TH2F("h2PhotonVs3PiMass", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);M_{#gamma#gamma} (GeV)", 600, 0.3, 2.3, 600, 0.1, 0.17);
	dHist_2PhotonVs2PiMass = new TH2F("h2PhotonVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);M_{#gamma#gamma} (GeV)", 600, 0., 1.5, 600, 0.1, 0.17);
	dHist_2PhotonVsOmegaPiPlusMass = new TH2F("h2PhotonVsOmegaPiPlusMass", ";M_{#omega#pi^{+}} (GeV);M_{#gamma#gamma} (GeV)", 600, 0.3, 2.3, 600, 0.1, 0.17);
	dHist_2PhotonVsOmegaPiMinusMass = new TH2F("h2PhotonVsOmegaPiMinusMass", ";M_{#omega#pi^{-}} (GeV);M_{#gamma#gamma} (GeV)", 600, 0.3, 2.3, 600, 0.1, 0.17);
	dHist_2PhotonVsOmega2PiMass = new TH2F("h2PhotonVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);M_{#gamma#gamma} (GeV)", 600, 1., 3., 600, 0.1, 0.17);


	dHist_5PiMass = new TH1F("5PiMass", ";M_{#pi^{+}#pi^{+}#pi^{-}#pi^{-}#pi^{0}} (GeV)", 600, 1., 3.);
	dHist_5PiMass_Measured = new TH1F("5PiMass_Measured", ";M_{#pi^{+}#pi^{+}#pi^{-}#pi^{-}#pi^{0}} (GeV)", 600, 1., 3.);
	dHist_3PiMass = new TH1F("3PiMass", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV)", 600, 0.3, 2.3);
	dHist_3PiMass_Measured = new TH1F("3PiMass_Measured", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV)", 600, 0.3, 2.3);
	dHist_Omega2PiMass = new TH1F("Omega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV)", 600, 1., 3.);
	dHist_3vs2PiMass = new TH2F("3vs2PiMass", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);M_{#pi^{+}#pi^{-}} (GeV)", 600, 0.3, 2.3, 600, 0, 1.5);
	dHist_OmegaPiPlusMass = new TH1F("OmegaPiPlusMass", ";M_{#omega#pi^{+}} (GeV)", 600, 0.3, 2.3);
	dHist_OmegaPiMinusMass = new TH1F("OmegaPiMinusMass", ";M_{#omega#pi^{-}} (GeV)", 600, 0.3, 2.3);
	dHist_OmegaPiPlusVsOmegaPiMinus = new TH2F("OmegaPiPlusVsOmegaPiMinus", ";M_{#omega#pi^{+}};M_{#omega#pi^{-}}", 600, 0.3, 2.3, 600, 0.3, 2.3);
	dHist_ProtonPiPlusMass = new TH1F("ProtonPiPlusMass", ";M_{p#pi^{+}} (GeV)", 600, 0.3, 3.);
	dHist_ProtonPiPlusVsOmegaPiMinus = new TH2F("ProtonPiPlusVsOmegaPiMinus", ";M_{#omega#pi^{-}} (GeV);M_{p#pi^{+}} (GeV)", 600, 0.3, 2.3, 600, 0.3, 3.);

      	dHist_ProtonPiMinusMass = new TH1F("ProtonPiMinusMass", ";M_{p#pi^{-}} (GeV)", 600, 0.3, 3.);
	dHist_ProtonPiMinusVsOmegaPiPlus = new TH2F("ProtonPiMinusVsOmegaPiPlus", ";M_{#omega#pi^{+}} (GeV);M_{p#pi^{-}} (GeV)", 600, 0.3, 2.3, 600, 0.5, 3.);
	dHist_3vs4PiMass_plus = new TH2F("3vs4PiMass_plus", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);M_{#pi^{+}#pi^{+}#pi^{-}#pi^{0}} (GeV)", 600, 0.3, 2.3, 600, 0.3, 2.3);
	dHist_3vs4PiMass_minus = new TH2F("3vs4PiMass_minus", ";M_{#pi^{+}#pi^{-}#pi^{0}} (GeV);M_{#pi^{+}#pi^{-}#pi^{-}#pi^{0}} (GeV)", 600, 0.3, 2.3, 600, 0.3, 2.3);

	dHist_OmegaPiPlusVsOmega2PiMass = new TH2F("OmegaPiPlusVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);M_{#omega#pi^{+}} (GeV)", 600, 1., 3., 600, 0.3, 2.3);
	dHist_OmegaPiMinusVsOmega2PiMass = new TH2F("OmegaPiMinusVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);M_{#omega#pi^{-}} (GeV)", 600, 1., 3., 600, 0.3, 2.3);

	dHist_t_proton_total = new TH1F("t_proton_total", ";-t_{p} (GeV)^{2}", 600, 0.0, 3.); 
	dHist_t_omega_total = new TH1F("t_omega_total", ";-t_{#omega} (GeV)^{2}", 600, 0., 3.);
	dHist_t_Deltaplusplus_total = new TH1F("t_Deltaplusplus_total", ";-t_{#Delta^{++}} (GeV)^{2}", 600, 0., 3.);
	dHist_t_Delta0_total = new TH1F("t_Delta0_total", ";-t_{#Delta^{0}} (GeV)^{2}", 600, 0., 3.);
	for(int i = 0; i < 10; i++){
	  dHist_t_proton[i] = new TH1F(Form("t_proton_%i", i+1), ";-t_{p} (GeV)^{2}", 600, 0.0, 3.0); 
	  dHist_t_omega[i] = new TH1F(Form("t_omega_%i", i+1), ";-t_{#omega} (GeV)^{2}", 600, 0., 3.);
	  dHist_t_Deltaplusplus[i] = new TH1F(Form("t_Deltaplusplus_%i", i+1), ";-t_{#Delta^{++}} (GeV)^{2}", 600, 0., 3.);
	  dHist_t_Delta0[i] = new TH1F(Form("t_Delta0_%i", i+1), ";-t_{#Delta^{0}} (GeV)^{2}", 600, 0., 3.);
	}

	dHist_t_pomega_correlation = new TH2F("t_pomega_correlation", ";-t_{p} (GeV)^{2};-t_{#omega} (GeV)^{2}", 600, 0., 3., 600, 0., 3.);
	dHist_t_Delta_correlation_total = new TH2F("t_Delta_correlation_total", ";-t_{#Delta^{++}} (GeV)^{2};-t_{#Delta^{0}} (GeV)^{2}", 600, 0., 3., 600, 0., 3.);
	for(int i = 0; i < 10; i++){
	  dHist_t_Delta_correlation[i] = new TH2F(Form("t_Delta_correlation_%i", i+1), ";-t_{#Delta^{++}} (GeV)^{2};-t_{#Delta^{0}} (GeV)^{2}", 600, 0., 3., 600, 0., 3.);
	}

	dHist_2PiMass = new TH1F("2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 1.5);
	dHist_2PiMass_plus = new TH1F("2PiMass_plus", ";M_{#pi^{+}#pi^{0}} (GeV)", 600, 0., 1.5);
	dHist_2PiMass_minus = new TH1F("2PiMass_minus", ";M_{#pi^{-}#pi^{0}} (GeV)", 600, 0., 1.5);
	dHist_2PiMass_plus_omegacut = new TH1F("2PiMass_plus_omegacut", ";M_{#pi^{+}#pi^{0}} (GeV)", 600, 0., 1.5);
	dHist_2PiMass_minus_omegacut = new TH1F("2PiMass_minus_omegacut", ";M_{#pi^{-}#pi^{0}} (GeV)", 600, 0., 1.5);
	dHist_2PiMass_plus_notomega = new TH1F("2PiMass_plus_notomega", ";M_{#pi^{+}#pi^{0}} (GeV)", 600, 0., 1.5);
	dHist_2PiMass_minus_notomega = new TH1F("2PiMass_minus_notomega", ";M_{#pi^{-}#pi^{0}} (GeV)", 600, 0., 1.5);
	dHist_3vs2PiMass_plus = new TH2F("3vs2PiMass_plus", ";M_{#pi^{+}#pi^{-}#pi^{0}};M_{#pi^{+}#pi^{0}}", 600, 0.3, 2.3, 600, 0, 1.5);
	dHist_3vs2PiMass_minus = new TH2F("3vs2PiMass_minus", ";M_{#pi^{+}#pi^{-}#pi^{0}};M_{#pi^{-}#pi^{0}}", 600, 0.3, 2.3, 600, 0, 1.5);

	dHist_2PiMass_t_proton_total = new TH2F("2PiMass_t_proton_total", ";-t_{p} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	dHist_2PiMass_t_omega_total = new TH2F("2PiMass_t_omega_total", ";-t_{#omega} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	dHist_2PiMass_t_Deltaplusplus_total = new TH2F("2PiMass_t_Deltaplusplus_total", ";-t_{#Delta^{++}} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	dHist_2PiMass_t_Delta0_total = new TH2F("2PiMass_t_Delta0_total", ";-t_{#Delta^{0}} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	for(int i = 0; i < 10; i++){
	  dHist_2PiMass_t_proton[i] = new TH2F(Form("2PiMass_t_proton_%i", i+1), ";-t_{p} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	  dHist_2PiMass_t_omega[i] = new TH2F(Form("2PiMass_t_omega_%i", i+1), ";-t_{#omega} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	  dHist_2PiMass_t_Deltaplusplus[i] = new TH2F(Form("2PiMass_t_Deltaplusplus_%i", i+1), ";-t_{#Delta^{++}} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	  dHist_2PiMass_t_Delta0[i] = new TH2F(Form("2PiMass_t_Delta0_%i", i+1), ";-t_{#Delta^{0}} (GeV)^{2};M_{#pi^{+}#pi^{-}} (GeV)", 600, 0., 3., 600, 0., 1.5);
	}

	dHist_omega2piMass_t_omega = new TH2F("omega2piMass_t_omega", ";-t_{#omega} (GeV)^{2};M_{#omega#pi#pi} (GeV)", 600, 0., 3., 600, 1., 3.);
	dHist_omegapiplusMass_t_omega = new TH2F("omegapiplusMass_t_omega", ";-t_{#omega} (GeV)^{2};M_{#omega#pi^{+}} (GeV)", 600, 0., 3., 600, .3, 2.3);
	dHist_omegapiminusMass_t_omega = new TH2F("omegapiminusMass_t_omega", ";-t_{#omega} (GeV)^{2};M_{#omega#pi^{-}} (GeV)", 600, 0., 3., 600, .3, 2.3);
	dHist_omega2piMass_t_Delta = new TH2F("omega2piMass_t_Delta", ";-t_{#Delta^{++}} (GeV)^{2};M_{#omega#pi#pi} (GeV)", 600, 0., 3., 600, 1., 3.);
	dHist_omegapiminusMass_t_Delta = new TH2F("omegapiminusMass_t_Delta", ";-t_{#Delta^{++}} (GeV)^{2};M_{#omega#pi^{-}} (GeV)", 600, 0., 3., 600, .3, 2.3);
	dHist_omega2piMass_t_Delta0 = new TH2F("omega2piMass_t_Delta0", ";-t_{#Delta} (GeV)^{2};M_{#omega#pi#pi} (GeV)", 600, 0., 3., 600, 1., 3.);
	dHist_omegapiplusMass_t_Delta0 = new TH2F("omegapiplusMass_t_Delta0", ";-t_{#Delta^{0}} (GeV)^{2};M_{#omega#pi^{+}} (GeV)", 600, 0., 3., 600, .3, 2.3);

	
	//Decay Angles
	dHist_CosThetaPlus = new TH1F("CosThetaPlus", ";cos(#theta_{b_{1}^{+}})", 600, -1., 1.);
	dHist_CosThetaMinus = new TH1F("CosThetaMinus", ";cos(#theta_{b_{1}^{-}})", 600, -1., 1.);
	dHist_PhiPlus = new TH1F("PhiPlus", ";#phi_{b_{1}^{+}} (deg)", 600, -180., 180.);
	dHist_PhiMinus = new TH1F("PhiMinus", ";#phi_{b_{1}^{-}} (deg)", 600, -180., 180.);

	dHist_CosTheta_omega_fromPlus = new TH1F("CosTheta_omega_fromPlus", ";cos(#theta_{#omega})", 600, -1., 1.);
	dHist_Phi_omega_fromPlus = new TH1F("Phi_omega_fromPlus", ";#phi_{#omega} (deg)", 600, -180., 180.);
	dHist_CosTheta_omega_fromMinus = new TH1F("CosTheta_omega_fromMinus", ";cos(#theta_{#omega})", 600, -1., 1.);
	dHist_Phi_omega_fromMinus = new TH1F("Phi_omega_fromMinus", ";#phi_{#omega} (deg)", 600, -180., 180.);

	dHist_CosTheta_H = new TH1F("CosTheta_H", ";cos(#theta_{H})", 600, -1., 1.);
	dHist_Phi_H = new TH1F("Phi_H", ";#phi_{H} (deg)", 600, -180., 180.);

	dHist_CosThetaPlusVsOmega2PiMass = new TH2F("CosThetaPlusVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);cos(#theta_{b_{1}^{+}})", 600, 1., 3., 600, -1., 1.);
	dHist_CosThetaMinusVsOmega2PiMass = new TH2F("CosThetaMinusVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);cos(#theta_{b_{1}^{-}})", 600, 1., 3., 600, -1., 1.);
	dHist_PhiPlusVsOmega2PiMass = new TH2F("PhiPlusVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);#phi_{b_{1}^{+}} (deg)", 600, 1., 3., 600, -180., 180.);
	dHist_PhiMinusVsOmega2PiMass = new TH2F("PhiMinusVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);#phi_{b_{1}^{-}} (deg)", 600, 1., 3., 600, -180., 180.);

	dHist_CosTheta_omegaVsOmegaPiPlusMass = new TH2F("CosTheta_omegaVsOmegaPiPlusMass", ";M_{#omega#pi^{+}} (GeV);cos(#theta_{#omega})", 600, .3, 2.3, 600, -1., 1.);
	dHist_CosTheta_omegaVsOmegaPiMinusMass = new TH2F("CosTheta_omegaVsOmegaPiMinusMass", ";M_{#omega#pi^{-}} (GeV);cos(#theta_{#omega})", 600, .3, 2.3, 600, -1., 1.);
	dHist_Phi_omegaVsOmegaPiPlusMass = new TH2F("Phi_omegaVsOmegaPiPlusMass", ";M_{#omega#pi^{+}} (GeV);#phi_{#omega} (deg)", 600, .3, 2.3, 600, -180., 180.);
	dHist_Phi_omegaVsOmegaPiMinusMass = new TH2F("Phi_omegaVsOmegaPiMinusMass", ";M_{#omega#pi^{-}} (GeV);#phi_{#omega} (deg)", 600, .3, 2.3, 600, -180., 180.);

	for(int i = 0; i < 10; i++){
	  dHist_PhiVsCosTheta_b1_Plus[i] = new TH2F(Form("PhiVsCosTheta_b1_Plus_%i", i+1), ";#phi_{b_{1}^{+}} (deg);cos(#theta_{b_{1}^{+}})", 600, -180., 180., 600, -1., 1);
	  dHist_PhiVsCosTheta_b1_Minus[i] = new TH2F(Form("PhiVsCosTheta_b1_Minus_%i", i+1), ";#phi_{b_{1}^{-}} (deg);cos(#theta_{b_{1}^{-}})", 600, -180., 180., 600, -1., 1);
	}
	dHist_PhiVsCosTheta_b1_Plus_total = new TH2F("PhiVsCosTheta_b1_Plus_total", ";#phi_{b_{1}^{+}} (deg);cos(#theta_{b_{1}^{+}})", 600, -180., 180., 600, -1., 1);
	dHist_PhiVsCosTheta_b1_Minus_total = new TH2F("PhiVsCosTheta_b1_Minus_total", ";#phi_{b_{1}^{-}} (deg);cos(#theta_{b_{1}^{-}})", 600, -180., 180., 600, -1., 1);


	dHist_PhiVsCosTheta_omega_Plus = new TH2F("PhiVsCosTheta_omega_Plus", ";#phi_{#omega} (deg) (from b_{1}^{+});cos(#theta_{#omega}) (from b_{1}^{+})", 600, -180., 180., 600, -1., 1);
	dHist_PhiVsCosTheta_omega_Minus = new TH2F("PhiVsCosTheta_omega_Minus", ";#phi_{#omega} (deg) (from b_{1}^{-});cos(#theta_{#omega}) (from b_{1}^{-})", 600, -180., 180., 600, -1., 1);

	dHist_PhiVsCosTheta_H = new TH2F("PhiVsCosTheta_H", ";#phi_{H} (deg);cos(#theta_{H})", 600, -180., 180., 600, -1., 1);

	dHist_M_CosTheta_omega_Plus = new TH2F("M_CosTheta_omega_Plus", ";M_{#omega#pi^{+}} (GeV);cos(#theta_{#omega}) (from b_{1}^{+})", 600, 0.3, 2.3, 600, -1., 1.);
	dHist_M_CosTheta_omega_Minus = new TH2F("M_CosTheta_omega_Minus", ";M_{#omega#pi^{-}} (GeV);cos(#theta_{#omega}) (from b_{1}^{-})", 600, 0.3, 2.3, 600, -1., 1.);

	dHist_M_CosThetaH_plus = new TH2F("M_CosThetaH_plus", ";M_{#omega#pi^{+}} (GeV);cos(#theta_{H})", 600, 0.3, 2.3, 600, -1., 1.);
	dHist_M_CosThetaH_minus = new TH2F("M_CosThetaH_minus", ";M_{#omega#pi^{-}} (GeV);cos(#theta_{H})", 600, 0.3, 2.3, 600, -1., 1.);

	dHist_M_Phi_omega_Plus = new TH2F("M_Phi_omega_Plus", ";M_{#omega#pi^{+}} (GeV);#phi_{#omega} (deg) (from b_{1}^{+})", 600, 0.3, 2.3, 600, -180., 180.);
	dHist_M_Phi_omega_Minus = new TH2F("M_Phi_omega_Minus", ";M_{#omega#pi^{-}} (GeV);#phi_{#omega} (deg) (from b_{1}^{-})", 600, 0.3, 2.3, 600, -180., 180.);

	dHist_M_PhiH_plus = new TH2F("M_PhiH_plus", ";M_{#omega#pi^{+}} (GeV);#phi_{H} (deg)", 600, 0.3, 2.3, 600, -180., 180.);
	dHist_M_PhiH_minus = new TH2F("M_PhiH_minus", ";M_{#omega#pi^{-}} (GeV);#phi_{H} (deg)", 600, 0.3, 2.3, 600, -180., 180.);

	dHist_CosThetaPlusVs2PiMass = new TH2F("CosThetaPlusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);cos(#theta_{b_{1}^{+}})", 600, 0., 1.5, 600, -1., 1.);
	dHist_PhiPlusVs2PiMass = new TH2F("PhiPlusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);#phi_{b_{1}^{+}} (deg)", 600, 0., 1.5, 600, -180., 180.);
	dHist_CosThetaMinusVs2PiMass = new TH2F("CosThetaMinusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);cos(#theta_{b_{1}^{-}})", 600, 0., 1.5, 600, -1., 1.);
	dHist_PhiMinusVs2PiMass = new TH2F("PhiMinusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);#phi_{b_{1}^{-}} (deg)", 600, 0., 1.5, 600, -180., 180.);

	dHist_CosTheta_omegaPlusVs2PiMass = new TH2F("CosTheta_omegaPlusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);cos(#theta_{#omega}) (from b_{1}^{+})", 600, 0., 1.5, 600, -1., 1.);
	dHist_Phi_omegaPlusVs2PiMass = new TH2F("Phi_omegaPlusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);#phi_{#omega} (deg) (from b_{1}^{+})", 600, 0., 1.5, 600, -180., 180.);
	dHist_CosTheta_omegaMinusVs2PiMass = new TH2F("CosTheta_omegaMinusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);cos(#theta_{#omega}) (from b_{1}^{-})", 600, 0., 1.5, 600, -1., 1.);
	dHist_Phi_omegaMinusVs2PiMass = new TH2F("Phi_omegaMinusVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);#phi_{#omega} (deg) (from b_{1}^{-})", 600, 0., 1.5, 600, -180., 180.);

	dHist_CosThetaPlusVsProtonPiMinusMass = new TH2F("CosThetaPlusVsProtonPiMinusMass", ";M_{p#pi^{-}} (GeV);cos(#theta_{b_{1}^{+}})", 600, 0.3, 3., 600, -1., 1.);
	dHist_PhiPlusVsProtonPiMinusMass = new TH2F("PhiPlusVsProtonPiMinusMass", ";M_{p#pi^{-}} (GeV);#phi_{b_{1}^{+}}", 600, 0.3, 3., 600, -180., 180.);
	dHist_CosThetaMinusVsProtonPiPlusMass = new TH2F("CosThetaMinusVsProtonPiPlusMass", ";M_{p#pi^{+}} (GeV);cos(#theta_{b_{1}^{-}})", 600, 0.3, 3., 600, -1., 1.);
	dHist_PhiMinusVsProtonPiPlusMass = new TH2F("PhiMinusVsProtonPiPlusMass", ";M_{p#pi^{+}} (GeV);#phi_{b_{1}^{-}}", 600, 0.3, 3., 600, -180., 180.);
	
	dHist_CosThetaOmega_OR = new TH1F("CosThetaOmega_OR", ";cos(#theta_{#omega})", 600, -1., 1.);
	dHist_PhiOmega_OR = new TH1F("PhiOmega_OR", ";#phi_{#omega} (deg)", 600, -180., 180.);
	dHist_CosThetaRho_OR = new TH1F("CosThetaRho_OR", ";cos(#theta_{#rho})", 600, -1., 1.);
	dHist_PhiRho_OR = new TH1F("PhiRho_OR", ";#phi_{#rho} (deg)", 600, -180., 180.);

	dHist_CosThetaOmega_ORVsOmega2PiMass = new TH2F("CosThetaOmega_ORVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);cos(#theta_{#omega})", 600, 1., 3., 600, -1., 1.);
	dHist_PhiOmega_ORVsOmega2PiMass = new TH2F("PhiOmega_ORVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);#phi_{#omega} (deg)", 600, 1., 3., 600, -180., 180.);
	dHist_CosThetaRho_ORVsOmega2PiMass = new TH2F("CosThetaRho_ORVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);cos(#theta_{#rho})", 600, 1., 3., 600, -1., 1.);
	dHist_PhiRho_ORVsOmega2PiMass = new TH2F("PhiRho_ORVsOmega2PiMass", ";M_{#omega#pi^{+}#pi^{-}} (GeV);#phi_{#rho} (deg)", 600, 1., 3., 600, -180., 180.);

	dHist_CosThetaOmega_ORVs3PiMass = new TH2F("CosThetaOmega_ORVs3PiMass", ";M_{#pi^{+}#pi^{-}#pi{0}} (GeV);cos(#theta_{#omega})", 600, 0.3, 2.3, 600, -1., 1.);
	dHist_PhiOmega_ORVs3PiMass = new TH2F("PhiOmega_ORVs3PiMass", ";M_{#pi^{+}#pi^{-}#pi{0}} (GeV);#phi_{#omega} (deg)", 600, 0.3, 2.3, 600, -180., 180.);
	//dHist_CosThetaOmega_ORVs3PiMass_omegacut = new TH2F("CosThetaOmega_ORVs3PiMass_omegacut", ";M_{#omega} (GeV);cos(#theta_{#omega})", 600, 0.3, 2.3, 600, -1., 1.);
	//dHist_PhiOmega_ORVs3PiMass_omegacut = new TH2F("PhiOmega_ORVs3PiMass_omegacut", ";M_{#omega} (GeV);#phi_{#omega} (deg)", 600, 0.3, 2.3, 600, -180., 180.);

	dHist_CosThetaRho_ORVs2PiMass = new TH2F("CosThetaRho_ORVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);cos(#theta_{#rho})", 600, 0., 1.5, 600, -1., 1.);
	dHist_PhiRho_ORVs2PiMass = new TH2F("PhiRho_ORVs2PiMass", ";M_{#pi^{+}#pi^{-}} (GeV);#phi_{#rho} (deg)", 600, 0., 1.5, 600, -180., 180.);
	dHist_CosThetaRho_ORVs2PiMass_omegacut = new TH2F("CosThetaRho_ORVs2PiMass_omegacut", ";M_{#pi^{+}#pi^{-}} (GeV);cos(#theta_{#rho})", 600, 0., 1.5, 600, -1., 1.);
	dHist_PhiRho_ORVs2PiMass_omegacut = new TH2F("PhiRho_ORVs2PiMass_omegacut", ";M_{#pi^{+}#pi^{-}} (GeV);#phi_{#rho} (deg)", 600, 0., 1.5, 600, -180., 180.);

	dHist_CosThetaOmegaVsPhiOmega_OR = new TH2F("CosThetaOmegaVsPhiOmega_OR", ";#phi_{#omega} (deg);cos(#theta_{#omega})", 600, -180., 180., 600, -1., 1.);
	dHist_CosThetaRhoVsPhiRho_OR = new TH2F("CosThetaRhoVsPhiRho_OR", ";#phi_{#rho} (deg);cos(#theta_{#rho})", 600, -180., 180., 600, -1., 1.);

	for(int i = 0; i < 10; i++){
	  dHist_Dalitz_OmegaPi[i] = new TH2F(Form("Dalitz_OmegaPi_%i", i+1), ";M^{2}_{#omega#pi^{+}} (GeV^{2});M^{2}_{#omega#pi^{-}} (GeV^{2})", 600, 0.5, 5.5, 600, 0.5, 5.5);
	}

	dHist_5PiMass_DeltaPlusPlusPeak = new TH1F("h5PiMass_DeltaPlusPlusPeak", ";M_{5#pi} (GeV)", 600, 1., 3.);
	dHist_Omega2PiMass_DeltaPlusPlusPeak = new TH1F("hOmega2PiMass_DeltaPlusPlusPeak", ";M_{#omega#pi^{+}#pi^{-}} (GeV)", 600, 1., 3.);
	dHist_OmegaPiMinusMass_DeltaPlusPlusPeak = new TH1F("hOmegaPiMinusMass_DeltaPlusPlusPeak", ";M_{#omega#pi^{-}}, (GeV)", 600, 0.3, 2.3);
	dHist_ProtonPiPlusMass_DeltaPlusPlusPeak = new TH1F("hProtonPiPlusMass_DeltaPlusPlusPeak", ";M_{p#pi^{+}} (GeV)", 600, .3, 3.);
	dHist_t_DeltaPlusPlus_DeltaPlusPlusPeak = new TH1F("ht_DeltaPlusPlus_DeltaPlusPlusPeak", ";-t_{#Delta^{++}} (GeV^{2})", 600, 0., 3.);
	dHist_t_proton_DeltaPlusPlusPeak = new TH1F("ht_proton_DeltaPlusPlusPeak", ";-t_{p} (GeV^{2})", 600, 0., 3.);
	dHist_t_omega_DeltaPlusPlusPeak = new TH1F("ht_omega_DeltaPlusPlusPeak", ";-t_{#omega} (GeV^{2})", 600, 0., 3.);
	dHist_CosTheta_DeltaPlusPlusPeak = new TH1F("hCosTheta_DeltaPlusPlusPeak", ";cos(#theta)", 600, -1., 1.);
	dHist_Phi_DeltaPlusPlusPeak = new TH1F("hPhi_DeltaPlusPlusPeak", ";#phi (deg)", 600, -180., 180.);
	dHist_CosThetaH_DeltaPlusPlusPeak = new TH1F("hCosThetaH_DeltaPlusPlusPeak", ";cos(#theta_{H})", 600, -1., 1.);
	dHist_PhiH_DeltaPlusPlusPeak = new TH1F("hPhiH_DeltaPlusPlusPeak", ";#phi_{H} (deg)", 600, -180., 180.);
	dHist_CosThetaVsPhi_DeltaPlusPlusPeak = new TH2F("hCosThetaVsPhi_DeltaPlusPlusPeak", ";#phi (deg);cos(#theta)", 600, -180., 180., 600, -1., 1.);
	dHist_CosThetaHVsPhiH_DeltaPlusPlusPeak = new TH2F("hCosThetaVsPhiH_DeltaPlusPlusPeak", ";#phi_{H} (deg);cos(#theta)", 600, -180., 180., 600, -1., 1.);

	dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak = new TH2F("hPhi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak", ";E_{#gamma} (GeV);#Phi_{lab,p#pi^{+}} (deg)", 600, 0., 12., 600, -180., 180.);
	dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak_high_t = new TH2F("hPhi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak_high_t", ";E_{#gamma} (GeV);#Phi_{lab,p#pi^{+}} (deg)", 600, 0., 12., 600, -180., 180.);

	dHist_CosTheta_M_omega2pi_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_omega2pi_DeltaPlusPlusPeak", ";M_{#omega#pi#pi} (GeV);cos(#theta)", 600, 1., 3., 600, -1., 1.);
	dHist_CosTheta_M_omegapiminus_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_omegapiminus_DeltaPlusPlusPeak", ";M_{#omega#pi^{-}} (GeV);cos(#theta)", 600, .3, 2.3, 600, -1., 1.);
	dHist_CosTheta_M_omegapiplus_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_omegapiplus_DeltaPlusPlusPeak", ";M_{#omega#pi^{+}} (GeV);cos(#theta)", 600, .3, 2.3, 600, -1., 1.);
	dHist_CosTheta_M_protonpiminus_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_protonpiminus_DeltaPlusPlusPeak", ";M_{p#pi^{-}} (GeV);cos(#theta)", 600, .3, 3., 600, -1., 1.);
	dHist_CosTheta_M_protonpiplus_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_protonpiplus_DeltaPlusPlusPeak", ";M_{p#pi^{+}} (GeV);cos(#theta)", 600, .3, 3., 600, -1., 1.);
	dHist_CosTheta_M_omega_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_omega_DeltaPlusPlusPeak", ";M_{#omega} (GeV);cos(#theta)", 600, .3, 2.3, 600, -1., 1.);
	dHist_CosTheta_M_proton2pi_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_proton2pi_DeltaPlusPlusPeak", ";M_{p#pi^{+}#pi^{-}} (GeV);cos(#theta)", 600, 1., 4., 600, -1., 1.);
	dHist_CosTheta_M_2pi_DeltaPlusPlusPeak = new TH2F("hCosTheta_M_2pi_DeltaPlusPlusPeak", ";M_{#pi^{+}#pi^{-}} (GeV);cos(#theta)", 600, 0., 1.5, 600, -1., 1.);

	dHist_M_proton2pi_protonpiplus_DeltaPlusPlusPeak = new TH2F("hM_proton2pi_protonpiplus_DeltaPlusPlusPeak", ";M_{p#pi^{+}} (GeV);M_{p#pi^{+}#pi^{-}} (GeV)", 600, .3, 2.3, 600, 1., 4.);
	dHist_M_proton2pi_protonpiminus_DeltaPlusPlusPeak = new TH2F("hM_proton2pi_protonpiminus_DeltaPlusPlusPeak", ";M_{p#pi^{-}} (GeV);M_{p#pi^{+}#pi^{-}} (GeV)", 600, .3, 2.3, 600, 1., 4.);
	dHist_M_proton2pi_pipluspiminus_DeltaPlusPlusPeak = new TH2F("hM_proton2pi_pipluspiminus_DeltaPlusPlusPeak", ";M_{#pi^{+}#pi^{-}} (GeV);M_{p#pi^{+}#pi^{-}} (GeV)", 600, 0., 1.5, 600, 1., 4.);

	dHist_Omega_Dalitz_xy_DeltaPlusPlusPeak = new TH2F("hOmega_Dalitz_xy_DeltaPlusPlusPeak", ";x_{Dalitz};y_{Dalitz}", 600, -2., 2., 600, -2., 2.);



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

	// fill generated from Thrown_Tree
	if(dOption.Contains("thrown") && dTreeInterface->Get_Branch("NumCombos") == NULL) {
                dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", 1.0);
                TLorentzVector locBeamP4 = dThrownBeam->Get_P4();
                FillAmpTools_FlatTree(locBeamP4, locFinalStateThrownP4);
                Fill_FlatTree();

                return kTRUE;
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
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_5PiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_5PiMass_Measured;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_3PiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_3PiMass_Measured;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Omega2PiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_3vs2PiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaPiPlusMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaPiMinusMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_ProtonPiPlusMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_ProtonPiPlusVsOmegaPiMinus;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Man_t;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_t_prime;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_2PiMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_2PiMass_charged;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_ProtonPiMinusMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_ProtonPiMinusVsOmegaPiPlus;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_3vs4PiMass_plus;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_3vs4PiMass_minus;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_DecayAngles;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_DecayAngles_OR;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Dalitz;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_DeltaPlusPlusPeak;

	/************************************************* FLAG WHICH CUTS TO USE *****************************************/
	bool cut_t_prime = false;
	bool veto_Deltaplusplus = false;
	bool veto_Delta0 = false;

	/************************************************* PARSE THROWN TOPOLOGY ***************************************/
	TString locThrownTopology = Get_ThrownTopologyString();

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
		TLorentzVector locDecayingPi0P4 = dPhoton1Wrapper->Get_P4() + dPhoton2Wrapper->Get_P4();//dDecayingPi0Wrapper->Get_P4();
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

		TLorentzVector loc5PiP4 = locPiPlus1P4 + locPiMinus1P4 + locPiPlus2P4 + locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc5PiP4_Measured = locPiPlus1P4_Measured + locPiMinus1P4_Measured + locPiPlus2P4_Measured + locPiMinus2P4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

		TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;

		//3Pi Momenta:
		TLorentzVector loc3PiP4_11 = locPiPlus1P4 + locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc3PiP4_12 = locPiPlus1P4 + locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc3PiP4_21 = locPiPlus2P4 + locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc3PiP4_22 = locPiPlus2P4 + locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;
		//Corresponding 2Pi Momenta (not from 3Pi combo):
		TLorentzVector loc2PiP4_11 = locPiPlus2P4 + locPiMinus2P4;
		TLorentzVector loc2PiP4_12 = locPiPlus2P4 + locPiMinus1P4;
		TLorentzVector loc2PiP4_21 = locPiPlus1P4 + locPiMinus2P4;
		TLorentzVector loc2PiP4_22 = locPiPlus1P4 + locPiMinus1P4;
		//Proton + "bachelor" 2pi momenta:
		TLorentzVector locProton2PiP4_11 = locProtonP4 + loc2PiP4_11;
		TLorentzVector locProton2PiP4_12 = locProtonP4 + loc2PiP4_12;
		TLorentzVector locProton2PiP4_21 = locProtonP4 + loc2PiP4_21;
		TLorentzVector locProton2PiP4_22 = locProtonP4 + loc2PiP4_22;

		double locProton2PiMass11 = locProton2PiP4_11.M();
		double locProton2PiMass12 = locProton2PiP4_12.M();
		double locProton2PiMass21 = locProton2PiP4_21.M();
		double locProton2PiMass22 = locProton2PiP4_22.M();

		//"Charged" 2pi Momenta (using these to look for charged rhos)
		TLorentzVector loc2PiP4_p01 = locPiPlus1P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc2PiP4_p02 = locPiPlus2P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc2PiP4_m01 = locPiMinus1P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector loc2PiP4_m02 = locPiMinus2P4 + locPhoton1P4 + locPhoton2P4;

		//3Pi Momenta (Measured):
		TLorentzVector loc3PiP4_11_Measured = locPiPlus1P4_Measured + locPiMinus1P4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector loc3PiP4_12_Measured = locPiPlus1P4_Measured + locPiMinus2P4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector loc3PiP4_21_Measured = locPiPlus2P4_Measured + locPiMinus1P4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector loc3PiP4_22_Measured = locPiPlus2P4_Measured + locPiMinus2P4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

		//proton pi+ momenta and masses
		TLorentzVector locProtonPiPlusP4_1 = locProtonP4 + locPiPlus1P4;
		TLorentzVector locProtonPiPlusP4_2 = locProtonP4 + locPiPlus2P4;

		double locProtonPiPlusMass1 = locProtonPiPlusP4_1.M();
		double locProtonPiPlusMass2 = locProtonPiPlusP4_2.M();

		//proton pi+ lab production angle Phi
		double locProtonPiPlus1_Phi_lab = locProtonPiPlusP4_1.Phi() * 180. / TMath::Pi();
		double locProtonPiPlus2_Phi_lab = locProtonPiPlusP4_2.Phi() * 180. / TMath::Pi();

		//proton pi- momenta and masses
		TLorentzVector locProtonPiMinusP4_1 = locProtonP4 + locPiMinus1P4;
		TLorentzVector locProtonPiMinusP4_2 = locProtonP4 + locPiMinus2P4;

		double locProtonPiMinusMass1 = locProtonPiMinusP4_1.M();
		double locProtonPiMinusMass2 = locProtonPiMinusP4_2.M();

		//Dalitz parameters for the omega
		double locDalitz_s_11 = (loc3PiP4_11 - locPi0P4).M2();
		double locDalitz_t_11 = (loc3PiP4_11 - locPiMinus1P4).M2();
		double locDalitz_u_11 = (loc3PiP4_11 - locPiPlus1P4).M2();
		double locDalitz_sc_11 = (1/3.0) * (loc3PiP4_11.M2() + locPiPlus1P4.M2() + locPiMinus1P4.M2() + locPi0P4.M2());
		double locDalitz_d_11 = 2.0 * loc3PiP4_11.M() * (loc3PiP4_11.M() - locPiPlus1P4.M() - locPiMinus1P4.M() - locPi0P4.M());

		double locDalitz_x_11 = TMath::Sqrt(3.0) * (locDalitz_t_11 - locDalitz_u_11) / locDalitz_d_11;
		double locDalitz_y_11 = 3.0 * (locDalitz_sc_11 - locDalitz_s_11) / locDalitz_d_11;

		double locDalitz_s_12 = (loc3PiP4_12 - locPi0P4).M2();
		double locDalitz_t_12 = (loc3PiP4_12 - locPiMinus2P4).M2();
		double locDalitz_u_12 = (loc3PiP4_12 - locPiPlus1P4).M2();
		double locDalitz_sc_12 = (1/3.0) * (loc3PiP4_12.M2() + locPiPlus1P4.M2() + locPiMinus2P4.M2() + locPi0P4.M2());
		double locDalitz_d_12 = 2.0 * loc3PiP4_12.M() * (loc3PiP4_12.M() - locPiPlus1P4.M() - locPiMinus2P4.M() - locPi0P4.M());

		double locDalitz_x_12 = TMath::Sqrt(3.0) * (locDalitz_t_12 - locDalitz_u_12) / locDalitz_d_12;
		double locDalitz_y_12 = 3.0 * (locDalitz_sc_12 - locDalitz_s_12) / locDalitz_d_12;

		double locDalitz_s_21 = (loc3PiP4_21 - locPi0P4).M2();
		double locDalitz_t_21 = (loc3PiP4_21 - locPiMinus1P4).M2();
		double locDalitz_u_21 = (loc3PiP4_21 - locPiPlus2P4).M2();
		double locDalitz_sc_21 = (1/3.0) * (loc3PiP4_21.M2() + locPiPlus2P4.M2() + locPiMinus1P4.M2() + locPi0P4.M2());
		double locDalitz_d_21 = 2.0 * loc3PiP4_21.M() * (loc3PiP4_21.M() - locPiPlus2P4.M() - locPiMinus1P4.M() - locPi0P4.M());

		double locDalitz_x_21 = TMath::Sqrt(3.0) * (locDalitz_t_21 - locDalitz_u_21) / locDalitz_d_21;
		double locDalitz_y_21 = 3.0 * (locDalitz_sc_21 - locDalitz_s_21) / locDalitz_d_21;

		double locDalitz_s_22 = (loc3PiP4_22 - locPi0P4).M2();
		double locDalitz_t_22 = (loc3PiP4_22 - locPiMinus2P4).M2();
		double locDalitz_u_22 = (loc3PiP4_22 - locPiPlus2P4).M2();
		double locDalitz_sc_22 = (1/3.0) * (loc3PiP4_22.M2() + locPiPlus2P4.M2() + locPiMinus2P4.M2() + locPi0P4.M2());
		double locDalitz_d_22 = 2.0 * loc3PiP4_22.M() * (loc3PiP4_22.M() - locPiPlus2P4.M() - locPiMinus2P4.M() - locPi0P4.M());

		double locDalitz_x_22 = TMath::Sqrt(3.0) * (locDalitz_t_22 - locDalitz_u_22) / locDalitz_d_22;
		double locDalitz_y_22 = 3.0 * (locDalitz_sc_22 - locDalitz_s_22) / locDalitz_d_22;





		/******************************************* DEFINE DECAY ANGLES ****************************************************/

		//Boost to gamma-p rest frame:
		TLorentzVector locGammapP4 = locBeamP4 + dTargetP4;
		TVector3 locGammapBoost = locGammapP4.BoostVector();  //create boost vector from gamma-p r.f. to lab frame - use negative of this to boost to gamma-p rest frame
		TLorentzVector locBeamP4_gpRest = locBeamP4;
		locBeamP4_gpRest.Boost(-1.0*locGammapBoost); //boost beam to gamma-p r.f.
		TLorentzVector loc5PiP4_gpRest = loc5PiP4;
		loc5PiP4_gpRest.Boost(-1.0*locGammapBoost);  //boost 5pi system to gamma-p r.f.
		//create the 4 different possible b1 combos:
		TLorentzVector loc4PiPlus1P4_gpRest = loc3PiP4_11 + locPiPlus2P4;
		TLorentzVector loc4PiPlus2P4_gpRest = loc3PiP4_12 + locPiPlus2P4;
		TLorentzVector loc4PiMinus1P4_gpRest = loc3PiP4_11 + locPiMinus2P4;
		TLorentzVector loc4PiMinus2P4_gpRest = loc3PiP4_21 + locPiMinus2P4;
		//Boost b1 combos to gamma-p r.f.
		loc4PiPlus1P4_gpRest.Boost(-1.0*locGammapBoost);
		loc4PiPlus2P4_gpRest.Boost(-1.0*locGammapBoost);
		loc4PiMinus1P4_gpRest.Boost(-1.0*locGammapBoost);
		loc4PiMinus2P4_gpRest.Boost(-1.0*locGammapBoost);

		//Define unit 3-vectors:
		TVector3 lock = locBeamP4_gpRest.Vect().Unit();
		TVector3 locz = loc5PiP4_gpRest.Vect().Unit();
		TVector3 locy = lock.Cross(locz).Unit();
		TVector3 locx = locy.Cross(locz);

		//Boost b1 from gamma-p to b1pi r.f.
		TVector3 loc5PiBoost = loc5PiP4_gpRest.BoostVector();  //create boost vector from b1pi r.f. to gamma-p r.f.
		TLorentzVector loc4PiPlus1P4_b1piRest = loc4PiPlus1P4_gpRest;
		loc4PiPlus1P4_b1piRest.Boost(-1.0*loc5PiBoost);
		TLorentzVector loc4PiPlus2P4_b1piRest = loc4PiPlus2P4_gpRest;
		loc4PiPlus2P4_b1piRest.Boost(-1.0*loc5PiBoost);
		TLorentzVector loc4PiMinus1P4_b1piRest = loc4PiMinus1P4_gpRest;
		loc4PiMinus1P4_b1piRest.Boost(-1.0*loc5PiBoost);
		TLorentzVector loc4PiMinus2P4_b1piRest = loc4PiMinus2P4_gpRest;
		loc4PiMinus2P4_b1piRest.Boost(-1.0*loc5PiBoost);
		TVector3 loc4PiPlus1P3 = loc4PiPlus1P4_b1piRest.Vect();
		TVector3 loc4PiPlus2P3 = loc4PiPlus2P4_b1piRest.Vect();
		TVector3 loc4PiMinus1P3 = loc4PiMinus1P4_b1piRest.Vect();
		TVector3 loc4PiMinus2P3 = loc4PiMinus2P4_b1piRest.Vect();

		double locCosTheta1Plus = loc4PiPlus1P3.Dot(locz) / loc4PiPlus1P3.Mag();
		double locCosTheta2Plus = loc4PiPlus2P3.Dot(locz) / loc4PiPlus2P3.Mag();
		double locCosTheta1Minus = loc4PiMinus1P3.Dot(locz) / loc4PiMinus1P3.Mag();
		double locCosTheta2Minus = loc4PiMinus2P3.Dot(locz) / loc4PiMinus2P3.Mag();

		double locPhi1Plus = 180. * TMath::ATan2(loc4PiPlus1P3.Dot(locy), loc4PiPlus1P3.Dot(locx)) / TMath::Pi();
		double locPhi2Plus = 180. * TMath::ATan2(loc4PiPlus2P3.Dot(locy), loc4PiPlus2P3.Dot(locx)) / TMath::Pi();
		double locPhi1Minus = 180. * TMath::ATan2(loc4PiMinus1P3.Dot(locy), loc4PiMinus1P3.Dot(locx)) / TMath::Pi();
		double locPhi2Minus = 180. * TMath::ATan2(loc4PiMinus2P3.Dot(locy), loc4PiMinus2P3.Dot(locx)) / TMath::Pi();

		//Define omega unit vectors
		TVector3 locz_omega_p1 = loc4PiPlus1P3.Unit(); //b1+_1 direction in b1pi r.f.
		TVector3 locy_omega_p1 = locz.Cross(locz_omega_p1).Unit();
		TVector3 locx_omega_p1 = locy_omega_p1.Cross(locz_omega_p1);
		TVector3 locz_omega_p2 = loc4PiPlus2P3.Unit(); //b1+_2 direction in b1pi r.f.
		TVector3 locy_omega_p2 = locz.Cross(locz_omega_p2).Unit();
		TVector3 locx_omega_p2 = locy_omega_p2.Cross(locz_omega_p2);
		TVector3 locz_omega_m1 = loc4PiMinus1P3.Unit(); //b1-_1 direction in b1pi r.f.
		TVector3 locy_omega_m1 = locz.Cross(locz_omega_m1).Unit();
		TVector3 locx_omega_m1 = locy_omega_m1.Cross(locz_omega_m1);
		TVector3 locz_omega_m2 = loc4PiMinus2P3.Unit(); //b1-_2 direction in b1pi r.f.
		TVector3 locy_omega_m2 = locz.Cross(locz_omega_m2).Unit();
		TVector3 locx_omega_m2 = locy_omega_m2.Cross(locz_omega_m2);

		//Boost omegas from lab to omegapi r.f.
		//1) b1+1 -> omega11 pi+2
		//2) b1+1 -> omega21 pi+1
		//need to get omega11 and omega21 from lab to b1+_1 r.f. 
		TVector3 locb1Plus1Boost = loc4PiPlus1P4_b1piRest.BoostVector(); //Create boost vector from b1 r.f. to b1pi r.f
		TLorentzVector loc3PiP4_11_b1p1Rest = loc3PiP4_11;
		loc3PiP4_11_b1p1Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_11_b1p1Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_11_b1p1Rest.Boost(-1.0*locb1Plus1Boost);
		TVector3 loc3PiP3_11_b1p1Rest = loc3PiP4_11_b1p1Rest.Vect();
		TLorentzVector loc3PiP4_21_b1p1Rest = loc3PiP4_21;
		loc3PiP4_21_b1p1Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_21_b1p1Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_21_b1p1Rest.Boost(-1.0*locb1Plus1Boost);
		TVector3 loc3PiP3_21_b1p1Rest = loc3PiP4_21_b1p1Rest.Vect();

		//3) b1+2 -> omega12 pi+2
		//4) b1+2 -> omega22 pi+1
		TVector3 locb1Plus2Boost = loc4PiPlus2P4_b1piRest.BoostVector();
		TLorentzVector loc3PiP4_12_b1p2Rest = loc3PiP4_12;
		loc3PiP4_12_b1p2Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_12_b1p2Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_12_b1p2Rest.Boost(-1.0*locb1Plus2Boost);
		TVector3 loc3PiP3_12_b1p2Rest = loc3PiP4_12_b1p2Rest.Vect();
		TLorentzVector loc3PiP4_22_b1p2Rest = loc3PiP4_22;
		loc3PiP4_22_b1p2Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_22_b1p2Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_22_b1p2Rest.Boost(-1.0*locb1Plus2Boost);
		TVector3 loc3PiP3_22_b1p2Rest = loc3PiP4_22_b1p2Rest.Vect();

		//5) b1-1 -> omega11 pi-2
		//6) b1-1 -> omega12 pi-1
		TVector3 locb1Minus1Boost = loc4PiMinus1P4_b1piRest.BoostVector();
		TLorentzVector loc3PiP4_11_b1m1Rest = loc3PiP4_11;
		loc3PiP4_11_b1m1Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_11_b1m1Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_11_b1m1Rest.Boost(-1.0*locb1Minus1Boost);
		TVector3 loc3PiP3_11_b1m1Rest = loc3PiP4_11_b1m1Rest.Vect();
		TLorentzVector loc3PiP4_12_b1m1Rest = loc3PiP4_12;
		loc3PiP4_12_b1m1Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_12_b1m1Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_12_b1m1Rest.Boost(-1.0*locb1Minus1Boost);
		TVector3 loc3PiP3_12_b1m1Rest = loc3PiP4_12_b1m1Rest.Vect();
		//7) b1-2 -> omega21 pi-2
		//8) b1-2 -> omega22 pi-1
		TVector3 locb1Minus2Boost = loc4PiMinus2P4_b1piRest.BoostVector();
		TLorentzVector loc3PiP4_21_b1m2Rest = loc3PiP4_21;
		loc3PiP4_21_b1m2Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_21_b1m2Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_21_b1m2Rest.Boost(-1.0*locb1Minus2Boost);
		TVector3 loc3PiP3_21_b1m2Rest = loc3PiP4_21_b1m2Rest.Vect();
		TLorentzVector loc3PiP4_22_b1m2Rest = loc3PiP4_22;
		loc3PiP4_22_b1m2Rest.Boost(-1.0*locGammapBoost);
		loc3PiP4_22_b1m2Rest.Boost(-1.0*loc5PiBoost);
		loc3PiP4_22_b1m2Rest.Boost(-1.0*locb1Minus2Boost);
		TVector3 loc3PiP3_22_b1m2Rest = loc3PiP4_22_b1m2Rest.Vect();

		//Define omega decay angles:
		double locCosTheta_omega11_b1p1 = loc3PiP3_11_b1p1Rest.Dot(locz_omega_p1) / loc3PiP3_11_b1p1Rest.Mag();
		double locCosTheta_omega21_b1p1 = loc3PiP3_21_b1p1Rest.Dot(locz_omega_p1) / loc3PiP3_21_b1p1Rest.Mag();
		double locCosTheta_omega12_b1p2 = loc3PiP3_12_b1p2Rest.Dot(locz_omega_p2) / loc3PiP3_12_b1p2Rest.Mag();
		double locCosTheta_omega22_b1p2 = loc3PiP3_22_b1p2Rest.Dot(locz_omega_p2) / loc3PiP3_22_b1p2Rest.Mag();
		double locCosTheta_omega11_b1m1 = loc3PiP3_11_b1m1Rest.Dot(locz_omega_m1) / loc3PiP3_11_b1m1Rest.Mag();
		double locCosTheta_omega12_b1m1 = loc3PiP3_12_b1m1Rest.Dot(locz_omega_m1) / loc3PiP3_12_b1m1Rest.Mag();
		double locCosTheta_omega21_b1m2 = loc3PiP3_21_b1m2Rest.Dot(locz_omega_m2) / loc3PiP3_21_b1m2Rest.Mag();
		double locCosTheta_omega22_b1m2 = loc3PiP3_22_b1m2Rest.Dot(locz_omega_m2) / loc3PiP3_22_b1m2Rest.Mag();

		double locPhi_omega11_b1p1 = 180. * TMath::ATan2(loc3PiP3_11_b1p1Rest.Dot(locy_omega_p1), loc3PiP3_11_b1p1Rest.Dot(locx_omega_p1)) / TMath::Pi();
		double locPhi_omega21_b1p1 = 180. * TMath::ATan2(loc3PiP3_21_b1p1Rest.Dot(locy_omega_p1), loc3PiP3_21_b1p1Rest.Dot(locx_omega_p1)) / TMath::Pi();
		double locPhi_omega12_b1p2 = 180. * TMath::ATan2(loc3PiP3_12_b1p2Rest.Dot(locy_omega_p2), loc3PiP3_12_b1p2Rest.Dot(locx_omega_p2)) / TMath::Pi();
		double locPhi_omega22_b1p2 = 180. * TMath::ATan2(loc3PiP3_22_b1p2Rest.Dot(locy_omega_p2), loc3PiP3_22_b1p2Rest.Dot(locx_omega_p2)) / TMath::Pi();
		double locPhi_omega11_b1m1 = 180. * TMath::ATan2(loc3PiP3_11_b1m1Rest.Dot(locy_omega_m1), loc3PiP3_11_b1m1Rest.Dot(locx_omega_m1)) / TMath::Pi();
		double locPhi_omega12_b1m1 = 180. * TMath::ATan2(loc3PiP3_12_b1m1Rest.Dot(locy_omega_m1), loc3PiP3_12_b1m1Rest.Dot(locx_omega_m1)) / TMath::Pi();
		double locPhi_omega21_b1m2 = 180. * TMath::ATan2(loc3PiP3_21_b1m2Rest.Dot(locy_omega_m2), loc3PiP3_21_b1m2Rest.Dot(locx_omega_m2)) / TMath::Pi();
		double locPhi_omega22_b1m2 = 180. * TMath::ATan2(loc3PiP3_22_b1m2Rest.Dot(locy_omega_m2), loc3PiP3_22_b1m2Rest.Dot(locx_omega_m2)) / TMath::Pi();


		//Define unit vectors in omega direction in omegapi r.f.
		TVector3 locz_H_11_b1p1 = loc3PiP3_11_b1p1Rest.Unit();
		TVector3 locy_H_11_b1p1 = locz_omega_p1.Cross(locz_H_11_b1p1).Unit();
		TVector3 locx_H_11_b1p1 = locy_H_11_b1p1.Cross(locz_H_11_b1p1);
		TVector3 locz_H_21_b1p1 = loc3PiP3_21_b1p1Rest.Unit();
		TVector3 locy_H_21_b1p1 = locz_omega_p1.Cross(locz_H_21_b1p1).Unit();
		TVector3 locx_H_21_b1p1 = locy_H_21_b1p1.Cross(locz_H_21_b1p1);

		TVector3 locz_H_12_b1p2 = loc3PiP3_12_b1p2Rest.Unit();
		TVector3 locy_H_12_b1p2 = locz_omega_p2.Cross(locz_H_12_b1p2).Unit();
		TVector3 locx_H_12_b1p2 = locy_H_12_b1p2.Cross(locz_H_12_b1p2);
		TVector3 locz_H_22_b1p2 = loc3PiP3_22_b1p2Rest.Unit();
		TVector3 locy_H_22_b1p2 = locz_omega_p2.Cross(locz_H_22_b1p2).Unit();
		TVector3 locx_H_22_b1p2 = locy_H_22_b1p2.Cross(locz_H_22_b1p2);

		TVector3 locz_H_11_b1m1 = loc3PiP3_11_b1m1Rest.Unit();
		TVector3 locy_H_11_b1m1 = locz_omega_m1.Cross(locz_H_11_b1m1).Unit();
		TVector3 locx_H_11_b1m1 = locy_H_11_b1m1.Cross(locz_H_11_b1m1);
		TVector3 locz_H_12_b1m1 = loc3PiP3_12_b1m1Rest.Unit();
		TVector3 locy_H_12_b1m1 = locz_omega_m1.Cross(locz_H_12_b1m1).Unit();
		TVector3 locx_H_12_b1m1 = locy_H_12_b1m1.Cross(locz_H_12_b1m1);

		TVector3 locz_H_21_b1m2 = loc3PiP3_21_b1m2Rest.Unit();
		TVector3 locy_H_21_b1m2 = locz_omega_m2.Cross(locz_H_21_b1m2).Unit();
		TVector3 locx_H_21_b1m2 = locy_H_21_b1m2.Cross(locz_H_21_b1m2);
		TVector3 locz_H_22_b1m2 = loc3PiP3_22_b1m2Rest.Unit();
		TVector3 locy_H_22_b1m2 = locz_omega_m2.Cross(locz_H_22_b1m2).Unit();
		TVector3 locx_H_22_b1m2 = locy_H_22_b1m2.Cross(locz_H_22_b1m2);

		//Boost charged pions from lab to omega rest frame:
		//Omega11_b1p1 r.f.
		TVector3 locOmega11_b1p1Boost = loc3PiP4_11_b1p1Rest.BoostVector();
		TLorentzVector locPiPlus1P4_o11_b1p1Rest = locPiPlus1P4;
		locPiPlus1P4_o11_b1p1Rest.Boost(-1.0*locGammapBoost);
		locPiPlus1P4_o11_b1p1Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus1P4_o11_b1p1Rest.Boost(-1.0*locb1Plus1Boost);
		locPiPlus1P4_o11_b1p1Rest.Boost(-1.0*locOmega11_b1p1Boost);
		TVector3 locPiPlus1P3_o11_b1p1Rest = locPiPlus1P4_o11_b1p1Rest.Vect();
		TLorentzVector locPiMinus1P4_o11_b1p1Rest = locPiMinus1P4;
		locPiMinus1P4_o11_b1p1Rest.Boost(-1.0*locGammapBoost);
		locPiMinus1P4_o11_b1p1Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus1P4_o11_b1p1Rest.Boost(-1.0*locb1Plus1Boost);
		locPiMinus1P4_o11_b1p1Rest.Boost(-1.0*locOmega11_b1p1Boost);
		TVector3 locPiMinus1P3_o11_b1p1Rest = locPiMinus1P4_o11_b1p1Rest.Vect();

		//Omega21_b1p1 r.f.
		TVector3 locOmega21_b1p1Boost = loc3PiP4_21_b1p1Rest.BoostVector();
		TLorentzVector locPiPlus2P4_o21_b1p1Rest = locPiPlus2P4;
		locPiPlus2P4_o21_b1p1Rest.Boost(-1.0*locGammapBoost);
		locPiPlus2P4_o21_b1p1Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus2P4_o21_b1p1Rest.Boost(-1.0*locb1Plus1Boost);
		locPiPlus2P4_o21_b1p1Rest.Boost(-1.0*locOmega21_b1p1Boost);
		TVector3 locPiPlus2P3_o21_b1p1Rest = locPiPlus2P4_o21_b1p1Rest.Vect();
		TLorentzVector locPiMinus1P4_o21_b1p1Rest = locPiMinus1P4;
		locPiMinus1P4_o21_b1p1Rest.Boost(-1.0*locGammapBoost);
		locPiMinus1P4_o21_b1p1Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus1P4_o21_b1p1Rest.Boost(-1.0*locb1Plus1Boost);
		locPiMinus1P4_o21_b1p1Rest.Boost(-1.0*locOmega21_b1p1Boost);
		TVector3 locPiMinus1P3_o21_b1p1Rest = locPiMinus1P4_o21_b1p1Rest.Vect();

		//Omega12_b1p2 r.f.
		TVector3 locOmega12_b1p2Boost = loc3PiP4_12_b1p2Rest.BoostVector();
		TLorentzVector locPiPlus1P4_o12_b1p2Rest = locPiPlus1P4;
		locPiPlus1P4_o12_b1p2Rest.Boost(-1.0*locGammapBoost);
		locPiPlus1P4_o12_b1p2Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus1P4_o12_b1p2Rest.Boost(-1.0*locb1Plus2Boost);
		locPiPlus1P4_o12_b1p2Rest.Boost(-1.0*locOmega12_b1p2Boost);
		TVector3 locPiPlus1P3_o12_b1p2Rest = locPiPlus1P4_o12_b1p2Rest.Vect();
		TLorentzVector locPiMinus2P4_o12_b1p2Rest = locPiMinus2P4;
		locPiMinus2P4_o12_b1p2Rest.Boost(-1.0*locGammapBoost);
		locPiMinus2P4_o12_b1p2Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus2P4_o12_b1p2Rest.Boost(-1.0*locb1Plus2Boost);
		locPiMinus2P4_o12_b1p2Rest.Boost(-1.0*locOmega12_b1p2Boost);
		TVector3 locPiMinus2P3_o12_b1p2Rest = locPiMinus2P4_o12_b1p2Rest.Vect();

		//Omega22_b1p2 r.f.
		TVector3 locOmega22_b1p2Boost = loc3PiP4_22_b1p2Rest.BoostVector();
		TLorentzVector locPiPlus2P4_o22_b1p2Rest = locPiPlus2P4;
		locPiPlus2P4_o22_b1p2Rest.Boost(-1.0*locGammapBoost);
		locPiPlus2P4_o22_b1p2Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus2P4_o22_b1p2Rest.Boost(-1.0*locb1Plus2Boost);
		locPiPlus2P4_o22_b1p2Rest.Boost(-1.0*locOmega22_b1p2Boost);
		TVector3 locPiPlus2P3_o22_b1p2Rest = locPiPlus2P4_o22_b1p2Rest.Vect();
		TLorentzVector locPiMinus2P4_o22_b1p2Rest = locPiMinus2P4;
		locPiMinus2P4_o22_b1p2Rest.Boost(-1.0*locGammapBoost);
		locPiMinus2P4_o22_b1p2Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus2P4_o22_b1p2Rest.Boost(-1.0*locb1Plus2Boost);
		locPiMinus2P4_o22_b1p2Rest.Boost(-1.0*locOmega22_b1p2Boost);
		TVector3 locPiMinus2P3_o22_b1p2Rest = locPiMinus2P4_o22_b1p2Rest.Vect();

		//Omega11_b1m1 r.f.
		TVector3 locOmega11_b1m1Boost = loc3PiP4_11_b1m1Rest.BoostVector();
		TLorentzVector locPiPlus1P4_o11_b1m1Rest = locPiPlus1P4;
		locPiPlus1P4_o11_b1m1Rest.Boost(-1.0*locGammapBoost);
		locPiPlus1P4_o11_b1m1Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus1P4_o11_b1m1Rest.Boost(-1.0*locb1Minus1Boost);
		locPiPlus1P4_o11_b1m1Rest.Boost(-1.0*locOmega11_b1m1Boost);
		TVector3 locPiPlus1P3_o11_b1m1Rest = locPiPlus1P4_o11_b1m1Rest.Vect();
		TLorentzVector locPiMinus1P4_o11_b1m1Rest = locPiMinus1P4;
		locPiMinus1P4_o11_b1m1Rest.Boost(-1.0*locGammapBoost);
		locPiMinus1P4_o11_b1m1Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus1P4_o11_b1m1Rest.Boost(-1.0*locb1Minus1Boost);
		locPiMinus1P4_o11_b1m1Rest.Boost(-1.0*locOmega11_b1m1Boost);
		TVector3 locPiMinus1P3_o11_b1m1Rest = locPiMinus1P4_o11_b1m1Rest.Vect();

		//Omega12_b1m1 r.f.
		TVector3 locOmega12_b1m1Boost = loc3PiP4_12_b1m1Rest.BoostVector();
		TLorentzVector locPiPlus1P4_o12_b1m1Rest = locPiPlus1P4;
		locPiPlus1P4_o12_b1m1Rest.Boost(-1.0*locGammapBoost);
		locPiPlus1P4_o12_b1m1Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus1P4_o12_b1m1Rest.Boost(-1.0*locb1Minus1Boost);
		locPiPlus1P4_o12_b1m1Rest.Boost(-1.0*locOmega12_b1m1Boost);
		TVector3 locPiPlus1P3_o12_b1m1Rest = locPiPlus1P4_o12_b1m1Rest.Vect();
		TLorentzVector locPiMinus2P4_o12_b1m1Rest = locPiMinus2P4;
		locPiMinus2P4_o12_b1m1Rest.Boost(-1.0*locGammapBoost);
		locPiMinus2P4_o12_b1m1Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus2P4_o12_b1m1Rest.Boost(-1.0*locb1Minus1Boost);
		locPiMinus2P4_o12_b1m1Rest.Boost(-1.0*locOmega12_b1m1Boost);
		TVector3 locPiMinus2P3_o12_b1m1Rest = locPiMinus2P4_o12_b1m1Rest.Vect();

		//Omega21_b1m2 r.f.
		TVector3 locOmega21_b1m2Boost = loc3PiP4_21_b1m2Rest.BoostVector();
		TLorentzVector locPiPlus2P4_o21_b1m2Rest = locPiPlus2P4;
		locPiPlus2P4_o21_b1m2Rest.Boost(-1.0*locGammapBoost);
		locPiPlus2P4_o21_b1m2Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus2P4_o21_b1m2Rest.Boost(-1.0*locb1Minus2Boost);
		locPiPlus2P4_o21_b1m2Rest.Boost(-1.0*locOmega21_b1m2Boost);
		TVector3 locPiPlus2P3_o21_b1m2Rest = locPiPlus2P4_o21_b1m2Rest.Vect();
		TLorentzVector locPiMinus1P4_o21_b1m2Rest = locPiMinus1P4;
		locPiMinus1P4_o21_b1m2Rest.Boost(-1.0*locGammapBoost);
		locPiMinus1P4_o21_b1m2Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus1P4_o21_b1m2Rest.Boost(-1.0*locb1Minus2Boost);
		locPiMinus1P4_o21_b1m2Rest.Boost(-1.0*locOmega21_b1m2Boost);
		TVector3 locPiMinus1P3_o21_b1m2Rest = locPiMinus1P4_o21_b1m2Rest.Vect();

		//Omega22_b1m2 r.f.
		TVector3 locOmega22_b1m2Boost = loc3PiP4_22_b1m2Rest.BoostVector();
		TLorentzVector locPiPlus2P4_o22_b1m2Rest = locPiPlus2P4;
		locPiPlus2P4_o22_b1m2Rest.Boost(-1.0*locGammapBoost);
		locPiPlus2P4_o22_b1m2Rest.Boost(-1.0*loc5PiBoost);
		locPiPlus2P4_o22_b1m2Rest.Boost(-1.0*locb1Minus2Boost);
		locPiPlus2P4_o22_b1m2Rest.Boost(-1.0*locOmega22_b1m2Boost);
		TVector3 locPiPlus2P3_o22_b1m2Rest = locPiPlus2P4_o22_b1m2Rest.Vect();
		TLorentzVector locPiMinus2P4_o22_b1m2Rest = locPiMinus2P4;
		locPiMinus2P4_o22_b1m2Rest.Boost(-1.0*locGammapBoost);
		locPiMinus2P4_o22_b1m2Rest.Boost(-1.0*loc5PiBoost);
		locPiMinus2P4_o22_b1m2Rest.Boost(-1.0*locb1Minus2Boost);
		locPiMinus2P4_o22_b1m2Rest.Boost(-1.0*locOmega22_b1m2Boost);
		TVector3 locPiMinus2P3_o22_b1m2Rest = locPiMinus2P4_o22_b1m2Rest.Vect();

		//Normal vectors to omega decay planes:
		TVector3 locNormal11_b1p1 = locPiPlus1P3_o11_b1p1Rest.Cross(locPiMinus1P3_o11_b1p1Rest).Unit();
		TVector3 locNormal21_b1p1 = locPiPlus2P3_o21_b1p1Rest.Cross(locPiMinus1P3_o21_b1p1Rest).Unit();
		TVector3 locNormal12_b1p2 = locPiPlus1P3_o12_b1p2Rest.Cross(locPiMinus2P3_o12_b1p2Rest).Unit();
		TVector3 locNormal22_b1p2 = locPiPlus2P3_o22_b1p2Rest.Cross(locPiMinus2P3_o22_b1p2Rest).Unit();
		TVector3 locNormal11_b1m1 = locPiPlus1P3_o11_b1m1Rest.Cross(locPiMinus1P3_o11_b1m1Rest).Unit();
		TVector3 locNormal12_b1m1 = locPiPlus1P3_o12_b1m1Rest.Cross(locPiMinus2P3_o12_b1m1Rest).Unit();
		TVector3 locNormal21_b1m2 = locPiPlus2P3_o21_b1m2Rest.Cross(locPiMinus1P3_o21_b1m2Rest).Unit();
		TVector3 locNormal22_b1m2 = locPiPlus2P3_o22_b1m2Rest.Cross(locPiMinus2P3_o22_b1m2Rest).Unit();

		//Get H decay angles:
		double locCosTheta_H_11_b1p1 = locNormal11_b1p1.Dot(locz_H_11_b1p1);
		double locCosTheta_H_21_b1p1 = locNormal21_b1p1.Dot(locz_H_21_b1p1);
		double locCosTheta_H_12_b1p2 = locNormal12_b1p2.Dot(locz_H_12_b1p2);
		double locCosTheta_H_22_b1p2 = locNormal22_b1p2.Dot(locz_H_22_b1p2);
		double locCosTheta_H_11_b1m1 = locNormal11_b1m1.Dot(locz_H_11_b1m1);
		double locCosTheta_H_12_b1m1 = locNormal12_b1m1.Dot(locz_H_12_b1m1);
		double locCosTheta_H_21_b1m2 = locNormal21_b1m2.Dot(locz_H_21_b1m2);
		double locCosTheta_H_22_b1m2 = locNormal22_b1m2.Dot(locz_H_22_b1m2);

		double locPhi_H_11_b1p1 = 180. * TMath::ATan2(locNormal11_b1p1.Dot(locy_H_11_b1p1), locNormal11_b1p1.Dot(locx_H_11_b1p1)) / TMath::Pi();
		double locPhi_H_21_b1p1 = 180. * TMath::ATan2(locNormal21_b1p1.Dot(locy_H_21_b1p1), locNormal21_b1p1.Dot(locx_H_21_b1p1)) / TMath::Pi();
		double locPhi_H_12_b1p2 = 180. * TMath::ATan2(locNormal12_b1p2.Dot(locy_H_12_b1p2), locNormal12_b1p2.Dot(locx_H_12_b1p2)) / TMath::Pi();
		double locPhi_H_22_b1p2 = 180. * TMath::ATan2(locNormal22_b1p2.Dot(locy_H_22_b1p2), locNormal22_b1p2.Dot(locx_H_22_b1p2)) / TMath::Pi();
		double locPhi_H_11_b1m1 = 180. * TMath::ATan2(locNormal11_b1m1.Dot(locy_H_11_b1m1), locNormal11_b1m1.Dot(locx_H_11_b1m1)) / TMath::Pi();
		double locPhi_H_12_b1m1 = 180. * TMath::ATan2(locNormal12_b1m1.Dot(locy_H_12_b1m1), locNormal12_b1m1.Dot(locx_H_12_b1m1)) / TMath::Pi();
		double locPhi_H_21_b1m2 = 180. * TMath::ATan2(locNormal21_b1m2.Dot(locy_H_21_b1m2), locNormal21_b1m2.Dot(locx_H_21_b1m2)) / TMath::Pi();
		double locPhi_H_22_b1m2 = 180. * TMath::ATan2(locNormal22_b1m2.Dot(locy_H_22_b1m2), locNormal22_b1m2.Dot(locx_H_22_b1m2)) / TMath::Pi();


		/***************************************** DECAY ANGLES (ASSUMING OMEGARHO EVENTS) ****************************************/
		//Boost vector from lab to gp r.f. is defined above (locGammapBoost)
		//Boost vector from gp to 5pi r.f. is defined above (loc5PiBoost)

		//Boost possible omega combos to 5pi r.f.
		TLorentzVector loc3Pi11P4_5piRest = loc3PiP4_11;
		TLorentzVector loc3Pi12P4_5piRest = loc3PiP4_12;
		TLorentzVector loc3Pi21P4_5piRest = loc3PiP4_21;
		TLorentzVector loc3Pi22P4_5piRest = loc3PiP4_22;
		loc3Pi11P4_5piRest.Boost(-1.0*locGammapBoost);
		loc3Pi12P4_5piRest.Boost(-1.0*locGammapBoost);
		loc3Pi21P4_5piRest.Boost(-1.0*locGammapBoost);
		loc3Pi22P4_5piRest.Boost(-1.0*locGammapBoost);
		loc3Pi11P4_5piRest.Boost(-1.0*loc5PiBoost);
		loc3Pi12P4_5piRest.Boost(-1.0*loc5PiBoost);
		loc3Pi21P4_5piRest.Boost(-1.0*loc5PiBoost);
		loc3Pi22P4_5piRest.Boost(-1.0*loc5PiBoost);
		//Boost possible rho combos to 5pi r.f.
		TLorentzVector loc2Pi11P4_5piRest = loc2PiP4_11;
		TLorentzVector loc2Pi12P4_5piRest = loc2PiP4_12;
		TLorentzVector loc2Pi21P4_5piRest = loc2PiP4_21;
		TLorentzVector loc2Pi22P4_5piRest = loc2PiP4_22;
		loc2Pi11P4_5piRest.Boost(-1.0*locGammapBoost);
		loc2Pi12P4_5piRest.Boost(-1.0*locGammapBoost);
		loc2Pi21P4_5piRest.Boost(-1.0*locGammapBoost);
		loc2Pi22P4_5piRest.Boost(-1.0*locGammapBoost);
		loc2Pi11P4_5piRest.Boost(-1.0*loc5PiBoost);
		loc2Pi12P4_5piRest.Boost(-1.0*loc5PiBoost);
		loc2Pi21P4_5piRest.Boost(-1.0*loc5PiBoost);
		loc2Pi22P4_5piRest.Boost(-1.0*loc5PiBoost);

		//Basis vectors:
		TVector3 locz_omegarho = loc5PiP4_gpRest.Vect().Unit();
		TVector3 locy_omegarho = (locBeamP4_gpRest.Vect().Unit()).Cross(locz_omegarho).Unit();
		TVector3 locx_omegarho = locy_omegarho.Cross(locz_omegarho);

		TVector3 locOmega11P3 = loc3Pi11P4_5piRest.Vect();
		TVector3 locOmega12P3 = loc3Pi12P4_5piRest.Vect();
		TVector3 locOmega21P3 = loc3Pi21P4_5piRest.Vect();
		TVector3 locOmega22P3 = loc3Pi22P4_5piRest.Vect();
		TVector3 locRho11P3 = loc2Pi11P4_5piRest.Vect();
		TVector3 locRho12P3 = loc2Pi12P4_5piRest.Vect();
		TVector3 locRho21P3 = loc2Pi21P4_5piRest.Vect();
		TVector3 locRho22P3 = loc2Pi22P4_5piRest.Vect();

		double locCosThetaOmega11 = locOmega11P3.Dot(locz_omegarho) / locOmega11P3.Mag();
		double locCosThetaOmega12 = locOmega12P3.Dot(locz_omegarho) / locOmega12P3.Mag();
		double locCosThetaOmega21 = locOmega21P3.Dot(locz_omegarho) / locOmega21P3.Mag();
		double locCosThetaOmega22 = locOmega22P3.Dot(locz_omegarho) / locOmega22P3.Mag();

		double locCosThetaRho11 = locRho11P3.Dot(locz_omegarho) / locRho11P3.Mag();
		double locCosThetaRho12 = locRho12P3.Dot(locz_omegarho) / locRho12P3.Mag();
		double locCosThetaRho21 = locRho21P3.Dot(locz_omegarho) / locRho21P3.Mag();
		double locCosThetaRho22 = locRho22P3.Dot(locz_omegarho) / locRho22P3.Mag();

		double locPhiOmega11 = TMath::ATan2(locOmega11P3.Dot(locy_omegarho), locOmega11P3.Dot(locx_omegarho)) * 180. / TMath::Pi();
		double locPhiOmega12 = TMath::ATan2(locOmega12P3.Dot(locy_omegarho), locOmega12P3.Dot(locx_omegarho)) * 180. / TMath::Pi();
		double locPhiOmega21 = TMath::ATan2(locOmega21P3.Dot(locy_omegarho), locOmega21P3.Dot(locx_omegarho)) * 180. / TMath::Pi();
		double locPhiOmega22 = TMath::ATan2(locOmega22P3.Dot(locy_omegarho), locOmega22P3.Dot(locx_omegarho)) * 180. / TMath::Pi();

		double locPhiRho11 = TMath::ATan2(locRho11P3.Dot(locy_omegarho), locRho11P3.Dot(locx_omegarho)) * 180. / TMath::Pi();
		double locPhiRho12 = TMath::ATan2(locRho12P3.Dot(locy_omegarho), locRho12P3.Dot(locx_omegarho)) * 180. / TMath::Pi();
		double locPhiRho21 = TMath::ATan2(locRho21P3.Dot(locy_omegarho), locRho21P3.Dot(locx_omegarho)) * 180. / TMath::Pi();
		double locPhiRho22 = TMath::ATan2(locRho22P3.Dot(locy_omegarho), locRho22P3.Dot(locx_omegarho)) * 180. / TMath::Pi();


		/********************************** DECAY ANGLES (ASSUMING DELTA++B1- SYSTEM) **********************************/

		//omegapi- with recoil proton pi+
		//boost omega from gamma p to omegapi- r.f.
		//omega11, second pi+ is with the proton: (also works for omega12)
		TVector3 locOmegaPiMinus1Boost = loc4PiMinus1P4_gpRest.BoostVector();
		//omega21, first pi+ is with the proton: (also works for omega22)
		TVector3 locOmegaPiMinus2Boost = loc4PiMinus2P4_gpRest.BoostVector();

		TLorentzVector locOmega11P4_OmegaPiRest = loc3PiP4_11;
		locOmega11P4_OmegaPiRest.Boost(-1.*locGammapBoost);
		locOmega11P4_OmegaPiRest.Boost(-1.*locOmegaPiMinus1Boost);
		TLorentzVector locOmega12P4_OmegaPiRest = loc3PiP4_12;
		locOmega12P4_OmegaPiRest.Boost(-1.*locGammapBoost);
		locOmega12P4_OmegaPiRest.Boost(-1.*locOmegaPiMinus1Boost);
		TLorentzVector locOmega21P4_OmegaPiRest = loc3PiP4_21;
		locOmega21P4_OmegaPiRest.Boost(-1.*locGammapBoost);
		locOmega21P4_OmegaPiRest.Boost(-1.*locOmegaPiMinus2Boost);
		TLorentzVector locOmega22P4_OmegaPiRest = loc3PiP4_22;
		locOmega22P4_OmegaPiRest.Boost(-1.*locGammapBoost);
		locOmega22P4_OmegaPiRest.Boost(-1.*locOmegaPiMinus2Boost);

		TVector3 locOmega11P3_OmegaPiRest = locOmega11P4_OmegaPiRest.Vect();
		TVector3 locOmega12P3_OmegaPiRest = locOmega12P4_OmegaPiRest.Vect();
		TVector3 locOmega21P3_OmegaPiRest = locOmega21P4_OmegaPiRest.Vect();
		TVector3 locOmega22P3_OmegaPiRest = locOmega22P4_OmegaPiRest.Vect();

		//Define unit 3-vectors
		//lock defined above, unit vector in beam direction, gammap rest frame
		//Set 1, assuming second pi+ is in the Delta++ (omega11 and 12)
		TVector3 locz1 = loc4PiMinus1P4_gpRest.Vect().Unit();
		TVector3 locy1 = lock.Cross(locz1).Unit();
		TVector3 locx1 = locy1.Cross(locz1);
		//Set 2, assuming first pi+ is in the Delta++ (omega21 and 22)
		TVector3 locz2 = loc4PiMinus2P4_gpRest.Vect().Unit();
		TVector3 locy2 = lock.Cross(locz2).Unit();
		TVector3 locx2 = locy2.Cross(locz2);

		double locCosTheta11 = locOmega11P3_OmegaPiRest.Dot(locz1) / locOmega11P3_OmegaPiRest.Mag();
		double locPhi11 = TMath::ATan2(locOmega11P3_OmegaPiRest.Dot(locy1), locOmega11P3_OmegaPiRest.Dot(locx1)) * 180. / TMath::Pi();
		double locCosTheta12 = locOmega12P3_OmegaPiRest.Dot(locz1) / locOmega12P3_OmegaPiRest.Mag();
		double locPhi12 = TMath::ATan2(locOmega12P3_OmegaPiRest.Dot(locy1), locOmega12P3_OmegaPiRest.Dot(locx1)) * 180. / TMath::Pi();
		double locCosTheta21 = locOmega21P3_OmegaPiRest.Dot(locz2) / locOmega21P3_OmegaPiRest.Mag();
		double locPhi21 = TMath::ATan2(locOmega21P3_OmegaPiRest.Dot(locy2), locOmega21P3_OmegaPiRest.Dot(locx2)) * 180. / TMath::Pi();
		double locCosTheta22 = locOmega22P3_OmegaPiRest.Dot(locz2) / locOmega22P3_OmegaPiRest.Mag();
		double locPhi22 = TMath::ATan2(locOmega22P3_OmegaPiRest.Dot(locy2), locOmega22P3_OmegaPiRest.Dot(locx2)) * 180. / TMath::Pi();

		//Define unit vectors for omega decay plane
		//omega11
		TVector3 loczH_11 = locOmega11P3_OmegaPiRest.Unit();
		TVector3 locyH_11 = locz1.Cross(loczH_11).Unit();
		TVector3 locxH_11 = locyH_11.Cross(loczH_11);
		//omega12
		TVector3 loczH_12 = locOmega12P3_OmegaPiRest.Unit();
		TVector3 locyH_12 = locz1.Cross(loczH_12).Unit();
		TVector3 locxH_12 = locyH_12.Cross(loczH_12);
		//omega21
		TVector3 loczH_21 = locOmega21P3_OmegaPiRest.Unit();
		TVector3 locyH_21 = locz2.Cross(loczH_21).Unit();
		TVector3 locxH_21 = locyH_21.Cross(loczH_21);
		//omega22
		TVector3 loczH_22 = locOmega22P3_OmegaPiRest.Unit();
		TVector3 locyH_22 = locz2.Cross(loczH_22).Unit();
		TVector3 locxH_22 = locyH_22.Cross(loczH_22);

		//Boost the charged pions (from the omega) into the omega r.f.
		//Omega11 - pi+1, pi-1
		TLorentzVector locOmega11PiPlusP4 = locPiPlus1P4;
		TLorentzVector locOmega11PiMinusP4 = locPiMinus1P4;
		//Boost from lab to gammap to omegapi r.f.
		locOmega11PiPlusP4.Boost(-1.*locGammapBoost);
		locOmega11PiMinusP4.Boost(-1.*locGammapBoost);
		locOmega11PiPlusP4.Boost(-1.*locOmegaPiMinus1Boost);
		locOmega11PiMinusP4.Boost(-1.*locOmegaPiMinus1Boost);
		//Boost from omegapi to omega r.f.
		TVector3 locOmega11Boost = locOmega11P4_OmegaPiRest.BoostVector();
		locOmega11PiPlusP4.Boost(-1.*locOmega11Boost);
		locOmega11PiMinusP4.Boost(-1.*locOmega11Boost);

		//Omega12 - pi+1, pi-2
		TLorentzVector locOmega12PiPlusP4 = locPiPlus1P4;
		TLorentzVector locOmega12PiMinusP4 = locPiMinus2P4;
		//Boost from lab to gammap to omegapi r.f.
		locOmega12PiPlusP4.Boost(-1.*locGammapBoost);
		locOmega12PiMinusP4.Boost(-1.*locGammapBoost);
		locOmega12PiPlusP4.Boost(-1.*locOmegaPiMinus1Boost);
		locOmega12PiMinusP4.Boost(-1.*locOmegaPiMinus1Boost);
		//Boost from omegapi to omega r.f.
		TVector3 locOmega12Boost = locOmega12P4_OmegaPiRest.BoostVector();
		locOmega12PiPlusP4.Boost(-1.*locOmega12Boost);
		locOmega12PiMinusP4.Boost(-1.*locOmega12Boost);

		//Omega21 - pi+2, pi-1
		TLorentzVector locOmega21PiPlusP4 = locPiPlus2P4;
		TLorentzVector locOmega21PiMinusP4 = locPiMinus1P4;
		//Boost from lab to gammap to omegapi r.f.
		locOmega21PiPlusP4.Boost(-1.*locGammapBoost);
		locOmega21PiMinusP4.Boost(-1.*locGammapBoost);
		locOmega21PiPlusP4.Boost(-1.*locOmegaPiMinus2Boost);
		locOmega21PiMinusP4.Boost(-1.*locOmegaPiMinus2Boost);
		//Boost from omegapi to omega r.f.
		TVector3 locOmega21Boost = locOmega21P4_OmegaPiRest.BoostVector();
		locOmega21PiPlusP4.Boost(-1.*locOmega21Boost);
		locOmega21PiMinusP4.Boost(-1.*locOmega21Boost);

		//Omega22 - pi+2, pi-2
		TLorentzVector locOmega22PiPlusP4 = locPiPlus2P4;
		TLorentzVector locOmega22PiMinusP4 = locPiMinus2P4;
		//Boost from lab to gammap to omegapi r.f.
		locOmega22PiPlusP4.Boost(-1.*locGammapBoost);
		locOmega22PiMinusP4.Boost(-1.*locGammapBoost);
		locOmega22PiPlusP4.Boost(-1.*locOmegaPiMinus2Boost);
		locOmega22PiMinusP4.Boost(-1.*locOmegaPiMinus2Boost);
		//Boost from omegapi to omega r.f.
		TVector3 locOmega22Boost = locOmega22P4_OmegaPiRest.BoostVector();
		locOmega22PiPlusP4.Boost(-1.*locOmega22Boost);
		locOmega22PiMinusP4.Boost(-1.*locOmega22Boost);

		//Define normal vectors and decay angles for the omega decay planes
		//Omega11
		TVector3 locOmega11PiPlusP3 = locOmega11PiPlusP4.Vect();
		TVector3 locOmega11PiMinusP3 = locOmega11PiMinusP4.Vect();
		TVector3 locnormal11 = locOmega11PiPlusP3.Cross(locOmega11PiMinusP3).Unit();
		double locCosThetaH_11 = locnormal11.Dot(loczH_11);
		double locPhiH_11 = TMath::ATan2(locnormal11.Dot(locyH_11), locnormal11.Dot(locxH_11)) * 180. / TMath::Pi();

		//Omega12
		TVector3 locOmega12PiPlusP3 = locOmega12PiPlusP4.Vect();
		TVector3 locOmega12PiMinusP3 = locOmega12PiMinusP4.Vect();
		TVector3 locnormal12 = locOmega12PiPlusP3.Cross(locOmega12PiMinusP3).Unit();
		double locCosThetaH_12 = locnormal12.Dot(loczH_12);
		double locPhiH_12 = TMath::ATan2(locnormal12.Dot(locyH_12), locnormal12.Dot(locxH_12)) * 180. / TMath::Pi();

		//Omega21
		TVector3 locOmega21PiPlusP3 = locOmega21PiPlusP4.Vect();
		TVector3 locOmega21PiMinusP3 = locOmega21PiMinusP4.Vect();
		TVector3 locnormal21 = locOmega21PiPlusP3.Cross(locOmega21PiMinusP3).Unit();
		double locCosThetaH_21 = locnormal21.Dot(loczH_21);
		double locPhiH_21 = TMath::ATan2(locnormal21.Dot(locyH_21), locnormal21.Dot(locxH_21)) * 180. / TMath::Pi();

		//Omega22
		TVector3 locOmega22PiPlusP3 = locOmega22PiPlusP4.Vect();
		TVector3 locOmega22PiMinusP3 = locOmega22PiMinusP4.Vect();
		TVector3 locnormal22 = locOmega22PiPlusP3.Cross(locOmega22PiMinusP3).Unit();
		double locCosThetaH_22 = locnormal22.Dot(loczH_22);
		double locPhiH_22 = TMath::ATan2(locnormal22.Dot(locyH_22), locnormal21.Dot(locxH_22)) * 180. / TMath::Pi();



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

		if(fabs(locBeamDeltaT) < 0.5*4.008) { // prompt signal receives a weight of 1
			locAccWeight = 1.;
			locAccid = false;
		}
                else { // accidentals recieve a weight of 1/# RF bunches included in TTree (8 in this case)
			locAccWeight = -1./8.;
			locAccid = true;
		}


		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);


		/***************************** DELTA++ VETO ***************************/
		if(veto_Deltaplusplus == true) {
		  if(locProtonPiPlusMass1 < 1.4){
		    dComboWrapper->Set_IsComboCut(true);
		    continue;
		  }
		  if(locProtonPiPlusMass2 < 1.4){
		    dComboWrapper->Set_IsComboCut(true);
		    continue;
		  }
		}

		/***************************** DELTA0 VETO ***************************/
		if(veto_Delta0 == true) {
		  if(locProtonPiMinusMass1 < 1.8){
		    dComboWrapper->Set_IsComboCut(true);
		    continue;
		  }
		  if(locProtonPiMinusMass2 < 1.8){
		   dComboWrapper->Set_IsComboCut(true);
		   continue;
		  }
		}


		/******************************** Fill histogram of thrown topologies *************************************/
		
		if(locThrownTopology != ""){
		  dHist_ThrownTopologies->Fill(locThrownTopology.Data(), locAccWeight);

		  if(dHist_InvariantMass_ThrownTopology.find(locThrownTopology) != dHist_InvariantMass_ThrownTopology.end())
		    dHist_InvariantMass_ThrownTopology[locThrownTopology]->Fill(loc5PiP4.M(), locAccWeight);
		}
		
	       
		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
		  dHist_BeamEnergy->Fill(locBeamP4.E(), locAccWeight);
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
		  dHist_MissingMassSquared->Fill(locMissingMassSquared, locAccWeight);
		  locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}

		/****************************************** Set up omega sideband weights ************************************/

		//5Pi Mass
		double loc5PiMass = loc5PiP4.M();

		double loc3PiMass11 = loc3PiP4_11.M();
		double loc3PiMass12 = loc3PiP4_12.M();
		double loc3PiMass21 = loc3PiP4_21.M();
		double loc3PiMass22 = loc3PiP4_22.M();

		double loc2PiMass11 = loc2PiP4_11.M();
		double loc2PiMass12 = loc2PiP4_12.M();
		double loc2PiMass21 = loc2PiP4_21.M();
		double loc2PiMass22 = loc2PiP4_22.M();

		TLorentzVector locOmegaPiPlusP4_11 = loc3PiP4_11 + locPiPlus2P4; 
		TLorentzVector locOmegaPiPlusP4_12 = loc3PiP4_12 + locPiPlus2P4; 
		TLorentzVector locOmegaPiPlusP4_21 = loc3PiP4_21 + locPiPlus1P4; 
		TLorentzVector locOmegaPiPlusP4_22 = loc3PiP4_22 + locPiPlus1P4; 

		double locOmegaPiPlusMass11 = locOmegaPiPlusP4_11.M();
		double locOmegaPiPlusMass12 = locOmegaPiPlusP4_12.M();
		double locOmegaPiPlusMass21 = locOmegaPiPlusP4_21.M();
		double locOmegaPiPlusMass22 = locOmegaPiPlusP4_22.M();

		TLorentzVector locOmegaPiMinusP4_11 = loc3PiP4_11 + locPiMinus2P4; 
		TLorentzVector locOmegaPiMinusP4_12 = loc3PiP4_12 + locPiMinus1P4; 
		TLorentzVector locOmegaPiMinusP4_21 = loc3PiP4_21 + locPiMinus2P4; 
		TLorentzVector locOmegaPiMinusP4_22 = loc3PiP4_22 + locPiMinus1P4; 

		double locOmegaPiMinusMass11 = locOmegaPiMinusP4_11.M();
		double locOmegaPiMinusMass12 = locOmegaPiMinusP4_12.M();
		double locOmegaPiMinusMass21 = locOmegaPiMinusP4_21.M();
		double locOmegaPiMinusMass22 = locOmegaPiMinusP4_22.M();


		//Here we'll define weights to select the omega events
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



		/****************************************** HISTOGRAM 4-MOMENTUM TRANSFER SQUARED (-t) ************************************************/
		TLorentzVector locSqrt_t_proton = loc5PiP4 - locBeamP4;
		TLorentzVector locSqrt_s = locBeamP4 + dTargetP4;

		double loct_proton = fabs(locSqrt_t_proton.Dot(locSqrt_t_proton));

		double loc_s = fabs(locSqrt_s.Dot(locSqrt_s));
		double locEstargamma = (loc_s - locProtonP4.M2())/(2*TMath::Sqrt(loc_s));
		double locEstarX = (loc_s + loc5PiP4.M2() - locProtonP4.M2())/(2*TMath::Sqrt(loc_s));
		double locpstarX = TMath::Sqrt(locEstarX*locEstarX - loc5PiP4.M2());
		double loct_proton_min = fabs(loc5PiP4.M2() - 2*locEstargamma*(locEstarX - locpstarX));

		double loct_proton_prime = loct_proton - loct_proton_min;

		TLorentzVector locSqrt_t_omega11 = locBeamP4 - loc3PiP4_11;
		TLorentzVector locSqrt_t_omega12 = locBeamP4 - loc3PiP4_12;
		TLorentzVector locSqrt_t_omega21 = locBeamP4 - loc3PiP4_21;
		TLorentzVector locSqrt_t_omega22 = locBeamP4 - loc3PiP4_22;

		double loct_omega11 = fabs(locSqrt_t_omega11.Dot(locSqrt_t_omega11));
		double loct_omega12 = fabs(locSqrt_t_omega12.Dot(locSqrt_t_omega12));
		double loct_omega21 = fabs(locSqrt_t_omega21.Dot(locSqrt_t_omega21));
		double loct_omega22 = fabs(locSqrt_t_omega22.Dot(locSqrt_t_omega22));

		TLorentzVector locSqrt_t_Deltaplusplus1 = dTargetP4 - locProtonPiPlusP4_1;
		TLorentzVector locSqrt_t_Deltaplusplus2 = dTargetP4 - locProtonPiPlusP4_2;

		double loct_Deltaplusplus1 = fabs(locSqrt_t_Deltaplusplus1.Dot(locSqrt_t_Deltaplusplus1));
		double loct_Deltaplusplus2 = fabs(locSqrt_t_Deltaplusplus2.Dot(locSqrt_t_Deltaplusplus2));

		TLorentzVector locSqrt_t_Delta0_1 = dTargetP4 - locProtonPiMinusP4_1;
		TLorentzVector locSqrt_t_Delta0_2 = dTargetP4 - locProtonPiMinusP4_2;

		double loct_Delta0_1 = fabs(locSqrt_t_Delta0_1.Dot(locSqrt_t_Delta0_1));
		double loct_Delta0_2 = fabs(locSqrt_t_Delta0_2.Dot(locSqrt_t_Delta0_2));

		//Uniqueness Tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_t_prime;
		locUsedThisCombo_t_prime[Unknown].insert(locBeamID);
		locUsedThisCombo_t_prime[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_t_prime[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_t_prime[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_t_prime[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_t_prime[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_t_prime[Gamma].insert(locPhoton2NeutralID);


		//compare to what's been used so far
		if(locUsedSoFar_t_prime.find(locUsedThisCombo_t_prime) == locUsedSoFar_t_prime.end())
		  {
		    dHist_t_proton_total->Fill(loct_proton_prime, weight11 * locAccWeight);
		    dHist_t_proton_total->Fill(loct_proton_prime, weight12 * locAccWeight);
		    dHist_t_proton_total->Fill(loct_proton_prime, weight21 * locAccWeight);
		    dHist_t_proton_total->Fill(loct_proton_prime, weight22 * locAccWeight);

		    if(loc5PiMass > 1.2 && loc5PiMass < 1.8){
		      dHist_t_pomega_correlation->Fill(loct_proton_prime, loct_omega11, weight11 * locAccWeight);
		      dHist_t_pomega_correlation->Fill(loct_proton_prime, loct_omega12, weight12 * locAccWeight);
		      dHist_t_pomega_correlation->Fill(loct_proton_prime, loct_omega21, weight21 * locAccWeight);
		      dHist_t_pomega_correlation->Fill(loct_proton_prime, loct_omega22, weight22 * locAccWeight);
		    }
		    
		    dHist_2PiMass_t_proton_total->Fill(loct_proton_prime, loc2PiMass11, weight11 * locAccWeight);
		    dHist_2PiMass_t_proton_total->Fill(loct_proton_prime, loc2PiMass12, weight12 * locAccWeight);
		    dHist_2PiMass_t_proton_total->Fill(loct_proton_prime, loc2PiMass21, weight21 * locAccWeight);
		    dHist_2PiMass_t_proton_total->Fill(loct_proton_prime, loc2PiMass22, weight22 * locAccWeight);

		    //bin the -t plots in omegapipi mass:
		    for(int i = 0; i < 10; i++){
		      double massmin = 1.4 + 0.1*i;
		      double massmax = 1.5 + 0.1*i;

		      if(loc5PiMass > massmin && loc5PiMass < massmax){
			dHist_t_proton[i]->Fill(loct_proton_prime, weight11 * locAccWeight);
			dHist_t_proton[i]->Fill(loct_proton_prime, weight12 * locAccWeight);
			dHist_t_proton[i]->Fill(loct_proton_prime, weight21 * locAccWeight);
			dHist_t_proton[i]->Fill(loct_proton_prime, weight22 * locAccWeight);
			
			dHist_2PiMass_t_proton[i]->Fill(loct_proton_prime, loc2PiMass11, weight11 * locAccWeight);
			dHist_2PiMass_t_proton[i]->Fill(loct_proton_prime, loc2PiMass12, weight12 * locAccWeight);
			dHist_2PiMass_t_proton[i]->Fill(loct_proton_prime, loc2PiMass21, weight21 * locAccWeight);
			dHist_2PiMass_t_proton[i]->Fill(loct_proton_prime, loc2PiMass22, weight22 * locAccWeight);
		      }
		    }

		    locUsedSoFar_t_prime.insert(locUsedThisCombo_t_prime);
		  }

		/************************************ CUT ON -t' *****************************************************/
		if(loct_proton_prime > 0.5 && cut_t_prime == true){
		  dComboWrapper->Set_IsComboCut(true);
		  continue;
		}


		//Uniqueness Tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_Man_t;
		locUsedThisCombo_Man_t[Unknown].insert(locBeamID);
		locUsedThisCombo_Man_t[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_Man_t[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_Man_t[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_Man_t[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_Man_t[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_Man_t[Gamma].insert(locPhoton2NeutralID);


		//compare to what's been used so far
		if(locUsedSoFar_Man_t.find(locUsedThisCombo_Man_t) == locUsedSoFar_Man_t.end())
		  {
		    //unique combo: histogram it, and register this combo of particles
		    dHist_t_omega_total->Fill(loct_omega11, weight11 * locAccWeight);
		    dHist_t_omega_total->Fill(loct_omega12, weight12 * locAccWeight);
		    dHist_t_omega_total->Fill(loct_omega21, weight21 * locAccWeight);
		    dHist_t_omega_total->Fill(loct_omega22, weight22 * locAccWeight);

		    dHist_t_Deltaplusplus_total->Fill(loct_Deltaplusplus2, weight11 * locAccWeight);
		    dHist_t_Deltaplusplus_total->Fill(loct_Deltaplusplus2, weight12 * locAccWeight);
		    dHist_t_Deltaplusplus_total->Fill(loct_Deltaplusplus1, weight21 * locAccWeight);
		    dHist_t_Deltaplusplus_total->Fill(loct_Deltaplusplus1, weight22 * locAccWeight);

		    dHist_t_Delta0_total->Fill(loct_Delta0_2, weight11 * locAccWeight);
		    dHist_t_Delta0_total->Fill(loct_Delta0_1, weight12 * locAccWeight);
		    dHist_t_Delta0_total->Fill(loct_Delta0_2, weight21 * locAccWeight);
		    dHist_t_Delta0_total->Fill(loct_Delta0_1, weight22 * locAccWeight);
		 
		    dHist_2PiMass_t_omega_total->Fill(loct_omega11, loc2PiMass11, weight11 * locAccWeight);
		    dHist_2PiMass_t_omega_total->Fill(loct_omega12, loc2PiMass12, weight12 * locAccWeight);
		    dHist_2PiMass_t_omega_total->Fill(loct_omega21, loc2PiMass21, weight21 * locAccWeight);
		    dHist_2PiMass_t_omega_total->Fill(loct_omega22, loc2PiMass22, weight22 * locAccWeight);

		    dHist_2PiMass_t_Deltaplusplus_total->Fill(loct_Deltaplusplus2, loc2PiMass11, weight11 * locAccWeight);
		    dHist_2PiMass_t_Deltaplusplus_total->Fill(loct_Deltaplusplus2, loc2PiMass12, weight12 * locAccWeight);
		    dHist_2PiMass_t_Deltaplusplus_total->Fill(loct_Deltaplusplus1, loc2PiMass21, weight21 * locAccWeight);
		    dHist_2PiMass_t_Deltaplusplus_total->Fill(loct_Deltaplusplus1, loc2PiMass22, weight22 * locAccWeight);

		    dHist_2PiMass_t_Delta0_total->Fill(loct_Delta0_2, loc2PiMass11, weight11 * locAccWeight);
		    dHist_2PiMass_t_Delta0_total->Fill(loct_Delta0_1, loc2PiMass12, weight12 * locAccWeight);
		    dHist_2PiMass_t_Delta0_total->Fill(loct_Delta0_2, loc2PiMass21, weight21 * locAccWeight);
		    dHist_2PiMass_t_Delta0_total->Fill(loct_Delta0_1, loc2PiMass22, weight22 * locAccWeight);

		    dHist_t_Delta_correlation_total->Fill(loct_Deltaplusplus2, loct_Delta0_2, weight11 * locAccWeight);
		    dHist_t_Delta_correlation_total->Fill(loct_Deltaplusplus2, loct_Delta0_1, weight12 * locAccWeight);
		    dHist_t_Delta_correlation_total->Fill(loct_Deltaplusplus1, loct_Delta0_2, weight21 * locAccWeight);
		    dHist_t_Delta_correlation_total->Fill(loct_Deltaplusplus1, loct_Delta0_1, weight22 * locAccWeight);

		    //bin the -t plots in omegapipi mass:
		    for(int i = 0; i < 10; i++){
		      double massmin = 1.4 + 0.1*i;
		      double massmax = 1.5 + 0.1*i;

		      if(loc5PiMass > massmin && loc5PiMass < massmax){
			dHist_t_omega[i]->Fill(loct_omega11, weight11 * locAccWeight);
			dHist_t_omega[i]->Fill(loct_omega12, weight12 * locAccWeight);
			dHist_t_omega[i]->Fill(loct_omega21, weight21 * locAccWeight);
			dHist_t_omega[i]->Fill(loct_omega22, weight22 * locAccWeight);

			dHist_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus2, weight11 * locAccWeight);
			dHist_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus2, weight12 * locAccWeight);
			dHist_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus1, weight21 * locAccWeight);
			dHist_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus1, weight22 * locAccWeight);

			dHist_t_Delta0[i]->Fill(loct_Delta0_2, weight11 * locAccWeight);
			dHist_t_Delta0[i]->Fill(loct_Delta0_1, weight12 * locAccWeight);
			dHist_t_Delta0[i]->Fill(loct_Delta0_2, weight21 * locAccWeight);
			dHist_t_Delta0[i]->Fill(loct_Delta0_1, weight22 * locAccWeight);

			dHist_t_Delta_correlation[i]->Fill(loct_Deltaplusplus2, loct_Delta0_2, weight11 * locAccWeight);
			dHist_t_Delta_correlation[i]->Fill(loct_Deltaplusplus2, loct_Delta0_1, weight12 * locAccWeight);
			dHist_t_Delta_correlation[i]->Fill(loct_Deltaplusplus1, loct_Delta0_2, weight21 * locAccWeight);
			dHist_t_Delta_correlation[i]->Fill(loct_Deltaplusplus1, loct_Delta0_1, weight22 * locAccWeight);

			dHist_2PiMass_t_omega[i]->Fill(loct_omega11, loc2PiMass11, weight11 * locAccWeight);
			dHist_2PiMass_t_omega[i]->Fill(loct_omega12, loc2PiMass12, weight12 * locAccWeight);
			dHist_2PiMass_t_omega[i]->Fill(loct_omega21, loc2PiMass21, weight21 * locAccWeight);
			dHist_2PiMass_t_omega[i]->Fill(loct_omega22, loc2PiMass22, weight22 * locAccWeight);

			dHist_2PiMass_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus2, loc2PiMass11, weight11 * locAccWeight);
			dHist_2PiMass_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus2, loc2PiMass12, weight12 * locAccWeight);
			dHist_2PiMass_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus1, loc2PiMass21, weight21 * locAccWeight);
			dHist_2PiMass_t_Deltaplusplus[i]->Fill(loct_Deltaplusplus1, loc2PiMass22, weight22 * locAccWeight);

			dHist_2PiMass_t_Delta0[i]->Fill(loct_Delta0_2, loc2PiMass11, weight11 * locAccWeight);
			dHist_2PiMass_t_Delta0[i]->Fill(loct_Delta0_1, loc2PiMass12, weight12 * locAccWeight);
			dHist_2PiMass_t_Delta0[i]->Fill(loct_Delta0_2, loc2PiMass21, weight21 * locAccWeight);
			dHist_2PiMass_t_Delta0[i]->Fill(loct_Delta0_1, loc2PiMass22, weight22 * locAccWeight);
		      }
		    }

		    dHist_omega2piMass_t_omega->Fill(loct_omega11, loc5PiMass, weight11 * locAccWeight);
		    dHist_omega2piMass_t_omega->Fill(loct_omega12, loc5PiMass, weight12 * locAccWeight);
		    dHist_omega2piMass_t_omega->Fill(loct_omega21, loc5PiMass, weight21 * locAccWeight);
		    dHist_omega2piMass_t_omega->Fill(loct_omega22, loc5PiMass, weight22 * locAccWeight);

		    dHist_omegapiplusMass_t_omega->Fill(loct_omega11, locOmegaPiPlusMass11, weight11 * locAccWeight);
		    dHist_omegapiplusMass_t_omega->Fill(loct_omega12, locOmegaPiPlusMass12, weight12 * locAccWeight);
		    dHist_omegapiplusMass_t_omega->Fill(loct_omega21, locOmegaPiPlusMass21, weight21 * locAccWeight);
		    dHist_omegapiplusMass_t_omega->Fill(loct_omega22, locOmegaPiPlusMass22, weight22 * locAccWeight);

		    dHist_omegapiminusMass_t_omega->Fill(loct_omega11, locOmegaPiMinusMass11, weight11 * locAccWeight);
		    dHist_omegapiminusMass_t_omega->Fill(loct_omega12, locOmegaPiMinusMass12, weight12 * locAccWeight);
		    dHist_omegapiminusMass_t_omega->Fill(loct_omega21, locOmegaPiMinusMass21, weight21 * locAccWeight);
		    dHist_omegapiminusMass_t_omega->Fill(loct_omega22, locOmegaPiMinusMass22, weight22 * locAccWeight);

		    dHist_omega2piMass_t_Delta->Fill(loct_Deltaplusplus2, loc5PiMass, weight11 * locAccWeight);
		    dHist_omega2piMass_t_Delta->Fill(loct_Deltaplusplus2, loc5PiMass, weight12 * locAccWeight);
		    dHist_omega2piMass_t_Delta->Fill(loct_Deltaplusplus1, loc5PiMass, weight21 * locAccWeight);
		    dHist_omega2piMass_t_Delta->Fill(loct_Deltaplusplus1, loc5PiMass, weight22 * locAccWeight);

		    dHist_omegapiminusMass_t_Delta->Fill(loct_Deltaplusplus2, locOmegaPiMinusMass11, weight11 * locAccWeight);
		    dHist_omegapiminusMass_t_Delta->Fill(loct_Deltaplusplus2, locOmegaPiMinusMass12, weight12 * locAccWeight);
		    dHist_omegapiminusMass_t_Delta->Fill(loct_Deltaplusplus1, locOmegaPiMinusMass21, weight21 * locAccWeight);
		    dHist_omegapiminusMass_t_Delta->Fill(loct_Deltaplusplus1, locOmegaPiMinusMass22, weight22 * locAccWeight);

		    dHist_omega2piMass_t_Delta0->Fill(loct_Delta0_2, loc5PiMass, weight11 * locAccWeight);
		    dHist_omega2piMass_t_Delta0->Fill(loct_Delta0_1, loc5PiMass, weight12 * locAccWeight);
		    dHist_omega2piMass_t_Delta0->Fill(loct_Delta0_2, loc5PiMass, weight21 * locAccWeight);
		    dHist_omega2piMass_t_Delta0->Fill(loct_Delta0_1, loc5PiMass, weight22 * locAccWeight);

		    dHist_omegapiplusMass_t_Delta0->Fill(loct_Delta0_2, locOmegaPiPlusMass11, weight11 * locAccWeight);
		    dHist_omegapiplusMass_t_Delta0->Fill(loct_Delta0_1, locOmegaPiPlusMass12, weight12 * locAccWeight);
		    dHist_omegapiplusMass_t_Delta0->Fill(loct_Delta0_2, locOmegaPiPlusMass21, weight21 * locAccWeight);
		    dHist_omegapiplusMass_t_Delta0->Fill(loct_Delta0_1, locOmegaPiPlusMass22, weight22 * locAccWeight);

		    locUsedSoFar_Man_t.insert(locUsedThisCombo_Man_t);
		  }


		/************************************* HISTOGRAM 5PI MASS ********************************************************/
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_5PiMass;
		locUsedThisCombo_5PiMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_5PiMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_5PiMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_5PiMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_5PiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_5PiMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_5PiMass.find(locUsedThisCombo_5PiMass) == locUsedSoFar_5PiMass.end())
		  {
		    dHist_5PiMass->Fill(loc5PiMass, locAccWeight);

		    dHist_2PhotonMass->Fill(locPi0P4.M(), locAccWeight);

		    dHist_2PhotonVs3PiMass->Fill(loc3PiMass11, locPi0P4.M(), locAccWeight);
		    dHist_2PhotonVs3PiMass->Fill(loc3PiMass12, locPi0P4.M(), locAccWeight);
		    dHist_2PhotonVs3PiMass->Fill(loc3PiMass21, locPi0P4.M(), locAccWeight);
		    dHist_2PhotonVs3PiMass->Fill(loc3PiMass22, locPi0P4.M(), locAccWeight);

		    dHist_2PhotonVs2PiMass->Fill(loc2PiMass11, locPi0P4.M(), weight11 * locAccWeight);
		    dHist_2PhotonVs2PiMass->Fill(loc2PiMass12, locPi0P4.M(), weight12 * locAccWeight);
		    dHist_2PhotonVs2PiMass->Fill(loc2PiMass21, locPi0P4.M(), weight21 * locAccWeight);
		    dHist_2PhotonVs2PiMass->Fill(loc2PiMass22, locPi0P4.M(), weight22 * locAccWeight);

		    dHist_2PhotonVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass11, locPi0P4.M(), weight11 * locAccWeight);
		    dHist_2PhotonVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass12, locPi0P4.M(), weight12 * locAccWeight);
		    dHist_2PhotonVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass21, locPi0P4.M(), weight21 * locAccWeight);
		    dHist_2PhotonVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass22, locPi0P4.M(), weight22 * locAccWeight);

		    dHist_2PhotonVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass11, locPi0P4.M(), weight11 * locAccWeight);
		    dHist_2PhotonVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass12, locPi0P4.M(), weight12 * locAccWeight);
		    dHist_2PhotonVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass21, locPi0P4.M(), weight21 * locAccWeight);
		    dHist_2PhotonVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass22, locPi0P4.M(), weight22 * locAccWeight);

		    dHist_2PhotonVsOmega2PiMass->Fill(loc5PiMass, locPi0P4.M(), weight11 * locAccWeight);
		    dHist_2PhotonVsOmega2PiMass->Fill(loc5PiMass, locPi0P4.M(), weight12 * locAccWeight);
		    dHist_2PhotonVsOmega2PiMass->Fill(loc5PiMass, locPi0P4.M(), weight21 * locAccWeight);
		    dHist_2PhotonVsOmega2PiMass->Fill(loc5PiMass, locPi0P4.M(), weight22 * locAccWeight);

		    dHist_Phi_lab_ProtonPiPlus_vs_Ebeam->Fill(locBeamP4.E(), locProtonPiPlus1_Phi_lab, locAccWeight);
		    dHist_Phi_lab_ProtonPiPlus_vs_Ebeam->Fill(locBeamP4.E(), locProtonPiPlus2_Phi_lab, locAccWeight);

		    locUsedSoFar_5PiMass.insert(locUsedThisCombo_5PiMass);
		  }


		/************************************* HISTOGRAM 5PI MASS (MEASURED) ********************************************************/

		//5Pi Mass (Measured)
		double loc5PiMass_Measured = loc5PiP4_Measured.M();
		
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_5PiMass_Measured;
		locUsedThisCombo_5PiMass_Measured[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_5PiMass_Measured[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_5PiMass_Measured[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_5PiMass_Measured[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_5PiMass_Measured[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_5PiMass_Measured[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_5PiMass_Measured.find(locUsedThisCombo_5PiMass_Measured) == locUsedSoFar_5PiMass_Measured.end())
		  {
		    dHist_5PiMass_Measured->Fill(loc5PiMass_Measured, locAccWeight);
		    locUsedSoFar_5PiMass_Measured.insert(locUsedThisCombo_5PiMass_Measured);
		  }

		/********************************** HISTOGRAM 3PI MASS *************************************************************/

		//3Pi Mass

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_3PiMass;
		locUsedThisCombo_3PiMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_3PiMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_3PiMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_3PiMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_3PiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_3PiMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_3PiMass.find(locUsedThisCombo_3PiMass) == locUsedSoFar_3PiMass.end())
		  {
		    dHist_3PiMass->Fill(loc3PiMass11, locAccWeight);
		    dHist_3PiMass->Fill(loc3PiMass12, locAccWeight);
		    dHist_3PiMass->Fill(loc3PiMass21, locAccWeight);
		    dHist_3PiMass->Fill(loc3PiMass22, locAccWeight);
		    locUsedSoFar_3PiMass.insert(locUsedThisCombo_3PiMass);
		  }

		/********************************** HISTOGRAM 3PI MASS (MEASURED) *************************************************************/

		//3Pi Mass
		double loc3PiMass11_Measured = loc3PiP4_11_Measured.M();
		double loc3PiMass12_Measured = loc3PiP4_12_Measured.M();
		double loc3PiMass21_Measured = loc3PiP4_21_Measured.M();
		double loc3PiMass22_Measured = loc3PiP4_22_Measured.M();

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_3PiMass_Measured;
		locUsedThisCombo_3PiMass_Measured[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_3PiMass_Measured[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_3PiMass_Measured[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_3PiMass_Measured[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_3PiMass_Measured[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_3PiMass_Measured[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_3PiMass_Measured.find(locUsedThisCombo_3PiMass_Measured) == locUsedSoFar_3PiMass_Measured.end())
		  {
		    dHist_3PiMass_Measured->Fill(loc3PiMass11_Measured, locAccWeight);
		    dHist_3PiMass_Measured->Fill(loc3PiMass12_Measured, locAccWeight);
		    dHist_3PiMass_Measured->Fill(loc3PiMass21_Measured, locAccWeight);
		    dHist_3PiMass_Measured->Fill(loc3PiMass22_Measured, locAccWeight);
		    locUsedSoFar_3PiMass_Measured.insert(locUsedThisCombo_3PiMass_Measured);
		  }

		/************************************ HISTOGRAM OMEGAPI+PI- MASS *********************************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_Omega2PiMass;
		locUsedThisCombo_Omega2PiMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_Omega2PiMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_Omega2PiMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_Omega2PiMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_Omega2PiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_Omega2PiMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_Omega2PiMass.find(locUsedThisCombo_Omega2PiMass) == locUsedSoFar_Omega2PiMass.end())
		  {
		    dHist_Omega2PiMass->Fill(loc5PiMass, weight11 * locAccWeight);
		    dHist_Omega2PiMass->Fill(loc5PiMass, weight12 * locAccWeight);
		    dHist_Omega2PiMass->Fill(loc5PiMass, weight21 * locAccWeight);
		    dHist_Omega2PiMass->Fill(loc5PiMass, weight22 * locAccWeight);
		    locUsedSoFar_Omega2PiMass.insert(locUsedThisCombo_Omega2PiMass);
		  }

		/******************************** Fill histogram of thrown topologies with omega cut *************************************/
		if(locThrownTopology != ""){
		dHist_ThrownTopologies_omegacut->Fill(locThrownTopology.Data(),weight11 * locAccWeight);
		dHist_ThrownTopologies_omegacut->Fill(locThrownTopology.Data(),weight12 * locAccWeight);
		dHist_ThrownTopologies_omegacut->Fill(locThrownTopology.Data(),weight21 * locAccWeight);
		dHist_ThrownTopologies_omegacut->Fill(locThrownTopology.Data(),weight22 * locAccWeight);

		if(dHist_InvariantMass_ThrownTopology_omegacut.find(locThrownTopology) != dHist_InvariantMass_ThrownTopology_omegacut.end()){
		  dHist_InvariantMass_ThrownTopology_omegacut[locThrownTopology]->Fill(loc5PiP4.M(), weight11 * locAccWeight);
		  dHist_InvariantMass_ThrownTopology_omegacut[locThrownTopology]->Fill(loc5PiP4.M(), weight12 * locAccWeight);
		  dHist_InvariantMass_ThrownTopology_omegacut[locThrownTopology]->Fill(loc5PiP4.M(), weight21 * locAccWeight);
		  dHist_InvariantMass_ThrownTopology_omegacut[locThrownTopology]->Fill(loc5PiP4.M(), weight22 * locAccWeight);
		}

		//weighting for omega peak only (no sideband subtraction)
		double weight_nosideband11;
		if (loc3PiMass11 > 0.7479 && loc3PiMass11 < 0.8169) weight_nosideband11 = +1;
		else weight_nosideband11 = 0;
		double weight_nosideband12;
		if (loc3PiMass12 > 0.7479 && loc3PiMass12 < 0.8169) weight_nosideband12 = +1;
		else weight_nosideband12 = 0;
		double weight_nosideband21;
		if (loc3PiMass21 > 0.7479 && loc3PiMass21 < 0.8169) weight_nosideband21 = +1;
		else weight_nosideband21 = 0;
		double weight_nosideband22;
		if (loc3PiMass22 > 0.7479 && loc3PiMass22 < 0.8169) weight_nosideband22 = +1;
		else weight_nosideband22 = 0;

		//no sideband 
		dHist_ThrownTopologies_omegacut_nosideband->Fill(locThrownTopology.Data(),weight_nosideband11 * locAccWeight);
		dHist_ThrownTopologies_omegacut_nosideband->Fill(locThrownTopology.Data(),weight_nosideband12 * locAccWeight);
		dHist_ThrownTopologies_omegacut_nosideband->Fill(locThrownTopology.Data(),weight_nosideband21 * locAccWeight);
		dHist_ThrownTopologies_omegacut_nosideband->Fill(locThrownTopology.Data(),weight_nosideband22 * locAccWeight);

		if(dHist_InvariantMass_ThrownTopology_omegacut_nosideband.find(locThrownTopology) != dHist_InvariantMass_ThrownTopology_omegacut_nosideband.end()){
		  dHist_InvariantMass_ThrownTopology_omegacut_nosideband[locThrownTopology]->Fill(loc5PiP4.M(), weight_nosideband11 * locAccWeight);
		  dHist_InvariantMass_ThrownTopology_omegacut_nosideband[locThrownTopology]->Fill(loc5PiP4.M(), weight_nosideband12 * locAccWeight);
		  dHist_InvariantMass_ThrownTopology_omegacut_nosideband[locThrownTopology]->Fill(loc5PiP4.M(), weight_nosideband21 * locAccWeight);
		  dHist_InvariantMass_ThrownTopology_omegacut_nosideband[locThrownTopology]->Fill(loc5PiP4.M(), weight_nosideband22 * locAccWeight);
		}
		}

		/******************************************* HISTOGRAM 3PI MASS VS 2PI MASS ***************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_3vs2PiMass;
		locUsedThisCombo_3vs2PiMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_3vs2PiMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_3vs2PiMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_3vs2PiMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_3vs2PiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_3vs2PiMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_3vs2PiMass.find(locUsedThisCombo_3vs2PiMass) == locUsedSoFar_3vs2PiMass.end())
		  {
		    dHist_3vs2PiMass->Fill(loc3PiMass11, loc2PiMass11, locAccWeight);
		    dHist_3vs2PiMass->Fill(loc3PiMass12, loc2PiMass12, locAccWeight);
		    dHist_3vs2PiMass->Fill(loc3PiMass21, loc2PiMass21, locAccWeight);
		    dHist_3vs2PiMass->Fill(loc3PiMass22, loc2PiMass22, locAccWeight);
		    locUsedSoFar_3vs2PiMass.insert(locUsedThisCombo_3vs2PiMass);
		  }

		/***************************************** HISTOGRAM OMEGAPI+ MASS ************************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_OmegaPiPlusMass;
		locUsedThisCombo_OmegaPiPlusMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_OmegaPiPlusMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_OmegaPiPlusMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_OmegaPiPlusMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_OmegaPiPlusMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_OmegaPiPlusMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_OmegaPiPlusMass.find(locUsedThisCombo_OmegaPiPlusMass) == locUsedSoFar_OmegaPiPlusMass.end())
		  {
		    dHist_OmegaPiPlusMass->Fill(locOmegaPiPlusMass11, weight11 * locAccWeight);
		    dHist_OmegaPiPlusMass->Fill(locOmegaPiPlusMass12, weight12 * locAccWeight);
		    dHist_OmegaPiPlusMass->Fill(locOmegaPiPlusMass21, weight21 * locAccWeight);
		    dHist_OmegaPiPlusMass->Fill(locOmegaPiPlusMass22, weight22 * locAccWeight);
		    locUsedSoFar_OmegaPiPlusMass.insert(locUsedThisCombo_OmegaPiPlusMass);
		  }

		/***************************************** HISTOGRAM OMEGAPI- MASS ************************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_OmegaPiMinusMass;
		locUsedThisCombo_OmegaPiMinusMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_OmegaPiMinusMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_OmegaPiMinusMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_OmegaPiMinusMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_OmegaPiMinusMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_OmegaPiMinusMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_OmegaPiMinusMass.find(locUsedThisCombo_OmegaPiMinusMass) == locUsedSoFar_OmegaPiMinusMass.end())
		  {
		    dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusMass11, weight11 * locAccWeight);
		    dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusMass12, weight12 * locAccWeight);
		    dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusMass21, weight21 * locAccWeight);
		    dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusMass22, weight22 * locAccWeight);
		    dHist_OmegaPiPlusVsOmegaPiMinus->Fill(locOmegaPiPlusMass11, locOmegaPiMinusMass11, weight11 * locAccWeight);
		    dHist_OmegaPiPlusVsOmegaPiMinus->Fill(locOmegaPiPlusMass12, locOmegaPiMinusMass12, weight12 * locAccWeight);
		    dHist_OmegaPiPlusVsOmegaPiMinus->Fill(locOmegaPiPlusMass21, locOmegaPiMinusMass21, weight21 * locAccWeight);
		    dHist_OmegaPiPlusVsOmegaPiMinus->Fill(locOmegaPiPlusMass22, locOmegaPiMinusMass22, weight22 * locAccWeight);

		    dHist_OmegaPiPlusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiPlusMass11, weight11 * locAccWeight);
		    dHist_OmegaPiPlusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiPlusMass12, weight12 * locAccWeight);
		    dHist_OmegaPiPlusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiPlusMass21, weight21 * locAccWeight);
		    dHist_OmegaPiPlusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiPlusMass22, weight22 * locAccWeight);

		    dHist_OmegaPiMinusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiMinusMass11, weight11 * locAccWeight);
		    dHist_OmegaPiMinusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiMinusMass12, weight12 * locAccWeight);
		    dHist_OmegaPiMinusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiMinusMass21, weight21 * locAccWeight);
		    dHist_OmegaPiMinusVsOmega2PiMass->Fill(loc5PiMass, locOmegaPiMinusMass22, weight22 * locAccWeight);

		    locUsedSoFar_OmegaPiMinusMass.insert(locUsedThisCombo_OmegaPiMinusMass);
		  }

		/***************************************** HISTOGRAM PROTONPI+ MASS ***************************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_ProtonPiPlusMass;
		locUsedThisCombo_ProtonPiPlusMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_ProtonPiPlusMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_ProtonPiPlusMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_ProtonPiPlusMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_ProtonPiPlusMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_ProtonPiPlusMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_ProtonPiPlusMass.find(locUsedThisCombo_ProtonPiPlusMass) == locUsedSoFar_ProtonPiPlusMass.end())
		  {
		    dHist_ProtonPiPlusMass->Fill(locProtonPiPlusMass1, locAccWeight);
		    dHist_ProtonPiPlusMass->Fill(locProtonPiPlusMass2, locAccWeight);
		    locUsedSoFar_ProtonPiPlusMass.insert(locUsedThisCombo_ProtonPiPlusMass);
		  }

		/***************************************** HISTOGRAM PROTONPI+ MASS VS OMEGAPI- MASS ***********************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus;
		locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_ProtonPiPlusVsOmegaPiMinus.find(locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus) == locUsedSoFar_ProtonPiPlusVsOmegaPiMinus.end())
		  {
		    dHist_ProtonPiPlusVsOmegaPiMinus->Fill(locOmegaPiMinusMass11, locProtonPiPlusMass2, weight11 * locAccWeight);
		    dHist_ProtonPiPlusVsOmegaPiMinus->Fill(locOmegaPiMinusMass12, locProtonPiPlusMass2, weight12 * locAccWeight);
		    dHist_ProtonPiPlusVsOmegaPiMinus->Fill(locOmegaPiMinusMass21, locProtonPiPlusMass1, weight21 * locAccWeight);
		    dHist_ProtonPiPlusVsOmegaPiMinus->Fill(locOmegaPiMinusMass22, locProtonPiPlusMass1, weight22 * locAccWeight);
		    locUsedSoFar_ProtonPiPlusVsOmegaPiMinus.insert(locUsedThisCombo_ProtonPiPlusVsOmegaPiMinus);
		  }



		/*************************************** HISTOGRAM MASS OF "BACHELOR" PIONS ********************************************************/
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_2PiMass;
		locUsedThisCombo_2PiMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_2PiMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_2PiMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_2PiMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_2PiMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_2PiMass[Gamma].insert(locPhoton2NeutralID);


		//compare to what's been used so far
		if(locUsedSoFar_2PiMass.find(locUsedThisCombo_2PiMass) == locUsedSoFar_2PiMass.end())
		  {
		    dHist_2PiMass->Fill(loc2PiMass11, weight11 * locAccWeight);
		    dHist_2PiMass->Fill(loc2PiMass12, weight12 * locAccWeight);
		    dHist_2PiMass->Fill(loc2PiMass21, weight21 * locAccWeight);
		    dHist_2PiMass->Fill(loc2PiMass22, weight22 * locAccWeight);
		    locUsedSoFar_2PiMass.insert(locUsedThisCombo_2PiMass);
		  }

		/************************************* HISTOGRAM 2PI MASSES TO LOOK FOR CHARGED RHO CONTRIBUTIONS ************************/

		double loc2PiMass_p01 = loc2PiP4_p01.M();
		double loc2PiMass_p02 = loc2PiP4_p02.M();
		double loc2PiMass_m01 = loc2PiP4_m01.M();
		double loc2PiMass_m02 = loc2PiP4_m02.M();

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_2PiMass_charged;
		locUsedThisCombo_2PiMass_charged[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_2PiMass_charged[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_2PiMass_charged[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_2PiMass_charged[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_2PiMass_charged[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_2PiMass_charged[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_2PiMass_charged.find(locUsedThisCombo_2PiMass_charged) == locUsedSoFar_2PiMass_charged.end())
		  {
		    dHist_2PiMass_plus->Fill(loc2PiMass_p01, locAccWeight);
		    dHist_2PiMass_plus->Fill(loc2PiMass_p02, locAccWeight);

		    dHist_2PiMass_minus->Fill(loc2PiMass_m01, locAccWeight);
		    dHist_2PiMass_minus->Fill(loc2PiMass_m02, locAccWeight);

		    dHist_2PiMass_plus_omegacut->Fill(loc2PiMass_p01, weight11 * locAccWeight);
		    dHist_2PiMass_plus_omegacut->Fill(loc2PiMass_p01, weight12 * locAccWeight);
		    dHist_2PiMass_plus_omegacut->Fill(loc2PiMass_p02, weight21 * locAccWeight);
		    dHist_2PiMass_plus_omegacut->Fill(loc2PiMass_p02, weight22 * locAccWeight);

		    dHist_2PiMass_minus_omegacut->Fill(loc2PiMass_m01, weight11 * locAccWeight);
		    dHist_2PiMass_minus_omegacut->Fill(loc2PiMass_m01, weight21 * locAccWeight);
		    dHist_2PiMass_minus_omegacut->Fill(loc2PiMass_m02, weight12 * locAccWeight);
		    dHist_2PiMass_minus_omegacut->Fill(loc2PiMass_m02, weight22 * locAccWeight);

		    dHist_2PiMass_plus_notomega->Fill(loc2PiMass_p01, weight21 * locAccWeight);
		    dHist_2PiMass_plus_notomega->Fill(loc2PiMass_p01, weight22 * locAccWeight);
		    dHist_2PiMass_plus_notomega->Fill(loc2PiMass_p02, weight11 * locAccWeight);
		    dHist_2PiMass_plus_notomega->Fill(loc2PiMass_p02, weight12 * locAccWeight);

		    dHist_2PiMass_minus_notomega->Fill(loc2PiMass_m01, weight12 * locAccWeight);
		    dHist_2PiMass_minus_notomega->Fill(loc2PiMass_m01, weight22 * locAccWeight);
		    dHist_2PiMass_minus_notomega->Fill(loc2PiMass_m02, weight11 * locAccWeight);
		    dHist_2PiMass_minus_notomega->Fill(loc2PiMass_m02, weight21 * locAccWeight);

		    dHist_3vs2PiMass_plus->Fill(loc3PiMass11, loc2PiMass_p01, locAccWeight);
		    dHist_3vs2PiMass_plus->Fill(loc3PiMass12, loc2PiMass_p01, locAccWeight);
		    dHist_3vs2PiMass_plus->Fill(loc3PiMass21, loc2PiMass_p02, locAccWeight);
		    dHist_3vs2PiMass_plus->Fill(loc3PiMass22, loc2PiMass_p02, locAccWeight);

		    dHist_3vs2PiMass_minus->Fill(loc3PiMass11, loc2PiMass_m01, locAccWeight);
		    dHist_3vs2PiMass_minus->Fill(loc3PiMass21, loc2PiMass_m01, locAccWeight);
		    dHist_3vs2PiMass_minus->Fill(loc3PiMass12, loc2PiMass_m02, locAccWeight);
		    dHist_3vs2PiMass_minus->Fill(loc3PiMass22, loc2PiMass_m02, locAccWeight);

		    locUsedSoFar_2PiMass_charged.insert(locUsedThisCombo_2PiMass_charged);
		  }

		/************************************* HISTOGRAM QUANTITIES RELATED TO DELTA++B1- ENHANCEMENT ******************/
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_DeltaPlusPlusPeak;
		locUsedThisCombo_DeltaPlusPlusPeak[Proton].insert(locProtonTrackID);
		locUsedThisCombo_DeltaPlusPlusPeak[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_DeltaPlusPlusPeak[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_DeltaPlusPlusPeak[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_DeltaPlusPlusPeak[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_DeltaPlusPlusPeak[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_DeltaPlusPlusPeak[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_DeltaPlusPlusPeak.find(locUsedThisCombo_DeltaPlusPlusPeak) == locUsedSoFar_DeltaPlusPlusPeak.end())
		  {
		    if(locProtonPiPlusMass1 < 1.4 && loct_Deltaplusplus1 < 0.5) //proton and piplus1 combo in the Delta peak, also cut on momentum transfer off Delta
		      {
			//fill hists involving omega21 and omega22 combos
			dHist_5PiMass_DeltaPlusPlusPeak->Fill(loc5PiMass, locAccWeight);
			dHist_ProtonPiPlusMass_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass1, locAccWeight); //with or without omega cut?

			if(locOmegaPiMinusMass21 > 1.1 && locOmegaPiMinusMass21 < 1.35){
			  dHist_Omega2PiMass_DeltaPlusPlusPeak->Fill(loc5PiMass, weight21 * locAccWeight);
			  dHist_OmegaPiMinusMass_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass21, weight21 * locAccWeight);
			  dHist_t_DeltaPlusPlus_DeltaPlusPlusPeak->Fill(loct_Deltaplusplus1, weight21 * locAccWeight);
			  dHist_t_proton_DeltaPlusPlusPeak->Fill(loct_proton_prime, weight21 * locAccWeight);
			  dHist_t_omega_DeltaPlusPlusPeak->Fill(loct_omega21, weight21 * locAccWeight);
			  dHist_CosTheta_DeltaPlusPlusPeak->Fill(locCosTheta21, weight21 * locAccWeight);
			  dHist_Phi_DeltaPlusPlusPeak->Fill(locPhi21, weight21 * locAccWeight);
			  dHist_CosThetaH_DeltaPlusPlusPeak->Fill(locCosThetaH_21, weight21 * locAccWeight);
			  dHist_PhiH_DeltaPlusPlusPeak->Fill(locPhiH_21, weight21 * locAccWeight);
			  dHist_CosThetaVsPhi_DeltaPlusPlusPeak->Fill(locPhi21, locCosTheta21, weight21 * locAccWeight);
			  dHist_CosThetaHVsPhiH_DeltaPlusPlusPeak->Fill(locPhiH_21, locCosThetaH_21, weight21 * locAccWeight);
			  dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak->Fill(locBeamP4.E(), locProtonPiPlus1_Phi_lab, weight21 * locAccWeight);
			  dHist_CosTheta_M_omega2pi_DeltaPlusPlusPeak->Fill(loc5PiMass, locCosTheta21, weight21 * locAccWeight);
			  dHist_CosTheta_M_omegapiminus_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass21, locCosTheta21, weight21 * locAccWeight);
			  dHist_CosTheta_M_omegapiplus_DeltaPlusPlusPeak->Fill(locOmegaPiPlusMass21, locCosTheta21, weight21 * locAccWeight);
			  dHist_CosTheta_M_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass2, locCosTheta21, weight21 * locAccWeight);
			  dHist_CosTheta_M_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass1, locCosTheta21, weight21 * locAccWeight);
			  dHist_CosTheta_M_omega_DeltaPlusPlusPeak->Fill(loc3PiMass21, locCosTheta21, locAccWeight);
			  dHist_CosTheta_M_proton2pi_DeltaPlusPlusPeak->Fill(locProton2PiMass21, locCosTheta21, weight21 * locAccWeight);
			  dHist_CosTheta_M_2pi_DeltaPlusPlusPeak->Fill(loc2PiMass21, locCosTheta21, weight21 * locAccWeight);
			  dHist_M_proton2pi_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass1, locProton2PiMass21, weight21 * locAccWeight);
			  dHist_M_proton2pi_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass2, locProton2PiMass21, weight21 * locAccWeight);
			  dHist_M_proton2pi_pipluspiminus_DeltaPlusPlusPeak->Fill(loc2PiMass21, locProton2PiMass21, weight21 * locAccWeight);
			  dHist_Omega_Dalitz_xy_DeltaPlusPlusPeak->Fill(locDalitz_x_21, locDalitz_y_21, weight21 * locAccWeight);
			}
			if(locOmegaPiMinusMass22 > 1.1 && locOmegaPiMinusMass22 < 1.35){
			  dHist_Omega2PiMass_DeltaPlusPlusPeak->Fill(loc5PiMass, weight22 * locAccWeight);
			  dHist_OmegaPiMinusMass_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass22, weight22 * locAccWeight);
			  dHist_t_DeltaPlusPlus_DeltaPlusPlusPeak->Fill(loct_Deltaplusplus1, weight22 * locAccWeight);
			  dHist_t_proton_DeltaPlusPlusPeak->Fill(loct_proton_prime, weight22 * locAccWeight);
			  dHist_t_omega_DeltaPlusPlusPeak->Fill(loct_omega22, weight22 * locAccWeight);
			  dHist_CosTheta_DeltaPlusPlusPeak->Fill(locCosTheta22, weight22 * locAccWeight);
			  dHist_Phi_DeltaPlusPlusPeak->Fill(locPhi22, weight22 * locAccWeight);
			  dHist_CosThetaH_DeltaPlusPlusPeak->Fill(locCosThetaH_22, weight22 * locAccWeight);
			  dHist_PhiH_DeltaPlusPlusPeak->Fill(locPhiH_22, weight22 * locAccWeight);
			  dHist_CosThetaVsPhi_DeltaPlusPlusPeak->Fill(locPhi22, locCosTheta22, weight22 * locAccWeight);
			  dHist_CosThetaHVsPhiH_DeltaPlusPlusPeak->Fill(locPhiH_22, locCosThetaH_22, weight22 * locAccWeight);
			  dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak->Fill(locBeamP4.E(), locProtonPiPlus1_Phi_lab, weight21 * locAccWeight);
			  dHist_CosTheta_M_omega2pi_DeltaPlusPlusPeak->Fill(loc5PiMass, locCosTheta22, weight22 * locAccWeight);
			  dHist_CosTheta_M_omegapiminus_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass22, locCosTheta22, weight22 * locAccWeight);
			  dHist_CosTheta_M_omegapiplus_DeltaPlusPlusPeak->Fill(locOmegaPiPlusMass22, locCosTheta22, weight22 * locAccWeight);
			  dHist_CosTheta_M_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass1, locCosTheta22, weight22 * locAccWeight);
			  dHist_CosTheta_M_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass1, locCosTheta22, weight22 * locAccWeight);
			  dHist_CosTheta_M_omega_DeltaPlusPlusPeak->Fill(loc3PiMass22, locCosTheta22, locAccWeight);
			  dHist_CosTheta_M_proton2pi_DeltaPlusPlusPeak->Fill(locProton2PiMass22, locCosTheta22, weight22 * locAccWeight);
			  dHist_CosTheta_M_2pi_DeltaPlusPlusPeak->Fill(loc2PiMass22, locCosTheta22, weight22 * locAccWeight);
			  dHist_M_proton2pi_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass1, locProton2PiMass22, weight22 * locAccWeight);
			  dHist_M_proton2pi_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass1, locProton2PiMass22, weight22 * locAccWeight);
			  dHist_M_proton2pi_pipluspiminus_DeltaPlusPlusPeak->Fill(loc2PiMass22, locProton2PiMass22, weight22 * locAccWeight);
			  dHist_Omega_Dalitz_xy_DeltaPlusPlusPeak->Fill(locDalitz_x_22, locDalitz_y_22, weight22 * locAccWeight);
			}

		      }
		    else if(locProtonPiPlusMass2 < 1.4 && loct_Deltaplusplus2 < 0.5) //proton and piplus2 combo in the Delta peak
		      {
			//fill hists involving omega11 and omega12 combos
			dHist_5PiMass_DeltaPlusPlusPeak->Fill(loc5PiMass, locAccWeight);
			dHist_ProtonPiPlusMass_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass2, locAccWeight); //with or without omega cut?

			if(locOmegaPiMinusMass11 > 1.1 && locOmegaPiMinusMass11 < 1.35){
			  dHist_Omega2PiMass_DeltaPlusPlusPeak->Fill(loc5PiMass, weight11 * locAccWeight);
			  dHist_OmegaPiMinusMass_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass11, weight11 * locAccWeight);
			  dHist_t_DeltaPlusPlus_DeltaPlusPlusPeak->Fill(loct_Deltaplusplus2, weight11 * locAccWeight);
			  dHist_t_proton_DeltaPlusPlusPeak->Fill(loct_proton_prime, weight11 * locAccWeight);
			  dHist_t_omega_DeltaPlusPlusPeak->Fill(loct_omega11, weight11 * locAccWeight);
			  dHist_CosTheta_DeltaPlusPlusPeak->Fill(locCosTheta11, weight11 * locAccWeight);
			  dHist_Phi_DeltaPlusPlusPeak->Fill(locPhi11, weight11 * locAccWeight);
			  dHist_CosThetaH_DeltaPlusPlusPeak->Fill(locCosThetaH_11, weight11 * locAccWeight);
			  dHist_PhiH_DeltaPlusPlusPeak->Fill(locPhiH_11, weight11 * locAccWeight);
			  dHist_CosThetaVsPhi_DeltaPlusPlusPeak->Fill(locPhi11, locCosTheta11, weight11 * locAccWeight);
			  dHist_CosThetaHVsPhiH_DeltaPlusPlusPeak->Fill(locPhiH_11, locCosThetaH_11, weight11 * locAccWeight);
			  dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak->Fill(locBeamP4.E(), locProtonPiPlus2_Phi_lab, weight11 * locAccWeight);
			  dHist_CosTheta_M_omega2pi_DeltaPlusPlusPeak->Fill(loc5PiMass, locCosTheta11, weight11 * locAccWeight);
			  dHist_CosTheta_M_omegapiminus_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass11, locCosTheta11, weight11 * locAccWeight);
			  dHist_CosTheta_M_omegapiplus_DeltaPlusPlusPeak->Fill(locOmegaPiPlusMass11, locCosTheta11, weight11 * locAccWeight);
			  dHist_CosTheta_M_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass2, locCosTheta11, weight11 * locAccWeight);
			  dHist_CosTheta_M_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass2, locCosTheta11, weight11 * locAccWeight);
			  dHist_CosTheta_M_omega_DeltaPlusPlusPeak->Fill(loc3PiMass11, locCosTheta11, locAccWeight);
			  dHist_CosTheta_M_proton2pi_DeltaPlusPlusPeak->Fill(locProton2PiMass11, locCosTheta11, weight11 * locAccWeight);
			  dHist_CosTheta_M_2pi_DeltaPlusPlusPeak->Fill(loc2PiMass11, locCosTheta11, weight11 * locAccWeight);
			  dHist_M_proton2pi_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass2, locProton2PiMass11, weight11 * locAccWeight);
			  dHist_M_proton2pi_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass2, locProton2PiMass11, weight11 * locAccWeight);
			  dHist_M_proton2pi_pipluspiminus_DeltaPlusPlusPeak->Fill(loc2PiMass11, locProton2PiMass11, weight11 * locAccWeight);
			  dHist_Omega_Dalitz_xy_DeltaPlusPlusPeak->Fill(locDalitz_x_11, locDalitz_y_11, weight11 * locAccWeight);
			}
			if(locOmegaPiMinusMass12 > 1.1 && locOmegaPiMinusMass12 < 1.35){
			  dHist_Omega2PiMass_DeltaPlusPlusPeak->Fill(loc5PiMass, weight12 * locAccWeight);
			  dHist_OmegaPiMinusMass_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass12, weight12 * locAccWeight);
			  dHist_t_DeltaPlusPlus_DeltaPlusPlusPeak->Fill(loct_Deltaplusplus2, weight12 * locAccWeight);
			  dHist_t_proton_DeltaPlusPlusPeak->Fill(loct_proton_prime, weight12 * locAccWeight);
			  dHist_t_omega_DeltaPlusPlusPeak->Fill(loct_omega12, weight12 * locAccWeight);
			  dHist_CosTheta_DeltaPlusPlusPeak->Fill(locCosTheta12, weight12 * locAccWeight);
			  dHist_Phi_DeltaPlusPlusPeak->Fill(locPhi12, weight12 * locAccWeight);
			  dHist_CosThetaH_DeltaPlusPlusPeak->Fill(locCosThetaH_12, weight12 * locAccWeight);
			  dHist_PhiH_DeltaPlusPlusPeak->Fill(locPhiH_12, weight12 * locAccWeight);
			  dHist_CosThetaVsPhi_DeltaPlusPlusPeak->Fill(locPhi12, locCosTheta12, weight12 * locAccWeight);
			  dHist_CosThetaHVsPhiH_DeltaPlusPlusPeak->Fill(locPhiH_12, locCosThetaH_12, weight12 * locAccWeight);
			  dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak->Fill(locBeamP4.E(), locProtonPiPlus2_Phi_lab, weight12 * locAccWeight);
			  dHist_CosTheta_M_omega2pi_DeltaPlusPlusPeak->Fill(loc5PiMass, locCosTheta12, weight12 * locAccWeight);
			  dHist_CosTheta_M_omegapiminus_DeltaPlusPlusPeak->Fill(locOmegaPiMinusMass12, locCosTheta12, weight12 * locAccWeight);
			  dHist_CosTheta_M_omegapiplus_DeltaPlusPlusPeak->Fill(locOmegaPiPlusMass12, locCosTheta12, weight12 * locAccWeight);
			  dHist_CosTheta_M_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass1, locCosTheta12, weight12 * locAccWeight);
			  dHist_CosTheta_M_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass2, locCosTheta12, weight12 * locAccWeight);
			  dHist_CosTheta_M_omega_DeltaPlusPlusPeak->Fill(loc3PiMass12, locCosTheta12, locAccWeight);
			  dHist_CosTheta_M_proton2pi_DeltaPlusPlusPeak->Fill(locProton2PiMass12, locCosTheta12, weight12 * locAccWeight);
			  dHist_CosTheta_M_2pi_DeltaPlusPlusPeak->Fill(loc2PiMass12, locCosTheta12, weight12 * locAccWeight);
			  dHist_M_proton2pi_protonpiplus_DeltaPlusPlusPeak->Fill(locProtonPiPlusMass2, locProton2PiMass12, weight12 * locAccWeight);
			  dHist_M_proton2pi_protonpiminus_DeltaPlusPlusPeak->Fill(locProtonPiMinusMass1, locProton2PiMass12, weight12 * locAccWeight);
			  dHist_M_proton2pi_pipluspiminus_DeltaPlusPlusPeak->Fill(loc2PiMass12, locProton2PiMass12, weight12 * locAccWeight);
			  dHist_Omega_Dalitz_xy_DeltaPlusPlusPeak->Fill(locDalitz_x_12, locDalitz_y_12, weight12 * locAccWeight);
			}
		      }
		    else if(locProtonPiPlusMass1 < 1.4 && loct_Deltaplusplus1 > 0.5 && loct_Deltaplusplus1 < 1.5)
		      {
			if(locOmegaPiMinusMass21 > 1.1 && locOmegaPiMinusMass21 < 1.35)
			  {
			    dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak_high_t->Fill(locBeamP4.E(), locProtonPiPlus1_Phi_lab, weight21 * locAccWeight);
			  }
			if(locOmegaPiMinusMass22 > 1.1 && locOmegaPiMinusMass22 < 1.35)
			  {
			    dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak_high_t->Fill(locBeamP4.E(), locProtonPiPlus1_Phi_lab, weight22 * locAccWeight);
			  }
		      }
		    else if(locProtonPiPlusMass2 < 1.4 && loct_Deltaplusplus2 > 0.5 && loct_Deltaplusplus2 < 1.5)
		      {
			if(locOmegaPiMinusMass11 > 1.1 && locOmegaPiMinusMass11 < 1.35)
			  {
			    dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak_high_t->Fill(locBeamP4.E(), locProtonPiPlus2_Phi_lab, weight11 * locAccWeight);
			  }
			if(locOmegaPiMinusMass12 > 1.1 && locOmegaPiMinusMass12 < 1.35)
			  {
			    dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak_high_t->Fill(locBeamP4.E(), locProtonPiPlus2_Phi_lab, weight12 * locAccWeight);
			  }
		      }
		    locUsedSoFar_DeltaPlusPlusPeak.insert(locUsedThisCombo_DeltaPlusPlusPeak);
		  }


		/************************************* HISTOGRAM PROTONPI- MASS **************************************************************/
		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_ProtonPiMinusMass;
		locUsedThisCombo_ProtonPiMinusMass[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_ProtonPiMinusMass[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_ProtonPiMinusMass[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_ProtonPiMinusMass[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_ProtonPiMinusMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_ProtonPiMinusMass[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_ProtonPiMinusMass.find(locUsedThisCombo_ProtonPiMinusMass) == locUsedSoFar_ProtonPiMinusMass.end())
		  {
		    dHist_ProtonPiMinusMass->Fill(locProtonPiMinusMass1, locAccWeight);
		    dHist_ProtonPiMinusMass->Fill(locProtonPiMinusMass2, locAccWeight);
		    locUsedSoFar_ProtonPiMinusMass.insert(locUsedThisCombo_ProtonPiMinusMass);
		  }

		/******************************************** HISTOGRAM MASS OF PROTONPI- VS OMEGAPI+ ***********************************************************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus;
		locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_ProtonPiMinusVsOmegaPiPlus.find(locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus) == locUsedSoFar_ProtonPiMinusVsOmegaPiPlus.end())
		  {
		    dHist_ProtonPiMinusVsOmegaPiPlus->Fill(locOmegaPiPlusMass11, locProtonPiMinusMass2, weight11 * locAccWeight);
		    dHist_ProtonPiMinusVsOmegaPiPlus->Fill(locOmegaPiPlusMass12, locProtonPiMinusMass1, weight12 * locAccWeight);
		    dHist_ProtonPiMinusVsOmegaPiPlus->Fill(locOmegaPiPlusMass21, locProtonPiMinusMass2, weight21 * locAccWeight);
		    dHist_ProtonPiMinusVsOmegaPiPlus->Fill(locOmegaPiPlusMass22, locProtonPiMinusMass1, weight22 * locAccWeight);
		    locUsedSoFar_ProtonPiMinusVsOmegaPiPlus.insert(locUsedThisCombo_ProtonPiMinusVsOmegaPiPlus);
		  }

		/*********************************************** HISTOGRAM 4PI(B1+) MASS VS 3PI MASS *********************************************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_3vs4PiMass_plus;
		locUsedThisCombo_3vs4PiMass_plus[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_3vs4PiMass_plus[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_3vs4PiMass_plus[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_3vs4PiMass_plus[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_3vs4PiMass_plus[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_3vs4PiMass_plus[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_3vs4PiMass_plus.find(locUsedThisCombo_3vs4PiMass_plus) == locUsedSoFar_3vs4PiMass_plus.end())
		  {
		    dHist_3vs4PiMass_plus->Fill(loc3PiMass11, locOmegaPiPlusMass11, locAccWeight);
		    dHist_3vs4PiMass_plus->Fill(loc3PiMass12, locOmegaPiPlusMass12, locAccWeight);
		    dHist_3vs4PiMass_plus->Fill(loc3PiMass21, locOmegaPiPlusMass21, locAccWeight);
		    dHist_3vs4PiMass_plus->Fill(loc3PiMass22, locOmegaPiPlusMass22, locAccWeight);
		    locUsedSoFar_3vs4PiMass_plus.insert(locUsedThisCombo_3vs4PiMass_plus);
		  }

		/*********************************************** HISTOGRAM 4PI(B1-) MASS VS 3PI MASS *********************************************************************/

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_3vs4PiMass_minus;
		locUsedThisCombo_3vs4PiMass_minus[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_3vs4PiMass_minus[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_3vs4PiMass_minus[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_3vs4PiMass_minus[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_3vs4PiMass_minus[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_3vs4PiMass_minus[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_3vs4PiMass_minus.find(locUsedThisCombo_3vs4PiMass_minus) == locUsedSoFar_3vs4PiMass_minus.end())
		  {
		    dHist_3vs4PiMass_minus->Fill(loc3PiMass11, locOmegaPiMinusMass11, locAccWeight);
		    dHist_3vs4PiMass_minus->Fill(loc3PiMass12, locOmegaPiMinusMass12, locAccWeight);
		    dHist_3vs4PiMass_minus->Fill(loc3PiMass21, locOmegaPiMinusMass21, locAccWeight);
		    dHist_3vs4PiMass_minus->Fill(loc3PiMass22, locOmegaPiMinusMass22, locAccWeight);
		    locUsedSoFar_3vs4PiMass_minus.insert(locUsedThisCombo_3vs4PiMass_minus);
		  }

		/********************************** HISTOGRAM DECAY ANGLES ***************************************************/
		//Uniqueness Tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_DecayAngles;
		locUsedThisCombo_DecayAngles[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_DecayAngles[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_DecayAngles[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_DecayAngles[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_DecayAngles[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_DecayAngles[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_DecayAngles.find(locUsedThisCombo_DecayAngles) == locUsedSoFar_DecayAngles.end())
		  {
		    dHist_CosThetaPlus->Fill(locCosTheta1Plus, weight11 * locAccWeight);
		    dHist_CosThetaPlus->Fill(locCosTheta1Plus, weight21 * locAccWeight);
		    dHist_CosThetaPlus->Fill(locCosTheta2Plus, weight12 * locAccWeight);
		    dHist_CosThetaPlus->Fill(locCosTheta2Plus, weight22 * locAccWeight);
		    dHist_CosThetaMinus->Fill(locCosTheta1Minus, weight11 * locAccWeight);
		    dHist_CosThetaMinus->Fill(locCosTheta1Minus, weight12 * locAccWeight);
		    dHist_CosThetaMinus->Fill(locCosTheta2Minus, weight21 * locAccWeight);
		    dHist_CosThetaMinus->Fill(locCosTheta2Minus, weight22 * locAccWeight);

		    dHist_PhiPlus->Fill(locPhi1Plus, weight11 * locAccWeight);
		    dHist_PhiPlus->Fill(locPhi1Plus, weight21 * locAccWeight);
		    dHist_PhiPlus->Fill(locPhi2Plus, weight12 * locAccWeight);
		    dHist_PhiPlus->Fill(locPhi2Plus, weight22 * locAccWeight);
		    dHist_PhiMinus->Fill(locPhi1Minus, weight11 * locAccWeight);
		    dHist_PhiMinus->Fill(locPhi1Minus, weight12 * locAccWeight);
		    dHist_PhiMinus->Fill(locPhi2Minus, weight21 * locAccWeight);
		    dHist_PhiMinus->Fill(locPhi2Minus, weight22 * locAccWeight);

		    dHist_CosTheta_omega_fromPlus->Fill(locCosTheta_omega11_b1p1, weight11 * locAccWeight);
		    dHist_CosTheta_omega_fromPlus->Fill(locCosTheta_omega21_b1p1, weight21 * locAccWeight);
		    dHist_CosTheta_omega_fromPlus->Fill(locCosTheta_omega12_b1p2, weight12 * locAccWeight);
		    dHist_CosTheta_omega_fromPlus->Fill(locCosTheta_omega22_b1p2, weight22 * locAccWeight);
		    dHist_CosTheta_omega_fromMinus->Fill(locCosTheta_omega11_b1m1, weight11 * locAccWeight);
		    dHist_CosTheta_omega_fromMinus->Fill(locCosTheta_omega12_b1m1, weight12 * locAccWeight);
		    dHist_CosTheta_omega_fromMinus->Fill(locCosTheta_omega21_b1m2, weight21 * locAccWeight);
		    dHist_CosTheta_omega_fromMinus->Fill(locCosTheta_omega22_b1m2, weight22 * locAccWeight);

		    dHist_Phi_omega_fromPlus->Fill(locPhi_omega11_b1p1, weight11 * locAccWeight);
		    dHist_Phi_omega_fromPlus->Fill(locPhi_omega21_b1p1, weight21 * locAccWeight);
		    dHist_Phi_omega_fromPlus->Fill(locPhi_omega12_b1p2, weight12 * locAccWeight);
		    dHist_Phi_omega_fromPlus->Fill(locPhi_omega22_b1p2, weight22 * locAccWeight);
		    dHist_Phi_omega_fromMinus->Fill(locPhi_omega11_b1m1, weight11 * locAccWeight);
		    dHist_Phi_omega_fromMinus->Fill(locPhi_omega12_b1m1, weight12 * locAccWeight);
		    dHist_Phi_omega_fromMinus->Fill(locPhi_omega21_b1m2, weight21 * locAccWeight);
		    dHist_Phi_omega_fromMinus->Fill(locPhi_omega22_b1m2, weight22 * locAccWeight);

		    dHist_CosTheta_H->Fill(locCosTheta_H_11_b1p1, weight11 * locAccWeight);
		    dHist_CosTheta_H->Fill(locCosTheta_H_21_b1p1, weight21 * locAccWeight);
		    dHist_CosTheta_H->Fill(locCosTheta_H_12_b1p2, weight12 * locAccWeight);
		    dHist_CosTheta_H->Fill(locCosTheta_H_22_b1p2, weight22 * locAccWeight);
		    dHist_CosTheta_H->Fill(locCosTheta_H_11_b1m1, weight11 * locAccWeight);
		    dHist_CosTheta_H->Fill(locCosTheta_H_12_b1m1, weight12 * locAccWeight);
		    dHist_CosTheta_H->Fill(locCosTheta_H_21_b1m2, weight21 * locAccWeight);
		    dHist_CosTheta_H->Fill(locCosTheta_H_22_b1m2, weight22 * locAccWeight);

		    dHist_Phi_H->Fill(locPhi_H_11_b1p1, weight11 * locAccWeight);
		    dHist_Phi_H->Fill(locPhi_H_21_b1p1, weight21 * locAccWeight);
		    dHist_Phi_H->Fill(locPhi_H_12_b1p2, weight12 * locAccWeight);
		    dHist_Phi_H->Fill(locPhi_H_22_b1p2, weight22 * locAccWeight);
		    dHist_Phi_H->Fill(locPhi_H_11_b1m1, weight11 * locAccWeight);
		    dHist_Phi_H->Fill(locPhi_H_12_b1m1, weight12 * locAccWeight);
		    dHist_Phi_H->Fill(locPhi_H_21_b1m2, weight21 * locAccWeight);
		    dHist_Phi_H->Fill(locPhi_H_22_b1m2, weight22 * locAccWeight);

		    dHist_CosThetaPlusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta1Plus, weight11 * locAccWeight);
		    dHist_CosThetaPlusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta1Plus, weight21 * locAccWeight);
		    dHist_CosThetaPlusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta2Plus, weight12 * locAccWeight);
		    dHist_CosThetaPlusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta2Plus, weight22 * locAccWeight);
		    dHist_CosThetaMinusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta1Minus, weight11 * locAccWeight);
		    dHist_CosThetaMinusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta1Minus, weight12 * locAccWeight);
		    dHist_CosThetaMinusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta2Minus, weight21 * locAccWeight);
		    dHist_CosThetaMinusVsOmega2PiMass->Fill(loc5PiMass, locCosTheta2Minus, weight22 * locAccWeight);

		    dHist_PhiPlusVsOmega2PiMass->Fill(loc5PiMass, locPhi1Plus, weight11 * locAccWeight);
		    dHist_PhiPlusVsOmega2PiMass->Fill(loc5PiMass, locPhi1Plus, weight21 * locAccWeight);
		    dHist_PhiPlusVsOmega2PiMass->Fill(loc5PiMass, locPhi2Plus, weight12 * locAccWeight);
		    dHist_PhiPlusVsOmega2PiMass->Fill(loc5PiMass, locPhi2Plus, weight22 * locAccWeight);
		    dHist_PhiMinusVsOmega2PiMass->Fill(loc5PiMass, locPhi1Minus, weight11 * locAccWeight);
		    dHist_PhiMinusVsOmega2PiMass->Fill(loc5PiMass, locPhi1Minus, weight12 * locAccWeight);
		    dHist_PhiMinusVsOmega2PiMass->Fill(loc5PiMass, locPhi2Minus, weight21 * locAccWeight);
		    dHist_PhiMinusVsOmega2PiMass->Fill(loc5PiMass, locPhi2Minus, weight22 * locAccWeight);

		    dHist_CosTheta_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass11, locCosTheta_omega11_b1p1, weight11 * locAccWeight);
		    dHist_CosTheta_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass21, locCosTheta_omega21_b1p1, weight21 * locAccWeight);
		    dHist_CosTheta_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass12, locCosTheta_omega12_b1p2, weight12 * locAccWeight);
		    dHist_CosTheta_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass22, locCosTheta_omega22_b1p2, weight22 * locAccWeight);
		    dHist_CosTheta_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass11, locCosTheta_omega11_b1m1, weight11 * locAccWeight);
		    dHist_CosTheta_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass12, locCosTheta_omega12_b1m1, weight12 * locAccWeight);
		    dHist_CosTheta_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass21, locCosTheta_omega21_b1m2, weight21 * locAccWeight);
		    dHist_CosTheta_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass22, locCosTheta_omega22_b1m2, weight22 * locAccWeight);

		    dHist_Phi_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass11, locPhi_omega11_b1p1, weight11 * locAccWeight);
		    dHist_Phi_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass21, locPhi_omega21_b1p1, weight21 * locAccWeight);
		    dHist_Phi_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass12, locPhi_omega12_b1p2, weight12 * locAccWeight);
		    dHist_Phi_omegaVsOmegaPiPlusMass->Fill(locOmegaPiPlusMass22, locPhi_omega22_b1p2, weight22 * locAccWeight);
		    dHist_Phi_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass11, locPhi_omega11_b1m1, weight11 * locAccWeight);
		    dHist_Phi_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass12, locPhi_omega12_b1m1, weight12 * locAccWeight);
		    dHist_Phi_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass21, locPhi_omega21_b1m2, weight21 * locAccWeight);
		    dHist_Phi_omegaVsOmegaPiMinusMass->Fill(locOmegaPiMinusMass22, locPhi_omega22_b1m2, weight22 * locAccWeight);

		    for(int i = 0; i < 10; i++){
		      double massmin = 1. + 0.2*i;
		      double massmax = 1.2 + 0.2*i;
		      if(loc5PiMass > massmin && loc5PiMass < massmax){
			dHist_PhiVsCosTheta_b1_Plus[i]->Fill(locPhi1Plus, locCosTheta1Plus, weight11 * locAccWeight);
			dHist_PhiVsCosTheta_b1_Plus[i]->Fill(locPhi1Plus, locCosTheta1Plus, weight21 * locAccWeight);
			dHist_PhiVsCosTheta_b1_Plus[i]->Fill(locPhi2Plus, locCosTheta2Plus, weight12 * locAccWeight);
			dHist_PhiVsCosTheta_b1_Plus[i]->Fill(locPhi2Plus, locCosTheta2Plus, weight22 * locAccWeight);
			dHist_PhiVsCosTheta_b1_Minus[i]->Fill(locPhi1Minus, locCosTheta1Minus, weight11 * locAccWeight);
			dHist_PhiVsCosTheta_b1_Minus[i]->Fill(locPhi1Minus, locCosTheta1Minus, weight12 * locAccWeight);
			dHist_PhiVsCosTheta_b1_Minus[i]->Fill(locPhi2Minus, locCosTheta2Minus, weight21 * locAccWeight);
			dHist_PhiVsCosTheta_b1_Minus[i]->Fill(locPhi2Minus, locCosTheta2Minus, weight22 * locAccWeight);
		      }
		    }

		    dHist_PhiVsCosTheta_b1_Plus_total->Fill(locPhi1Plus, locCosTheta1Plus, weight11 * locAccWeight);
		    dHist_PhiVsCosTheta_b1_Plus_total->Fill(locPhi1Plus, locCosTheta1Plus, weight21 * locAccWeight);
		    dHist_PhiVsCosTheta_b1_Plus_total->Fill(locPhi2Plus, locCosTheta2Plus, weight12 * locAccWeight);
		    dHist_PhiVsCosTheta_b1_Plus_total->Fill(locPhi2Plus, locCosTheta2Plus, weight22 * locAccWeight);
		    dHist_PhiVsCosTheta_b1_Minus_total->Fill(locPhi1Minus, locCosTheta1Minus, weight11 * locAccWeight);
		    dHist_PhiVsCosTheta_b1_Minus_total->Fill(locPhi1Minus, locCosTheta1Minus, weight12 * locAccWeight);
		    dHist_PhiVsCosTheta_b1_Minus_total->Fill(locPhi2Minus, locCosTheta2Minus, weight21 * locAccWeight);
		    dHist_PhiVsCosTheta_b1_Minus_total->Fill(locPhi2Minus, locCosTheta2Minus, weight22 * locAccWeight);


		    dHist_CosThetaPlusVs2PiMass->Fill(loc2PiMass11, locCosTheta1Plus, weight11 * locAccWeight);
		    dHist_CosThetaPlusVs2PiMass->Fill(loc2PiMass21, locCosTheta1Plus, weight21 * locAccWeight);
		    dHist_CosThetaPlusVs2PiMass->Fill(loc2PiMass12, locCosTheta2Plus, weight12 * locAccWeight);
		    dHist_CosThetaPlusVs2PiMass->Fill(loc2PiMass22, locCosTheta2Plus, weight22 * locAccWeight);
		    dHist_CosThetaMinusVs2PiMass->Fill(loc2PiMass11, locCosTheta1Minus, weight11 * locAccWeight);
		    dHist_CosThetaMinusVs2PiMass->Fill(loc2PiMass12, locCosTheta1Minus, weight12 * locAccWeight);
		    dHist_CosThetaMinusVs2PiMass->Fill(loc2PiMass21, locCosTheta2Minus, weight21 * locAccWeight);
		    dHist_CosThetaMinusVs2PiMass->Fill(loc2PiMass22, locCosTheta2Minus, weight22 * locAccWeight);

		    dHist_PhiPlusVs2PiMass->Fill(loc2PiMass11, locPhi1Plus, weight11 * locAccWeight);
		    dHist_PhiPlusVs2PiMass->Fill(loc2PiMass21, locPhi1Plus, weight21 * locAccWeight);
		    dHist_PhiPlusVs2PiMass->Fill(loc2PiMass12, locPhi2Plus, weight12 * locAccWeight);
		    dHist_PhiPlusVs2PiMass->Fill(loc2PiMass22, locPhi2Plus, weight22 * locAccWeight);
		    dHist_PhiMinusVs2PiMass->Fill(loc2PiMass11, locPhi1Minus, weight11 * locAccWeight);
		    dHist_PhiMinusVs2PiMass->Fill(loc2PiMass12, locPhi1Minus, weight12 * locAccWeight);
		    dHist_PhiMinusVs2PiMass->Fill(loc2PiMass21, locPhi2Minus, weight21 * locAccWeight);
		    dHist_PhiMinusVs2PiMass->Fill(loc2PiMass22, locPhi2Minus, weight22 * locAccWeight);

		    dHist_CosTheta_omegaPlusVs2PiMass->Fill(loc2PiMass11, locCosTheta_omega11_b1p1, weight11 * locAccWeight);
		    dHist_CosTheta_omegaPlusVs2PiMass->Fill(loc2PiMass21, locCosTheta_omega21_b1p1, weight21 * locAccWeight);
		    dHist_CosTheta_omegaPlusVs2PiMass->Fill(loc2PiMass12, locCosTheta_omega12_b1p2, weight12 * locAccWeight);
		    dHist_CosTheta_omegaPlusVs2PiMass->Fill(loc2PiMass22, locCosTheta_omega22_b1p2, weight22 * locAccWeight);
		    dHist_CosTheta_omegaMinusVs2PiMass->Fill(loc2PiMass11, locCosTheta_omega11_b1m1, weight11 * locAccWeight);
		    dHist_CosTheta_omegaMinusVs2PiMass->Fill(loc2PiMass12, locCosTheta_omega12_b1m1, weight12 * locAccWeight);
		    dHist_CosTheta_omegaMinusVs2PiMass->Fill(loc2PiMass21, locCosTheta_omega21_b1m2, weight21 * locAccWeight);
		    dHist_CosTheta_omegaMinusVs2PiMass->Fill(loc2PiMass22, locCosTheta_omega22_b1m2, weight22 * locAccWeight);

		    dHist_Phi_omegaPlusVs2PiMass->Fill(loc2PiMass11, locPhi_omega11_b1p1, weight11 * locAccWeight);
		    dHist_Phi_omegaPlusVs2PiMass->Fill(loc2PiMass21, locPhi_omega21_b1p1, weight21 * locAccWeight);
		    dHist_Phi_omegaPlusVs2PiMass->Fill(loc2PiMass12, locPhi_omega12_b1p2, weight12 * locAccWeight);
		    dHist_Phi_omegaPlusVs2PiMass->Fill(loc2PiMass22, locPhi_omega22_b1p2, weight22 * locAccWeight);
		    dHist_Phi_omegaMinusVs2PiMass->Fill(loc2PiMass11, locPhi_omega11_b1m1, weight11 * locAccWeight);
		    dHist_Phi_omegaMinusVs2PiMass->Fill(loc2PiMass12, locPhi_omega12_b1m1, weight12 * locAccWeight);
		    dHist_Phi_omegaMinusVs2PiMass->Fill(loc2PiMass21, locPhi_omega21_b1m2, weight21 * locAccWeight);
		    dHist_Phi_omegaMinusVs2PiMass->Fill(loc2PiMass22, locPhi_omega22_b1m2, weight22 * locAccWeight);

		    dHist_CosThetaPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass2, locCosTheta1Plus, weight11 * locAccWeight);
		    dHist_CosThetaPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass2, locCosTheta1Plus, weight21 * locAccWeight);
		    dHist_CosThetaPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass1, locCosTheta2Plus, weight12 * locAccWeight);
		    dHist_CosThetaPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass1, locCosTheta2Plus, weight22 * locAccWeight);
		    dHist_CosThetaMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass2, locCosTheta1Minus, weight11 * locAccWeight);
		    dHist_CosThetaMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass2, locCosTheta1Minus, weight12 * locAccWeight);
		    dHist_CosThetaMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass1, locCosTheta2Minus, weight21 * locAccWeight);
		    dHist_CosThetaMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass1, locCosTheta2Minus, weight22 * locAccWeight);

		    dHist_PhiPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass2, locPhi1Plus, weight11 * locAccWeight);
		    dHist_PhiPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass2, locPhi1Plus, weight21 * locAccWeight);
		    dHist_PhiPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass1, locPhi2Plus, weight12 * locAccWeight);
		    dHist_PhiPlusVsProtonPiMinusMass->Fill(locProtonPiMinusMass1, locPhi2Plus, weight22 * locAccWeight);
		    dHist_PhiMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass2, locPhi1Minus, weight11 * locAccWeight);
		    dHist_PhiMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass2, locPhi1Minus, weight12 * locAccWeight);
		    dHist_PhiMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass1, locPhi2Minus, weight21 * locAccWeight);
		    dHist_PhiMinusVsProtonPiPlusMass->Fill(locProtonPiPlusMass1, locPhi2Minus, weight22 * locAccWeight);

		    dHist_PhiVsCosTheta_omega_Plus->Fill(locPhi_omega11_b1p1, locCosTheta_omega11_b1p1, weight11 * locAccWeight);
		    dHist_PhiVsCosTheta_omega_Plus->Fill(locPhi_omega12_b1p2, locCosTheta_omega12_b1p2, weight12 * locAccWeight);
		    dHist_PhiVsCosTheta_omega_Plus->Fill(locPhi_omega21_b1p1, locCosTheta_omega21_b1p1, weight21 * locAccWeight);
		    dHist_PhiVsCosTheta_omega_Plus->Fill(locPhi_omega22_b1p2, locCosTheta_omega22_b1p2, weight22 * locAccWeight);
		    dHist_PhiVsCosTheta_omega_Minus->Fill(locPhi_omega11_b1m1, locCosTheta_omega11_b1m1, weight11 * locAccWeight);
		    dHist_PhiVsCosTheta_omega_Minus->Fill(locPhi_omega12_b1m1, locCosTheta_omega12_b1m1, weight12 * locAccWeight);
		    dHist_PhiVsCosTheta_omega_Minus->Fill(locPhi_omega21_b1m2, locCosTheta_omega21_b1m2, weight21 * locAccWeight);
		    dHist_PhiVsCosTheta_omega_Minus->Fill(locPhi_omega22_b1m2, locCosTheta_omega22_b1m2, weight22 * locAccWeight);

		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_11_b1p1, locCosTheta_H_11_b1p1, weight11 * locAccWeight);
		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_12_b1p2, locCosTheta_H_12_b1p2, weight12 * locAccWeight);
		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_21_b1p1, locCosTheta_H_21_b1p1, weight21 * locAccWeight);
		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_22_b1p2, locCosTheta_H_22_b1p2, weight22 * locAccWeight);
		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_11_b1m1, locCosTheta_H_11_b1m1, weight11 * locAccWeight);
		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_21_b1m2, locCosTheta_H_21_b1m2, weight21 * locAccWeight);
		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_12_b1m1, locCosTheta_H_12_b1m1, weight12 * locAccWeight);
		    dHist_PhiVsCosTheta_H->Fill(locPhi_H_22_b1m2, locCosTheta_H_22_b1m2, weight22 * locAccWeight);

		    dHist_M_CosTheta_omega_Plus->Fill(locOmegaPiPlusMass11, locCosTheta_omega11_b1p1, weight11 * locAccWeight);
		    dHist_M_CosTheta_omega_Plus->Fill(locOmegaPiPlusMass12, locCosTheta_omega12_b1p2, weight12 * locAccWeight);
		    dHist_M_CosTheta_omega_Plus->Fill(locOmegaPiPlusMass21, locCosTheta_omega21_b1p1, weight21 * locAccWeight);
		    dHist_M_CosTheta_omega_Plus->Fill(locOmegaPiPlusMass22, locCosTheta_omega22_b1p2, weight22 * locAccWeight);
		    dHist_M_CosTheta_omega_Minus->Fill(locOmegaPiMinusMass11, locCosTheta_omega11_b1m1, weight11 * locAccWeight);
		    dHist_M_CosTheta_omega_Minus->Fill(locOmegaPiMinusMass12, locCosTheta_omega12_b1m1, weight12 * locAccWeight);
		    dHist_M_CosTheta_omega_Minus->Fill(locOmegaPiMinusMass21, locCosTheta_omega21_b1m2, weight21 * locAccWeight);
		    dHist_M_CosTheta_omega_Minus->Fill(locOmegaPiMinusMass22, locCosTheta_omega22_b1m2, weight22 * locAccWeight);

		    dHist_M_Phi_omega_Plus->Fill(locOmegaPiPlusMass11, locPhi_omega11_b1p1, weight11 * locAccWeight);
		    dHist_M_Phi_omega_Plus->Fill(locOmegaPiPlusMass12, locPhi_omega12_b1p2, weight12 * locAccWeight);
		    dHist_M_Phi_omega_Plus->Fill(locOmegaPiPlusMass21, locPhi_omega21_b1p1, weight21 * locAccWeight);
		    dHist_M_Phi_omega_Plus->Fill(locOmegaPiPlusMass22, locPhi_omega22_b1p2, weight22 * locAccWeight);
		    dHist_M_Phi_omega_Minus->Fill(locOmegaPiMinusMass11, locPhi_omega11_b1m1, weight11 * locAccWeight);
		    dHist_M_Phi_omega_Minus->Fill(locOmegaPiMinusMass12, locPhi_omega12_b1m1, weight12 * locAccWeight);
		    dHist_M_Phi_omega_Minus->Fill(locOmegaPiMinusMass21, locPhi_omega21_b1m2, weight21 * locAccWeight);
		    dHist_M_Phi_omega_Minus->Fill(locOmegaPiMinusMass22, locPhi_omega22_b1m2, weight22 * locAccWeight);

		    dHist_M_CosThetaH_plus->Fill(locOmegaPiPlusMass11, locCosTheta_H_11_b1p1, weight11 * locAccWeight);
		    dHist_M_CosThetaH_plus->Fill(locOmegaPiPlusMass12, locCosTheta_H_12_b1p2, weight12 * locAccWeight);
		    dHist_M_CosThetaH_plus->Fill(locOmegaPiPlusMass21, locCosTheta_H_21_b1p1, weight21 * locAccWeight);
		    dHist_M_CosThetaH_plus->Fill(locOmegaPiPlusMass22, locCosTheta_H_22_b1p2, weight22 * locAccWeight);
		    dHist_M_CosThetaH_minus->Fill(locOmegaPiMinusMass11, locCosTheta_H_11_b1m1, weight11 * locAccWeight);
		    dHist_M_CosThetaH_minus->Fill(locOmegaPiMinusMass21, locCosTheta_H_21_b1m2, weight21 * locAccWeight);
		    dHist_M_CosThetaH_minus->Fill(locOmegaPiMinusMass12, locCosTheta_H_12_b1m1, weight12 * locAccWeight);
		    dHist_M_CosThetaH_minus->Fill(locOmegaPiMinusMass22, locCosTheta_H_22_b1m2, weight22 * locAccWeight);

		    dHist_M_PhiH_plus->Fill(locOmegaPiPlusMass11, locPhi_H_11_b1p1, weight11 * locAccWeight);
		    dHist_M_PhiH_plus->Fill(locOmegaPiPlusMass12, locPhi_H_12_b1p2, weight12 * locAccWeight);
		    dHist_M_PhiH_plus->Fill(locOmegaPiPlusMass21, locPhi_H_21_b1p1, weight21 * locAccWeight);
		    dHist_M_PhiH_plus->Fill(locOmegaPiPlusMass22, locPhi_H_22_b1p2, weight22 * locAccWeight);
		    dHist_M_PhiH_minus->Fill(locOmegaPiMinusMass11, locPhi_H_11_b1m1, weight11 * locAccWeight);
		    dHist_M_PhiH_minus->Fill(locOmegaPiMinusMass21, locPhi_H_21_b1m2, weight21 * locAccWeight);
		    dHist_M_PhiH_minus->Fill(locOmegaPiMinusMass12, locPhi_H_12_b1m1, weight12 * locAccWeight);
		    dHist_M_PhiH_minus->Fill(locOmegaPiMinusMass22, locPhi_H_22_b1m2, weight22 * locAccWeight);



		    locUsedSoFar_DecayAngles.insert(locUsedThisCombo_DecayAngles);
		  }

		/********************************** HISTOGRAM DECAY ANGLES (OMEGARHO) ***************************************************/
		//Uniqueness Tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_DecayAngles_OR;
		locUsedThisCombo_DecayAngles_OR[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_DecayAngles_OR[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_DecayAngles_OR[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_DecayAngles_OR[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_DecayAngles_OR[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_DecayAngles_OR[Gamma].insert(locPhoton2NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_DecayAngles_OR.find(locUsedThisCombo_DecayAngles_OR) == locUsedSoFar_DecayAngles_OR.end())
		  {
		    dHist_CosThetaOmega_OR->Fill(locCosThetaOmega11, weight11 * locAccWeight);
		    dHist_CosThetaOmega_OR->Fill(locCosThetaOmega12, weight12 * locAccWeight);
		    dHist_CosThetaOmega_OR->Fill(locCosThetaOmega21, weight21 * locAccWeight);
		    dHist_CosThetaOmega_OR->Fill(locCosThetaOmega22, weight22 * locAccWeight);

		    dHist_PhiOmega_OR->Fill(locPhiOmega11, weight11 * locAccWeight);
		    dHist_PhiOmega_OR->Fill(locPhiOmega12, weight12 * locAccWeight);
		    dHist_PhiOmega_OR->Fill(locPhiOmega21, weight21 * locAccWeight);
		    dHist_PhiOmega_OR->Fill(locPhiOmega22, weight22 * locAccWeight);

		    dHist_CosThetaRho_OR->Fill(locCosThetaRho11, weight11 * locAccWeight);
		    dHist_CosThetaRho_OR->Fill(locCosThetaRho12, weight12 * locAccWeight);
		    dHist_CosThetaRho_OR->Fill(locCosThetaRho21, weight21 * locAccWeight);
		    dHist_CosThetaRho_OR->Fill(locCosThetaRho22, weight22 * locAccWeight);

		    dHist_PhiRho_OR->Fill(locPhiRho11, weight11 * locAccWeight);
		    dHist_PhiRho_OR->Fill(locPhiRho12, weight12 * locAccWeight);
		    dHist_PhiRho_OR->Fill(locPhiRho21, weight21 * locAccWeight);
		    dHist_PhiRho_OR->Fill(locPhiRho22, weight22 * locAccWeight);

		    dHist_CosThetaOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaOmega11, weight11 * locAccWeight);
		    dHist_CosThetaOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaOmega12, weight12 * locAccWeight);
		    dHist_CosThetaOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaOmega21, weight21 * locAccWeight);
		    dHist_CosThetaOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaOmega22, weight22 * locAccWeight);

		    dHist_PhiOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiOmega11, weight11 * locAccWeight);
		    dHist_PhiOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiOmega12, weight12 * locAccWeight);
		    dHist_PhiOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiOmega21, weight21 * locAccWeight);
		    dHist_PhiOmega_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiOmega22, weight22 * locAccWeight);

		    dHist_CosThetaRho_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaRho11, weight11 * locAccWeight);
		    dHist_CosThetaRho_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaRho12, weight12 * locAccWeight);
		    dHist_CosThetaRho_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaRho21, weight21 * locAccWeight);
		    dHist_CosThetaRho_ORVsOmega2PiMass->Fill(loc5PiMass, locCosThetaRho22, weight22 * locAccWeight);

		    dHist_PhiRho_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiRho11, weight11 * locAccWeight);
		    dHist_PhiRho_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiRho12, weight12 * locAccWeight);
		    dHist_PhiRho_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiRho21, weight21 * locAccWeight);
		    dHist_PhiRho_ORVsOmega2PiMass->Fill(loc5PiMass, locPhiRho22, weight22 * locAccWeight);

		    dHist_CosThetaOmega_ORVs3PiMass->Fill(loc3PiMass11, locCosThetaOmega11, locAccWeight);
		    dHist_CosThetaOmega_ORVs3PiMass->Fill(loc3PiMass12, locCosThetaOmega12, locAccWeight);
		    dHist_CosThetaOmega_ORVs3PiMass->Fill(loc3PiMass21, locCosThetaOmega21, locAccWeight);
		    dHist_CosThetaOmega_ORVs3PiMass->Fill(loc3PiMass22, locCosThetaOmega22, locAccWeight);

		    dHist_PhiOmega_ORVs3PiMass->Fill(loc3PiMass11, locPhiOmega11, locAccWeight);
		    dHist_PhiOmega_ORVs3PiMass->Fill(loc3PiMass12, locPhiOmega12, locAccWeight);
		    dHist_PhiOmega_ORVs3PiMass->Fill(loc3PiMass21, locPhiOmega21, locAccWeight);
		    dHist_PhiOmega_ORVs3PiMass->Fill(loc3PiMass22, locPhiOmega22, locAccWeight);

		    dHist_CosThetaRho_ORVs2PiMass->Fill(loc2PiMass11, locCosThetaRho11, locAccWeight);
		    dHist_CosThetaRho_ORVs2PiMass->Fill(loc2PiMass12, locCosThetaRho12, locAccWeight);
		    dHist_CosThetaRho_ORVs2PiMass->Fill(loc2PiMass21, locCosThetaRho21, locAccWeight);
		    dHist_CosThetaRho_ORVs2PiMass->Fill(loc2PiMass22, locCosThetaRho22, locAccWeight);

		    dHist_PhiRho_ORVs2PiMass->Fill(loc2PiMass11, locPhiRho11, locAccWeight);
		    dHist_PhiRho_ORVs2PiMass->Fill(loc2PiMass12, locPhiRho12, locAccWeight);
		    dHist_PhiRho_ORVs2PiMass->Fill(loc2PiMass21, locPhiRho21, locAccWeight);
		    dHist_PhiRho_ORVs2PiMass->Fill(loc2PiMass22, locPhiRho22, locAccWeight);

		    dHist_CosThetaRho_ORVs2PiMass_omegacut->Fill(loc2PiMass11, locCosThetaRho11, weight11 * locAccWeight);
		    dHist_CosThetaRho_ORVs2PiMass_omegacut->Fill(loc2PiMass12, locCosThetaRho12, weight12 * locAccWeight);
		    dHist_CosThetaRho_ORVs2PiMass_omegacut->Fill(loc2PiMass21, locCosThetaRho21, weight21 * locAccWeight);
		    dHist_CosThetaRho_ORVs2PiMass_omegacut->Fill(loc2PiMass22, locCosThetaRho22, weight22 * locAccWeight);

		    dHist_PhiRho_ORVs2PiMass_omegacut->Fill(loc2PiMass11, locPhiRho11, weight11 * locAccWeight);
		    dHist_PhiRho_ORVs2PiMass_omegacut->Fill(loc2PiMass12, locPhiRho12, weight12 * locAccWeight);
		    dHist_PhiRho_ORVs2PiMass_omegacut->Fill(loc2PiMass21, locPhiRho21, weight21 * locAccWeight);
		    dHist_PhiRho_ORVs2PiMass_omegacut->Fill(loc2PiMass22, locPhiRho22, weight22 * locAccWeight);

		    dHist_CosThetaOmegaVsPhiOmega_OR->Fill(locPhiOmega11, locCosThetaOmega11, weight11 * locAccWeight);
		    dHist_CosThetaOmegaVsPhiOmega_OR->Fill(locPhiOmega12, locCosThetaOmega12, weight12 * locAccWeight);
		    dHist_CosThetaOmegaVsPhiOmega_OR->Fill(locPhiOmega21, locCosThetaOmega21, weight21 * locAccWeight);
		    dHist_CosThetaOmegaVsPhiOmega_OR->Fill(locPhiOmega22, locCosThetaOmega22, weight22 * locAccWeight);

		    dHist_CosThetaRhoVsPhiRho_OR->Fill(locPhiRho11, locCosThetaRho11, weight11 * locAccWeight);
		    dHist_CosThetaRhoVsPhiRho_OR->Fill(locPhiRho12, locCosThetaRho12, weight12 * locAccWeight);
		    dHist_CosThetaRhoVsPhiRho_OR->Fill(locPhiRho21, locCosThetaRho21, weight21 * locAccWeight);
		    dHist_CosThetaRhoVsPhiRho_OR->Fill(locPhiRho22, locCosThetaRho22, weight22 * locAccWeight);

		    locUsedSoFar_DecayAngles_OR.insert(locUsedThisCombo_DecayAngles_OR);
		  }

		/*************************** HISTOGRAM DALITZ PLOTS ******************************************************/
		double locOmegaPiPlusMass2_11 = locOmegaPiPlusP4_11.M2();
		double locOmegaPiPlusMass2_12 = locOmegaPiPlusP4_12.M2();
		double locOmegaPiPlusMass2_21 = locOmegaPiPlusP4_21.M2();
		double locOmegaPiPlusMass2_22 = locOmegaPiPlusP4_22.M2();

		double locOmegaPiMinusMass2_11 = locOmegaPiMinusP4_11.M2();
		double locOmegaPiMinusMass2_12 = locOmegaPiMinusP4_12.M2();
		double locOmegaPiMinusMass2_21 = locOmegaPiMinusP4_21.M2();
		double locOmegaPiMinusMass2_22 = locOmegaPiMinusP4_22.M2();

		//Uniqueness tracking
		map<Particle_t, set<Int_t> > locUsedThisCombo_Dalitz;
		locUsedThisCombo_Dalitz[PiPlus].insert(locPiPlus1TrackID);
		locUsedThisCombo_Dalitz[PiMinus].insert(locPiMinus1TrackID);
		locUsedThisCombo_Dalitz[PiPlus].insert(locPiPlus2TrackID);
		locUsedThisCombo_Dalitz[PiMinus].insert(locPiMinus2TrackID);
		locUsedThisCombo_Dalitz[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_Dalitz[Gamma].insert(locPhoton2NeutralID);

		if(locUsedSoFar_Dalitz.find(locUsedThisCombo_Dalitz) == locUsedSoFar_Dalitz.end())
		  {
		    for(int i = 0; i < 10; i++){
		      double DalitzMassMin = 1.4 + 0.1*i;
		      double DalitzMassMax = 1.5 + 0.1*i;
		      if(loc5PiMass > DalitzMassMin && loc5PiMass < DalitzMassMax){
			dHist_Dalitz_OmegaPi[i]->Fill(locOmegaPiPlusMass2_11, locOmegaPiMinusMass2_11, weight11 * locAccWeight);
			dHist_Dalitz_OmegaPi[i]->Fill(locOmegaPiPlusMass2_12, locOmegaPiMinusMass2_12, weight12 * locAccWeight);
			dHist_Dalitz_OmegaPi[i]->Fill(locOmegaPiPlusMass2_21, locOmegaPiMinusMass2_21, weight21 * locAccWeight);
			dHist_Dalitz_OmegaPi[i]->Fill(locOmegaPiPlusMass2_22, locOmegaPiMinusMass2_22, weight22 * locAccWeight);
		      }
		    }
		    locUsedSoFar_Dalitz.insert(locUsedThisCombo_Dalitz);
		  }


		/***************************** SET UP ISCOMBOCUT FLAG *************************************/
		/*
		if(locProtonPiPlusMass1 < 1.4 && loct_Deltaplusplus1 < 0.5){
		  if(locOmegaPiMinusMass21 > 1.1 && locOmegaPiMinusMass21 < 1.35){
		    dComboWrapper->Set_IsComboCut(false);
		  }
		  if(locOmegaPiMinusMass22 > 1.1 && locOmegaPiMinusMass22 < 1.35){
		    dComboWrapper->Set_IsComboCut(false);
		  }
		}
		else if(locProtonPiPlusMass2 < 1.4 && loct_Deltaplusplus2 < 0.5){
		  if(locOmegaPiMinusMass11 > 1.1 && locOmegaPiMinusMass11 < 1.35){
		    dComboWrapper->Set_IsComboCut(false);
		  }
		  if(locOmegaPiMinusMass12 > 1.1 && locOmegaPiMinusMass12 < 1.35){
		    dComboWrapper->Set_IsComboCut(false);
		  }
		}
		else{
		  dComboWrapper->Set_IsComboCut(true);
		}
		*/


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

  				// loose pi0 mass
 				if(fabs(locPi0P4.M() - 0.135) > 0.015) continue; 

  				// loose -t cut for now
 				double loct = -1. * (locRecoilP4 - dTargetP4).M2();
 				if(loct > 1.0) continue;

  				// set weight from 2D 3pi mass correlation (see Chung et. al. 1975)...
				double loc2Dweight;
				double locLmin = 0.690;
				double locLmax = 0.735;
				double locomegamin = 0.760;
				double locomegamax = 0.805;
				double locHmin = 0.830;
				double locHmax = 0.875;

				if((loc3PiMass > locomegamin && loc3PiMass < locomegamax) || (loc3PiMass_alt > locomegamin && loc3PiMass_alt < locomegamax)) loc2Dweight = 1.0;
				else if((loc3PiMass > locLmin && loc3PiMass < locLmax) && (loc3PiMass_alt > locLmin && loc3PiMass_alt < locLmax)) loc2Dweight = -0.625;
				else if((loc3PiMass > locLmin && loc3PiMass < locLmax) && (loc3PiMass_alt > locHmin && loc3PiMass_alt < locHmax)) loc2Dweight = -0.625;
				else if((loc3PiMass > locHmin && loc3PiMass < locHmax) && (loc3PiMass_alt > locLmin && loc3PiMass_alt < locLmax)) loc2Dweight = -0.625;
				else if((loc3PiMass > locHmin && loc3PiMass < locHmax) && (loc3PiMass_alt > locHmin && loc3PiMass_alt < locHmax)) loc2Dweight = -0.625;
				else if((loc3PiMass > locLmin && loc3PiMass < locLmax) || (loc3PiMass > locHmin && loc3PiMass < locHmax)) loc2Dweight = -0.5;
				else if((loc3PiMass_alt > locLmin && loc3PiMass_alt < locLmax) || (loc3PiMass_alt > locHmin && loc3PiMass_alt < locHmax)) loc2Dweight = -0.5;
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
 				dFlatTreeInterface->Fill_Fundamental<Float_t>("Weight", loc2Dweight*locAccWeight);

				//remove random trigger background from phasespace MC
				if(dOption.Contains("phasespace")) {
				  Bool_t locIsGeneratorFlag = (dThrownBeam->Get_P4().E() == dComboBeamWrapper->Get_P4().E() && fabs(dThrownBeam->Get_X4().T() - dComboBeamWrapper->Get_X4().T()) < 2.004) ? kTRUE : kFALSE;
				  if( !(locIsGeneratorFlag || dComboBeamWrapper->Get_IsGenerator()) ) {
				    dComboWrapper->Set_IsComboCut(true);
				    continue;
				  }
				}

  				// set ordered final state P4 for filling flat tree
 				//FillAmpTools_FlatTree(locBeamP4, locFinalStateP4); 
 				FillAmpTools_FlatTree(locBeamP4, locFinalStateThrownP4);

  				Fill_FlatTree(); //for the active combo
 			}
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
