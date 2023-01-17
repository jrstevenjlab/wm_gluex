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
	  
	if(dTreeInterface->Get_Branch("NumCombos") != NULL) {

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

	  //CUT ON SHOWER QUALITY (commented for now to get a plot of full spectrum
	  //dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5)); 

	  //BEAM ENERGY
	  dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	  dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.0, 12.0));

	  //KINEMATICS
	  dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));
	
	  //INITIALIZE ACTIONS
	  //If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	  Initialize_Actions();
	}
	/******************************** USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/	
	// NORMAL HISTOGRAMS
	// angular distributions
	dHist_CosTheta   = new TH1F("CosTheta",   "; cos#theta",        400,   -1.,   1.);
	dHist_Phi 	 = new TH1F("Phi",        "; #phi (#circ)",     720, -180., 180.);
	dHist_CosTheta_H = new TH1F("CosTheta_H", "; cos#theta_{H}",    400,   -1.,   1.);
	dHist_Phi_H      = new TH1F("Phi_H",      "; #phi_{H} (#circ)", 720, -180., 180.);
	
	dHist_CosThetaVsPhi	= new TH2F("CosThetaVsPhi",     "; cos#theta; #phi (#circ)",         400, -1.0, 1.0, 720, -180.0, 180.0);
	dHist_CosTheta_HVsPhi_H = new TH2F("CosTheta_HVsPhi_H", "; cos#theta_{H}; #phi_{H} (#circ)", 400, -1.0, 1.0, 720, -180.0, 180.0);

	dHist_TVsCosTheta	= new TH2F("TVsCosTheta", "; -t (GeV)^{2}; cos#theta", 100, 0.0, 5.0, 400, -1.0, 1.0);
	dHist_TVsPhi		= new TH2F("TVsPhi", "; -t (GeV)^{2}; #phi (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);
	dHist_TVsCosTheta_H	= new TH2F("TVsCosTheta_H", "; -t (GeV)^{2}; cos#theta_{H}", 100, 0.0, 5.0, 400, -1.0, 1.0);
	dHist_TVsPhi_H		= new TH2F("TVsPhi_H", "; -t (GeV)^{2}; #phi_{H} (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);

	// other
	dHist_MM2  		  = new TH1F("MM2", 		    ";MM^{2} (GeV)^{2}", 	  200,-0.1,  0.1);
	dHist_MM2NoPhoton	  = new TH1F("MM2NoPhoton", 	    "No Photon;MM^{2} (GeV)^{2}", 200,-0.1,  0.1);
	dHist_KinFitChiSq         = new TH1F("KinFitChiSq",         ";#chi^{2}/NDF",              100, 0.0, 20.0);
	dHist_T 		  = new TH1F("T", 		    ";-t (GeV)^{2}",              100, 0.0,  5.0);
	dHist_EnergyUnusedShowers = new TH1F("EnergyUnusedShowers", ";E_{unused Showers} (GeV)",   20, 0.0,  2.0);
	dHist_ShowerQuality 	  = new TH1F("ShowerQuality", 	    ";Shower Quality",            100, 0.0,  1.0);
	// pi0 combos
	dHist_Pi0Mass	= new TH1F("Pi0Mass",  ";#gamma_{2}#gamma_{3} Mass (GeV)", 700, 0.06, 0.2);
	dHist_Pi02Mass	= new TH1F("Pi02Mass", ";#gamma_{1}#gamma_{2} Mass (GeV)", 500, 0.0,  1.0);
	dHist_Pi03Mass	= new TH1F("Pi03Mass", ";#gamma_{1}#gamma_{3} Mass (GeV)", 500, 0.0,  1.0);
	// 3pi combos
	dHist_3PiMass 	= new TH1F("3PiMass",  ";#pi^{+}#pi^{-}#gamma_{2}#gamma_{3} Mass (GeV)", 600, 0.1, 1.3);
	dHist_3Pi2Mass	= new TH1F("3Pi2Mass", ";#pi^{+}#pi^{-}#gamma_{1}#gamma_{2} Mass (GeV)", 600, 0.1, 1.3);
	dHist_3Pi3Mass	= new TH1F("3Pi3Mass", ";#pi^{+}#pi^{-}#gamma_{1}#gamma_{3} Mass (GeV)", 600, 0.1, 1.3);
	// other pion plots
	dHist_ProtonPiPlusMass	= new TH1F("ProtonPiPlusMass", ";p#pi^{+} Mass (GeV)", 		          750, 1.0, 2.5);	
	dHist_Pi0PiMinusMass	= new TH1F("Pi0PiMinusMass",   ";#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 600, 0.1, 1.3);

	dHist_ProtonPiPlusMassVsPi0PiMinusMass 		= new TH2F("ProtonPiPlusMassVsPi0PiMinusMass",    ";p#pi^{+} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 
								   750, 1.0, 2.5, 600, 0.1, 1.3);	
	dHist_ProtonPiPlusMassVsPi0PiMinusMass_PK 	= new TH2F("ProtonPiPlusMassVsPi0PiMinusMass_PK", "#omega Peak Only ;p#pi^{+} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 
								   750, 1.0, 2.5, 600, 0.1, 1.3);	
	dHist_ProtonPiPlusMassVsPi0PiMinusMass_SB 	= new TH2F("ProtonPiPlusMassVsPi0PiMinusMass_SB", "Sidebands Only;p#pi^{+} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 
								   750, 1.0, 2.5, 600, 0.1, 1.3);	

	dHist_Pi0GammaMass 			= new TH1F("Pi0GammaMass", ";#gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV)", 600, 0.1, 1.3);
	// Pi0GammaVsVariable distributions	
	dHist_Pi0GammaMassVsMM2			= new TH2F("Pi0GammaMassVsMM2", 		"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); MM^{2} (GeV)^{2}",
							   600, 0.1, 1.3, 200,   -0.1,   0.1);
	dHist_Pi0GammaMassVsMM2NoPhoton		= new TH2F("Pi0GammaMassVsMM2NoPhoton", 	"No Photon; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); MM^{2} (GeV)^{2}",
							   600, 0.1, 1.3, 200,   -0.1,   0.1);
	dHist_Pi0GammaMassVsCosTheta 		= new TH2F("Pi0GammaMassVsCosTheta", 		"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); cos#theta",
							   600, 0.1, 1.3, 400,   -1.0,   1.0);
	dHist_Pi0GammaMassVsPhi 		= new TH2F("Pi0GammaMassVsPhi",		        "; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #phi (#circ)",
							   600, 0.1, 1.3, 720, -180.0, 180.0);
	dHist_Pi0GammaMassVsCosTheta_H		= new TH2F("Pi0GammaMassVsCosTheta_H",          "; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); cos#theta_{H}",
							   600, 0.1, 1.3, 400,   -1.0,   1.0);
	dHist_Pi0GammaMassVsPhi_H 		= new TH2F("Pi0GammaMassVsPhi_H",		"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #phi_{H} (#circ)",
							   600, 0.1, 1.3, 720, -180.0, 180.0);
	dHist_Pi0GammaMassVsKinFitChiSq 	= new TH2F("Pi0GammaMassVsKinFitChiSq",	        "; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); KinFit #chi^{2}",
							   600, 0.1, 1.3, 100,    0.0,  20.0);
	dHist_Pi0GammaMassVsT 			= new TH2F("Pi0GammaMassVsT", 			"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); -t (GeV)^{2}", 
							   600, 0.1, 1.3, 100,    0.0,   5.0);
	dHist_Pi0GammaMassVsEnergyUnusedShowers = new TH2F("Pi0GammaMassVsEnergyUnusedShowers", "; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); E_{unused Showers} (GeV)", 
							   600, 0.1, 1.3, 20,     0.0,   2.0);
	dHist_Pi0GammaMassVsShowerQuality 	= new TH2F("Pi0GammaMassVsShowerQuality",       "; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); ShowerQuality (GeV)", 
							   600, 0.1, 1.3, 100,    0.0,   1.0);
	dHist_Pi0GammaMassVsPi0Mass 		= new TH2F("Pi0GammaMassVsPi0Mass", 		"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #gamma_{2}#gamma_{3} Mass (GeV)",
							   600, 0.1, 1.3, 700,    0.06,  0.2);
	dHist_Pi0GammaMassVs3PiMass 		= new TH2F("Pi0GammaMassVs3PiMass", 		"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #pi^{+}#pi^{-}#gamma_{2}#gamma_{3} Mass (GeV)",
							   600, 0.1, 1.3, 600,    0.1,   1.3);
	dHist_Pi0GammaMassVsProtonPiPlusMass    = new TH2F("Pi0GammaMassVsProtonPiPlusMass", 	"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); p#pi^{+} Mass (GeV); ",
							   750, 0.1, 1.3, 600,    1.0,   2.5);
	dHist_Pi0GammaMassVsPi0PiMinusMass      = new TH2F("Pi0GammaMassVsPi0PiMinusMass", 	"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)",
							   600, 0.1, 1.3, 600,    0.1,   1.3);
	dHist_Pi0GammaMassVsOmegaPiMinusMass    = new TH2F("Pi0GammaMassVsOmegaPiMinusMass", 	"; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)",
							   600, 0.1, 1.3, 600,    0.9,   2.1);	
		
	dHist_OmegaPiMinusMass 			= new TH1F("OmegaPiMinusMass", "; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 600, 0.9, 2.1);
	// OmegaPiMinusVsVariable distributions
	dHist_OmegaPiMinusMassVsPi0PiMinusMass  = new TH2F("OmegaPiMinusMassVsPi0PiMinusMass", 	"; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)",
							   600, 0.9, 2.1, 600,    0.1,   2.1);
	dHist_OmegaPiMinusMassVsCosTheta 	= new TH2F("OmegaPiMinusMassVsCosTheta",	"; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta_{H}",
							   600, 0.9, 2.1, 400,   -1.0,   1.0);
	dHist_OmegaPiMinusMassVsPhi 		= new TH2F("OmegaPiMinusMassVsPhi", 		"; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi (#circ)",
							   600, 0.9, 2.1, 720, -180.0, 180.0);
	dHist_OmegaPiMinusMassVsCosTheta_H 	= new TH2F("OmegaPiMinusMassVsCosTheta_H",	"; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta_{H}",
							   600, 0.9, 2.1, 400,   -1.0,   1.0);
	dHist_OmegaPiMinusMassVsPhi_H 		= new TH2F("OmegaPiMinusMassVsPhi_H", 		"; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi_{H} (#circ)",
							   600, 0.9, 2.1, 720, -180.0, 180.0);

	/* 
	   | *********************
	   | * THROWN TOPOLOGIES *
	   | *********************
	*/ 		
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);	// determine if analyzing simulated data       
	if(dIsMC) {
	  // histograms for counting the # of reactions themselves
	  dHistTop_PercentStart 		= new TH1F("PercentStart",        "PercentStart",        10, -0.5, 9.5);
	  dHistTop_PercentAngular 		= new TH1F("PercentAngular",      "PercentAngular",      10, -0.5, 9.5);
	  dHistTop_PercentChi2 			= new TH1F("PercentChi2",         "PercentChi2",         10, -0.5, 9.5);
	  dHistTop_PercentT	 		= new TH1F("PercentT",            "PercentT",            10, -0.5, 9.5);
	  dHistTop_PercentPhoton 		= new TH1F("PercentPhoton",       "PercentPhoton",       10, -0.5, 9.5);
	  dHistTop_PercentPi0 			= new TH1F("PercentPi0",          "PercentPi0",          10, -0.5, 9.5);
	  dHistTop_PercentPi0Combo 		= new TH1F("PercentPi0Combo",     "PercentPi0Combo",     10, -0.5, 9.5);
	  dHistTop_Percent3PiCombo 		= new TH1F("Percent3PiCombo",     "Percent3PiCombo",     10, -0.5, 9.5);
	  dHistTop_PercentProtonPiPlus 		= new TH1F("PercentProtonPiPlus", "PercentProtonPiPlus", 10, -0.5, 9.5);
	  dHistTop_PercentPi0PiMinus 		= new TH1F("PercentPi0PiMinus",   "PercentPi0PiMinus",   10, -0.5, 9.5);
	  dHistTop_PercentPi0Gamma 		= new TH1F("PercentPi0Gamma",     "PercentPi0Gamma",     10, -0.5, 9.5);
	  dHistTop_PercentOmega 		= new TH1F("PercentOmega",        "PercentOmega",        10, -0.5, 9.5);

	  // TOPOLOGY-MASS HISTOGRAMS
	  // Consider making a vector with the signal and bkgTT strings in them and then for-looping through the vector to assign the variables (avoiding double typing each hist)
	  // 	also make a custom string to write "sig" and "bkg" for each hist, then use Form("HistName_%s", mcType.Data())
	  TString signalTT	= "3#gamma#pi^{#plus}#pi^{#minus}p[#pi^{0},#omega]"; // signal Thrown Topology string
	  TString bkgTT 	= "background"; // custom string to just congregate all other topologies into "background"

	  dHistThrown_BeamEnergy    		= new TH1F("ThrownBeamEnergy",	  	";E_{beam} (GeV)", 100, 0.0,13.0);
	  dHistThrown_BeamEnergySig 		= new TH1F("ThrownBeamEnergy_McSignal", ";E_{beam} (Gev)", 100, 0.0, 13.0); // for tracking number of thrown signal events	  
	  
	  dHistTop_MM2[signalTT] 		= new TH1F("MM2_Sig",         "Signal; MM^{2} (GeV)^{2}",              200, -0.1, 0.1);
	  dHistTop_MM2[bkgTT]	 		= new TH1F("MM2_Bkg",         "Background; MM^{2} (GeV)^{2}",          200, -0.1, 0.1);
	  dHistTop_MM2NoPhoton[signalTT]	= new TH1F("MM2NoPhoton_Sig", "NoPhoton Signal; MM^{2} (GeV)^{2}",     200, -0.1, 0.1);
	  dHistTop_MM2NoPhoton[bkgTT]		= new TH1F("MM2NoPhoton_Bkg", "NoPhoton Background; MM^{2} (GeV)^{2}", 200, -0.1, 0.1);
	  
	  dHistTop_CosTheta[signalTT] 		= new TH1F("CosTheta_Sig",   "Signal; cos#theta",            400,   -1.,   1.);
	  dHistTop_CosTheta[bkgTT]    		= new TH1F("CosTheta_Bkg",   "Background; cos#theta",        400,   -1.,   1.);
	  dHistTop_Phi[signalTT] 		= new TH1F("Phi_Sig",        "Signal; #phi (#circ)",         720, -180., 180.);
	  dHistTop_Phi[bkgTT] 			= new TH1F("Phi_Bkg",        "Background; #phi (#circ)",     720, -180., 180.);
	  dHistTop_CosTheta_H[signalTT] 	= new TH1F("CosTheta_H_Sig", "Signal; cos#theta_{H}",        400,   -1.,   1.);
	  dHistTop_CosTheta_H[bkgTT]    	= new TH1F("CosTheta_H_Bkg", "Background; cos#theta_{H}",    400,   -1.,   1.);
	  dHistTop_Phi_H[signalTT] 		= new TH1F("Phi_H_Sig",      "Signal; #phi_{H} (#circ)",     720, -180., 180.);
	  dHistTop_Phi_H[bkgTT] 		= new TH1F("Phi_H_Bkg",      "Background; #phi_{H} (#circ)", 720, -180., 180.);

	  dHistThrown_CosTheta[signalTT]        = new TH1F("CosTheta_ThrownSig",   "ThrownSignal; cos#theta",            400,   -1.,   1.);
	  dHistThrown_Phi[signalTT] 		= new TH1F("Phi_ThrownSig",        "ThrownSignal; #phi (#circ)",         720, -180., 180.);
	  dHistThrown_CosTheta_H[signalTT] 	= new TH1F("CosTheta_H_ThrownSig", "ThrownSignal; cos#theta_{H}",        400,   -1.,   1.);
	  dHistThrown_Phi_H[signalTT] 		= new TH1F("Phi_H_ThrownSig",      "ThrownSignal; #phi_{H} (#circ)",     720, -180., 180.);
	  
	  dHistTop_CosThetaVsPhi[signalTT]	= new TH2F("CosThetaVsPhi_Sig",
							   "Signal; cos#theta; #phi (#circ)",
							   400, -1.0, 1.0, 720, -180.0, 180.0);
	  dHistTop_CosThetaVsPhi[bkgTT]		= new TH2F("CosThetaVsPhi_Bkg",
							   "Background; cos#theta; #phi (#circ)",
							   400, -1.0, 1.0, 720, -180.0, 180.0);
	  dHistTop_CosTheta_HVsPhi_H[signalTT]	= new TH2F("CosTheta_HVsPhi_H_Sig",
							   "Signal; cos#theta_{H}; #phi_{H} (#circ)",
							   400, -1.0, 1.0, 720, -180.0, 180.0);
	  dHistTop_CosTheta_HVsPhi_H[bkgTT]	= new TH2F("CosTheta_HVsPhi_H_Bkg",
							   "Background; cos#theta_{H}; #phi_{H} (#circ)",
							   400, -1.0, 1.0, 720, -180.0, 180.0);
	  dHistThrown_CosThetaVsPhi[signalTT]	= new TH2F("CosThetaVsPhi_ThrownSig",
							   "ThrownSignal; cos#theta; #phi (#circ)",
							   400, -1.0, 1.0, 720, -180.0, 180.0);
	  dHistThrown_CosTheta_HVsPhi_H[signalTT]  = new TH2F("CosTheta_HVsPhi_H_ThrownSig",
							      "ThrownSignal; cos#theta_{H}; #phi_{H} (#circ)",
							      400, -1.0, 1.0, 720, -180.0, 180.0);

	  dHistTop_TVsCosTheta[signalTT]	= new TH2F("TVsCosTheta_Sig", 
							   "Signal; -t (GeV)^{2}; cos#theta", 100, 0.0, 5.0, 400, -1.0, 1.0);
	  dHistTop_TVsCosTheta[bkgTT]		= new TH2F("TVsCosTheta_Bkg", 
							   "Background; -t (GeV)^{2}; cos#theta", 100, 0.0, 5.0, 400, -1.0, 1.0);
	  dHistTop_TVsPhi[signalTT]		= new TH2F("TVsPhi_Sig", 
							   "Signal; -t (GeV)^{2}; #phi (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);
	  dHistTop_TVsPhi[bkglTT]		= new TH2F("TVsPhi_Bkg", 
							   "Background; -t (GeV)^{2}; #phi (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);
	  dHistTop_TVsCosTheta_H[signalTT]	= new TH2F("TVsCosTheta_H_Sig", 
							   "Signal; -t (GeV)^{2}; cos#theta_{H}", 100, 0.0, 5.0, 400, -1.0, 1.0);
	  dHistTop_TVsCosTheta_H[bkgTT]	 	= new TH2F("TVsCosTheta_H_Bkg", 
							   "Background; -t (GeV)^{2}; cos#theta_{H}", 100, 0.0, 5.0, 400, -1.0, 1.0);
	  dHistTop_TVsPhi_H[signalTT]		= new TH2F("TVsPhi_H_Sig", 
							   "Signal; -t (GeV)^{2}; #phi_{H} (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);
	  dHistTop_TVsPhi_H[bkgTT]		= new TH2F("TVsPhi_H_Bkg", 
							   "Background; -t (GeV)^{2}; #phi_{H} (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);

	  dHistThrown_TVsCosTheta[signalTT]	= new TH2F("TVsCosTheta_ThrownSig", 
							   "ThrownSignal; -t (GeV)^{2}; cos#theta", 100, 0.0, 5.0, 400, -1.0, 1.0);
	  dHistThrown_TVsPhi[signalTT]		= new TH2F("TVsPhi_ThrownSig", 
							   "ThrownSignal; -t (GeV)^{2}; #phi (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);
	  dHistThrown_TVsCosTheta_H[signalTT]	= new TH2F("TVsCosTheta_H_ThrownSig", 
							   "ThrownSignal; -t (GeV)^{2}; cos#theta_{H}", 100, 0.0, 5.0, 400, -1.0, 1.0);
	  dHistThrown_TVsPhi_H[signalTT]        = new TH2F("TVsPhi_H_ThrownSig", 
							   "ThrownSignal; -t (GeV)^{2}; #phi_{H} (#circ)", 100, 0,0, 5.0, 720, -180.0, 180.0);
	  
	  dHistTop_KinFitChiSq[signalTT] 	= new TH1F("KinFitChiSq_Sig",         "Signal; #chi^{2}/NDF",                  100, 0.0, 20.0);
	  dHistTop_KinFitChiSq[bkgTT] 		= new TH1F("KinFitChiSq_Bkg",         "Background; #chi^{2}/NDF",              100, 0.0, 20.0);
	  dHistTop_T[signalTT] 			= new TH1F("T_Sig",                   "Signal; -t (GeV)^{2}",                  100, 0.0,  5.0);
	  dHistTop_T[bkgTT] 			= new TH1F("T_Bkg",                   "Background; -t (GeV)^{2}",              100, 0.0,  5.0);
	  dHistTop_EnergyUnusedShowers[signalTT]= new TH1F("EnergyUnusedShowers_Sig", "Signal; E_{unused Showers} (GeV)",      20, 0.0,  2.0);
	  dHistTop_EnergyUnusedShowers[bkgTT] 	= new TH1F("EnergyUnusedShowers_Bkg", "Background; E_{unused Showers} (GeV)",  20, 0.0,  2.0);
	  dHistTop_ShowerQuality[signalTT] 	= new TH1F("ShowerQuality_Sig",       "Signal; Shower Quality",                100, 0.0,  1.0);
	  dHistTop_ShowerQuality[bkgTT] 	= new TH1F("ShowerQuality_Bkg",       "Background; Shower Quality",            100, 0.0,  1.0);

	  dHistTop_Pi0Mass[signalTT] 		= new TH1F("Pi0Mass_Sig",  "Signal; #gamma_{2}#gamma_{3} Mass (GeV)",     700, 0.06, 0.20);
	  dHistTop_Pi0Mass[bkgTT]		= new TH1F("Pi0Mass_Bkg",  "Background; #gamma_{2}#gamma_{3} Mass (GeV)", 700, 0.06, 0.20);
	  dHistTop_Pi02Mass[signalTT] 		= new TH1F("Pi02Mass_Sig", "Signal; #gamma_{1}#gamma_{2} Mass (GeV)",     500, 0.0,  1.0);
	  dHistTop_Pi02Mass[bkgTT]		= new TH1F("Pi02Mass_Bkg", "Background; #gamma_{1}#gamma_{2} Mass (GeV)", 500, 0.0,  1.0);
	  dHistTop_Pi03Mass[signalTT] 		= new TH1F("Pi03Mass_Sig", "Signal; #gamma_{1}#gamma_{3} Mass (GeV)",     500, 0.0,  1.0);
	  dHistTop_Pi03Mass[bkgTT]		= new TH1F("Pi03Mass_Bkg", "Background; #gamma_{1}#gamma_{3} Mass (GeV)", 500, 0.0,  1.0);
	    
	  dHistTop_3PiMass[signalTT] 		= new TH1F("3PiMass_Sig",  "Signal; #pi^{+}#pi^{-}#gamma_{2}#gamma_{3} Mass (GeV)",     600, 0.1, 1.3);
	  dHistTop_3PiMass[bkgTT] 		= new TH1F("3PiMass_Bkg",  "Background; #pi^{+}#pi^{-}#gamma_{2}#gamma_{3} Mass (GeV)", 600, 0.1, 1.3);
	  dHistTop_3Pi2Mass[signalTT] 		= new TH1F("3Pi2Mass_Sig", "Signal; #pi^{+}#pi^{-}#gamma_{1}#gamma_{2} Mass (GeV)",     600, 0.1, 1.3);
	  dHistTop_3Pi2Mass[bkgTT] 		= new TH1F("3Pi2Mass_Bkg", "Background; #pi^{+}#pi^{-}#gamma_{1}#gamma_{2} Mass (GeV)", 600, 0.1, 1.3);
	  dHistTop_3Pi3Mass[signalTT] 		= new TH1F("3Pi3Mass_Sig", "Signal; #pi^{+}#pi^{-}#gamma_{1}#gamma_{3} Mass (GeV)",     600, 0.1, 1.3);
	  dHistTop_3Pi3Mass[bkgTT] 		= new TH1F("3Pi3Mass_Bkg", "Background; #pi^{+}#pi^{-}#gamma_{1}#gamma_{3} Mass (GeV)", 600, 0.1, 1.3);	  
	  
	  dHistTop_ProtonPiPlusMass[signalTT] 	= new TH1F("ProtonPiPlusMass_Sig", "Signal; p#pi^{+} Mass (GeV)",     750, 1.0, 2.5);	  
	  dHistTop_ProtonPiPlusMass[bkgTT] 	= new TH1F("ProtonPiPlusMass_Bkg", "Background; p#pi^{+} Mass (GeV)", 750, 1.0, 2.5);	  
	  dHistTop_Pi0PiMinusMass[signalTT] 	= new TH1F("Pi0PiMinusMass_Sig",   "Signal; #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)",     600, 0.1, 1.3);
	  dHistTop_Pi0PiMinusMass[bkgTT] 	= new TH1F("Pi0PiMinusMass_Bkg",   "Background; #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 600, 0.1, 1.3);

	  dHistTop_ProtonPiPlusMassVsPi0PiMinusMass[signalTT] 	= new TH2F("ProtonPiPlusMassVsPi0PiMinusMass_Sig", "Signal; p#pi^{+} Mass (GeV) ;#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 
									   750, 1.0, 2.5, 600, 0.1, 1.3);
	  dHistTop_ProtonPiPlusMassVsPi0PiMinusMass[bkgTT] 	= new TH2F("ProtonPiPlusMassVsPi0PiMinusMass_Bkg", "Background; p#pi^{+} Mass (GeV) ;#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 
									   750, 1.0, 2.5, 600, 0.1, 1.3);
	  
	  dHistTop_OmegaPiMinusMass[signalTT]			= new TH1F("OmegaPiMinusMass_Sig", "Signal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)",     600, 0.9, 2.1);
	  dHistTop_OmegaPiMinusMass[bkgTT]			= new TH1F("OmegaPiMinusMass_Bkg", "Background; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", 600, 0.9, 2.1);	  

	  dHistTop_OmegaPiMinusMassVsPi0PiMinusMass[signalTT] 	= new TH2F("OmegaPiMinusMassVsPi0PiMinusMass_Sig", 
									   "Signal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)",
									   600, 0.9, 2.1, 600,    0.1,  2.1);
	  dHistTop_OmegaPiMinusMassVsPi0PiMinusMass[bkgTT] 	= new TH2F("OmegaPiMinusMassVsPi0PiMinusMass_Bkg", 
									   "Background; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)",
									   600, 0.9, 2.1, 600,    0.1,  2.1);
	  dHistTop_OmegaPiMinusMassVsCosTheta[signalTT] 	= new TH2F("OmegaPiMinusMassVsCosTheta_Sig",
									   "Signal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta",
									   600, 0.9, 2.1, 400,   -1.0,  1.0);
	  dHistTop_OmegaPiMinusMassVsCosTheta[bkgTT] 		= new TH2F("OmegaPiMinusMassVsCosTheta_Bkg", 
									   "Background; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta",
									   600, 0.9, 2.1, 400,   -1.0,  1.0);
	  dHistTop_OmegaPiMinusMassVsPhi[signalTT] 		= new TH2F("OmegaPiMinusMassVsPhi_Sig",
									   "Signal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi (#circ)",
									   600, 0.9, 2.1, 720, -180.0, 180.0);
	  dHistTop_OmegaPiMinusMassVsPhi[bkgTT] 		= new TH2F("OmegaPiMinusMassVsPhi_Bkg",
									   "Background; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi (#circ)",
									   600, 0.9, 2.1, 720, -180.0, 180.0);
	  dHistTop_OmegaPiMinusMassVsCosTheta_H[signalTT] 	= new TH2F("OmegaPiMinusMassVsCosTheta_H_Sig",
									   "Signal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta_{H}",
									   600, 0.9, 2.1, 400,   -1.0,  1.0);
	  dHistTop_OmegaPiMinusMassVsCosTheta_H[bkgTT] 		= new TH2F("OmegaPiMinusMassVsCosTheta_H_Bkg", 
									   "Background; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta_{H}",
									   600, 0.9, 2.1, 400,   -1.0,  1.0);
	  dHistTop_OmegaPiMinusMassVsPhi_H[signalTT] 		= new TH2F("OmegaPiMinusMassVsPhi_H_Sig",
									   "Signal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi_{H} (#circ)",
									   600, 0.9, 2.1, 720, -180.0, 180.0);
	  dHistTop_OmegaPiMinusMassVsPhi_H[bkgTT] 		= new TH2F("OmegaPiMinusMassVsPhi_H_Bkg",
									   "Background; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi_{H} (#circ)",
									   600, 0.9, 2.1, 720, -180.0, 180.0);
	  dHistThrown_OmegaPiMinusMassVsCosTheta[signalTT] 	= new TH2F("OmegaPiMinusMassVsCosTheta_ThrownSig",
									   "ThrownSignal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta",
									   600, 0.9, 2.1, 400,   -1.0,  1.0);
	  dHistThrown_OmegaPiMinusMassVsPhi[signalTT] 		= new TH2F("OmegaPiMinusMassVsPhi_ThrownSig",
									   "ThrownSignal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi (#circ)",
									   600, 0.9, 2.1, 720, -180.0, 180.0);
	  dHistThrown_OmegaPiMinusMassVsCosTheta_H[signalTT] 	= new TH2F("OmegaPiMinusMassVsCosTheta_H_ThrownSig",
									   "ThrownSignal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); cos#theta_{H}",
									   600, 0.9, 2.1, 400,   -1.0,  1.0);
	  dHistThrown_OmegaPiMinusMassVsPhi_H[signalTT] 	= new TH2F("OmegaPiMinusMassVsPhi_H_ThrownSig",
									   "ThrownSignal; #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV); #phi_{H} (#circ)",
									   600, 0.9, 2.1, 720, -180.0, 180.0);
	  	 		    
	  // Pi0Gamma plots track what happens before sideband subtraction. Here I can see distribtuions of specific variables better
	  dHistThrown_Pi0GammaMass[signalTT] 	= new TH1F("Pi0GammaMass_ThrownSig", "ThrownSignal; #pi^{0}#gamma Mass (GeV)", 600, 0.1, 1.3);
	  vector<TString> locVecTopologies;
	  locVecTopologies.insert(locVecTopologies.end(), {
	      "4#gamma#pi^{#plus}#pi^{#minus}p[2#pi^{0}]", 	    // rho^- decay
		"2#gamma#pi^{#plus}#pi^{#minus}p[#pi^{0}]",	    // extra "photon" is either stray or misidentified pion/kaon
		"3#gamma#pi^{#plus}#pi^{#minus}p[#pi^{0},#omega]",  // signal
		"4#gamma#pi^{#plus}#pi^{#minus}p[2#pi^{0},#omega]", // neutral b1 decay
		"background" 					    // the "background" plot for these will just be all other topologies not listed here (should be ~1% of total)
		});
	  for(uint i=0; i<locVecTopologies.size(); i++) {	    	    	    
	    dHistTop_Pi0GammaMass[locVecTopologies[i]] 				= new TH1F(Form("Pi0GammaMass_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV)", locVecTopologies[i].Data()),
											   600, 0.1, 1.3);
	    dHistTop_Pi0GammaMassVsMM2[locVecTopologies[i]]			= new TH2F(Form("Pi0GammaMassVsMM2_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); MM^{2} (GeV)^{2}", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 200,   -0.1,   0.1);
	    dHistTop_Pi0GammaMassVsMM2NoPhoton[locVecTopologies[i]]		= new TH2F(Form("Pi0GammaMassVsMM2NoPhoton_TT%d", i),
											   Form("No Photon TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); MM^{2} (GeV)^{2}", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 200,   -0.1,   0.1);
	    dHistTop_Pi0GammaMassVsCosTheta[locVecTopologies[i]]        	= new TH2F(Form("Pi0GammaMassVsCosTheta_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); cos#theta", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 400,   -1.0,   1.0);
	    dHistTop_Pi0GammaMassVsPhi[locVecTopologies[i]] 			= new TH2F(Form("Pi0GammaMassVsPhi_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #phi (#circ)", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 720, -180.0, 180.0);
	    dHistTop_Pi0GammaMassVsCosTheta_H[locVecTopologies[i]]        	= new TH2F(Form("Pi0GammaMassVsCosTheta_H_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); cos#theta_{H}", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 400,   -1.0,   1.0);
	    dHistTop_Pi0GammaMassVsPhi_H[locVecTopologies[i]] 			= new TH2F(Form("Pi0GammaMassVsPhi_H_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #phi_{H} (#circ)", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 720, -180.0, 180.0);
	    dHistTop_Pi0GammaMassVsKinFitChiSq[locVecTopologies[i]] 		= new TH2F(Form("Pi0GammaMassVsKinFitChiSq_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); KinFit #chi^{2}", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 100,    0.0,  20.0);
	    dHistTop_Pi0GammaMassVsT[locVecTopologies[i]] 			= new TH2F(Form("Pi0GammaMassVsT_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); -t (GeV)^{2}",  locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 100,    0.0,   5.0);
	    dHistTop_Pi0GammaMassVsEnergyUnusedShowers[locVecTopologies[i]] 	= new TH2F(Form("Pi0GammaMassVsEnergyUnusedShowers_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); E_{unused Showers} (GeV)",  locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 20,     0.0,   2.0);
	    dHistTop_Pi0GammaMassVsShowerQuality[locVecTopologies[i]] 		= new TH2F(Form("Pi0GammaMassVsShowerQuality_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); ShowerQuality (GeV)",  locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 100,    0.0,   1.0);
	    dHistTop_Pi0GammaMassVsPi0Mass[locVecTopologies[i]] 		= new TH2F(Form("Pi0GammaMassVsPi0Mass_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #gamma_{2}#gamma_{3} Mass (GeV)", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 700,    0.06,  0.2);
	    dHistTop_Pi0GammaMassVs3PiMass[locVecTopologies[i]] 		= new TH2F(Form("Pi0GammaMassVs3PiMass_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #pi^{+}#pi^{-}#gamma_{2}#gamma_{3} Mass (GeV)", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 600,    0.1,   1.3);
	    dHistTop_Pi0GammaMassVsProtonPiPlusMass[locVecTopologies[i]]    	= new TH2F(Form("Pi0GammaMassVsProtonPiPlusMass_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); p#pi^{+} Mass (GeV); ", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 750,    1.0,   2.5);
	    dHistTop_Pi0GammaMassVsPi0PiMinusMass[locVecTopologies[i]]      	= new TH2F(Form("Pi0GammaMassVsPi0PiMinusMass_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", locVecTopologies[i].Data()),
											   600, 0.1, 1.3, 600,    0.1,   1.3);
	    dHistTop_Pi0GammaMassVsOmegaPiMinusMass[locVecTopologies[i]] 	= new TH2F(Form("Pi0GammaMassVsOmegaPiMinusMass_TT%d", i),
											   Form("TT: %s; #gamma_{1}#gamma_{2}#gamma_{3} Mass (GeV); #gamma_{1}#gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", locVecTopologies[i].Data()), 
											   600, 0.1, 1.3, 600,    0.9,   2.1);
	    
	    dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_PK[locVecTopologies[i]]   = new TH2F(Form("ProtonPiPlusMassVsPi0PiMinusMass_PK_TT%d", i),
											   Form("#omega Peak TT: %s; p#pi^{+} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", locVecTopologies[i].Data()),
											   750, 1.0, 2.5, 600,    0.1,   1.3);
	    dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_SB[locVecTopologies[i]]	= new TH2F(Form("ProtonPiPlusMassVsPi0PiMinusMass_SB_TT%d", i),
											   Form("Sideband TT: %s; p#pi^{+} Mass (GeV); #gamma_{2}#gamma_{3}#pi^{-} Mass (GeV)", locVecTopologies[i].Data()),
											   750, 1.0, 2.5, 600,    0.1,   1.3);
	  }

	  if(isVarStudy) { // when changing this, make sure it will still fit with varStudy.C code
	    dHistTop_OmegaPiMinusMassVsPi0Gamma[signalTT] 		= new TH2F("hTop_OmegaPiMinusMassVsPi0GammaWidth_Sig",       "Signal",     600, 0.9, 2.1, 1200, 0.1,  1.3);
	    dHistTop_OmegaPiMinusMassVsPi0Gamma[bkgTT] 			= new TH2F("hTop_OmegaPiMinusMassVsPi0GammaWidth_Bkg",       "Background", 600, 0.9, 2.1, 1200, 0.1,  1.3);
	    dHistTop_OmegaPiMinusMassVsKinFitChiSq[signalTT] 		= new TH2F("hTop_OmegaPiMinusMassVsKinFitChiSq_Sig",         "Signal",     600, 0.9, 2.1, 1000, 0.0, 10.0);
	    dHistTop_OmegaPiMinusMassVsKinFitChiSq[bkgTT] 		= new TH2F("hTop_OmegaPiMinusMassVsKinFitChiSq_Bkg",         "Background", 600, 0.9, 2.1, 1000, 0.0, 10.0);
	    dHistTop_OmegaPiMinusMassVsEnergyUnusedShowers[signalTT]    = new TH2F("hTop_OmegaPiMinusMassVsEnergyUnusedShowers_Sig", "Signal",     600, 0.9, 2.1, 2000, 0.0,  2.0);
	    dHistTop_OmegaPiMinusMassVsEnergyUnusedShowers[bkgTT] 	= new TH2F("hTop_OmegaPiMinusMassVsEnergyUnusedShowers_Bkg", "Background", 600, 0.9, 2.1, 2000, 0.0,  2.0);
	    dHistTop_OmegaPiMinusMassVsPi0Width[signalTT] 		= new TH2F("hTop_OmegaPiMinusMassVsPi0Width_Sig", 	     "Signal",     600, 0.9, 2.1, 1300, 0.07, 0.20);
	    dHistTop_OmegaPiMinusMassVsPi0Width[bkgTT] 			= new TH2F("hTop_OmegaPiMinusMassVsPi0Width_Bkg", 	     "Background", 600, 0.9, 2.1, 1300, 0.07, 0.20);
	    dHistTop_OmegaPiMinusMassVsShowerQuality[signalTT] 		= new TH2F("hTop_OmegaPiMinusMassVsShowerQuality_Sig", 	     "Signal",     600, 0.9, 2.1, 1000, 0.0,  1.0);
	    dHistTop_OmegaPiMinusMassVsShowerQuality[bkgTT] 		= new TH2F("hTop_OmegaPiMinusMassVsShowerQuality_Bkg",       "Background", 600, 0.9, 2.1, 1000, 0.0,  1.0);	      
	  }
	} // ends if(dIsMC) statement	
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
	if(locRunNumber != dPreviousRunNumber) {
	  dIsPolarizedFlag 	= dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
	  dPreviousRunNumber 	= locRunNumber;
	}
	
	/************************************************* PARSE THROWN TOPOLOGY ***************************************/
	TString locThrownTopology 	= (dIsMC) ? Get_ThrownTopologyString() : "";	// (MC) Get the contributing reaction for this tree entry
	TString signalTT      		= "3#gamma#pi^{#plus}#pi^{#minus}p[#pi^{0},#omega]";
	TString bkgTT 			= "background";

	/******************************************* THROWN EVENTS ********************************************/

	if(Get_NumThrown() > 0) { // get thrown particles for this one event
	  TLorentzVector thrownBeamP4   = dThrownBeam->Get_P4();
	  TLorentzVector thrownTargetP4 = dTargetP4;
	  dHistThrown_BeamEnergy->Fill(thrownBeamP4.E());
	  
	  if(locThrownTopology.EqualTo(signalTT) && dIsMC) {
	    dHistThrown_BeamEnergySig->Fill(thrownBeamP4.E()); // integrate this and other beam hist from 8 to 12 to get true comparison to MC and Signal	    
	    TLorentzVector thrownPhoton1P4; 
	    TLorentzVector thrownPhoton2P4; 
	    TLorentzVector thrownPhoton3P4; 

	    TLorentzVector thrownPi0P4;
	    TLorentzVector thrownPiPlusP4;
	    TLorentzVector thrownPiMinusP4;
	    TLorentzVector thrownProtonP4;
	    TLorentzVector thrownOmegaP4;
	    TLorentzVector emptyVec; // empty vec for comparison uses

	    for(UInt_t i=0; i<Get_NumThrown(); i++) { // get momenta of each thrown particle (for our signal)	      	      
	      dThrownWrapper->Set_ArrayIndex(i);
	      Int_t thrownPID = dThrownWrapper->Get_PID();
	      TLorentzVector thrownWrapperP4 = dThrownWrapper->Get_P4();
	      
	      if(thrownPID ==  1 && thrownPhoton1P4 == emptyVec)      {thrownPhoton1P4  = thrownWrapperP4;}
	      else if(thrownPID ==  1 && thrownPhoton2P4 == emptyVec) {thrownPhoton2P4  = thrownWrapperP4;}
	      else if(thrownPID ==  1 && thrownPhoton3P4 == emptyVec) {thrownPhoton3P4  = thrownWrapperP4;}
	      
	      if(thrownPID ==  7) {thrownPi0P4     = thrownWrapperP4;}
	      if(thrownPID ==  8) {thrownPiPlusP4  = thrownWrapperP4;}
	      if(thrownPID ==  9) {thrownPiMinusP4 = thrownWrapperP4;}
	      if(thrownPID == 14) {thrownProtonP4  = thrownWrapperP4;}	      	      
	      if(thrownPID == 33) {thrownOmegaP4   = thrownWrapperP4;}
	    }
	    // identify the radiated photon by checking if other two photons add to pi0 mass
	    // 	photon1 is pretty ensured to be the radiated one, but its a nice cross-check
	    TLorentzVector thrownPhotonRadiatedP4;
	    if(fabs((thrownPhoton1P4 + thrownPhoton2P4).M() -.135) < .015)      {thrownPhotonRadiatedP4 = thrownPhoton3P4;}
	    else if(fabs((thrownPhoton1P4 + thrownPhoton3P4).M() -.135) < .015) {thrownPhotonRadiatedP4 = thrownPhoton2P4;}
	    else if(fabs((thrownPhoton2P4 + thrownPhoton3P4).M() -.135) < .015) {thrownPhotonRadiatedP4 = thrownPhoton1P4;}

	    // this handles reconstructed bggen, where only photons are "detected" and a pi0 or omega is not generated
	    if(thrownOmegaP4==emptyVec) 
	      thrownOmegaP4 = thrownPhoton1P4 + thrownPhoton2P4 + thrownPhoton3P4;
	    
	    // setup variables for cuts/plots
	    double thrownT 			= ((thrownProtonP4+thrownPiPlusP4) - thrownTargetP4).M2();
	    TLorentzVector thrownDeltaPlusPlus 	= thrownPiPlusP4 + thrownProtonP4;	    
	    TLorentzVector thrownOmegaPiMinusP4 = thrownOmegaP4 + thrownPiMinusP4;

	    // Most cuts are ID-ing particles, but these restrict us to the same kinematics
	    if( fabs(thrownT) < 0.5 && fabs(thrownBeamP4.E() - 10.0) < 2.0 && thrownDeltaPlusPlus.M() < 1.3) { 
	      /* **********************************
	       | * THROWN DECAY ANGLES  AND PLOTS *
	       | **********************************
	      */	      
	      TVector3 thrownBeamProtonBoost                       = (thrownBeamP4 + thrownTargetP4).BoostVector();
	      TLorentzVector thrownBeamP4_beamProtonRest           = thrownBeamP4; thrownBeamP4_beamProtonRest.Boost(-1.0*thrownBeamProtonBoost);
	      TLorentzVector thrownOmegaPiMinusP4_beamProtonRest   = thrownOmegaPiMinusP4; thrownOmegaPiMinusP4_beamProtonRest.Boost(-1.0*thrownBeamProtonBoost);
	      TLorentzVector thrownOmegaP4_beamProtonRest          = thrownOmegaP4; thrownOmegaP4_beamProtonRest.Boost(-1.0*thrownBeamProtonBoost);
	      TLorentzVector thrownPhotonRadiatedP4_beamProtonRest = thrownPhotonRadiatedP4; thrownPhotonRadiatedP4_beamProtonRest.Boost(-1.0*thrownBeamProtonBoost);

	      // define production plane unit vectors
	      TVector3 thrownUnitK = thrownBeamP4_beamProtonRest.Vect().Unit();
	      TVector3 thrownUnitZ = thrownOmegaPiMinusP4_beamProtonRest.Vect().Unit();
	      TVector3 thrownUnitY = thrownUnitK.Cross(thrownUnitZ).Unit();
	      TVector3 thrownUnitX = thrownUnitY.Cross(thrownUnitZ).Unit();
	      
	      // Get b1 boost vector (from production plane) and boost momenta into it
	      TVector3 thrownOmegaPiMinusBoost             = thrownOmegaPiMinusP4_beamProtonRest.BoostVector();
	      TLorentzVector thrownOmegaP4_b1Rest          = thrownOmegaP4_beamProtonRest; thrownOmegaP4_b1Rest.Boost(-1.0*thrownOmegaPiMinusBoost);
	      TLorentzVector thrownPhotonRadiatedP4_b1Rest = thrownPhotonRadiatedP4_beamProtonRest; thrownPhotonRadiatedP4_b1Rest.Boost(-1.0*thrownOmegaPiMinusBoost);
	      
	      // Define omega decay angles
	      TVector3 thrownOmegaP3_b1Rest = thrownOmegaP4_b1Rest.Vect();
	      double thrownCosTheta         = thrownOmegaP3_b1Rest.Dot(thrownUnitZ) / thrownOmegaP3_b1Rest.Mag();
	      double thrownPhi              = 180. * TMath::ATan2(thrownOmegaP3_b1Rest.Dot(thrownUnitY), thrownOmegaP3_b1Rest.Dot(thrownUnitX)) / TMath::Pi();
	      
	      // Define decay plane unit vectors for helicity frame
	      TVector3 thrownUnitZ_H = thrownOmegaP3_b1Rest.Unit();
	      TVector3 thrownUnitY_H = thrownUnitZ.Cross(thrownUnitZ_H).Unit();
	      TVector3 thrownUnitX_H = thrownUnitY_H.Cross(thrownUnitZ_H);

	      // Get Omega boost vector (to b1 rf) and boost radiated photon
	      TVector3 thrownOmegaBoost                       = thrownOmegaP4_b1Rest.BoostVector();
	      TLorentzVector thrownPhotonRadiatedP4_omegaRest = thrownPhotonRadiatedP4_b1Rest; thrownPhotonRadiatedP4_omegaRest.Boost(-1.0*thrownOmegaBoost);
	      
	      // Get helicity angles wrt radiated photon (in omega rf)
	      TVector3 thrownPhotonRadiatedP3_omegaRest = thrownPhotonRadiatedP4_omegaRest.Vect();
	      double thrownCosTheta_H                   = thrownPhotonRadiatedP3_omegaRest.Dot(thrownUnitZ_H) / thrownPhotonRadiatedP3_omegaRest.Mag();
	      double thrownPhi_H                        = 180. * TMath::ATan2(thrownPhotonRadiatedP3_omegaRest.Dot(thrownUnitY_H) , thrownPhotonRadiatedP3_omegaRest.Dot(thrownUnitX_H)) / TMath::Pi();

	      // fill normal hists
	      dHistThrown_Pi0GammaMass[locThrownTopology]->Fill(thrownOmegaP4.M());

	      // Fill histograms with angles
	      dHistThrown_CosTheta[locThrownTopology]->Fill(thrownCosTheta);
	      dHistThrown_Phi[locThrownTopology]->Fill(thrownPhi);
	      dHistThrown_CosTheta_H[locThrownTopology]->Fill(thrownCosTheta_H);
	      dHistThrown_Phi_H[locThrownTopology]->Fill(thrownPhi_H);

	      // fill 2d plots
	      dHistThrown_CosThetaVsPhi[locThrownTopology]->Fill(thrownCosTheta, thrownPhi);
	      dHistThrown_CosTheta_HVsPhi_H[locThrownTopology]->Fill(thrownCosTheta_H, thrownPhi_H);
	      dHistThrown_OmegaPiMinusMassVsCosTheta[locThrownTopology]->Fill(thrownOmegaPiMinusP4.M(), thrownCosTheta);
	      dHistThrown_OmegaPiMinusMassVsPhi[locThrownTopology]->Fill(thrownOmegaPiMinusP4.M(), thrownPhi);
	      dHistThrown_OmegaPiMinusMassVsCosTheta_H[locThrownTopology]->Fill(thrownOmegaPiMinusP4.M(), thrownCosTheta_H);
	      dHistThrown_OmegaPiMinusMassVsPhi_H[locThrownTopology]->Fill(thrownOmegaPiMinusP4.M(), thrownPhi_H);
	    }
	  }
	}
	if(dTreeInterface->Get_Branch("NumCombos") == NULL)
	  return kTRUE;

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/


	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();	
	
	/************************************************* LOOP OVER COMBOS *************************************************/
	
	
	
	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
	    	//Set branch array indices for combo and all combo particles
	  	dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID 		= dComboBeamWrapper->Get_BeamID();
		Int_t locPhoton1NeutralID 	= dPhoton1Wrapper->Get_NeutralID();
		Int_t locPiPlusTrackID 		= dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID 	= dPiMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID 		= dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton2NeutralID 	= dPhoton2Wrapper->Get_NeutralID();
		Int_t locPhoton3NeutralID 	= dPhoton3Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/
	
		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 	= dComboBeamWrapper->Get_P4();
		TLorentzVector locPhoton1P4 	= dPhoton1Wrapper->Get_P4();
		TLorentzVector locPiPlusP4 	= dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4 	= dPiMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 	= dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locPhoton2P4 	= dPhoton2Wrapper->Get_P4();
		TLorentzVector locPhoton3P4 	= dPhoton3Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured 	= dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPhoton1P4_Measured 	= dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured 	= dPiPlusWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured 	= dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured 	= dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton2P4_Measured 	= dPhoton2Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton3P4_Measured 	= dPhoton3Wrapper->Get_P4_Measured();

		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		TLorentzVector locBeamX4_Measured 	= dComboBeamWrapper->Get_X4_Measured();
		Double_t locBunchPeriod 		= dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t locDeltaT_RF 			= dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		Int_t locRelBeamBucket 			= dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t locNumOutOfTimeBunchesInTree 	= 4; // YOU need to specify this number. Is number of out-of-time beam bunches in tree 
								// (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 
		Bool_t locSkipNearestOutOfTimeBunch 	= true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t locNumOutOfTimeBunchesToUse 	= locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		Double_t locAccidentalScalingFactor 	= dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		Double_t locHistAccidWeightFactor 	= locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time

		if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		  dComboWrapper->Set_IsComboCut(true); 
		  continue; 
		} 

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured 		= locBeamP4_Measured + dTargetP4;
		TLorentzVector locMissingP4NoPhoton_Measured    = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured 	       -= locPhoton1P4_Measured + locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured;
		locMissingP4NoPhoton_Measured  -= locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured;
		
		TLorentzVector locPi0P4 	= locPhoton2P4 + locPhoton3P4;
		TLorentzVector locPi02P4	= locPhoton1P4 + locPhoton2P4;
		TLorentzVector locPi03P4	= locPhoton1P4 + locPhoton3P4;
		TLorentzVector loc3PiP4 	= locPi0P4 + locPiPlusP4 + locPiMinusP4;
		TLorentzVector loc3Pi2P4 	= locPhoton1P4 + locPhoton2P4 + locPiPlusP4 + locPiMinusP4;
		TLorentzVector loc3Pi3P4 	= locPhoton1P4 + locPhoton3P4 + locPiPlusP4 + locPiMinusP4;
		TLorentzVector locPi0GammaP4 	= locPi0P4 + locPhoton1P4;
		TLorentzVector locPi0PiMinusP4  = locPi0P4 + locPiMinusP4;

		TLorentzVector locOmegaPiMinusP4= locPi0GammaP4 + locPiMinusP4;
		TLorentzVector locDeltaPlusPlus = locProtonP4 + locPiPlusP4;

		/*********************************************** DEFINE DECAY ANGLES *********************************************/
		
		// boosting into production plane (beamGamma-proton rest frame)
		TLorentzVector locBeamProtonP4 			= locBeamP4 + dTargetP4; // create 4vec of gamma-p system
		TVector3 locBeamProtonBoost    			= locBeamProtonP4.BoostVector(); // make boost vector (beta vector) from gamma-p rest frame to lab. Negative of this gives lab->gamma-p boost
		TLorentzVector locBeamP4_beamProtonRest 	= locBeamP4; 		locBeamP4_beamProtonRest.Boost(-1.0*locBeamProtonBoost); // boost beam gamma into gamma-p rf
		TLorentzVector locOmegaPiMinusP4_beamProtonRest = locOmegaPiMinusP4; 	locOmegaPiMinusP4_beamProtonRest.Boost(-1.0*locBeamProtonBoost); // boost b1 into ""
		TLorentzVector locOmegaP4_beamProtonRest 	= locPi0GammaP4; 	locOmegaP4_beamProtonRest.Boost(-1.0*locBeamProtonBoost); // boost Omega into ""
		TLorentzVector locGammaP4_beamProtonRest	= locPhoton1P4;		locGammaP4_beamProtonRest.Boost(-1.0*locBeamProtonBoost); // boost radiated photon into ""
		
		// define production plane unit vectors
		TVector3 locUnitK = locBeamP4_beamProtonRest.Vect().Unit();
		TVector3 locUnitZ = locOmegaPiMinusP4_beamProtonRest.Vect().Unit();
		TVector3 locUnitY = locUnitK.Cross(locUnitZ).Unit();
		TVector3 locUnitX = locUnitY.Cross(locUnitZ);
		
		// Get b1 boost vector (from production plane) and boost momentums into it
		TVector3 locOmegaPiMinusBoost		= locOmegaPiMinusP4_beamProtonRest.BoostVector(); 
		TLorentzVector locOmegaP4_b1Rest  	= locOmegaP4_beamProtonRest; locOmegaP4_b1Rest.Boost(-1.0*locOmegaPiMinusBoost);
		TLorentzVector locGammaP4_b1Rest	= locGammaP4_beamProtonRest; locGammaP4_b1Rest.Boost(-1.0*locOmegaPiMinusBoost);
		
		// Define omega decay angles
		TVector3 locOmegaP3_b1Rest 	= locOmegaP4_b1Rest.Vect(); // get just the momentum vector for the dot products
		double locCosTheta 		= locOmegaP3_b1Rest.Dot(locUnitZ) / locOmegaP3_b1Rest.Mag();
		double locPhi 			= 180. * TMath::ATan2(locOmegaP3_b1Rest.Dot(locUnitY) , locOmegaP3_b1Rest.Dot(locUnitX)) / TMath::Pi();

		// Define decay plane unit vectors (helicity frame)
		TVector3 locUnitZ_H = locOmegaP4_b1Rest.Vect().Unit();	 // z_h is direction of omega in b_1 rest frame
		TVector3 locUnitY_H = locUnitZ.Cross(locUnitZ_H).Unit(); // y_h is cross product of b_1 frame velocity and z_h
		TVector3 locUnitX_H = locUnitY_H.Cross(locUnitZ_H);

		// Get Omega boost vector (to b1 rf) and boost radiated photon into it
		TVector3 locOmegaBoost 			= locOmegaP4_b1Rest.BoostVector();
		TLorentzVector locGammaP4_omegaRest 	= locGammaP4_b1Rest; locGammaP4_omegaRest.Boost(-1.0*locOmegaBoost);

		// Measure helicity angles wrt radiated photon direction (in omega rest frame)
		TVector3 locGammaP3_omegaRest 	= locGammaP4_omegaRest.Vect();
		double locCosTheta_H		= locGammaP3_omegaRest.Dot(locUnitZ_H) / locGammaP3_omegaRest.Mag();
		double locPhi_H			= 180. * TMath::ATan2(locGammaP3_omegaRest.Dot(locUnitY_H) , locGammaP3_omegaRest.Dot(locUnitX_H)) / TMath::Pi();

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/
		if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;
	        
		if(dIsMC) {dHistTop_PercentStart->Fill(locThrownTopology.Data(), locHistAccidWeightFactor);}

		// cut variables defined here
		double kinFitChiSq 		= dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit();
		double t 			= ( (locProtonP4+locPiPlusP4) - dTargetP4 ).M2();
		double numUnusedShowers 	= dComboWrapper->Get_NumUnusedShowers();
		double energyUnusedShowers 	= dComboWrapper->Get_Energy_UnusedShowers();		
		double showerQuality 		= dPhoton1Wrapper->Get_Shower_Quality();
		
		// setup cut values
		inTRange 		= (fabs(t) < 0.5)         ? true : false;	// limit momentum of t-channel exchange
		inKinFitChiSqRange 	= (kinFitChiSq < 3.0)     ? true : false;	// limit event ChiSquare from KinFit results
		inNumUnusedShowersRange = (numUnusedShowers <= 0) ? true : false;	// Essentially cut EnergyUnusedShowers < 100 MeV
		inReject3PiRange 	= (loc3PiP4.M() > 0.9)    ? true : false;	// veto omega->3pi decay
		inRejectWrongPi0Range 	= ( (fabs(locPi02P4.M() - 0.135) > 0.035) ||	
					    (fabs(locPi03P4.M() - 0.135) > 0.035) ) ? true : false; 	// pi0 combo cut: reject when photon 1+2 or 1+3 are the result of a pi0
		inDeltaPlusPlusRange 	= (locDeltaPlusPlus.M() < 1.3)              ? true : false; 	// select DeltaPlusPlus(1232)
		inPi0Range 		= (fabs(locPi0P4.M() - 0.135) < 0.015)      ? true : false; 	// select Pi0 from "right" photons		
		inPi0GammaRange 	= (fabs(locPi0GammaP4.M() - 0.78) < 0.019)  ? true : false;	// select omega
		inSidebandRange 	= (fabs(locPi0GammaP4.M() - (0.78 - 3*0.019)) < (0.019 / 2) ||  // Sideband Subtraction
					   fabs(locPi0GammaP4.M() - (0.78 + 3*0.019)) < (0.019 / 2)) ? true : false;		

		double sbWeight = (inSidebandRange) ? -1.0 : 1.0; // Multiply event weight  by -1.0 if in the sideband range
				
		/*  
		  | *****************************************************************************************************
		  | * EACH BLOCK OF HISTS BELOW TYPICALLY HAVE ALL CUTS APPLIED EXCEPT THE ONE OF INTEREST BEING HIST'D *
		  | *****************************************************************************************************
		 */
		//_______________ Missing Mass Squared_________________
		// all other blocks follow this setup, so refer to comments in this one
		if(inTRange && inKinFitChiSqRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  // since omega signal is in Pi0Gamma spectrum, and we sideband subtract on that spectrum, having this plot will allow us to extract what the variable of interest
		  //	looks like in signal region (inPi0GammaRange), in the sidebands, or the whole spectrum, all before sideband subtraction has been applied. This is done using ProjectionY function
		  dHist_Pi0GammaMassVsMM2->Fill(locPi0GammaP4.M(), locMissingP4_Measured.M2(), locHistAccidWeightFactor);
		  dHist_Pi0GammaMassVsMM2NoPhoton->Fill(locPi0GammaP4.M(), locMissingP4NoPhoton_Measured.M2(), locHistAccidWeightFactor);
		  if(dIsMC) { // bggen plots of above
		    if(dHistTop_Pi0GammaMassVsMM2.find(locThrownTopology) != dHistTop_Pi0GammaMassVsMM2.end()) { // fill each plot according to its respective thrown topology
		      dHistTop_Pi0GammaMassVsMM2[locThrownTopology]->Fill(locPi0GammaP4.M(), locMissingP4_Measured.M2(), locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsMM2NoPhoton[locThrownTopology]->Fill(locPi0GammaP4.M(), locMissingP4NoPhoton_Measured.M2(), locHistAccidWeightFactor);
		    }
		    else {
		      dHistTop_Pi0GammaMassVsMM2[bkgTT]->Fill(locPi0GammaP4.M(), locMissingP4_Measured.M2(), locHistAccidWeightFactor); // all other topologies not recorded above are tossed into "background"
		      dHistTop_Pi0GammaMassVsMM2NoPhoton[bkgTT]->Fill(locPi0GammaP4.M(), locMissingP4NoPhoton_Measured.M2(), locHistAccidWeightFactor);
		    }
		  }
		  // if events are in peak or sideaband than fill hists. The var "sbWeight" handles the appropriate weighting for whichever range its in
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_MM2->Fill(locMissingP4_Measured.M2(), sbWeight*locHistAccidWeightFactor);
		    dHist_MM2NoPhoton->Fill(locMissingP4NoPhoton_Measured.M2(), sbWeight*locHistAccidWeightFactor);
		    if(dIsMC) {
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_MM2[locThrownTopology]->Fill(locMissingP4_Measured.M2(), sbWeight*locHistAccidWeightFactor);
			dHistTop_MM2NoPhoton[locThrownTopology]->Fill(locMissingP4NoPhoton_Measured.M2(), sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_MM2[bkgTT]->Fill(locMissingP4_Measured.M2(), sbWeight*locHistAccidWeightFactor);
			dHistTop_MM2NoPhoton[bkgTT]->Fill(locMissingP4NoPhoton_Measured.M2(), sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  } 
		}
		//________________Angular Distributions________________		
		if(inTRange && inKinFitChiSqRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  dHist_Pi0GammaMassVsCosTheta->Fill(locPi0GammaP4.M(), locCosTheta, locHistAccidWeightFactor); 
		  dHist_Pi0GammaMassVsCosTheta_H->Fill(locPi0GammaP4.M(), locCosTheta_H, locHistAccidWeightFactor); 
		  dHist_Pi0GammaMassVsPhi->Fill(locPi0GammaP4.M(), locPhi, locHistAccidWeightFactor);
		  dHist_Pi0GammaMassVsPhi_H->Fill(locPi0GammaP4.M(), locPhi_H, locHistAccidWeightFactor);
		  if(dIsMC) { 
		    if(dHistTop_Pi0GammaMassVsCosTheta.find(locThrownTopology) != dHistTop_Pi0GammaMassVsCosTheta.end()) {
		      dHistTop_Pi0GammaMassVsCosTheta[locThrownTopology]->Fill(locPi0GammaP4.M(), locCosTheta, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsCosTheta_H[locThrownTopology]->Fill(locPi0GammaP4.M(), locCosTheta_H, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsPhi[locThrownTopology]->Fill(locPi0GammaP4.M(), locPhi, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsPhi_H[locThrownTopology]->Fill(locPi0GammaP4.M(), locPhi_H, locHistAccidWeightFactor);
		    }
		    else {
		      dHistTop_Pi0GammaMassVsCosTheta[bkgTT]->Fill(locPi0GammaP4.M(), locCosTheta, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsCosTheta_H[bkgTT]->Fill(locPi0GammaP4.M(), locCosTheta_H, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsPhi[bkgTT]->Fill(locPi0GammaP4.M(), locPhi, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsPhi_H[bkgTT]->Fill(locPi0GammaP4.M(), locPhi_H, locHistAccidWeightFactor);
		    }
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_CosTheta->Fill(locCosTheta, sbWeight*locHistAccidWeightFactor);
		    dHist_CosTheta_H->Fill(locCosTheta_H, sbWeight*locHistAccidWeightFactor);
		    dHist_Phi->Fill(locPhi, sbWeight*locHistAccidWeightFactor);
		    dHist_Phi_H->Fill(locPhi_H, sbWeight*locHistAccidWeightFactor);
		    
		    dHist_CosThetaVsPhi->Fill(locCosTheta, locPhi, sbWeight*locHistAccidWeightFactor);
		    dHist_CosTheta_HVsPhi_H->Fill(locCosTheta_H, locPhi_H, sbWeight*locHistAccidWeightFactor);

		    dHist_OmegaPiMinusMassVsCosTheta->Fill(locOmegaPiMinusP4.M(), locCosTheta, sbWeight*locHistAccidWeightFactor);
		    dHist_OmegaPiMinusMassVsCosTheta_H->Fill(locOmegaPiMinusP4.M(), locCosTheta_H, sbWeight*locHistAccidWeightFactor);
		    dHist_OmegaPiMinusMassVsPhi->Fill(locOmegaPiMinusP4.M(), locPhi, sbWeight*locHistAccidWeightFactor);
		    dHist_OmegaPiMinusMassVsPhi_H->Fill(locOmegaPiMinusP4.M(), locPhi_H, sbWeight*locHistAccidWeightFactor);

		    if(dIsMC) {
		      dHistTop_PercentAngular->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_CosTheta[locThrownTopology]->Fill(locCosTheta, sbWeight*locHistAccidWeightFactor);
			dHistTop_CosTheta_H[locThrownTopology]->Fill(locCosTheta_H, sbWeight*locHistAccidWeightFactor);
			dHistTop_Phi[locThrownTopology]->Fill(locPhi, sbWeight*locHistAccidWeightFactor);
			dHistTop_Phi_H[locThrownTopology]->Fill(locPhi_H, sbWeight*locHistAccidWeightFactor);
			
			dHistTop_CosThetaVsPhi[locThrownTopology]->Fill(locCosTheta, locPhi, sbWeight*locHistAccidWeightFactor);
			dHistTop_CosTheta_HVsPhi_H[locThrownTopology]->Fill(locCosTheta_H, locPhi_H, sbWeight*locHistAccidWeightFactor);

			dHistTop_OmegaPiMinusMassVsCosTheta[locThrownTopology]->Fill(locOmegaPiMinusP4.M(), locCosTheta, sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsCosTheta_H[locThrownTopology]->Fill(locOmegaPiMinusP4.M(), locCosTheta_H, sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsPhi[locThrownTopology]->Fill(locOmegaPiMinusP4.M(), locPhi, sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsPhi_H[locThrownTopology]->Fill(locOmegaPiMinusP4.M(), locPhi_H, sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_CosTheta[bkgTT]->Fill(locCosTheta, sbWeight*locHistAccidWeightFactor);
			dHistTop_CosTheta_H[bkgTT]->Fill(locCosTheta_H, sbWeight*locHistAccidWeightFactor);
			dHistTop_Phi[bkgTT]->Fill(locPhi, sbWeight*locHistAccidWeightFactor);
			dHistTop_Phi_H[bkgTT]->Fill(locPhi_H, sbWeight*locHistAccidWeightFactor);
			
			dHistTop_CosThetaVsPhi[bkgTT]->Fill(locCosTheta, locPhi, sbWeight*locHistAccidWeightFactor);
			dHistTop_CosTheta_HVsPhi_H[bkgTT]->Fill(locCosTheta_H, locPhi_H, sbWeight*locHistAccidWeightFactor);
			
			dHistTop_OmegaPiMinusMassVsCosTheta[bkgTT]->Fill(locOmegaPiMinusP4.M(), locCosTheta, sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsCosTheta_H[bkgTT]->Fill(locOmegaPiMinusP4.M(), locCosTheta_H, sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsPhi[bkgTT]->Fill(locOmegaPiMinusP4.M(), locPhi, sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsPhi_H[bkgTT]->Fill(locOmegaPiMinusP4.M(), locPhi_H, sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  }
		}
		//________ChiSquare__________
		if(inTRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  dHist_Pi0GammaMassVsKinFitChiSq->Fill(locPi0GammaP4.M(), kinFitChiSq, locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVsKinFitChiSq.find(locThrownTopology) != dHistTop_Pi0GammaMassVsKinFitChiSq.end()) {
		      dHistTop_Pi0GammaMassVsKinFitChiSq[locThrownTopology]->Fill(locPi0GammaP4.M(), kinFitChiSq, locHistAccidWeightFactor);		      
		    }
		    else {
		      dHistTop_Pi0GammaMassVsKinFitChiSq[bkgTT]->Fill(locPi0GammaP4.M(), kinFitChiSq, locHistAccidWeightFactor);
		    }
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_KinFitChiSq->Fill(kinFitChiSq, sbWeight*locHistAccidWeightFactor);
		    if(dIsMC) {
		      dHistTop_PercentChi2->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_KinFitChiSq[locThrownTopology]->Fill(kinFitChiSq, sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_KinFitChiSq[bkgTT]->Fill(kinFitChiSq, sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  }
		}
		//__________T_________		
		if(inKinFitChiSqRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  dHist_Pi0GammaMassVsT->Fill(locPi0GammaP4.M(), -1*t, locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVsT.find(locThrownTopology) != dHistTop_Pi0GammaMassVsT.end()) {
		      dHistTop_Pi0GammaMassVsT[locThrownTopology]->Fill(locPi0GammaP4.M(), -1*t, locHistAccidWeightFactor);		      
		    }
		    else {
		      dHistTop_Pi0GammaMassVsT[bkgTT]->Fill(locPi0GammaP4.M(), -1*t, locHistAccidWeightFactor);		      
		    }
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_T->Fill(-1*t, sbWeight*locHistAccidWeightFactor);
		    if(dIsMC){
		      dHistTop_PercentT->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_T[locThrownTopology]->Fill(-1*t, sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_T[bkgTT]->Fill(-1*t, sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  }
		}
		//__________Photons__________		
		// Excluded or unused cuts to be tested ?
		//if(dPhoton1Wrapper->Get_Shower_Quality() < 0.75 ) continue;
		//if(dIsMC) dHistTop_Cut_ShowerQuality->Fill(locThrownTopology.Data(), locHistAccidWeightFactor);
		//if( dPhoton1Wrapper->Get_Energy_BCAL() <= 0.0 && dPhoton1Wrapper->Get_Energy_FCAL() < 0.5 ) continue;
		//if( dPhoton2Wrapper->Get_Energy_BCAL() <= 0.0 && dPhoton2Wrapper->Get_Energy_FCAL() < 0.25 ) continue;
		//if( dPhoton3Wrapper->Get_Energy_BCAL() <= 0.0 && dPhoton3Wrapper->Get_Energy_FCAL() < 0.25 ) continue;

		if(inKinFitChiSqRange && inTRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  dHist_Pi0GammaMassVsShowerQuality->Fill(locPi0GammaP4.M(), showerQuality, locHistAccidWeightFactor);
		  dHist_Pi0GammaMassVsEnergyUnusedShowers->Fill(locPi0GammaP4.M(), energyUnusedShowers, locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVsShowerQuality.find(locThrownTopology) != dHistTop_Pi0GammaMassVsShowerQuality.end()) {
		      dHistTop_Pi0GammaMassVsShowerQuality[locThrownTopology]->Fill(locPi0GammaP4.M(), showerQuality, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsEnergyUnusedShowers[locThrownTopology]->Fill(locPi0GammaP4.M(), energyUnusedShowers, locHistAccidWeightFactor);
		    }
		    else {
		      dHistTop_Pi0GammaMassVsShowerQuality[bkgTT]->Fill(locPi0GammaP4.M(), showerQuality, locHistAccidWeightFactor);
		      dHistTop_Pi0GammaMassVsEnergyUnusedShowers[bkgTT]->Fill(locPi0GammaP4.M(), energyUnusedShowers, locHistAccidWeightFactor);
		    }
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_ShowerQuality->Fill(showerQuality, sbWeight*locHistAccidWeightFactor);
		    dHist_EnergyUnusedShowers->Fill(energyUnusedShowers, sbWeight*locHistAccidWeightFactor);
		    if(dIsMC){
		      dHistTop_PercentPhoton->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_ShowerQuality[locThrownTopology]->Fill(showerQuality, sbWeight*locHistAccidWeightFactor);
			dHistTop_EnergyUnusedShowers[locThrownTopology]->Fill(energyUnusedShowers, sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_ShowerQuality[bkgTT]->Fill(showerQuality, sbWeight*locHistAccidWeightFactor);
			dHistTop_EnergyUnusedShowers[bkgTT]->Fill(energyUnusedShowers, sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  }
		}
		
		//__________Pi0 Mass__________
		if(inKinFitChiSqRange && inTRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange) {
		  dHist_Pi0GammaMassVsPi0Mass->Fill(locPi0GammaP4.M(), locPi0P4.M(), locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVsPi0Mass.find(locThrownTopology) != dHistTop_Pi0GammaMassVsPi0Mass.end()) {
		      dHistTop_Pi0GammaMassVsPi0Mass[locThrownTopology]->Fill(locPi0GammaP4.M(), locPi0P4.M(), locHistAccidWeightFactor);
		    }
		    else {
		      dHistTop_Pi0GammaMassVsPi0Mass[bkgTT]->Fill(locPi0GammaP4.M(), locPi0P4.M(), locHistAccidWeightFactor);
		    }
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_Pi0Mass->Fill(locPi0P4.M(), sbWeight*locHistAccidWeightFactor);
		    if(dIsMC){
		      dHistTop_PercentPi0->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);		  		      
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_Pi0Mass[locThrownTopology]->Fill(locPi0P4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_Pi0Mass[bkgTT]->Fill(locPi0P4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  }
		}
		//__________Pi0 Combos__________ (includes "wrong" photon combinations being a pi0)
		if(inKinFitChiSqRange && inTRange && inNumUnusedShowersRange  && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_Pi02Mass->Fill(locPi02P4.M(), sbWeight*locHistAccidWeightFactor);
		    dHist_Pi03Mass->Fill(locPi03P4.M(), sbWeight*locHistAccidWeightFactor);
		    if(dIsMC) {
		      dHistTop_PercentPi0Combo->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_Pi02Mass[locThrownTopology]->Fill(locPi02P4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_Pi03Mass[locThrownTopology]->Fill(locPi03P4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_Pi02Mass[bkgTT]->Fill(locPi02P4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_Pi03Mass[bkgTT]->Fill(locPi03P4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  }
		}
		//__________3Pi Mass Combos__________ (includes the omega->3pi decay)
		if(inKinFitChiSqRange && inTRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inDeltaPlusPlusRange && inPi0Range) {
		  dHist_Pi0GammaMassVs3PiMass->Fill(locPi0GammaP4.M(), loc3PiP4.M(), locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVs3PiMass.find(locThrownTopology) != dHistTop_Pi0GammaMassVs3PiMass.end())
		      dHistTop_Pi0GammaMassVs3PiMass[locThrownTopology]->Fill(locPi0GammaP4.M(), loc3PiP4.M(), locHistAccidWeightFactor);
		    else
		      dHistTop_Pi0GammaMassVs3PiMass[bkgTT]->Fill(locPi0GammaP4.M(), loc3PiP4.M(), locHistAccidWeightFactor);
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_3PiMass->Fill(loc3PiP4.M(), sbWeight*locHistAccidWeightFactor);
		    dHist_3Pi2Mass->Fill(loc3Pi2P4.M(), sbWeight*locHistAccidWeightFactor);
		    dHist_3Pi3Mass->Fill(loc3Pi3P4.M(), sbWeight*locHistAccidWeightFactor);
		    if(dIsMC){
		      dHistTop_Percent3PiCombo->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_3PiMass[locThrownTopology]->Fill(loc3PiP4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_3Pi2Mass[locThrownTopology]->Fill(loc3Pi2P4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_3Pi3Mass[locThrownTopology]->Fill(loc3Pi3P4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_3PiMass[bkgTT]->Fill(loc3PiP4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_3Pi2Mass[bkgTT]->Fill(loc3Pi2P4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_3Pi3Mass[bkgTT]->Fill(loc3Pi3P4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		    }
		  }
		}
		
		//__________ProtonPiPlus__________ (look for Delta++ events)
		if(inKinFitChiSqRange && inTRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inPi0Range) {
		  dHist_Pi0GammaMassVsProtonPiPlusMass->Fill(locPi0GammaP4.M(), locDeltaPlusPlus.M(), locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVsProtonPiPlusMass.find(locThrownTopology) != dHistTop_Pi0GammaMassVsProtonPiPlusMass.end())
		      dHistTop_Pi0GammaMassVsProtonPiPlusMass[locThrownTopology]->Fill(locPi0GammaP4.M(), locDeltaPlusPlus.M(), locHistAccidWeightFactor);
		    else 
		      dHistTop_Pi0GammaMassVsProtonPiPlusMass[bkgTT]->Fill(locPi0GammaP4.M(), locDeltaPlusPlus.M(), locHistAccidWeightFactor);
		  }	  
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_ProtonPiPlusMass->Fill(locDeltaPlusPlus.M(), sbWeight*locHistAccidWeightFactor);		 		  
		    dHist_ProtonPiPlusMassVsPi0PiMinusMass->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		    
		    if(sbWeight > 0.) {dHist_ProtonPiPlusMassVsPi0PiMinusMass_PK->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);}
		    if(sbWeight < 0.) {dHist_ProtonPiPlusMassVsPi0PiMinusMass_SB->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);}
		    if(dIsMC) {
		      dHistTop_PercentProtonPiPlus->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);		      
		      if(dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_SB.find(locThrownTopology) != dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_SB.end()) {
			if(sbWeight > 0.) {dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_PK[locThrownTopology]->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);}
			if(sbWeight < 0.) {dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_SB[locThrownTopology]->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);}
		      }
		      if(locThrownTopology.EqualTo(signalTT)) {
			  dHistTop_ProtonPiPlusMass[locThrownTopology]->Fill(locDeltaPlusPlus.M(), sbWeight*locHistAccidWeightFactor);
			  dHistTop_ProtonPiPlusMassVsPi0PiMinusMass[locThrownTopology]->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_ProtonPiPlusMass[bkgTT]->Fill(locDeltaPlusPlus.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_ProtonPiPlusMassVsPi0PiMinusMass[bkgTT]->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
			if(sbWeight > 0.) {dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_PK[bkgTT]->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);}
			if(sbWeight < 0.) {dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_SB[bkgTT]->Fill(locDeltaPlusPlus.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);}
		      }
		    }
		  }
		}
		//__________Pi0PiMinus___________
		if(inKinFitChiSqRange && inTRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  dHist_Pi0GammaMassVsPi0PiMinusMass->Fill(locPi0GammaP4.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVsPi0PiMinusMass.find(locThrownTopology) != dHistTop_Pi0GammaMassVsPi0PiMinusMass.end())
		      dHistTop_Pi0GammaMassVsPi0PiMinusMass[locThrownTopology]->Fill(locPi0GammaP4.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);
		    else
		      dHistTop_Pi0GammaMassVsPi0PiMinusMass[bkgTT]->Fill(locPi0GammaP4.M(), locPi0PiMinusP4.M(), locHistAccidWeightFactor);
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_Pi0PiMinusMass->Fill(locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		    dHist_OmegaPiMinusMassVsPi0PiMinusMass->Fill(locOmegaPiMinusP4.M(), locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		    if(dIsMC) {
		      dHistTop_PercentPi0PiMinus->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) {
			dHistTop_Pi0PiMinusMass[locThrownTopology]->Fill(locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsPi0PiMinusMass[locThrownTopology]->Fill(locOmegaPiMinusP4.M(), locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		      }
		      else {
			dHistTop_Pi0PiMinusMass[bkgTT]->Fill(locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
			dHistTop_OmegaPiMinusMassVsPi0PiMinusMass[bkgTT]->Fill(locOmegaPiMinusP4.M(), locPi0PiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		      } 
		    }
		  }
		}
		//__________Pi0Gamma__________
		if(inKinFitChiSqRange && inTRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {		  
		  dHist_Pi0GammaMass->Fill(locPi0GammaP4.M(), locHistAccidWeightFactor);
		  if(dIsMC){
		    dHistTop_PercentPi0Gamma->Fill(locThrownTopology.Data(), locHistAccidWeightFactor);
		    if(dHistTop_Pi0GammaMass.find(locThrownTopology) != dHistTop_Pi0GammaMass.end())
		      dHistTop_Pi0GammaMass[locThrownTopology]->Fill(locPi0GammaP4.M(), locHistAccidWeightFactor);
		    else
		      dHistTop_Pi0GammaMass[bkgTT]->Fill(locPi0GammaP4.M(),locHistAccidWeightFactor);
		  }
		}
		//_________OmegaPiMinus________
		if(inKinFitChiSqRange && inTRange && inNumUnusedShowersRange && inRejectWrongPi0Range && inReject3PiRange && inDeltaPlusPlusRange && inPi0Range) {
		  dHist_Pi0GammaMassVsOmegaPiMinusMass->Fill(locPi0GammaP4.M(), locOmegaPiMinusP4.M(), locHistAccidWeightFactor);
		  if(dIsMC) {
		    if(dHistTop_Pi0GammaMassVsOmegaPiMinusMass.find(locThrownTopology) != dHistTop_Pi0GammaMassVsOmegaPiMinusMass.end())
		      dHistTop_Pi0GammaMassVsOmegaPiMinusMass[locThrownTopology]->Fill(locPi0GammaP4.M(), locOmegaPiMinusP4.M(), locHistAccidWeightFactor);
		    else
		      dHistTop_Pi0GammaMassVsOmegaPiMinusMass[bkgTT]->Fill(locPi0GammaP4.M(), locOmegaPiMinusP4.M(), locHistAccidWeightFactor);
		  }
		  if(inPi0GammaRange || inSidebandRange) {
		    dHist_OmegaPiMinusMass->Fill(locOmegaPiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		    if(dIsMC) {
		      dHistTop_PercentOmega->Fill(locThrownTopology.Data(), sbWeight*locHistAccidWeightFactor);
		      if(locThrownTopology.EqualTo(signalTT)) 
			dHistTop_OmegaPiMinusMass[locThrownTopology]->Fill(locOmegaPiMinusP4.M(), sbWeight*locHistAccidWeightFactor);			
		      else
			dHistTop_OmegaPiMinusMass[bkgTT]->Fill(locOmegaPiMinusP4.M(), sbWeight*locHistAccidWeightFactor);
		    }
		  }
		}
		
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
