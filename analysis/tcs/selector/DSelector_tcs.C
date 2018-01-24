#include "DSelector_tcs.h"

void DSelector_tcs::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "tcs.root"; //"" for none
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

	//PID HISTOGRAMS
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "NoCut"));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, true, 2.0, Proton, SYS_BCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, true, 1.0, Proton, SYS_TOF));

	//MASSES and KINEMATICS HISTOGRAMS
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, Jpsi, 400, 0.0, 4.0, "JpsiNoCut"));
	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.25, 0.25, "NoCut"));

	//KINFIT RESULTS HISTOGRAMS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS and KINFIT CL
	dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.25, 0.25));
	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 5e-7));
	
	//MASSES and KINEMATICS (after KINFIT CL cut)
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "KinCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, Jpsi, 400, 0.0, 4.0, "JpsiKinCut"));
	dAnalysisActions.push_back(new DHistogramAction_MissingP4(dComboWrapper, false, "KinCut"));

	//CUT E/P for e+/-
	//dAnalysisActions.push_back(new DCutAction_TrackShowerEOverP(dComboWrapper, true, SYS_BCAL, Positron, 0.5, "PositronEOverP"));
	//dAnalysisActions.push_back(new DCutAction_TrackShowerEOverP(dComboWrapper, true, SYS_BCAL, Electron, 0.5, "ElectronEOverP"));
	//dAnalysisActions.push_back(new DCutAction_TrackShowerEOverP(dComboWrapper, true, SYS_FCAL, Positron, 0.5, "PositronEOverP"));
	//dAnalysisActions.push_back(new DCutAction_TrackShowerEOverP(dComboWrapper, true, SYS_FCAL, Electron, 0.5, "ElectronEOverP"));
	
	//MASSES and KINEMATICS (after E/P cut)
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "EOverPCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, Jpsi, 400, 0.0, 4.0, "JpsiEOverPCut"));
	dAnalysisActions.push_back(new DHistogramAction_MissingP4(dComboWrapper, false, "EOverPCut"));
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false, "EOverPCut"));

	//CUT BEAM ENERGY 
	dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.0, 13.0));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	/*
	//CUT BCAL Preshower e+/-
	dAnalysisActions.push_back(new DCutAction_TrackBCALPreshowerFraction(dComboWrapper, true, Positron, 0.15, "PositronPreshower"));
	dAnalysisActions.push_back(new DCutAction_TrackBCALPreshowerFraction(dComboWrapper, true, Electron, 0.15, "ElectronPreshower"));

	//MASSES and KINEMATICS (after preshower cut)
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "PreshowerCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, Jpsi, 400, 0.0, 4.0, "JpsiPreshowerCut"));
	dAnalysisActions.push_back(new DHistogramAction_MissingP4(dComboWrapper, false, "PreshowerCut"));
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false, "PreshowerCut"));

	//CUT BEAM ENERGY 
	dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.0, 12.0));

	//MASSES and KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "BeamCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, Jpsi, 400, 0.0, 4.0, "JpsiBeamCut"));
	dAnalysisActions.push_back(new DHistogramAction_MissingP4(dComboWrapper, false, "BeamCut"));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));
	*/

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	// Some directory structure for histograms
	TDirectory *locFileDir = gDirectory;
        TDirectory *locDir = locFileDir->mkdir("CustomHistograms");

	TDirectory *locParticleDir[2];
	locParticleDir[0] = locDir->mkdir("Electron");
	locParticleDir[1] = locDir->mkdir("Positron");
	locDir->cd();

	double locMinMee = 0.0, locMaxMee = 4.0;
	double locMinEOverP = 0.0, locMaxEOverP = 1.5;
	double locMinP = 0.0, locMaxP = 12.0;
	double locMinTheta = 0.0, locMaxTheta = 120.;
	double locMinPreshowerEoverShowerE = 0.0, locMaxPreshowerEoverShowerE = 0.5;
	double locMinSigTrans = 0.0, locMaxSigTrans = 0.02;
	double locMinSigTheta = 0.0, locMaxSigTheta = 0.05;
	double locMinSigLong = 0.0, locMaxSigLong = 10.;
	double locMindEdx = 0.001, locMaxdEdx = 5.;
        double locMinDeltaT = -2.0, locMaxDeltaT = 2.;

	//TString locRFBinLabel[2] = {"", "Acci_"};
        //TString locCutCounterLabel[5] = {"Fiducial", "UE", "DeltaTheta", "DeltaPhiSignal", "DeltaPhiSideband"};
	TString locParticleLabel[2] = {"Electron", "Positron"};

	//MANUAL HISTOGRAMS:
	dHist_Mee = new TH1I("Mee", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_CDCdEdx = new TH1I("Mee_CDCdEdx", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_ShowerCut = new TH1I("Mee_ShowerCut", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_KinFit1 = new TH1I("Mee_KinFit1", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_KinFit2 = new TH1I("Mee_KinFit2", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_KinFit3 = new TH1I("Mee_KinFit3", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_KinFit4 = new TH1I("Mee_KinFit4", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);

	dHist_EleEOverP_PosEOverP_Jpsi = new TH2I("EleEOverP_PosEOverP_Jpsi", ";Positron E/p; Electron E/p", 150, locMinEOverP, locMaxEOverP, 150, locMinEOverP, locMaxEOverP);
	dHist_EleEOverP_PosEOverP_Bkgd = new TH2I("EleEOverP_PosEOverP_Bkgd", ";Positron E/p; Electron E/p", 150, locMinEOverP, locMaxEOverP, 150, locMinEOverP, locMaxEOverP);

	// Separate histograms for e+ and e-
	for(int loc_i = 0; loc_i < 2; loc_i++) {
		locParticleDir[loc_i]->cd(); // store in e+/e- directories
		dHist_EOverP_Mee[loc_i] = new TH2I(Form("%sEOverP_Mee",locParticleLabel[loc_i].Data()), Form(";Mass e^{+}e^{-} (GeV/c^{2}); %s E/p", locParticleLabel[loc_i].Data()), 400, locMinMee, locMaxMee, 150, locMinEOverP, locMaxEOverP);
		dHist_P_Mee[loc_i] = new TH2I(Form("%sP_Mee",locParticleLabel[loc_i].Data()), Form(";Mass e^{+}e^{-} (GeV/c^{2}); %s p (GeV)", locParticleLabel[loc_i].Data()), 400, locMinMee, locMaxMee, 120, locMinP, locMaxP);
		dHist_Theta_Mee[loc_i] = new TH2I(Form("%sTheta_Mee",locParticleLabel[loc_i].Data()), Form(";Mass e^{+}e^{-} (GeV/c^{2}); %s #theta (degrees)", locParticleLabel[loc_i].Data()), 400, locMinMee, locMaxMee, 120, locMinTheta, locMaxTheta);
		dHist_P_Mee_CDCdEdx[loc_i] = new TH2I(Form("%sP_Mee_CDCdEdx",locParticleLabel[loc_i].Data()), Form(";Mass e^{+}e^{-} (GeV/c^{2}); %s p (GeV)", locParticleLabel[loc_i].Data()), 400, locMinMee, locMaxMee, 120, locMinP, locMaxP);
                dHist_Theta_Mee_CDCdEdx[loc_i] = new TH2I(Form("%sTheta_Mee_CDCdEdx",locParticleLabel[loc_i].Data()), Form(";Mass e^{+}e^{-} (GeV/c^{2}); %s #theta (degrees)", locParticleLabel[loc_i].Data()), 400, locMinMee, locMaxMee, 120, locMinTheta, locMaxTheta);


		dHist_EOverP_Preshower_Jpsi[loc_i] = new TH2I(Form("%sBCALEOverP_Preshower_Jpsi",locParticleLabel[loc_i].Data()), Form(";%s E(Preshower)/E(Shower); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 200, locMinPreshowerEoverShowerE, locMaxPreshowerEoverShowerE, 150, locMinEOverP, locMaxEOverP);
		dHist_EOverP_Preshower_Bkgd[loc_i] = new TH2I(Form("%sBCALEOverP_Preshower_Bkgd",locParticleLabel[loc_i].Data()), Form(";%s E(Preshower)/E(Shower); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 200, locMinPreshowerEoverShowerE, locMaxPreshowerEoverShowerE, 150, locMinEOverP, locMaxEOverP);

		dHistEOverP_P_BCAL_Jpsi[loc_i] = new TH2I(Form("%sBCALEOverP_P_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 150, locMinEOverP, locMaxEOverP);
                dHistEOverP_P_FCAL_Jpsi[loc_i] = new TH2I(Form("%sFCALEOverP_P_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 10., 150, locMinEOverP, locMaxEOverP);
                dHistEOverP_Theta_BCAL_Jpsi[loc_i] = new TH2I(Form("%sBCALEOverP_Theta_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 10., 110., 150, locMinEOverP, locMaxEOverP);
                dHistEOverP_Theta_FCAL_Jpsi[loc_i] = new TH2I(Form("%sFCALEOverP_Theta_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 120, 0., 12., 150, locMinEOverP, locMaxEOverP);
		dHistEOverP_P_BCAL_Bkgd[loc_i] = new TH2I(Form("%sBCALEOverP_P_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 150, locMinEOverP, locMaxEOverP);
                dHistEOverP_P_FCAL_Bkgd[loc_i] = new TH2I(Form("%sFCALEOverP_P_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 10., 150, locMinEOverP, locMaxEOverP);
                dHistEOverP_Theta_BCAL_Bkgd[loc_i] = new TH2I(Form("%sBCALEOverP_Theta_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 10., 110., 150, locMinEOverP, locMaxEOverP);
                dHistEOverP_Theta_FCAL_Bkgd[loc_i] = new TH2I(Form("%sFCALEOverP_Theta_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s E/p", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 120, 0., 12., 150, locMinEOverP, locMaxEOverP);

		dHist_SigTrans_Preshower_Jpsi[loc_i] = new TH2I(Form("%sSigTrans_Preshower_Jpsi",locParticleLabel[loc_i].Data()), Form(";%s E(Preshower)/E(Shower); %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 200, locMinPreshowerEoverShowerE, locMaxPreshowerEoverShowerE, 100, locMinSigTrans, locMaxSigTrans);
		dHist_SigTheta_Preshower_Jpsi[loc_i] = new TH2I(Form("%sSigTheta_Preshower_Jpsi",locParticleLabel[loc_i].Data()), Form(";%s E(Preshower)/E(Shower); %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 200, locMinPreshowerEoverShowerE, locMaxPreshowerEoverShowerE, 100, locMinSigTheta, locMaxSigTheta);
		dHist_SigTrans_SigTheta_Jpsi[loc_i] = new TH2I(Form("%sSigTrans_SigTheta_Jpsi",locParticleLabel[loc_i].Data()), Form(";%s #sigma_{theta}; %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, locMinSigTheta, locMaxSigTheta, 100, locMinSigTrans, locMaxSigTrans);
		dHist_SigTrans_SigLong_Jpsi[loc_i] = new TH2I(Form("%sSigTrans_SigLong_Jpsi",locParticleLabel[loc_i].Data()), Form(";%s #sigma_{long}; %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, locMinSigLong, locMaxSigLong, 100, locMinSigTrans, locMaxSigTrans);
		dHist_SigLong_SigTheta_Jpsi[loc_i] = new TH2I(Form("%sSigLong_SigTheta_Jpsi",locParticleLabel[loc_i].Data()), Form(";%s #sigma_{theta}; %s #sigma_{long}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, locMinSigTheta, locMaxSigTheta, 100, locMinSigLong, locMaxSigLong);
		
		
		dHist_SigTrans_Preshower_Bkgd[loc_i] = new TH2I(Form("%sSigTrans_Preshower_Bkgd",locParticleLabel[loc_i].Data()), Form(";%s E(Preshower)/E(Shower); %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 200, locMinPreshowerEoverShowerE, locMaxPreshowerEoverShowerE, 100, locMinSigTrans, locMaxSigTrans);
		dHist_SigTheta_Preshower_Bkgd[loc_i] = new TH2I(Form("%sSigTheta_Preshower_Bkgd",locParticleLabel[loc_i].Data()), Form(";%s E(Preshower)/E(Shower); %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 200, locMinPreshowerEoverShowerE, locMaxPreshowerEoverShowerE, 100, locMinSigTheta, locMaxSigTheta);
		dHist_SigTrans_SigTheta_Bkgd[loc_i] = new TH2I(Form("%sSigTrans_SigTheta_Bkgd",locParticleLabel[loc_i].Data()), Form(";%s #sigma_{theta}; %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, locMinSigTheta, locMaxSigTheta, 100, locMinSigTrans, locMaxSigTrans);
		dHist_SigTrans_SigLong_Bkgd[loc_i] = new TH2I(Form("%sSigTrans_SigLong_Bkgd",locParticleLabel[loc_i].Data()), Form(";%s #sigma_{long}; %s #sigma_{trans}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, locMinSigLong, locMaxSigLong, 100, locMinSigTrans, locMaxSigTrans);
		dHist_SigLong_SigTheta_Bkgd[loc_i] = new TH2I(Form("%sSigLong_SigTheta_Bkgd",locParticleLabel[loc_i].Data()), Form(";%s #sigma_{theta}; %s #sigma_{long}", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, locMinSigTheta, locMaxSigTheta, 100, locMinSigLong, locMaxSigLong);

		dHistFDCdEdx_P_Jpsi[loc_i] = new TH2I(Form("%sFDC_dEdx_P_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s FDC dE/dx (keV/cm)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 120, 0., 12.0, 250, locMindEdx, locMaxdEdx);
                dHistCDCdEdx_P_Jpsi[loc_i] = new TH2I(Form("%sCDC_dEdx_P_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s CDC dE/dx (keV/cm)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 250, locMindEdx, locMaxdEdx);
		dHistCDCdEdx_Theta_Jpsi[loc_i] = new TH2I(Form("%sCDC_dEdx_Theta_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s CDC dE/dx (keV/cm)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 10., 110.0, 250, locMindEdx, locMaxdEdx);
		dHistCDCP_Theta_Jpsi[loc_i] = new TH2I(Form("%sCDC_P_Theta_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s p (GeV)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 10., 110.0, 100, 0., 5.);
                dHistBCALDeltaT_P_Jpsi[loc_i] = new TH2I(Form("%sBCAL_DeltaT_P_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s BCAL #Delta T (ns)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 100, locMinDeltaT, locMaxDeltaT);
                dHistTOFDeltaT_P_Jpsi[loc_i] = new TH2I(Form("%sTOF_DeltaT_P_Jpsi",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s TOF #Delta T (ns)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 100, locMinDeltaT, locMaxDeltaT);
		dHistFDCdEdx_P_Bkgd[loc_i] = new TH2I(Form("%sFDC_dEdx_P_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s FDC dE/dx (keV/cm)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 120, 0., 12.0, 250, locMindEdx, locMaxdEdx);
                dHistCDCdEdx_P_Bkgd[loc_i] = new TH2I(Form("%sCDC_dEdx_P_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s CDC dE/dx (keV/cm)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 250, locMindEdx, locMaxdEdx);
		dHistCDCdEdx_Theta_Bkgd[loc_i] = new TH2I(Form("%sCDC_dEdx_Theta_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s CDC dE/dx (keV/cm)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 10., 110.0, 250, locMindEdx, locMaxdEdx);
		dHistCDCP_Theta_Bkgd[loc_i] = new TH2I(Form("%sCDC_P_Theta_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s #theta (degrees); %s p (GeV)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 10., 110.0, 100, 0., 5.);
                dHistBCALDeltaT_P_Bkgd[loc_i] = new TH2I(Form("%sBCAL_DeltaT_P_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s BCAL #Delta T (ns)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 100, locMinDeltaT, locMaxDeltaT);
                dHistTOFDeltaT_P_Bkgd[loc_i] = new TH2I(Form("%sTOF_DeltaT_P_Bkgd",locParticleLabel[loc_i].Data()), Form("; %s p (GeV); %s TOF #Delta T (ns)", locParticleLabel[loc_i].Data(), locParticleLabel[loc_i].Data()), 100, 0., 5.0, 100, locMinDeltaT, locMaxDeltaT);
	}

	locDir->cd();
	dHist_EleBDT_PosBDT = new TH2I("EleBDT_PosBDT", "; Positron BDT; Electron BDT", 100, -1., 1., 100, -1., 1.);
	dHist_Mee_BDT = new TH2I("Mee_BDT", "; e^{#pm} BDT; Mass e^{+}e^{-} (GeV/c^{2})", 100, -1., 1., 400, locMinMee, locMaxMee);
	dHist_Mee_1BCAL = new TH1I("Mee_1BCAL", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_2BCAL = new TH1I("Mee_2BCAL", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_1BDT = new TH1I("Mee_1BDT", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);
	dHist_Mee_2BDT = new TH1I("Mee_2BDT", ";Mass e^{+}e^{-} (GeV/c^{2})", 400, locMinMee, locMaxMee);

	/***************************************** ADVANCED: CHOOSE BRANCHES TO READ ****************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	reader = new TMVA::Reader( "!Color:!Silent" );

	reader->AddVariable( "pos_E2p := pos_E/pos_p", &E2p);
	reader->AddVariable( "pos_Epre", &Epre );
	reader->AddVariable( "pos_p", &p );
	reader->AddVariable( "pos_theta", &theta );
	reader->AddVariable( "pos_dEdx", &dEdx );
	reader->AddVariable( "pos_sigLong", &sigLong );
	reader->AddVariable( "pos_sigTrans", &sigTrans );
	reader->AddVariable( "pos_sigTheta", &sigTheta );
	reader->AddVariable( "pos_deltaT",  &deltaT );

	reader->BookMVA( "BDT method", "/sciclone/home10/jrstevens01/2016-spring-ana/analysis/tcs/selector/dataset/weights/TMVAClassification_BDT.weights.xml" );

}

Bool_t DSelector_tcs::Process(Long64_t locEntry)
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
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locElectronTrackID = dElectronWrapper->Get_TrackID();
		Int_t locPositronTrackID = dPositronWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locElectronP4 = dElectronWrapper->Get_P4();
		TLorentzVector locPositronP4 = dPositronWrapper->Get_P4();
		TLorentzVector locDecayingJpsiP4 = locElectronP4 + locPositronP4;

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locElectronP4_Measured = dElectronWrapper->Get_P4_Measured();
		TLorentzVector locPositronP4_Measured = dPositronWrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/
		
		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locProtonP4_Measured + locElectronP4_Measured + locPositronP4_Measured;

		double locMee = (locElectronP4 + locPositronP4).M();

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/******************************************** DO YOUR STUFF HERE *******************************************/
		//////////////////////////////////////////////
		// E/p cuts (should be FCAL/BCAL dependent) //
		//////////////////////////////////////////////
		double locElectronFCALE = dElectronWrapper->Get_Energy_FCAL();
		double locElectronBCALE = dElectronWrapper->Get_Energy_BCAL();
		double locElectronBCALPreshowerE = dElectronWrapper->Get_Energy_BCALPreshower()/locElectronBCALE*sin(locElectronP4.Theta());
		double locElectronSigLongBCALShower = dElectronWrapper->Get_SigLong_BCAL();
		double locElectronSigThetaBCALShower = dElectronWrapper->Get_SigTheta_BCAL();
		double locElectronSigTransBCALShower = dElectronWrapper->Get_SigTrans_BCAL();

		double locPositronFCALE = dPositronWrapper->Get_Energy_FCAL();
		double locPositronBCALE = dPositronWrapper->Get_Energy_BCAL();
		double locPositronBCALPreshowerE = dPositronWrapper->Get_Energy_BCALPreshower()/locPositronBCALE*sin(locPositronP4.Theta());
		double locPositronSigLongBCALShower = dPositronWrapper->Get_SigLong_BCAL();
		double locPositronSigThetaBCALShower = dPositronWrapper->Get_SigTheta_BCAL();
		double locPositronSigTransBCALShower = dPositronWrapper->Get_SigTrans_BCAL();

		double locElectronEOverP = locElectronBCALE > 0. ? locElectronBCALE/locElectronP4.P() : locElectronFCALE/locElectronP4.P();
		double locPositronEOverP = locPositronBCALE > 0. ? locPositronBCALE/locPositronP4.P() : locPositronFCALE/locPositronP4.P();

		// PID plots after conversion cuts:
                //    dE/dx, DeltaT, E/p, Preshower vs p and theta
                double locElectronP = locElectronP4.Vect().Mag();
                double locPositronP = locPositronP4.Vect().Mag();
                double locElectronTheta = locElectronP4.Theta()*180./TMath::Pi();
                double locPositronTheta = locPositronP4.Theta()*180./TMath::Pi();
                double locElectronFDCdEdx = dElectronWrapper->Get_dEdx_FDC()*1e6;
                double locPositronFDCdEdx = dPositronWrapper->Get_dEdx_FDC()*1e6;
                double locElectronCDCdEdx = dElectronWrapper->Get_dEdx_CDC()*1e6;
                double locPositronCDCdEdx = dPositronWrapper->Get_dEdx_CDC()*1e6;
                double locRFTime = dComboWrapper->Get_RFTime_Measured();
                TVector3 locVertex = dPositronWrapper->Get_X4().Vect();
                double locTargetCenterZ = dComboWrapper->Get_TargetCenter().Z();
                double locPropagatedRFTime = locRFTime + (locVertex.Z() - locTargetCenterZ)/29.9792458;
		double locElectronDeltaT = dElectronWrapper->Get_X4().T() - locPropagatedRFTime;
                double locPositronDeltaT = dPositronWrapper->Get_X4().T() - locPropagatedRFTime;
		
		// require both tracks have showers in calorimeter
		if(!(locElectronBCALE>0. || locElectronFCALE>0.) || !(locPositronBCALE>0. || locPositronFCALE>0.))
			continue;

		// CDC dE/dx cut (temporary)
		double locElectrondEdxCut = 2.2 - 0.00015*locElectronTheta*locElectronTheta;
		double locPositrondEdxCut = 2.2 - 0.00015*locPositronTheta*locPositronTheta;
		//if(!((locElectronCDCdEdx > locElectrondEdxCut && locElectronTheta < 90.) || (locPositronCDCdEdx > locPositrondEdxCut && locPositronTheta < 90.) || (locElectronCDCdEdx < 0.001 && locPositronCDCdEdx < 0.001)))
		//	continue;
		
		// BCAL shower shape cut (temporary)
		double locElectronShowerShapeCut = 0.018 - 0.035*pow((locElectronBCALPreshowerE - 0.5), 2);
		double locPositronShowerShapeCut = 0.018 - 0.035*pow((locPositronBCALPreshowerE - 0.5), 2);
		//if(!((locElectronSigTransBCALShower < locElectronShowerShapeCut) ||  (locPositronSigTransBCALShower < locPositronShowerShapeCut) || (locElectronBCALE < 0.001 && locPositronBCALE < 0.001)))
		//	continue;

		double locMinBkgdMee = 2.0;
		double locMaxBkgdMee = 3.0;
		double locMinJpsiMee = 3.07;
		double locMaxJpsiMee = 3.12;

		// evaluate MVA
		E2p=locElectronEOverP, Epre=locElectronBCALPreshowerE, p=locElectronP, theta=locElectronTheta, dEdx=locElectronCDCdEdx, sigLong=locElectronSigLongBCALShower, sigTrans=locElectronSigTransBCALShower, sigTheta=locElectronSigThetaBCALShower, deltaT=locElectronDeltaT;
		double locElectronBDT = 1.0; //reader->EvaluateMVA("BDT method");

		E2p=locPositronEOverP, Epre=locPositronBCALPreshowerE, p=locPositronP, theta=locPositronTheta, dEdx=locPositronCDCdEdx, sigLong=locPositronSigLongBCALShower, sigTrans=locPositronSigTransBCALShower, sigTheta=locPositronSigThetaBCALShower, deltaT=locPositronDeltaT;
		double locPositronBDT = 1.0; //reader->EvaluateMVA("BDT method");

		if(locPositronBCALE > 0. || locElectronBCALE > 0.) {
			dHist_Mee_1BCAL->Fill(locMee);
			if(locPositronBCALE > 0. && locElectronEOverP > 0.8) 
				dHist_Mee_BDT->Fill(locPositronBDT, locMee);
			if(locElectronBCALE > 0. && locPositronEOverP > 0.8) 
				dHist_Mee_BDT->Fill(locElectronBDT, locMee);
		}
		if(locPositronBCALE > 0. && locElectronBCALE > 0.) {
			dHist_Mee_2BCAL->Fill(locMee);
			dHist_EleBDT_PosBDT->Fill(locPositronBDT, locElectronBDT);
		}

		// mass distributions with both tagged
		if((locPositronBDT > 0.1 && locElectronEOverP > 0.8) || (locElectronBDT > 0.1 && locPositronEOverP > 0.8))
			dHist_Mee_1BDT->Fill(locMee);
		if(locPositronBDT > 0.1 && locElectronBDT > 0.1)
			dHist_Mee_2BDT->Fill(locMee);

		// add J/psi and BH mass regions to check on selection...

		// fill histograms
		dHist_EOverP_Mee[0]->Fill(locMee, locElectronEOverP);
		dHist_EOverP_Mee[1]->Fill(locMee, locPositronEOverP);
		if(locMee > locMinJpsiMee && locMee < locMaxJpsiMee) {
			dHist_EleEOverP_PosEOverP_Jpsi->Fill(locPositronEOverP, locElectronEOverP);
		}
		else if(locMee > locMinBkgdMee && locMee < locMaxBkgdMee) {
                	dHist_EleEOverP_PosEOverP_Bkgd->Fill(locPositronEOverP, locElectronEOverP);           
		}

		if(locElectronBCALE > 0. && locPositronEOverP > 0.8) {
			if(locMee > locMinJpsiMee && locMee < locMaxJpsiMee) {
				dHist_EOverP_Preshower_Jpsi[0]->Fill(locElectronBCALPreshowerE, locElectronEOverP);
				dHistEOverP_P_BCAL_Jpsi[0]->Fill(locElectronP, locElectronEOverP);
				dHistEOverP_Theta_BCAL_Jpsi[0]->Fill(locElectronTheta, locElectronEOverP);
			}
			else if(locMee > locMinBkgdMee && locMee < locMaxBkgdMee) {
				dHist_EOverP_Preshower_Bkgd[0]->Fill(locElectronBCALPreshowerE, locElectronEOverP);
				dHistEOverP_P_BCAL_Bkgd[0]->Fill(locElectronP, locElectronEOverP);
				dHistEOverP_Theta_BCAL_Bkgd[0]->Fill(locElectronTheta, locElectronEOverP);
			}
		}
		if(locPositronBCALE > 0. && locElectronEOverP > 0.8) {
			if(locMee > locMinJpsiMee && locMee < locMaxJpsiMee) {
				dHist_EOverP_Preshower_Jpsi[1]->Fill(locPositronBCALPreshowerE, locPositronEOverP);
				dHistEOverP_P_BCAL_Jpsi[1]->Fill(locPositronP, locPositronEOverP);
				dHistEOverP_Theta_BCAL_Jpsi[1]->Fill(locPositronTheta, locPositronEOverP);
			}
			else if(locMee > locMinBkgdMee && locMee < locMaxBkgdMee) {
				dHist_EOverP_Preshower_Bkgd[1]->Fill(locPositronBCALPreshowerE, locPositronEOverP);
				dHistEOverP_P_BCAL_Bkgd[1]->Fill(locPositronP, locPositronEOverP);
				dHistEOverP_Theta_BCAL_Bkgd[1]->Fill(locPositronTheta, locPositronEOverP);
			}
		}		
		if(locElectronFCALE > 0. && locPositronEOverP > 0.8) {
			if(locMee > locMinJpsiMee && locMee < locMaxJpsiMee) {
				dHistEOverP_P_FCAL_Jpsi[0]->Fill(locElectronP, locElectronEOverP);
				dHistEOverP_Theta_FCAL_Jpsi[0]->Fill(locElectronTheta, locElectronEOverP);
				if(locElectronEOverP > 0.8 && locElectronEOverP < 1.1) {
					dHistFDCdEdx_P_Jpsi[0]->Fill(locElectronP, locElectronFDCdEdx);
					dHistTOFDeltaT_P_Jpsi[0]->Fill(locElectronP, locElectronDeltaT);
				}
			}
			else if(locMee > locMinBkgdMee && locMee < locMaxBkgdMee) {
				dHistEOverP_P_FCAL_Bkgd[0]->Fill(locElectronP, locElectronEOverP);
				dHistEOverP_Theta_FCAL_Bkgd[0]->Fill(locElectronTheta, locElectronEOverP);
				if(locElectronEOverP > 0.8 && locElectronEOverP < 1.1) {
					dHistFDCdEdx_P_Bkgd[0]->Fill(locElectronP, locElectronFDCdEdx);
					dHistTOFDeltaT_P_Bkgd[0]->Fill(locElectronP, locElectronDeltaT);
				}
			}
		}
		if(locPositronFCALE > 0. && locElectronEOverP > 0.8) {
			if(locMee > locMinJpsiMee && locMee < locMaxJpsiMee) {
				dHistEOverP_P_FCAL_Jpsi[1]->Fill(locPositronP, locPositronEOverP);
				dHistEOverP_Theta_FCAL_Jpsi[1]->Fill(locPositronTheta, locPositronEOverP);
				if(locPositronEOverP > 0.8 && locPositronEOverP < 1.1) {
					dHistFDCdEdx_P_Jpsi[1]->Fill(locPositronP, locPositronFDCdEdx);
					dHistTOFDeltaT_P_Jpsi[1]->Fill(locPositronP, locPositronDeltaT);
				}
			}
			else if(locMee > locMinBkgdMee && locMee < locMaxBkgdMee) {
				dHistEOverP_P_FCAL_Bkgd[1]->Fill(locPositronP, locPositronEOverP);
				dHistEOverP_Theta_FCAL_Bkgd[1]->Fill(locPositronTheta, locPositronEOverP);
				if(locPositronEOverP > 0.8 && locPositronEOverP < 1.1) {
					dHistFDCdEdx_P_Bkgd[0]->Fill(locPositronP, locPositronFDCdEdx);
					dHistTOFDeltaT_P_Bkgd[0]->Fill(locPositronP, locPositronDeltaT);
				}
			}
		}

		double minEOverP = 0.8;
		if(locElectronEOverP < minEOverP || locPositronEOverP < minEOverP) {
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
	
		dHist_P_Mee[0]->Fill(locMee, locElectronP);
		dHist_P_Mee[1]->Fill(locMee, locPositronP);
		dHist_Theta_Mee[0]->Fill(locMee, locElectronTheta);
		dHist_Theta_Mee[1]->Fill(locMee, locPositronTheta);	

		//////////////////////////////////////////
		// BCAL preshower and shower width cuts //
		//////////////////////////////////////////
		double locMinPreshowerE = 0.05;
		double locMaxSigTrans = 0.0065;
		double locMaxSigTheta = 0.0025;

		if(locElectronBCALE > 0.) {
			if(locMee > locMinJpsiMee && locMee < locMaxJpsiMee) {
				dHist_SigTrans_Preshower_Jpsi[0]->Fill(locElectronBCALPreshowerE, locElectronSigTransBCALShower);
				dHist_SigTheta_Preshower_Jpsi[0]->Fill(locElectronBCALPreshowerE, locElectronSigThetaBCALShower);
				dHist_SigTrans_SigTheta_Jpsi[0]->Fill(locElectronSigThetaBCALShower, locElectronSigTransBCALShower);
				dHist_SigTrans_SigLong_Jpsi[0]->Fill(locElectronSigLongBCALShower, locElectronSigTransBCALShower);
				dHist_SigLong_SigTheta_Jpsi[0]->Fill(locElectronSigThetaBCALShower, locElectronSigLongBCALShower);
				dHistFDCdEdx_P_Jpsi[0]->Fill(locElectronP, locElectronFDCdEdx);
				dHistCDCdEdx_P_Jpsi[0]->Fill(locElectronP, locElectronCDCdEdx);
				dHistCDCdEdx_Theta_Jpsi[0]->Fill(locElectronTheta, locElectronCDCdEdx);
				dHistCDCP_Theta_Jpsi[0]->Fill(locElectronTheta, locElectronP);
				dHistBCALDeltaT_P_Jpsi[0]->Fill(locElectronP, locElectronDeltaT);
			}
			else if(locMee > locMinBkgdMee && locMee < locMaxBkgdMee) {
				dHist_SigTrans_Preshower_Bkgd[0]->Fill(locElectronBCALPreshowerE, locElectronSigTransBCALShower);
				dHist_SigTheta_Preshower_Bkgd[0]->Fill(locElectronBCALPreshowerE, locElectronSigThetaBCALShower);
				dHist_SigTrans_SigTheta_Bkgd[0]->Fill(locElectronSigThetaBCALShower, locElectronSigTransBCALShower);
				dHist_SigTrans_SigLong_Bkgd[0]->Fill(locElectronSigLongBCALShower, locElectronSigTransBCALShower);
				dHist_SigLong_SigTheta_Bkgd[0]->Fill(locElectronSigThetaBCALShower, locElectronSigLongBCALShower);
				dHistFDCdEdx_P_Bkgd[0]->Fill(locElectronP, locElectronFDCdEdx);
				dHistCDCdEdx_P_Bkgd[0]->Fill(locElectronP, locElectronCDCdEdx);
				dHistCDCdEdx_Theta_Bkgd[0]->Fill(locElectronTheta, locElectronCDCdEdx);
				dHistCDCP_Theta_Bkgd[0]->Fill(locElectronTheta, locElectronP);
				dHistBCALDeltaT_P_Bkgd[0]->Fill(locElectronP, locElectronDeltaT);
			}

			//if(locElectronBCALPreshowerE < locMinPreshowerE) continue;
			//if(locElectronSigTransBCALShower > locMaxSigTrans && locElectronSigThetaBCALShower > locMaxSigTheta) continue;
		}
		if(locPositronBCALE > 0.) {
			if(locMee > locMinJpsiMee && locMee < locMaxJpsiMee) {
				dHist_SigTrans_Preshower_Jpsi[1]->Fill(locPositronBCALPreshowerE, locPositronSigTransBCALShower);
				dHist_SigTheta_Preshower_Jpsi[1]->Fill(locPositronBCALPreshowerE, locPositronSigThetaBCALShower);
				dHist_SigTrans_SigTheta_Jpsi[1]->Fill(locPositronSigThetaBCALShower, locPositronSigTransBCALShower);
				dHist_SigTrans_SigLong_Jpsi[1]->Fill(locPositronSigLongBCALShower, locPositronSigTransBCALShower);
				dHist_SigLong_SigTheta_Jpsi[1]->Fill(locPositronSigThetaBCALShower, locPositronSigLongBCALShower);
				dHistFDCdEdx_P_Jpsi[1]->Fill(locPositronP, locPositronFDCdEdx);
				dHistCDCdEdx_P_Jpsi[1]->Fill(locPositronP, locPositronCDCdEdx);
				dHistCDCdEdx_Theta_Jpsi[1]->Fill(locPositronTheta, locPositronCDCdEdx);
				dHistCDCP_Theta_Jpsi[1]->Fill(locPositronTheta, locPositronP);
				dHistBCALDeltaT_P_Jpsi[1]->Fill(locPositronP, locPositronDeltaT);
			}
			else if(locMee > locMinBkgdMee && locMee < locMaxBkgdMee) {
				dHist_SigTrans_Preshower_Bkgd[1]->Fill(locPositronBCALPreshowerE, locPositronSigTransBCALShower);
				dHist_SigTheta_Preshower_Bkgd[1]->Fill(locPositronBCALPreshowerE, locPositronSigThetaBCALShower);
				dHist_SigTrans_SigTheta_Bkgd[1]->Fill(locPositronSigThetaBCALShower, locPositronSigTransBCALShower);
				dHist_SigTrans_SigLong_Bkgd[1]->Fill(locPositronSigLongBCALShower, locPositronSigTransBCALShower);
				dHist_SigLong_SigTheta_Bkgd[1]->Fill(locPositronSigThetaBCALShower, locPositronSigLongBCALShower);
				dHistFDCdEdx_P_Bkgd[1]->Fill(locPositronP, locPositronFDCdEdx);
				dHistCDCdEdx_P_Bkgd[1]->Fill(locPositronP, locPositronCDCdEdx);
				dHistCDCdEdx_Theta_Bkgd[1]->Fill(locPositronTheta, locPositronCDCdEdx);
				dHistCDCP_Theta_Bkgd[1]->Fill(locPositronTheta, locPositronP);
				dHistBCALDeltaT_P_Bkgd[1]->Fill(locPositronP, locPositronDeltaT);
			}

			//if(locPositronBCALPreshowerE < locMinPreshowerE) continue;
			//if(locPositronSigTransBCALShower > locMaxSigTrans && locPositronSigThetaBCALShower > locMaxSigTheta) continue;
		}

		//if(locElectronBCALE > 0. || locPositronBCALE > 0.)
		dHist_Mee->Fill(locMee);

		// CDC dE/dx cut
		//double locElectrondEdxCut = 2.2 - 0.00015*locElectronTheta*locElectronTheta;
		//double locPositrondEdxCut = 2.2 - 0.00015*locPositronTheta*locPositronTheta;
		if((locElectronCDCdEdx > locElectrondEdxCut && locElectronTheta < 90.) || (locPositronCDCdEdx > locPositrondEdxCut && locPositronTheta < 90.) || (locElectronCDCdEdx < 0.001 && locPositronCDCdEdx < 0.001)) {
			dHist_Mee_CDCdEdx->Fill(locMee);
			dHist_P_Mee_CDCdEdx[0]->Fill(locMee, locElectronP);
	                dHist_P_Mee_CDCdEdx[1]->Fill(locMee, locPositronP);
        	        dHist_Theta_Mee_CDCdEdx[0]->Fill(locMee, locElectronTheta);
                	dHist_Theta_Mee_CDCdEdx[1]->Fill(locMee, locPositronTheta);
		}
		else continue;

		if(((locElectronSigTransBCALShower < locElectronShowerShapeCut) ||  (locPositronSigTransBCALShower < locPositronShowerShapeCut) || (locElectronBCALE < 0.001 && locPositronBCALE < 0.001)))
			dHist_Mee_ShowerCut->Fill(locMee);
		else continue;

		double locKinFitCL = dComboWrapper->Get_ConfidenceLevel_KinFit();
		if(locKinFitCL > 1e-4) dHist_Mee_KinFit4->Fill(locMee);
		if(locKinFitCL > 1e-3) dHist_Mee_KinFit3->Fill(locMee);
		if(locKinFitCL > 1e-2) dHist_Mee_KinFit2->Fill(locMee);
		if(locKinFitCL > 1e-1) dHist_Mee_KinFit1->Fill(locMee);

	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	return kTRUE;
}

void DSelector_tcs::Finalize(void)
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
