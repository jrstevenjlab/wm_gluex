#include "DSelector_npi_missing.h"

void DSelector_npi_missing::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "npi_missing.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
        dFlatTreeName = "flat_npi_missing"; //if blank, default name will be chosen

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
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.4, PiPlus, SYS_TOF));
	dAnalysisActions.push_back(new DCutAction_NoPIDHit(dComboWrapper, PiPlus, SYS_TOF));

	//MASSES
	dAnalysisActions.push_back(new DHistogramAction_MissingMass(dComboWrapper, false, 300, 0.0, 3.0, "MissingNeutron"));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.1));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	dHist_BeamDeltaT = new TH1I("BeamDeltaT", ";#Delta t (ns)", 200, -10.02, 10.02);
	
	double locPiPlusThetaMin = 0.; double locPiPlusThetaMax = 60.;
        double locPiPlusPMin = 0.; double locPiPlusPMax = 12.;
        double locNeutronThetaMin = 0.; double locNeutronThetaMax = 180.;
        double locNeutronPMin = 0.; double locNeutronPMax = 6.;
	double locDeltaPhiMin = 120.; double locDeltaPhiMax = 240.;

        TString locRFBinLabel[2] = {"", "Acci_"};
        //TString locCutCounterLabel[5] = {"Fiducial", "MM", "UE", "DeltaPhiSignal", "DeltaPhiSideband"};
        //int locDeltaPhiOffset = 3;
        for(int locRFBin=0; locRFBin<2; locRFBin++) {
                dHist_BeamEnergy[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
                dHist_PiPlusP_Theta_Init[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"PiPlusP_Theta_Init", "; #theta_{#pi^{+}} (degrees); p_{#pi^{+}} (GeV/c)", 360., locPiPlusThetaMin, locPiPlusThetaMax, 400, locPiPlusPMin, locPiPlusPMax);
		
		dHist_NeutronP_Theta_Init[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronP_Theta_Init", "; #theta_{n} (degrees); p_{n} (GeV/c)", 360., locNeutronThetaMin, locNeutronThetaMax, 400, locNeutronPMin, locNeutronPMax);
		dHist_NeutronBeta_P_Init[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronBeta_P_Init", "; p_{n} (GeV/c); Neutron #beta", 400, locNeutronPMin, locNeutronPMax, 120, 0., 1.2);
		dHist_NeutronDeltaBeta_P_Init[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronDeltaBeta_P_Init", ";p_{n} (GeV/c); Neutron #Delta#beta", 400, locNeutronPMin, locNeutronPMax, 200, -1., 1.);
		dHist_NeutronBeta_P[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronBeta_P", "; p_{n} (GeV/c); Neutron #beta", 400, locNeutronPMin, locNeutronPMax, 120, 0., 1.2);
		dHist_NeutronDeltaBeta_P[locRFBin]= new TH2I(locRFBinLabel[locRFBin]+"NeutronDeltaBeta_P", ";p_{n} (GeV/c); Neutron #Delta#beta", 400, locNeutronPMin, locNeutronPMax, 200, -1., 1.);
		
                dHist_VertexZ[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"VertexZ", "; z_{vertex} (cm)", 300., 0., 300.);
                dHist_VertexR[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"VertexR", "; r_{vertex} (cm)", 200., 0., 2.);
		
		dHist_MissingMass[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingMass", ";Missing Mass (GeV/c^{2})", 300, 0.0, 3.0);
		dHist_MissingMassSquared[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 100, -1.0, 1.0);
		dHist_MissingMass_Final[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingMass_Final", ";Missing Mass (GeV/c^{2})", 300, 0.0, 3.0);
		dHist_MissingMassSquared_Final[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingMassSquared_Final", ";Missing Mass Squared (GeV/c^{2})^{2}", 100, -1.0, 1.0);
		dHist_MissingMass_Final_noUnused[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingMass_Final_noUnused", ";Missing Mass (GeV/c^{2})", 300, 0.0, 3.0);
		dHist_MissingMassSquared_Final_noUnused[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingMassSquared_Final_noUnused", ";Missing Mass Squared (GeV/c^{2})^{2}", 100, -1.0, 1.0);

		dHist_UnusedTrack[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"UnusedTrack", ";# Unused Tracks", 10, 0., 10.);
		dHist_UnusedEnergyBCAL_t[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"UnusedEnergyBCAL_t", "; -t (GeV^2); BCAL Unused Energy (GeV)", 300, 0., 3., 200, 0.0, 2.0);
                dHist_UnusedEnergyFCAL_t[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"UnusedEnergyFCAL_t", "; -t (GeV^2); FCAL Unused Energy (GeV)", 300, 0., 3., 200, 0.0, 2.0);
                dHist_UnusedEnergyTotal_t[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"UnusedEnergyTotal_t", "; -t (GeV^2); BCAL+FCAL Unused Energy (GeV)", 300, 0., 3., 200, 0.0, 2.0);
		dHist_UnusedEnergyBCAL_FCAL[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"UnusedEnergyBCAL_FCAL", "; FCAL Unused Energy (GeV); BCAL Unused Energy (GeV)", 300, 0., 3., 300, 0.0, 3.0);

		dHist_DeltaPhi_DeltaTheta[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_DeltaTheta", "; #Delta#theta (degrees); #Delta#phi (degrees)", 300, -50., 100., 400, locDeltaPhiMin, locDeltaPhiMax);
		dHist_ShowerTheta_DeltaTheta[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"ShowerTheta_DeltaTheta", "; #Delta#theta (degrees); Neutron Shower Theta (degrees)", 300, -50., 100., 180, 0., 180.);
		dHist_ShowerTheta_DeltaTheta_noUnused[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"ShowerTheta_DeltaTheta_noUnused", "; #Delta#theta (degrees); Neutron Shower Theta (degrees)", 300, -50., 100., 180, 0., 180.);
		dHist_ShowerE_DeltaTheta[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"ShowerE_DeltaTheta", "; #Delta#theta (degrees); Neutron Shower Energy (GeV)", 300, -50., 100., 200, 0., 4.);
			

		dHist_NeutronShowerTheta_DeltaPhi[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronShowerTheta_DeltaPhi", ";#Delta#Phi (degrees); Neutron Shower #theta (deg)", 400, locDeltaPhiMin, locDeltaPhiMax, 180, 0., 180.);
		dHist_NeutronPreshowerE_DeltaPhi[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPreshowerE_DeltaPhi", ";#Delta#Phi (degrees); Neutron Preshower Energy (Gev)", 400, locDeltaPhiMin, locDeltaPhiMax, 100, 0., 0.2);

		dHist_NeutronPreshowerE_ShowerTheta_DeltaPhiCut[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPreshowerE_ShowerTheta_DeltaPhiCut", "; Neutron Shower #theta (deg); Neutron Preshower Energy (Gev)", 180, 0., 180., 100, 0., 0.2);
		
		dHist_NeutronShowerE_ShowerTheta[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronShowerE_ShowerTheta", "; Neutron Shower #theta (deg); Neutron Shower E", 180, 0., 180., 200, 0., 4.);
		dHist_NeutronMissingTheta_ShowerTheta[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronMissingTheta_ShowerTheta", "; Neutron Shower #theta (deg); Missing Neutron #theta (deg)", 180, 0., 180., 180, 0., 180.);
		dHist_NeutronPreshowerE_ShowerE[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPreshowerE_ShowerE", "; Neutron Shower Energy (GeV); Neutron Preshower Energy (Gev)", 100, 0., 1.0, 100, 0., 0.2);

		dHist_PiPlusMomentum_BeamE[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"PiPlusMomentum_BeamE", ";Beam Energy (GeV); #pi^{+} momentum (GeV/c)", 120, 0., 12., 120, 0., 12.);
		dHist_PiPlusMomentum_DeltaPhi[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"PiPlusMomentum_DeltaPhi", ";#Delta#Phi (degrees); #pi^{+} momentum (GeV/c)", 400, locDeltaPhiMin, locDeltaPhiMax, 120, 0., 12.);

		dHist_DeltaPhi_t[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_t", "; -t (GeV^2); #Delta#phi (degrees)", 300, 0., 3., 400, locDeltaPhiMin, locDeltaPhiMax);
		dHist_DeltaPhi_t_Init[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_t_Init", "; -t (GeV^2); #Delta#phi (degrees)", 300, 0., 3., 400, locDeltaPhiMin, locDeltaPhiMax);
		dHist_DeltaPhi_t_DeltaBetaCut[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_t_DeltaBetaCut", "; -t (GeV^2); #Delta#phi (degrees)", 300, 0., 3., 400, locDeltaPhiMin, locDeltaPhiMax);
		dHist_DeltaPhi_t_noUnused[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_t_noUnused", "; -t (GeV^2); #Delta#phi (degrees)", 300, 0., 3., 400, locDeltaPhiMin, locDeltaPhiMax);
		dHist_DeltaPhi_Phi_DeltaBetaCut[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_Phi_DeltaBetaCut", "; #phi_{#pi+}; #Delta#phi (degrees)", 360, -180., 180., 400, locDeltaPhiMin, locDeltaPhiMax);

		dHist_NeutronPhi_t[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPhi_t", ";-t (GeV^{2}); Neutron #phi (deg)", 300, 0., 3., 360, -180., 180.0);
		dHist_NeutronPhi_t_noUnused[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPhi_t_noUnused", ";-t (GeV^{2}); Neutron #phi (deg)", 300, 0., 3., 360, -180., 180.0);
		dHist_NeutronPhi_t_noUnused_tightDeltaTheta[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPhi_t_noUnused_tightDeltaTheta", ";-t (GeV^{2}); Neutron #phi (deg)", 300, 0., 3., 360, -180., 180.0);
		dHist_NeutronPhi_t_noUnused_tightMM[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPhi_t_noUnused_tightMM", ";-t (GeV^{2}); Neutron #phi (deg)", 300, 0., 3., 360, -180., 180.0);
		dHist_NeutronPhi_t_tightMM[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"NeutronPhi_t_tightMM", ";-t (GeV^{2}); Neutron #phi (deg)", 300, 0., 3., 360, -180., 180.0);

		dHist_Sideband1NeutronPhi_t[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"Sideband1NeutronPhi_t", ";-t (GeV^{2}); Neutron #phi (deg)", 300, 0., 3., 360, -180., 180.0);
		dHist_Sideband2NeutronPhi_t[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"Sideband2NeutronPhi_t", ";-t (GeV^{2}); Neutron #phi (deg)", 300, 0., 3., 360, -180., 180.0);

		dHist_Thrownt_Recot[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"Thrownt_Recot", ";Reco -t (GeV^{2}); Thrown -t (GeV^{2})", 300, 0., 3., 300, 0., 3.);
		dHist_PiPlusDeltaPOverP_Theta[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"PiPlusDeltaPOverP_Theta", "; Thrown #theta_{#pi^+} (degrees); #Delta p/p", 100, 0., 10., 100, -0.5, 0.5);
		dHist_ThrownDeltaThetaShower_DeltaThetaMissing[locRFBin] = new TH2I(locRFBinLabel[locRFBin]+"ThrownDeltaThetaShower_DeltaThetaMissing", "; #Delta#theta Missing Neutron (degrees); #Delta#theta Neutron Shower (degrees)", 200, -50., 50., 250, -50., 50.);
	}

	dHist_Thrownt = new TH1I("Thrownt", ";Thrown -t (GeV^{2})", 100, 0., 5.);
	dHist_ThrownNeutronP_t = new TH2I("ThrownNeutronP_t", ";Thrown -t (GeV^{2}); Thrown Neutron Momentum (GeV/c)", 100, 0., 5., 100, 0., 5.);
	dHist_ThrownPiPlusTheta_t = new TH2I("ThrownPiPlusP_t", ";Thrown -t (GeV^{2}); Thrown #pi^{+} #theta (degrees)", 100, 0., 5., 200, 0., 40.);
	dHist_ThrownNeutronP_Theta = new TH2I("ThrownNeutronP_Theta", ";Thrown Neutron #theta (degrees); Thrown Neutron Momentum (GeV/c)", 180, 0., 180., 100, 0., 5.);
	dHist_ThrownPiPlusP_Theta = new TH2I("ThrownPiPlusP_Theta", ";Thrown #pi^{+} #theta (degrees); Thrown #pi^{+} Momentum (GeV/c)", 200, 0., 40., 120, 0., 12.);

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

        //EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("t");
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Egamma");
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("PhiN");
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("DeltaPhi");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("UnusedE");
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Weight");

	// EXAMPLE CUT PARAMETERS:
	dNeutronDeltaBeta = new TF1("NeutronDeltaBeta", "(x/sqrt(x*x + [0]*[0]) - x/sqrt(x*x)) +  [1]", 0., 10.);
	dNeutronDeltaBeta->SetParameter(0, 0.9396);
	dNeutronDeltaBeta->SetParameter(1, 0.2);

	dDeltaPhiCut = 2.0;
	dDeltaThetaMinCut = -20.;
	dDeltaThetaMaxCut = 40;
	dNeutronShowerThetaCut = 40.;
	dBeamELowCut = 8.4;
	dBeamEHighCut = 9.0;

	/***************************************** ADVANCED: CHOOSE BRANCHES TO READ ****************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_npi_missing::Process(Long64_t locEntry)
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

	// Set beam energy range
	if(locRunNumber > 30000 && locRunNumber < 40000) {
		dBeamELowCut = 8.2;
		dBeamEHighCut = 8.8;
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

	double locThrownBeamE = 0; 
	if(Get_NumThrown()>0)
		dThrownBeam->Get_P4().E();

	// loop over thrown particles for resolution studies
	double locThrownNeutronTheta = 0;
	double locThrownNeutronP = 0;
	double locThrownPiPlusTheta = 0;
	double locThrownPiPlusP = 0;
	double locThrownt = 0;
	TLorentzVector locThrownPiPlusP4;
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		dThrownWrapper->Set_ArrayIndex(loc_i);
		if(dThrownWrapper->Get_P4().M() > 0.939) {
			locThrownNeutronTheta = dThrownWrapper->Get_P4().Theta();
			locThrownNeutronP = dThrownWrapper->Get_P4().P();
			locThrownt = -1.*(dTargetP4 - dThrownWrapper->Get_P4()).M2();
		}
		if(dThrownWrapper->Get_P4().M() < 0.142 &&  dThrownWrapper->Get_P4().M() > 0.137) {
			locThrownPiPlusTheta = dThrownWrapper->Get_P4().Theta();
			locThrownPiPlusP = dThrownWrapper->Get_P4().P();
			locThrownPiPlusP4 = dThrownWrapper->Get_P4();
		}
	}
	dHist_Thrownt->Fill(locThrownt);
	dHist_ThrownNeutronP_t->Fill(locThrownt, locThrownNeutronP);
	dHist_ThrownPiPlusTheta_t->Fill(locThrownt, locThrownPiPlusTheta*180./TMath::Pi());
	dHist_ThrownNeutronP_Theta->Fill(locThrownNeutronTheta*180./TMath::Pi(), locThrownNeutronP);
	dHist_ThrownPiPlusP_Theta->Fill(locThrownPiPlusTheta*180./TMath::Pi(), locThrownPiPlusP);

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

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		if(Get_NumThrown()>0)
			locPiPlusP4_Measured = locThrownPiPlusP4;

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		TLorentzVector locPiPlusX4_Measured = dPiPlusWrapper->Get_X4_Measured();

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/********************************************* ACCIDENTALS & RF TIMING ********************************************/

                double locRFTime = dComboWrapper->Get_RFTime_Measured();
                double locBeamDeltaT = locBeamX4_Measured.T() - (locRFTime + (locPiPlusX4_Measured.Z() - dTargetCenter.Z())/29.9792458);
		if(locBeamP4.E() > dBeamELowCut && locBeamP4.E() < dBeamEHighCut)
			dHist_BeamDeltaT->Fill(locBeamDeltaT);

                int locRFBin = -1;
                if(fabs(locBeamDeltaT) < 0.5*4.008) locRFBin = 0;
                else if (fabs(locBeamDeltaT) > 0.5*4.008 && fabs(locBeamDeltaT) < 3.5*4.008) locRFBin = 1;
                else continue;

		/********************************************* KINEMATIC VARIABLES **********************************************/

		// Momentum transfer t
		double loct = -1.*(locBeamP4 - locPiPlusP4).M2();

		/**************************************** HISTOGRAM BEAM ENERGY *****************************************/

                //Histogram beam energy (if haven't already)
                if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
                {
                        dHist_BeamEnergy[locRFBin]->Fill(locBeamP4.E());
                        locUsedSoFar_BeamEnergy.insert(locBeamID);
                }

                // Coherent peak cut
                if(locBeamP4.E() < dBeamELowCut || locBeamP4.E() > dBeamEHighCut)
                        continue;

		/********************************************* FIDUCIAL AND VERTEX CUTS ********************************************/

                double locPiPlusP = locPiPlusP4_Measured.P();
                double locPiPlusTheta = locPiPlusP4_Measured.Theta()*180./TMath::Pi();
                dHist_PiPlusP_Theta_Init[locRFBin]->Fill(locPiPlusTheta, locPiPlusP);
                if(locPiPlusTheta < 1.0) continue;

		double locVertexZ = locPiPlusX4_Measured.Z();
                dHist_VertexZ[locRFBin]->Fill(locVertexZ);
                if(locVertexZ < 45. || locVertexZ > 81.)
                        continue;

                double locVertexR = locPiPlusX4_Measured.Perp();
                dHist_VertexR[locRFBin]->Fill(locVertexR);
                if(locVertexR > 1.)
                        continue;
               
		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPiPlusP4_Measured;

		//Missing Mass
		double locMissingMass = locMissingP4_Measured.M();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMass[locRFBin]->Fill(locMissingMass);
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		if(locMissingMass < 0.3 || locMissingMass > 2.1)
			continue;

		// unused track veto
		int locUnusedTrack = 0;
		for(UInt_t loc_j = 0; loc_j < Get_NumChargedHypos(); ++loc_j) {
			//Set branch array indices corresponding to this particle
			dChargedHypoWrapper->Set_ArrayIndex(loc_j);

			// skip event if there is another physical track
			if(dChargedHypoWrapper->Get_TrackID() != locPiPlusTrackID) {
				locUnusedTrack++;
			}
		}
		dHist_UnusedTrack[locRFBin]->Fill(locUnusedTrack);
		if(locUnusedTrack > 0) 
			continue;

		// correlation of pi+ with neutron candidates
		for(UInt_t loc_k = 0; loc_k < Get_NumNeutralHypos(); ++loc_k) {
			
			//Set branch array indices corresponding to this particle
			dNeutralHypoWrapper->Set_ArrayIndex(loc_k);
			UInt_t locNeutralPID = dNeutralHypoWrapper->Get_PID();
			Int_t locNeutronNeutralID = dNeutralHypoWrapper->Get_NeutralID();
			if(locNeutralPID != Neutron)
				continue;

			TLorentzVector locNeutronShowerX4 = dNeutralHypoWrapper->Get_X4_Shower();
			TLorentzVector locNeutronShowerP4 = dNeutralHypoWrapper->Get_P4();
			double locNeutronShowerE = dNeutralHypoWrapper->Get_Energy_BCAL();
			double locNeutronPreshowerE = dNeutralHypoWrapper->Get_Energy_BCALPreshower();
			if(locNeutronShowerE <= 0.) {
				locNeutronShowerE = dNeutralHypoWrapper->Get_Energy_FCAL();
				locNeutronPreshowerE = 0.;
				continue; // skip FCAL showers for now
			}
			dHist_NeutronP_Theta_Init[locRFBin]->Fill(locNeutronShowerP4.Theta()*180./TMath::Pi(), locMissingP4_Measured.P());
			
			double locNeutronMass = 0.9396;
			double locNeutronDeltaT = locNeutronShowerX4.T() - (locRFTime + (locVertexZ - dTargetCenter.Z())/29.9792458);
			double locNeutronPathLength = (locNeutronShowerX4.Vect() - dTargetCenter).Mag();
			double locNeutronBeta = locNeutronPathLength/(locNeutronDeltaT*29.9792458);
			double locNeutronDeltaBeta = locMissingP4_Measured.P() / sqrt(pow(locMissingP4_Measured.P(),2) + locNeutronMass*locNeutronMass) - locNeutronBeta;
			dHist_NeutronBeta_P_Init[locRFBin]->Fill(locMissingP4_Measured.P(), locNeutronBeta);
			dHist_NeutronDeltaBeta_P_Init[locRFBin]->Fill(locMissingP4_Measured.P(), locNeutronDeltaBeta);
			
			// calculated neutron momentum from beta and then MM2
			double locNeutronMomentum_BetaMeasured = sqrt(0.9396*0.9396/(1-locNeutronBeta*locNeutronBeta));
			TVector3 locNeutronP3_BetaMeasured = locNeutronShowerX4.Vect() - dTargetCenter;
			locNeutronP3_BetaMeasured.SetMag(locNeutronMomentum_BetaMeasured);
			TLorentzVector locNeutronP4_BetaMeasured; 
			locNeutronP4_BetaMeasured.SetVectM(locNeutronP3_BetaMeasured, 0.9396);
			TLorentzVector locMissingP4_BetaMeasured = locBeamP4_Measured + dTargetP4;
			locMissingP4_BetaMeasured -= (locPiPlusP4_Measured + locNeutronP4_BetaMeasured);
			double locMissingMassSquared = locMissingP4_BetaMeasured.M2();
			dHist_MissingMassSquared[locRFBin]->Fill(locMissingMassSquared);

			// delta phi between pi+ track and neutron shower
			double locDeltaPhi = (locNeutronShowerP4.Phi() - locPiPlusP4.Phi())*180./TMath::Pi();
			if(locDeltaPhi>360.) locDeltaPhi -= 360.;
			if(locDeltaPhi<0.) locDeltaPhi += 360.;
			double locDeltaTheta = (locNeutronShowerP4.Theta() - locMissingP4_Measured.Theta())*180./TMath::Pi();

			// correct DeltaPhi for alignment
			double p0 = 1.0;
			double p1 = 0.126;
			double p2 = 179.77;
			double locDeltaPhiCorrected = locDeltaPhi - p0*sin(locNeutronShowerP4.Phi() + p1);

			// unused shower energy sum
			double unusedEnergyBCAL = 0.;
			double unusedEnergyFCAL = 0.;
			for(UInt_t loc_j = 0; loc_j < Get_NumNeutralHypos(); ++loc_j) {
				
				//Set branch array indices corresponding to this particle
				dNeutralHypoWrapper->Set_ArrayIndex(loc_j);
				UInt_t locPhotonNeutralPID = dNeutralHypoWrapper->Get_PID();
				if(locPhotonNeutralPID != Gamma)
					continue;
				if(fabs(dNeutralHypoWrapper->Get_Energy_BCAL()-locNeutronShowerE) < 0.01 || fabs(dNeutralHypoWrapper->Get_Energy_FCAL()-locNeutronShowerE) < 0.01) {
					//cout<<"same unused as neutron candidate"<<endl;
					continue;
				}

				unusedEnergyBCAL += dNeutralHypoWrapper->Get_Energy_BCAL();
				unusedEnergyFCAL += dNeutralHypoWrapper->Get_Energy_FCAL();
			}
		       
			// plot unused energy after all other signal cuts
			if((locDeltaTheta > dDeltaThetaMinCut && locDeltaTheta < dDeltaThetaMaxCut) && fabs(locDeltaPhiCorrected - 180.) < dDeltaPhiCut && locNeutronShowerP4.Theta()*180./TMath::Pi() > dNeutronShowerThetaCut && locNeutronDeltaBeta > dNeutronDeltaBeta->Eval(locMissingP4_Measured.P())) {
				dHist_UnusedEnergyBCAL_t[locRFBin]->Fill(loct, unusedEnergyBCAL);
				dHist_UnusedEnergyFCAL_t[locRFBin]->Fill(loct, unusedEnergyFCAL);
				dHist_UnusedEnergyTotal_t[locRFBin]->Fill(loct, unusedEnergyBCAL+unusedEnergyFCAL);
				dHist_UnusedEnergyBCAL_FCAL[locRFBin]->Fill(unusedEnergyFCAL, unusedEnergyBCAL);
			}

			dHist_DeltaPhi_t_Init[locRFBin]->Fill(loct, locDeltaPhiCorrected);

			// compare to thrown neutron 
			if(Get_NumThrown() > 0) {
				double locDeltaThetaShowerThrown = (locNeutronShowerP4.Theta()-locThrownNeutronTheta)*180./TMath::Pi();
				double locDeltaThetaMissingThrown = (locMissingP4_Measured.Theta()-locThrownNeutronTheta)*180./TMath::Pi();
				dHist_ThrownDeltaThetaShower_DeltaThetaMissing[locRFBin]->Fill(locDeltaThetaMissingThrown, locDeltaThetaShowerThrown);
			}

			// neutron kinematic distribution with later cuts applied
			if(locNeutronShowerP4.Theta()*180./TMath::Pi() > dNeutronShowerThetaCut && locNeutronDeltaBeta > dNeutronDeltaBeta->Eval(locMissingP4_Measured.P()))
				dHist_DeltaPhi_DeltaTheta[locRFBin]->Fill(locDeltaTheta, locDeltaPhiCorrected);

			
			// neutron shower distribution with DeltaPhiCut applied
			if(fabs(locDeltaPhiCorrected - 180.) < dDeltaPhiCut && locNeutronDeltaBeta > dNeutronDeltaBeta->Eval(locMissingP4_Measured.P())) {
				//if(fabs(locThrownBeamE - locBeamP4.E()) < 0.020) { //dComboBeamWrapper->Get_IsGenerator()) {
				dHist_ShowerTheta_DeltaTheta[locRFBin]->Fill(locDeltaTheta, locNeutronShowerP4.Theta()*180./TMath::Pi());
				dHist_ShowerE_DeltaTheta[locRFBin]->Fill(locDeltaTheta, locNeutronShowerE);
				//}
				if(unusedEnergyBCAL < 0.01) 
					dHist_ShowerTheta_DeltaTheta_noUnused[locRFBin]->Fill(locDeltaTheta, locNeutronShowerP4.Theta()*180./TMath::Pi());
			}

			// may be possible to cut tighter and remove more accidentals?
			if(locDeltaTheta < dDeltaThetaMinCut || locDeltaTheta > dDeltaThetaMaxCut)
				continue;
			   
			dHist_NeutronShowerTheta_DeltaPhi[locRFBin]->Fill(locDeltaPhiCorrected, locNeutronShowerP4.Theta()*180./TMath::Pi());
			dHist_NeutronPreshowerE_DeltaPhi[locRFBin]->Fill(locDeltaPhiCorrected, locNeutronPreshowerE);
			
			// neutron shower selection cuts: true neutron yield is small beyond this angle 
			if(locNeutronShowerP4.Theta()*180./TMath::Pi() < dNeutronShowerThetaCut)
				continue;
			
			// some plots with DeltaPhi cut applied
			if(fabs(locDeltaPhiCorrected - 180.) < dDeltaPhiCut) {
				dHist_NeutronPreshowerE_ShowerTheta_DeltaPhiCut[locRFBin]->Fill(locNeutronShowerP4.Theta()*180./TMath::Pi(), locNeutronPreshowerE);
				dHist_NeutronBeta_P[locRFBin]->Fill(locMissingP4_Measured.P(), locNeutronBeta);
				dHist_NeutronDeltaBeta_P[locRFBin]->Fill(locMissingP4_Measured.P(), locNeutronDeltaBeta);
			}

			dHist_DeltaPhi_t[locRFBin]->Fill(loct, locDeltaPhiCorrected);

			// apply neutron timing cut
			if(locNeutronDeltaBeta < dNeutronDeltaBeta->Eval(locMissingP4_Measured.P()))
				continue;

			dHist_PiPlusMomentum_DeltaPhi[locRFBin]->Fill(locDeltaPhiCorrected, locPiPlusP4.E());
			dHist_DeltaPhi_t_DeltaBetaCut[locRFBin]->Fill(loct, locDeltaPhiCorrected);
			dHist_DeltaPhi_Phi_DeltaBetaCut[locRFBin]->Fill(locNeutronShowerP4.Phi()*180./TMath::Pi(), locDeltaPhi);
			if(unusedEnergyBCAL < 0.01) 
				dHist_DeltaPhi_t_noUnused[locRFBin]->Fill(loct, locDeltaPhiCorrected);

			/****************************************** FILL FLAT TREE ******************************************/
			
			//FILL ANY CUSTOM BRANCHES FIRST!!
			double weight = 1.;
			
			// weights for 2 deltaT sideband bins
			Bool_t Acci  = (fabs(locBeamDeltaT)<(3.5*4.008)&&fabs(locBeamDeltaT)>(0.5*4.008)) ? true : false;
			if(Acci) weight = -1/6.;
			
			// remove any further deltaT sidebands if they exist
			if(fabs(locBeamDeltaT) > 3.5*4.008)
				continue;
			
			dFlatTreeInterface->Fill_Fundamental<Double_t>("t", -1.*loct);
			dFlatTreeInterface->Fill_Fundamental<Double_t>("Egamma", locBeamP4.E());
			dFlatTreeInterface->Fill_Fundamental<Double_t>("PhiN", locNeutronShowerP4.Phi());
			dFlatTreeInterface->Fill_Fundamental<Double_t>("DeltaPhi", locDeltaPhiCorrected);
			dFlatTreeInterface->Fill_Fundamental<Double_t>("UnusedE", unusedEnergyBCAL);
			dFlatTreeInterface->Fill_Fundamental<Double_t>("Weight", weight);

			//FILL FLAT TREE
			Fill_FlatTree(); //for the active combo
			
			// signal deltaPhi region
			if(fabs(locDeltaPhiCorrected - 180.) < dDeltaPhiCut) {
				
				dHist_PiPlusMomentum_BeamE[locRFBin]->Fill(locBeamP4.E(), locPiPlusP4.E());
				dHist_NeutronShowerE_ShowerTheta[locRFBin]->Fill(locNeutronShowerP4.Theta()*180./TMath::Pi(), locNeutronShowerE);
				dHist_NeutronMissingTheta_ShowerTheta[locRFBin]->Fill(locNeutronShowerP4.Theta()*180./TMath::Pi(), locMissingP4_Measured.Theta()*180./TMath::Pi());
				dHist_NeutronPreshowerE_ShowerE[locRFBin]->Fill(locNeutronShowerE, locNeutronPreshowerE);
				dHist_MissingMass_Final[locRFBin]->Fill(locMissingMass);
				dHist_MissingMassSquared_Final[locRFBin]->Fill(locMissingMassSquared);

				dHist_NeutronPhi_t[locRFBin]->Fill(loct, locNeutronShowerP4.Phi()*180./TMath::Pi());

				// some additional asymmetries with tighter cuts
				if(locMissingMass < 1.4)
					dHist_NeutronPhi_t_tightMM[locRFBin]->Fill(loct, locNeutronShowerP4.Phi()*180./TMath::Pi());
				if(unusedEnergyBCAL < 0.01) {
					dHist_MissingMass_Final_noUnused[locRFBin]->Fill(locMissingMass);
					dHist_MissingMassSquared_Final_noUnused[locRFBin]->Fill(locMissingMassSquared);
					dHist_NeutronPhi_t_noUnused[locRFBin]->Fill(loct, locNeutronShowerP4.Phi()*180./TMath::Pi());
					if(locMissingMass < 1.4)
						dHist_NeutronPhi_t_noUnused_tightMM[locRFBin]->Fill(loct, locNeutronShowerP4.Phi()*180./TMath::Pi());
					if(locDeltaTheta > -20. || locDeltaTheta > 20.) 
						dHist_NeutronPhi_t_noUnused_tightDeltaTheta[locRFBin]->Fill(loct, locNeutronShowerP4.Phi()*180./TMath::Pi());
				}

				if(Get_NumThrown() > 0) { 
					double locDeltaPOverP = (locPiPlusP - locThrownPiPlusP)/locThrownPiPlusP;
					dHist_PiPlusDeltaPOverP_Theta[locRFBin]->Fill(locThrownPiPlusTheta*180./TMath::Pi(), locDeltaPOverP);
					dHist_Thrownt_Recot[locRFBin]->Fill(loct, locThrownt);
				}
			}
			else if(fabs(locDeltaPhi - 180.) > 20. && fabs(locDeltaPhi - 180.) < 40.) { // sideband deltaPhi region
				dHist_Sideband1NeutronPhi_t[locRFBin]->Fill(loct, locNeutronShowerP4.Phi()*180./TMath::Pi());	
			}
			else if(fabs(locDeltaPhi - 180.) > 40. && fabs(locDeltaPhi - 180.) < 60.) { // sideband deltaPhi region
				dHist_Sideband2NeutronPhi_t[locRFBin]->Fill(loct, locNeutronShowerP4.Phi()*180./TMath::Pi());	
			}
		}
		
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	return kTRUE;
}

void DSelector_npi_missing::Finalize(void)
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
