#include "DSelector_compton.h"

void DSelector_compton::Init(TTree *locTree)
{
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "compton.root"; //"" for none
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
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.4, Proton, SYS_TOF));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, Proton, SYS_BCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 2.0, Gamma, SYS_FCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 2.5, Gamma, SYS_BCAL));
	//dAnalysisActions.push_back(new DCutAction_NoPIDHit(dComboWrapper, Gamma, SYS_BCAL));

	//MASSES
	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 1.0e-10));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.1, 0.1));
	
	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();

	/******************************** USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	double locPhotonThetaMin = 0.; double locPhotonThetaMax = 60.;
	double locPhotonPMin = 0.; double locPhotonPMax = 12.;
	double locProtonThetaMin = 0.; double locProtonThetaMax = 180.;
	double locProtonPMin = 0.; double locProtonPMax = 6.;

	TString locRFBinLabel[2] = {"", "Acci_"};
	TString locCutCounterLabel[5] = {"Fiducial", "UE", "DeltaTheta", "DeltaPhiSignal", "DeltaPhiSideband"};
	TString locBeamEBinLabel[5] = {"Egamma5_6", "Egamma6_7", "Egamma7_8.4", "Egamma8.4_9", "Egamma9_12"};
	//int locDeltaPhiOffset = 3;
	
	for(int locRFBin=0; locRFBin<2; locRFBin++) {
		dHist_BeamEnergy[locRFBin] = new TH1I(locRFBinLabel[locRFBin]+"BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

		for(int locBeamEBin=0; locBeamEBin<5; locBeamEBin++) {
			dHist_PhotonP_Theta_Init[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"PhotonP_Theta_Init"+locBeamEBinLabel[locBeamEBin], "; #theta_{#gamma} (degrees); p_{#gamma} (GeV/c)", 360., locPhotonThetaMin, locPhotonThetaMax, 400, locPhotonPMin, locPhotonPMax);
			dHist_PhotonTheta_ThetaDet_Init[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"PhotonTheta_ThetaDet_Init"+locBeamEBinLabel[locBeamEBin], "; Detector #theta_{#gamma} (degrees); Measured #theta_{#gamma} (degrees)", 360., locPhotonThetaMin, locPhotonThetaMax, 360., locPhotonThetaMin, locPhotonThetaMax);
			
			dHist_ProtonP_Theta_Init[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"ProtonP_Theta_Init"+locBeamEBinLabel[locBeamEBin], "; #theta_{p} (degrees); p_{p} (GeV/c)", 360., locProtonThetaMin, locProtonThetaMax, 400, locProtonPMin, locProtonPMax);
			dHist_VertexZ[locRFBin][locBeamEBin] = new TH1I(locRFBinLabel[locRFBin]+"VertexZ"+locBeamEBinLabel[locBeamEBin], "; z_{vertex} (cm)", 300., 0., 300.);
			dHist_VertexR[locRFBin][locBeamEBin] = new TH1I(locRFBinLabel[locRFBin]+"VertexR"+locBeamEBinLabel[locBeamEBin], "; r_{vertex} (cm)", 200., 0., 2.);
			
			dHist_MissingMassSquared[locRFBin][locBeamEBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingMassSquared"+locBeamEBinLabel[locBeamEBin], ";Missing Mass Squared (GeV/c^{2})^{2}", 500, -0.1, 0.1);
			dHist_MissingEnergy[locRFBin][locBeamEBin] = new TH1I(locRFBinLabel[locRFBin]+"MissingEnergy"+locBeamEBinLabel[locBeamEBin], ";Missing Energy (GeV)", 400, -2.0, 2.0);
			
			dHist_UnusedEnergyBCAL_t[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"UnusedEnergyBCAL_t"+locBeamEBinLabel[locBeamEBin], "; -t (GeV^2); BCAL Unused Energy (GeV)", 250, 0., 5., 200, 0.0, 2.0);
			dHist_UnusedEnergyFCAL_t[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"UnusedEnergyFCAL_t"+locBeamEBinLabel[locBeamEBin], "; -t (GeV^2); FCAL Unused Energy (GeV)", 250, 0., 5., 200, 0.0, 2.0);
			dHist_UnusedEnergyTotal_t[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"UnusedEnergyTotal_t"+locBeamEBinLabel[locBeamEBin], "; -t (GeV^2); BCAL+FCAL Unused Energy (GeV)", 250, 0., 5., 200, 0.0, 2.0);
			
			dHist_PhotonThetaMeasure_PhotonThetaMissing[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"PhotonThetaMeasure_PhotonThetaMissing"+locBeamEBinLabel[locBeamEBin], "; #theta #gamma_{Missing} (degrees); #theta #gamma_{Measured} (degrees)", 240, 0., 120., 240, 0., 120.);
			dHist_DeltaE_DeltaTheta_BCAL[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaE_DeltaTheta_BCAL"+locBeamEBinLabel[locBeamEBin], "; #Delta#theta #gamma_{BCAL} (degrees); #Delta E (GeV)", 200, -20., 20., 500, -9.0, 1.0);
			dHist_DeltaPhi_DeltaTheta_BCAL[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_DeltaTheta_BCAL"+locBeamEBinLabel[locBeamEBin], "; #Delta#theta #gamma_{BCAL} (degrees); #Delta#phi (degrees)", 200, -10., 10., 200, 170., 190.);
			dHist_DeltaPhi_DeltaTheta_FCAL[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_DeltaTheta_FCAL"+locBeamEBinLabel[locBeamEBin], "; #Delta#theta #gamma_{FCAL} (degrees); #Delta#phi (degrees)", 200, -10., 10., 200, 170., 190.);
			
			dHist_DeltaPhi_t[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"DeltaPhi_t"+locBeamEBinLabel[locBeamEBin], "; -t (GeV^2); #Delta#phi (degrees)", 250, 0., 5., 200, 170., 190.);
			
			// final histograms
			dHist_PhotonP_Theta_Final[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"PhotonP_Theta_Final_"+locBeamEBinLabel[locBeamEBin], "; #theta_{#gamma} (degrees); p_{#gamma} (GeV/c)", 360, locPhotonThetaMin, locPhotonThetaMax, 400, locPhotonPMin, locPhotonPMax);
			dHist_PhotonTheta_ThetaDet_Final[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"PhotonTheta_ThetaDet_Final"+locBeamEBinLabel[locBeamEBin], "; Detector #theta_{#gamma} (degrees); Measured #theta_{#gamma} (degrees)", 360., locPhotonThetaMin, locPhotonThetaMax, 360., locPhotonThetaMin, locPhotonThetaMax);
			dHist_ProtonP_Theta_Final[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"ProtonP_Theta_Final_"+locBeamEBinLabel[locBeamEBin], "; #theta_{p} (degrees); p_{p} (GeV/c)", 360, locProtonThetaMin, locProtonThetaMax, 400, locProtonPMin, locProtonPMax);
			
			dHist_ThetaCM_PhotonTheta[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"ThetaCM_PhotonTheta_"+locBeamEBinLabel[locBeamEBin], "; #theta_{#gamma} (degrees); #theta_{CM} (degrees)", 360, locPhotonThetaMin, locPhotonThetaMax, 360, 0., 180.);
			dHist_ThetaCM_ProtonTheta[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"ThetaCM_ProtonTheta_"+locBeamEBinLabel[locBeamEBin], "; #theta_{p} (degrees); #theta_{CM} (degrees)", 360, locProtonThetaMin, locProtonThetaMax, 360, 0., 180.);
			dHist_ThetaCM_t[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"ThetaCM_t_"+locBeamEBinLabel[locBeamEBin], "; -t (GeV^{2}); #theta_{CM} (degrees)", 250, 0., 5.0, 360, 0., 180.);
			dHist_u_t[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"u_t_"+locBeamEBinLabel[locBeamEBin], "; -t (GeV^{2}); -u (GeV^{2})", 250, 0., 5.0, 250, 0., 25.0);

			dHist_BCALSigTrans_SigLong[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"SigTrans_SigLong_"+locBeamEBinLabel[locBeamEBin], "; BCAL Shower #sigma_{#rho}; BCAL Shower #sigma_{#phi}", 300, 0, 10, 200, 0., 0.02);
			dHist_BCALSigTrans_SigTheta[locRFBin][locBeamEBin] = new TH2I(locRFBinLabel[locRFBin]+"SigTrans_SigTheta_"+locBeamEBinLabel[locBeamEBin], "; BCAL Shower #sigma_{#theta}; BCAL Shower #sigma_{#phi}", 200, 0., 0.05, 200, 0., 0.02);
		}			

		for(int locCutCounter=0; locCutCounter<5; locCutCounter++) {

			dHist_ProtonPhi_t[locRFBin][locCutCounter] = new TH2I(locRFBinLabel[locRFBin]+"ProtonPhi_t_"+locCutCounterLabel[locCutCounter], "; -t (GeV^2); #phi_{p} (degrees)", 250, 0., 5., 360, -180., 180.);
			dHist_ProtonPhi_ThetaCM[locRFBin][locCutCounter]= new TH2I(locRFBinLabel[locRFBin]+"ProtonPhi_ThetaCM_"+locCutCounterLabel[locCutCounter], "; #theta_{CM} (degrees); #phi_{p} (degrees)", 360, 0., 180., 360, -180., 180.);

			/*
			if(locCutCounter>=locDeltaPhiOffset) {
				dHist_PhotonP_Theta_Final[locRFBin][locCutCounter-locDeltaPhiOffset] = new TH2I(locRFBinLabel[locRFBin]+"PhotonP_Theta_Final_"+locCutCounterLabel[locCutCounter], "; #theta_{#gamma} (degrees); p_{#gamma} (GeV/c)", 360, locPhotonThetaMin, locPhotonThetaMax, 400, locPhotonPMin, locPhotonPMax);
				dHist_PhotonTheta_ThetaDet_Final[locRFBin][locCutCounter-locDeltaPhiOffset] = new TH2I(locRFBinLabel[locRFBin]+"PhotonTheta_ThetaDet_Final"+locCutCounterLabel[locCutCounter], "; Detector #theta_{#gamma} (degrees); Measured #theta_{#gamma} (degrees)", 360., locPhotonThetaMin, locPhotonThetaMax, 360., locPhotonThetaMin, locPhotonThetaMax);
				dHist_ProtonP_Theta_Final[locRFBin][locCutCounter-locDeltaPhiOffset] = new TH2I(locRFBinLabel[locRFBin]+"ProtonP_Theta_Final_"+locCutCounterLabel[locCutCounter], "; #theta_{p} (degrees); p_{p} (GeV/c)", 360, locProtonThetaMin, locProtonThetaMax, 400, locProtonPMin, locProtonPMax);
				dHist_ProtonP_Theta_Final_Test[locRFBin][locCutCounter-locDeltaPhiOffset] = new TH2I(locRFBinLabel[locRFBin]+"ProtonP_Theta_Final_Test_"+locCutCounterLabel[locCutCounter], "; #theta_{p} (degrees); p_{p} (GeV/c)", 360, locProtonThetaMin, locProtonThetaMax, 400, locProtonPMin, locProtonPMax);
				
				dHist_ThetaCM_PhotonTheta[locRFBin][locCutCounter-locDeltaPhiOffset] = new TH2I(locRFBinLabel[locRFBin]+"ThetaCM_PhotonTheta_"+locCutCounterLabel[locCutCounter], "; #theta_{#gamma} (degrees); #theta_{CM} (degrees)", 360, locPhotonThetaMin, locPhotonThetaMax, 360, 0., 180.);
				dHist_ThetaCM_ProtonTheta[locRFBin][locCutCounter-locDeltaPhiOffset] = new TH2I(locRFBinLabel[locRFBin]+"ThetaCM_ProtonTheta_"+locCutCounterLabel[locCutCounter], "; #theta_{p} (degrees); #theta_{CM} (degrees)", 360, locProtonThetaMin, locProtonThetaMax, 360, 0., 180.);
				dHist_ThetaCM_t[locRFBin][locCutCounter-locDeltaPhiOffset] = new TH2I(locRFBinLabel[locRFBin]+"ThetaCM_t_"+locCutCounterLabel[locCutCounter], "; -t (GeV^{2}); #theta_{CM} (degrees)", 250, 0., 5.0, 360, 0., 180.);
			}
			*/

		}
	}

}

Bool_t DSelector_compton::Process(Long64_t locEntry)
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
		Int_t locPhotonNeutralID = dPhotonWrapper->Get_NeutralID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPhotonP4 = dPhotonWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPhotonP4_Measured = dPhotonWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		TLorentzVector locPhotonX4_Measured = dPhotonWrapper->Get_X4_Shower();
		TLorentzVector locProtonX4_Measured = dProtonWrapper->Get_X4_Measured();

		// Get Shower Energy:
		double locPhotonBCALPreShowerE = dPhotonWrapper->Get_Energy_BCALPreshower();
		double locPhotonBCALShowerE = dPhotonWrapper->Get_Energy_BCAL();
		double locPhotonFCALShowerE = dPhotonWrapper->Get_Energy_FCAL();
		double locPhotonSigLongBCALShower = dPhotonWrapper->Get_SigLong_BCAL();
		double locPhotonSigThetaBCALShower = dPhotonWrapper->Get_SigTheta_BCAL();
		double locPhotonSigTransBCALShower = dPhotonWrapper->Get_SigTrans_BCAL();
		
		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/********************************************* ACCIDENTALS & RF TIMING ********************************************/

		double locRFTime = dComboWrapper->Get_RFTime_Measured();
		double locBeamDeltaT = locBeamX4_Measured.T() - (locRFTime + (locProtonX4_Measured.Z() - dTargetCenter.Z())/29.9792458);

		int locRFBin = -1;
		if(fabs(locBeamDeltaT) < 0.5*4.008) locRFBin = 0;
		else if (fabs(locBeamDeltaT) > 4.008*0.5 || fabs(locBeamDeltaT) < 4.008*2.5) locRFBin = 1;
		else continue;

		/********************************************* KINEMATIC VARIABLES **********************************************/

		// Get COM P4's for scattering angle
		TLorentzVector locInitialStateP4 = locBeamP4_Measured + dTargetP4;
		TVector3 locBoostVectorCM = -1.0*(locInitialStateP4.BoostVector());
		TLorentzVector locProtonP4_MeasuredCM = locProtonP4_Measured; locProtonP4_MeasuredCM.Boost(locBoostVectorCM);
		TLorentzVector locPhotonP4_MeasuredCM = locPhotonP4_Measured; locPhotonP4_MeasuredCM.Boost(locBoostVectorCM);
		
		// should directly be photon CM angle, but can't use that since BCAL saturates
		double locThetaCM = -1.*(locProtonP4_MeasuredCM.Theta()*180./TMath::Pi() - 180.); 
		double loct = -1.*(dTargetP4 - locProtonP4_Measured).M2();
		double locu = -1.*(locBeamP4_Measured - locProtonP4_Measured).M2();

		/**************************************** HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy[locRFBin]->Fill(locBeamP4.E());
			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		// cut low energy beam photons and make bins for Egamma (in arrays below)
		if(locBeamP4.E() < 5.0)
			continue;

		int locBeamEBin = 0;
		if(locBeamP4.E() > 6.)  locBeamEBin = 1;
		if(locBeamP4.E() > 7.)  locBeamEBin = 2;
		if(locBeamP4.E() > 8.4) locBeamEBin = 3;
		if(locBeamP4.E() > 9.)  locBeamEBin = 4;
		
		/********************************************* FIDUCIAL AND VERTEX CUTS ********************************************/
		
		double locPhotonP = locPhotonP4_Measured.P();
		double locPhotonTheta = locPhotonP4_Measured.Theta()*180./TMath::Pi();
		double locPhotonThetaDet = (locPhotonX4_Measured.Vect()-dTargetCenter).Theta()*180./TMath::Pi();
		dHist_PhotonP_Theta_Init[locRFBin][locBeamEBin]->Fill(locPhotonThetaDet, locPhotonP);
		if(locPhotonThetaDet < 2.5 || (locPhotonThetaDet > 10.3 && locPhotonThetaDet < 11.5)) continue;
		if(locPhotonBCALShowerE > 0.) {
			dHist_PhotonTheta_ThetaDet_Init[locRFBin][locBeamEBin]->Fill(locPhotonThetaDet, locPhotonTheta);
			if(locPhotonTheta < 10.5)
				continue; // remove BCAL showers with wrong measured theta 
		}

		double locProtonP = locProtonP4_Measured.P();
		double locProtonTheta = locProtonP4_Measured.Theta()*180./TMath::Pi();
		double locProtonPhi = locProtonP4_Measured.Phi()*180./TMath::Pi();
		dHist_ProtonP_Theta_Init[locRFBin][locBeamEBin]->Fill(locProtonTheta, locProtonP);
		if(locProtonP < 0.25) continue;
		if(locProtonTheta < 20.) continue; // remove high -t region populated by accidentals

		double locVertexZ = locProtonX4_Measured.Z();
                dHist_VertexZ[locRFBin][locBeamEBin]->Fill(locVertexZ);
                if(locVertexZ < 51. || locVertexZ > 76.)
                        continue;

		double locVertexR = locProtonX4_Measured.Perp();
                dHist_VertexR[locRFBin][locBeamEBin]->Fill(locVertexR);
                if(locVertexR > 1.)
                        continue;

		int locCutCounter = 0;
		if(locBeamEBin == 3) {
			dHist_ProtonPhi_t[locRFBin][locCutCounter]->Fill(loct, locProtonPhi);
			dHist_ProtonPhi_ThetaCM[locRFBin][locCutCounter]->Fill(locThetaCM, locProtonPhi);
		}

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPhotonP4_Measured + locProtonP4_Measured;
		TLorentzVector locPhotonP4_Missing_Measured = locBeamP4_Measured + dTargetP4;
		locPhotonP4_Missing_Measured -= locProtonP4_Measured;

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();
		double locMissingEnergy = locMissingP4_Measured.E();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[Gamma].insert(locPhotonNeutralID);
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared[locRFBin][locBeamEBin]->Fill(locMissingMassSquared);
			dHist_MissingEnergy[locRFBin][locBeamEBin]->Fill(locMissingEnergy);
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}
		
		// cut missing mass squared
		//if(fabs(locMissingMassSquared) > 0.05)
		//	continue;

		// unused shower energy sum
		double unusedEnergyBCAL = 0.;
		double unusedEnergyFCAL = 0.;
		for(UInt_t loc_j = 0; loc_j < Get_NumNeutralHypos(); ++loc_j) {
			
			//Set branch array indices corresponding to this particle
			dNeutralHypoWrapper->Set_ArrayIndex(loc_j);
			UInt_t locPhotonNeutralPID = dNeutralHypoWrapper->Get_PID();
			if(locPhotonNeutralPID != Gamma || dNeutralHypoWrapper->Get_NeutralID() == locPhotonNeutralID)
				continue;
			
			unusedEnergyBCAL += dNeutralHypoWrapper->Get_Energy_BCAL();
			unusedEnergyFCAL += dNeutralHypoWrapper->Get_Energy_FCAL();
		}
		dHist_UnusedEnergyBCAL_t[locRFBin][locBeamEBin]->Fill(loct, unusedEnergyBCAL);
		dHist_UnusedEnergyFCAL_t[locRFBin][locBeamEBin]->Fill(loct, unusedEnergyFCAL);
		dHist_UnusedEnergyTotal_t[locRFBin][locBeamEBin]->Fill(loct, unusedEnergyBCAL+unusedEnergyFCAL);

		// cut unused energy
		if(unusedEnergyBCAL+unusedEnergyFCAL > 0.1)
			continue;

		// cut BCAL shower shape
		if(locPhotonBCALShowerE > 0.) {
			dHist_BCALSigTrans_SigLong[locRFBin][locBeamEBin]->Fill(locPhotonSigLongBCALShower, locPhotonSigTransBCALShower);
			dHist_BCALSigTrans_SigTheta[locRFBin][locBeamEBin]->Fill(locPhotonSigThetaBCALShower, locPhotonSigTransBCALShower);
			if(locPhotonSigTransBCALShower > 0.0065 || locPhotonSigThetaBCALShower > 0.0025)
				continue;
		}

		locCutCounter++; // = 1
		if(locBeamEBin == 3) {
			dHist_ProtonPhi_t[locRFBin][locCutCounter]->Fill(loct, locProtonPhi);
			dHist_ProtonPhi_ThetaCM[locRFBin][locCutCounter]->Fill(locThetaCM, locProtonPhi);
		}

		// DeltaTheta
		double locDeltaTheta = (locPhotonP4_Measured.Theta() - locPhotonP4_Missing_Measured.Theta())*180./TMath::Pi();
		// DeltaE
		double locDeltaE = locPhotonBCALShowerE - locPhotonP4_Missing_Measured.E();
		if(locPhotonFCALShowerE > 0.) 
			locDeltaE = locPhotonFCALShowerE - locPhotonP4_Missing_Measured.E();

		// DeltaPhi
		double locDeltaPhi = (locPhotonP4_Measured.Phi() - locProtonP4_Measured.Phi())*180./TMath::Pi();
		if(locDeltaPhi > 360.) locDeltaPhi -= 360.;
		if(locDeltaPhi < 0.) locDeltaPhi += 360.;
		
		if(locPhotonBCALShowerE > 0.) {
			dHist_PhotonThetaMeasure_PhotonThetaMissing[locRFBin][locBeamEBin]->Fill(locPhotonP4_Missing_Measured.Theta()*180./TMath::Pi(), locPhotonP4_Measured.Theta()*180./TMath::Pi());
			dHist_DeltaPhi_DeltaTheta_BCAL[locRFBin][locBeamEBin]->Fill(locDeltaTheta, locDeltaPhi);
			if(fabs(locDeltaPhi - 180) < 1.)
				dHist_DeltaE_DeltaTheta_BCAL[locRFBin][locBeamEBin]->Fill(locDeltaTheta, locDeltaE);
		}
		else
			dHist_DeltaPhi_DeltaTheta_FCAL[locRFBin][locBeamEBin]->Fill(locDeltaTheta, locDeltaPhi);

		locCutCounter++; // = 2
		if(locBeamEBin == 3) {
			dHist_ProtonPhi_t[locRFBin][locCutCounter]->Fill(loct, locProtonPhi);
			dHist_ProtonPhi_ThetaCM[locRFBin][locCutCounter]->Fill(locThetaCM, locProtonPhi);
		}
		
		if(fabs(locDeltaTheta) > 4. || locDeltaE < -2.) 
			continue;

		dHist_DeltaPhi_t[locRFBin][locBeamEBin]->Fill(loct, locDeltaPhi);

		// assign delta phi bins to to match locCutCounter
		if(fabs(locDeltaPhi - 180.) < 0.5) locCutCounter += 1; // = 3
		else if(fabs(locDeltaPhi - 180.) > 10. && fabs(locDeltaPhi - 180.) < 20) locCutCounter += 2; // = 4
		else continue;

		/********************************************* FINAL ASYMMETRIES ********************************************/

		if(locBeamEBin == 3) {
			dHist_ProtonPhi_t[locRFBin][locCutCounter]->Fill(loct, locProtonPhi);
			dHist_ProtonPhi_ThetaCM[locRFBin][locCutCounter]->Fill(locThetaCM, locProtonPhi);
		}

		/********************************************* FINAL KINEMATICS ********************************************/
		
		if(locCutCounter == 3) {
			dHist_ProtonP_Theta_Final[locRFBin][locBeamEBin]->Fill(locProtonTheta, locProtonP);
			dHist_PhotonP_Theta_Final[locRFBin][locBeamEBin]->Fill(locPhotonThetaDet, locPhotonP);
			if(locPhotonBCALShowerE > 0.) {
				dHist_PhotonTheta_ThetaDet_Final[locRFBin][locBeamEBin]->Fill(locPhotonThetaDet, locPhotonTheta);
			}
			dHist_ThetaCM_ProtonTheta[locRFBin][locBeamEBin]->Fill(locProtonTheta, locThetaCM);
			dHist_ThetaCM_PhotonTheta[locRFBin][locBeamEBin]->Fill(locPhotonTheta, locThetaCM);
			dHist_ThetaCM_t[locRFBin][locBeamEBin]->Fill(loct, locThetaCM);
			dHist_u_t[locRFBin][locBeamEBin]->Fill(loct, locu);
		}

	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	return kTRUE;
}

void DSelector_compton::Finalize(void)
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
