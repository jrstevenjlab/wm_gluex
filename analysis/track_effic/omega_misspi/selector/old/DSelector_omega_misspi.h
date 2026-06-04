#ifndef DSelector_omega_misspi_h
#define DSelector_omega_misspi_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"

class DSelector_omega_misspi : public DSelector
{
	public:

		DSelector_omega_misspi(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_omega_misspi(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPiRecoWrapper;
		DChargedTrackHypothesis* dProtonWrapper;
		DKinematicData* dMissingPiMissWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		//DKinematicData* dDecayingPi0Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1F* dHist_MissingMassSquared;
		TH1F* dHist_MissingMass;
		TH1F* dHist_BeamEnergy;
		TH1F* dHist_MissingEnergy;

		TH1F* dHist_BeamDeltaT;
		TH1F* dHist_UnusedShowerEnergy;
		TH1F* dHist_NumTracks;
		TH1F* dHist_KinFitCL;

		TH1F* dHist_3PiMass[10];
		TH1F* dHist_2GammaMass;
		TH1F* dHist_OmegaMass_missing;
		TH1F* dHist_OmegaMass_1track[10];
		TH2F* dHist_MassCorrelation[10];

		TH1F* dHist_DeltaPhiMissing;
		TH1F* dHist_DeltaThetaMissing;
		TH1F* dHist_DeltaPMissing;
		TH2F* dHist_DeltaPhiVsDeltaTheta_missing;
		TH2F* dHist_DeltaPVsDeltaPhi_missing;
		TH2F* dHist_DeltaPVsDeltaTheta_missing;

		TH1F* dHist_DeltaPhiTruth;
		TH1F* dHist_DeltaThetaTruth;
		TH1F* dHist_DeltaPTruth;

		TH2F* dHist_OmegaMassVsPhi_0track;
		TH2F* dHist_OmegaMassVsTheta_0track;
		TH2F* dHist_OmegaMassVsP_0track;

		TH2F* dHist_OmegaMassVsPhi_1track;
		TH2F* dHist_OmegaMassVsTheta_1track;
		TH2F* dHist_OmegaMassVsP_1track;
		TH2F* dHist_3PiMassVsPhi;
		TH2F* dHist_3PiMassVsTheta;
		TH2F* dHist_3PiMassVsP;

		TH2F* dHist_OmegaMassVsTheta_0track_p[9];
		TH2F* dHist_OmegaMassVsTheta_1track_p[9];
		TH2F* dHist_3PiMassVsTheta_p[9];

		TH3F* dHist_DeltaPhi3D_truth;
		TH3F* dHist_DeltaTheta3D_truth;
		TH3F* dHist_DeltaP3D_truth;

		TH3F* dHist_DeltaPhi3D_miss;
		TH3F* dHist_DeltaTheta3D_miss;
		TH3F* dHist_DeltaP3D_miss;

		TH2F* dHist_ThetaVsP_reco;
		TH2F* dHist_ThetaVsP_missing;

		TH1F* dHist_PValue;
		TH2F* dHist_PValueVsTheta;

		TH1F* dHist_3PiMass_pvcut[7];
		TH1F* dHist_OmegaMass_1track_pvcut[7];
		TH2F* dHist_MassCorr_pvcut[7];

		TH2F* dHist_OmegaMassVsPhi_1track_pvcut[7];
		TH2F* dHist_OmegaMassVsTheta_1track_pvcut[7];
		TH2F* dHist_OmegaMassVsP_1track_pvcut[7];
		TH2F* dHist_3PiMassVsPhi_pvcut[7];
		TH2F* dHist_3PiMassVsTheta_pvcut[7];
		TH2F* dHist_3PiMassVsP_pvcut[7];

		TH2F* dHist_Theta_recoVstruth_pvcut[7];

		TH2F* dHist_OmegaMassVsTheta_1track_pvcut_p[9][7];
		TH2F* dHist_3PiMassVsTheta_pvcut_p[9][7];

		TH1F* dHist_3PiMass_thetacut[7];
		TH1F* dHist_OmegaMass_1track_thetacut[7];
		TH2F* dHist_MassCorr_thetacut[7];

		TH2F* dHist_OmegaMassVsPhi_1track_thetacut[7];
		TH2F* dHist_OmegaMassVsTheta_1track_thetacut[7];
		TH2F* dHist_OmegaMassVsP_1track_thetacut[7];
		TH2F* dHist_3PiMassVsPhi_thetacut[7];
		TH2F* dHist_3PiMassVsTheta_thetacut[7];
		TH2F* dHist_3PiMassVsP_thetacut[7];

		TH2F* dHist_Theta_recoVstruth_thetacut[7];

		TH2F* dHist_OmegaMassVsTheta_1track_thetacut_p[9][7];
		TH2F* dHist_3PiMassVsTheta_thetacut_p[9][7];

		TH1F* dHist_3PiMass_pcut[7];
		TH1F* dHist_OmegaMass_1track_pcut[7];
		TH2F* dHist_MassCorr_pcut[7];

		TH2F* dHist_OmegaMassVsPhi_1track_pcut[7];
		TH2F* dHist_OmegaMassVsTheta_1track_pcut[7];
		TH2F* dHist_OmegaMassVsP_1track_pcut[7];
		TH2F* dHist_3PiMassVsPhi_pcut[7];
		TH2F* dHist_3PiMassVsTheta_pcut[7];
		TH2F* dHist_3PiMassVsP_pcut[7];

		TH2F* dHist_Theta_recoVstruth_pcut[7];

		TH2F* dHist_OmegaMassVsTheta_1track_pcut_p[9][7];
		TH2F* dHist_3PiMassVsTheta_pcut_p[9][7];


	ClassDef(DSelector_omega_misspi, 3);
};

void DSelector_omega_misspi::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiRecoWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dMissingPiMissWrapper = dStep0Wrapper->Get_FinalParticle(3);

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	//dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_omega_misspi_h
