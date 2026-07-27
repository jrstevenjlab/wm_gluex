#ifndef DSelector_omega_misspip_h
#define DSelector_omega_misspip_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"

class DSelector_omega_misspip : public DSelector
{
	public:

		DSelector_omega_misspip(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_omega_misspip(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;
		DKinematicData* dMissingPiPlusWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		//DKinematicData* dDecayingPi0Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1F* dHist_MissingMassSquared;
		TH1F* dHist_MissingMass;
		TH1I* dHist_BeamEnergy;
		TH1F* dHist_MissingEnergy;

		TH1F* dHist_BeamDeltaT;
		TH1F* dHist_UnusedShowerEnergy;
		TH1F* dHist_NumTracks;
		TH1F* dHist_KinFitCL;

		TH1F* dHist_PValue;
		TH2F* dHist_PValueVsTheta;
		TH2F* dHist_PValueVsP;
		TH1F* dHist_PValue_sideband;
		TH2F* dHist_PValueVsTheta_sideband;
		TH2F* dHist_PValueVsP_sideband;
		TH1F* dHist_logPValue;
		TH2F* dHist_logPValueVsTheta;
		TH2F* dHist_logPValueVsP;
		TH1F* dHist_logPValue_sideband;
		TH2F* dHist_logPValueVsTheta_sideband;
		TH2F* dHist_logPValueVsP_sideband;

		TH2F* dHist_PValueVsTheta_reco;
		TH2F* dHist_PValueVsP_reco;
		TH2F* dHist_PValueVsTheta_sideband_reco;
		TH2F* dHist_PValueVsP_sideband_reco;
		TH2F* dHist_logPValueVsTheta_reco;
		TH2F* dHist_logPValueVsP_reco;
		TH2F* dHist_logPValueVsTheta_sideband_reco;
		TH2F* dHist_logPValueVsP_sideband_reco;

		TH1F* dHist_Chi2;
		TH2F* dHist_Chi2vsTheta;
		TH2F* dHist_Chi2vsP;
		TH1F* dHist_Chi2_sideband;
		TH2F* dHist_Chi2vsTheta_sideband;
		TH2F* dHist_Chi2vsP_sideband;

		TH1F* dHist_DeltaPRatio;
		TH1F* dHist_DeltaPRatio_truth;
		TH2F* dHist_DeltaPVsTheta;
		TH2F* dHist_DeltaPVsTheta_truth;

		TH1F* dHist_3PiMass_missing;
		TH1F* dHist_3PiMass_reco[15];
		TH1F* dHist_2GammaMass;
		TH1F* dHist_OmegaMass_missing;
		TH1F* dHist_OmegaMass_reco[15];
		TH2F* dHist_MassCorrelation[15];

		TH1F* dHist_3PiMass_pvcut[7];
		TH1F* dHist_OmegaMass_pvcut[7];
		TH2F* dHist_MassCorr_pvcut[7];
		TH2F* dHist_Theta_recoVstruth_pvcut[7];

		TH1F* dHist_3PiMass_deltapcut[7];
		TH1F* dHist_OmegaMass_deltapcut[7];
		TH2F* dHist_MassCorr_deltapcut[7];
		TH2F* dHist_Theta_recoVstruth_deltapcut[7];

		TH1F* dHist_DeltaPhiMissing;
		TH1F* dHist_DeltaThetaMissing;
		TH1F* dHist_DeltaPMissing;
		TH2F* dHist_DeltaPhiVsDeltaTheta_missing;
		TH2F* dHist_DeltaPVsDeltaPhi_missing;
		TH2F* dHist_DeltaPVsDeltaTheta_missing;

		TH1F* dHist_Theta_reco;
		TH1F* dHist_Theta_truth;
		TH2F* dHist_Theta_recoVstruth;

		TH1F* dHist_DeltaPhiTruth;
		TH1F* dHist_DeltaThetaTruth;
		TH1F* dHist_DeltaPTruth;

		TH2F* dHist_OmegaMassVsPhi_missing;
		TH2F* dHist_OmegaMassVsTheta_missing;
		TH2F* dHist_OmegaMassVsP_missing;

		TH2F* dHist_OmegaMassVsPhi_missing_1or0;
		TH2F* dHist_OmegaMassVsTheta_missing_1or0;
		TH2F* dHist_OmegaMassVsP_missing_1or0;

		TH2F* dHist_OmegaMassVsPhi_reco[7];
		TH2F* dHist_OmegaMassVsTheta_reco[7];
		TH2F* dHist_OmegaMassVsP_reco[7];
		TH2F* dHist_3PiMassVsPhi_reco[7];
		TH2F* dHist_3PiMassVsTheta_reco[7];
		TH2F* dHist_3PiMassVsP_reco[7];

		TH2F* dHist_OmegaMassVsPhi_pvcut[7];
		TH2F* dHist_OmegaMassVsTheta_pvcut[7];
		TH2F* dHist_OmegaMassVsP_pvcut[7];
		TH2F* dHist_3PiMassVsPhi_pvcut[7];
		TH2F* dHist_3PiMassVsTheta_pvcut[7];
		TH2F* dHist_3PiMassVsP_pvcut[7];

		TH2F* dHist_OmegaMassVsTheta_missing_p[9];
		TH2F* dHist_OmegaMassVsTheta_reco_p[9][7];
		TH2F* dHist_3PiMassVsTheta_reco_p[9][7];

		TH2F* dHist_OmegaMassVsTheta_pvcut_p[9][7];
		TH2F* dHist_3PiMassVsTheta_pvcut_p[9][7];

		TH2F* dHist_OmegaMassVsPhi_deltapcut[7];
		TH2F* dHist_OmegaMassVsTheta_deltapcut[7];
		TH2F* dHist_OmegaMassVsP_deltapcut[7];
		TH2F* dHist_3PiMassVsPhi_deltapcut[7];
		TH2F* dHist_3PiMassVsTheta_deltapcut[7];
		TH2F* dHist_3PiMassVsP_deltapcut[7];

		TH2F* dHist_OmegaMassVsTheta_deltapcut_p[9][7];
		TH2F* dHist_3PiMassVsTheta_deltapcut_p[9][7];

		TH3F* dHist_DeltaPhi3D_truth;
		TH3F* dHist_DeltaTheta_truth;
		TH3F* dHist_DeltaP_truth;

		TH3F* dHist_DeltaPhi_miss;
		TH3F* dHist_DeltaTheta_miss;
		TH3F* dHist_DeltaP_miss;

		TH2F* dHist_ThetaVsP_reco;
		TH2F* dHist_ThetaVsP_missing;

		TH1F* dHist_3PiMass_thetacut[7];
		TH1F* dHist_OmegaMass_thetacut[7];
		TH2F* dHist_MassCorr_thetacut[7];
		TH2F* dHist_Theta_recoVstruth_thetacut[7];
		TH2F* dHist_OmegaMassVsPhi_thetacut[7];
		TH2F* dHist_OmegaMassVsTheta_thetacut[7];
		TH2F* dHist_OmegaMassVsP_thetacut[7];
		TH2F* dHist_3PiMassVsPhi_thetacut[7];
		TH2F* dHist_3PiMassVsTheta_thetacut[7];
		TH2F* dHist_3PiMassVsP_thetacut[7];
		TH2F* dHist_OmegaMassVsTheta_thetacut_p[9][7];
		TH2F* dHist_3PiMassVsTheta_thetacut_p[9][7];

		TH1F* dHist_3PiMass_pcut[7];
		TH1F* dHist_OmegaMass_pcut[7];
		TH2F* dHist_MassCorr_pcut[7];
		TH2F* dHist_Theta_recoVstruth_pcut[7];
		TH2F* dHist_OmegaMassVsPhi_pcut[7];
		TH2F* dHist_OmegaMassVsTheta_pcut[7];
		TH2F* dHist_OmegaMassVsP_pcut[7];
		TH2F* dHist_3PiMassVsPhi_pcut[7];
		TH2F* dHist_3PiMassVsTheta_pcut[7];
		TH2F* dHist_3PiMassVsP_pcut[7];
		TH2F* dHist_OmegaMassVsTheta_pcut_p[9][7];
		TH2F* dHist_3PiMassVsTheta_pcut_p[9][7];




	ClassDef(DSelector_omega_misspip, 0);
};

void DSelector_omega_misspip::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dMissingPiPlusWrapper = dStep0Wrapper->Get_FinalParticle(3);

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	//dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_omega_misspip_h
