#ifndef DSelector_pomega2pi_omega3pi_h
#define DSelector_pomega2pi_omega3pi_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_pomega2pi_omega3pi : public DSelector
{
	public:

		DSelector_pomega2pi_omega3pi(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_pomega2pi_omega3pi(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		bool dIsMC;

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPiPlus1Wrapper;
		DChargedTrackHypothesis* dPiPlus2Wrapper;
		DChargedTrackHypothesis* dPiMinus1Wrapper;
		DChargedTrackHypothesis* dPiMinus2Wrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DKinematicData* dDecayingPi0Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1F* dHist_MissingMassSquared;
		TH1F* dHist_BeamEnergy;
		TH1F* dHist_BeamDeltaT;

		TH1F* dHist_ThrownTopologies;
		map<TString, TH1F*> dHist_InvariantMass_ThrownTopology;

		TH1F* dHist_ThrownTopologies_omegacut;
		map<TString, TH1F*> dHist_InvariantMass_ThrownTopology_omegacut;

		TH1F* dHist_ThrownTopologies_omegacut_nosideband;
		map<TString, TH1F*> dHist_InvariantMass_ThrownTopology_omegacut_nosideband;

		TH1F* dHist_2PhotonMass;
		TH2F* dHist_2PhotonVs3PiMass;
		TH2F* dHist_2PhotonVs2PiMass;
		TH2F* dHist_2PhotonVsOmegaPiPlusMass;
		TH2F* dHist_2PhotonVsOmegaPiMinusMass;
		TH2F* dHist_2PhotonVsOmega2PiMass;

		TH1F* dHist_5PiMass;
		TH1F* dHist_5PiMass_Measured;
		TH1F* dHist_3PiMass;
		TH1F* dHist_3PiMass_Measured;
		TH1F* dHist_Omega2PiMass;
		TH2F* dHist_3vs2PiMass;
		TH1F* dHist_OmegaPiPlusMass;
		TH1F* dHist_OmegaPiMinusMass;
		TH2F* dHist_OmegaPiPlusVsOmegaPiMinus;
		TH1F* dHist_ProtonPiPlusMass;
		TH2F* dHist_ProtonPiPlusVsOmegaPiMinus;

		TH1F* dHist_ProtonPiMinusMass;
		TH2F* dHist_ProtonPiMinusVsOmegaPiPlus;
		TH2F* dHist_3vs4PiMass_plus;
		TH2F* dHist_3vs4PiMass_minus;

		TH2F* dHist_OmegaPiPlusVsOmega2PiMass;
		TH2F* dHist_OmegaPiMinusVsOmega2PiMass;

		TH1F* dHist_t_proton_total;
		TH1F* dHist_t_omega_total;
		TH1F* dHist_t_Deltaplusplus_total;
		TH1F* dHist_t_Delta0_total; 
		TH1F* dHist_t_proton[10];
		TH1F* dHist_t_omega[10];
		TH1F* dHist_t_Deltaplusplus[10];
		TH1F* dHist_t_Delta0[10]; 
		TH2F* dHist_t_pomega_correlation;
		TH2F* dHist_t_Delta_correlation_total; //new
		TH2F* dHist_t_Delta_correlation[10]; //new

		TH1F* dHist_2PiMass;
		TH1F* dHist_2PiMass_plus;
		TH1F* dHist_2PiMass_minus;
		TH1F* dHist_2PiMass_plus_omegacut;
		TH1F* dHist_2PiMass_minus_omegacut;
		TH1F* dHist_2PiMass_plus_notomega;
		TH1F* dHist_2PiMass_minus_notomega;
		TH2F* dHist_3vs2PiMass_plus;
		TH2F* dHist_3vs2PiMass_minus;

		TH2F* dHist_2PiMass_t_proton_total;
		TH2F* dHist_2PiMass_t_omega_total;
		TH2F* dHist_2PiMass_t_Deltaplusplus_total;
		TH2F* dHist_2PiMass_t_Delta0_total; 
		TH2F* dHist_2PiMass_t_proton[10];
		TH2F* dHist_2PiMass_t_omega[10];
		TH2F* dHist_2PiMass_t_Deltaplusplus[10];
		TH2F* dHist_2PiMass_t_Delta0[10]; 

		TH2F* dHist_omega2piMass_t_omega;
		TH2F* dHist_omegapiplusMass_t_omega;
		TH2F* dHist_omegapiminusMass_t_omega;
		TH2F* dHist_omega2piMass_t_Delta;
		TH2F* dHist_omegapiminusMass_t_Delta;
		TH2F* dHist_omega2piMass_t_Delta0; 
		TH2F* dHist_omegapiplusMass_t_Delta0; 

		TH1F* dHist_CosThetaPlus;
		TH1F* dHist_CosThetaMinus;
		TH1F* dHist_PhiPlus;
		TH1F* dHist_PhiMinus;

		TH1F* dHist_CosTheta_omega_fromPlus;
		TH1F* dHist_Phi_omega_fromPlus;
		TH1F* dHist_CosTheta_omega_fromMinus;
		TH1F* dHist_Phi_omega_fromMinus;
		TH1F* dHist_CosTheta_H;
		TH1F* dHist_Phi_H;

		TH2F* dHist_CosThetaPlusVsOmega2PiMass;
		TH2F* dHist_CosThetaMinusVsOmega2PiMass;
		TH2F* dHist_PhiPlusVsOmega2PiMass;
		TH2F* dHist_PhiMinusVsOmega2PiMass;

		TH2F* dHist_CosTheta_omegaVsOmegaPiPlusMass;
		TH2F* dHist_CosTheta_omegaVsOmegaPiMinusMass;
		TH2F* dHist_Phi_omegaVsOmegaPiPlusMass;
		TH2F* dHist_Phi_omegaVsOmegaPiMinusMass;

		TH2F* dHist_PhiVsCosTheta_b1_Plus[10];
		TH2F* dHist_PhiVsCosTheta_b1_Minus[10];
		TH2F* dHist_PhiVsCosTheta_b1_Plus_total;
		TH2F* dHist_PhiVsCosTheta_b1_Minus_total;

		TH2F* dHist_PhiVsCosTheta_omega_Plus;
		TH2F* dHist_PhiVsCosTheta_omega_Minus;

		TH2F* dHist_PhiVsCosTheta_H;

		TH2F* dHist_M_CosTheta_omega_Plus;
		TH2F* dHist_M_CosTheta_omega_Minus;
		TH2F* dHist_M_CosThetaH_plus;
		TH2F* dHist_M_CosThetaH_minus;

		TH2F* dHist_M_Phi_omega_Plus;
		TH2F* dHist_M_Phi_omega_Minus;
		TH2F* dHist_M_PhiH_plus;
		TH2F* dHist_M_PhiH_minus;

		TH2F* dHist_CosThetaPlusVs2PiMass;
		TH2F* dHist_CosThetaMinusVs2PiMass;
		TH2F* dHist_PhiPlusVs2PiMass;
		TH2F* dHist_PhiMinusVs2PiMass;

		TH2F* dHist_CosTheta_omegaPlusVs2PiMass;
		TH2F* dHist_CosTheta_omegaMinusVs2PiMass;
		TH2F* dHist_Phi_omegaPlusVs2PiMass;
		TH2F* dHist_Phi_omegaMinusVs2PiMass;

		TH2F* dHist_CosThetaPlusVsProtonPiMinusMass;
		TH2F* dHist_PhiPlusVsProtonPiMinusMass;
		TH2F* dHist_CosThetaMinusVsProtonPiPlusMass;
		TH2F* dHist_PhiMinusVsProtonPiPlusMass;

		TH1F* dHist_CosThetaOmega_OR;
		TH1F* dHist_PhiOmega_OR;
		TH1F* dHist_CosThetaRho_OR;
		TH1F* dHist_PhiRho_OR;

		TH2F* dHist_CosThetaOmega_ORVsOmega2PiMass;
		TH2F* dHist_PhiOmega_ORVsOmega2PiMass;
		TH2F* dHist_CosThetaRho_ORVsOmega2PiMass;
		TH2F* dHist_PhiRho_ORVsOmega2PiMass;

		TH2F* dHist_CosThetaOmega_ORVs3PiMass;
		TH2F* dHist_PhiOmega_ORVs3PiMass;
		//TH2F* dHist_CosThetaOmega_ORVs3PiMass_omegacut;
		//TH2F* dHist_PhiOmega_ORVs3PiMass_omegacut;

		TH2F* dHist_CosThetaRho_ORVs2PiMass;
		TH2F* dHist_PhiRho_ORVs2PiMass;
		TH2F* dHist_CosThetaRho_ORVs2PiMass_omegacut;
		TH2F* dHist_PhiRho_ORVs2PiMass_omegacut;

		TH2F* dHist_CosThetaOmegaVsPhiOmega_OR;
		TH2F* dHist_CosThetaRhoVsPhiRho_OR;

		TH2F* dHist_Dalitz_OmegaPi[10];

		TH1F* dHist_5PiMass_DeltaPlusPlusPeak;
		TH1F* dHist_Omega2PiMass_DeltaPlusPlusPeak;
		TH1F* dHist_OmegaPiMinusMass_DeltaPlusPlusPeak;
		TH1F* dHist_ProtonPiPlusMass_DeltaPlusPlusPeak;
		TH1F* dHist_t_DeltaPlusPlus_DeltaPlusPlusPeak;
		TH1F* dHist_t_proton_DeltaPlusPlusPeak;
		TH1F* dHist_t_omega_DeltaPlusPlusPeak;
		TH1F* dHist_CosTheta_DeltaPlusPlusPeak;
		TH1F* dHist_Phi_DeltaPlusPlusPeak;
		TH1F* dHist_CosThetaH_DeltaPlusPlusPeak;
		TH1F* dHist_PhiH_DeltaPlusPlusPeak;
		TH2F* dHist_CosThetaVsPhi_DeltaPlusPlusPeak;
		TH2F* dHist_CosThetaHVsPhiH_DeltaPlusPlusPeak;

		TH2F* dHist_CosTheta_M_omega2pi_DeltaPlusPlusPeak;
		TH2F* dHist_CosTheta_M_omegapiminus_DeltaPlusPlusPeak;
		TH2F* dHist_CosTheta_M_protonpiminus_DeltaPlusPlusPeak;
		TH2F* dHist_CosTheta_M_omegapiplus_DeltaPlusPlusPeak;
		TH2F* dHist_CosTheta_M_protonpiplus_DeltaPlusPlusPeak;
		TH2F* dHist_CosTheta_M_omega_DeltaPlusPlusPeak;
		TH2F* dHist_CosTheta_M_proton2pi_DeltaPlusPlusPeak;
		TH2F* dHist_CosTheta_M_2pi_DeltaPlusPlusPeak;

		TH2F* dHist_M_proton2pi_protonpiplus_DeltaPlusPlusPeak;
		TH2F* dHist_M_proton2pi_protonpiminus_DeltaPlusPlusPeak;
		TH2F* dHist_M_proton2pi_pipluspiminus_DeltaPlusPlusPeak;

		TH2F* dHist_Phi_lab_ProtonPiPlus_vs_Ebeam;
		TH2F* dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak;
		TH2F* dHist_Phi_lab_ProtonPiPlus_vs_Ebeam_DeltaPlusPlusPeak_high_t;
		TH2F* dHist_Omega_Dalitz_xy_DeltaPlusPlusPeak;

	ClassDef(DSelector_pomega2pi_omega3pi, 0);
};

void DSelector_pomega2pi_omega3pi::Get_ComboWrappers(void)
{
	//Step 0
        dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
        dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
        dPiPlus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
        dPiPlus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
        dPiMinus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
        dPiMinus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(4));
        dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(5));

        //Step 1
        dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
        dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
        dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
        dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

}

#endif // DSelector_pomega2pi_omega3pi_h
