#ifndef DSelector_npi_missing_h
#define DSelector_npi_missing_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_npi_missing : public DSelector
{
	public:

		DSelector_npi_missing(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_npi_missing(){}

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
		DChargedTrackHypothesis* dPiPlusWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_BeamDeltaT;
		TH1I* dHist_BeamEnergy[2];
		TH2I* dHist_PiPlusP_Theta_Init[2], *dHist_NeutronP_Theta_Init[2], *dHist_PiPlusTheta_ThetaDet_Init[2];
		TH2I* dHist_NeutronBeta_P_Init[2], *dHist_NeutronDeltaBeta_P_Init[2], *dHist_NeutronBeta_P[2], *dHist_NeutronDeltaBeta_P[2];
                TH1I* dHist_VertexZ[2], *dHist_VertexR[2];
		TH1I* dHist_MissingMass[2], *dHist_MissingMass_Final[2], *dHist_MissingMass_Final_noUnused[2];
		TH1I* dHist_MissingMassSquared[2], *dHist_MissingMassSquared_Final[2], *dHist_MissingMassSquared_Final_noUnused[2];
		TH1I* dHist_UnusedTrack[2];
		TH2I* dHist_UnusedEnergyBCAL_t[2], *dHist_UnusedEnergyFCAL_t[2], *dHist_UnusedEnergyTotal_t[2];
		TH2I* dHist_UnusedEnergyBCAL_FCAL[2];
		TH2I* dHist_DeltaPhi_DeltaTheta[2];
		TH2I* dHist_ShowerTheta_DeltaTheta[2], *dHist_ShowerTheta_DeltaTheta_noUnused[2], *dHist_ShowerE_DeltaTheta[2];
                TH2I* dHist_DeltaPhi_t_Init[2], *dHist_DeltaPhi_t[2], *dHist_DeltaPhi_t_DeltaBetaCut[2], *dHist_DeltaPhi_t_noUnused[2];
		TH2I* dHist_DeltaPhi_Phi_DeltaBetaCut[2];

		TH2I* dHist_NeutronPreshowerE_ShowerTheta_DeltaPhiCut[2];
		TH2I* dHist_NeutronShowerTheta_DeltaPhi[2], *dHist_NeutronPreshowerE_DeltaPhi[2];
		TH2I* dHist_NeutronShowerE_ShowerTheta[2], *dHist_NeutronMissingTheta_ShowerTheta[2], *dHist_NeutronPreshowerE_ShowerE[2];
		TH2I* dHist_PiPlusMomentum_BeamE[2], *dHist_PiPlusMomentum_DeltaPhi[2];
		TH2I* dHist_NeutronPhi_t[2], *dHist_Sideband1NeutronPhi_t[2], *dHist_Sideband2NeutronPhi_t[2];

		TH2I* dHist_NeutronPhi_t_noUnused[2], *dHist_NeutronPhi_t_noUnused_tightDeltaTheta[2], *dHist_NeutronPhi_t_tightMM[2], *dHist_NeutronPhi_t_noUnused_tightMM[2];

		TH1I* dHist_Thrownt;
		TH2I* dHist_ThrownNeutronP_t, *dHist_ThrownPiPlusTheta_t, *dHist_ThrownNeutronP_Theta, *dHist_ThrownPiPlusP_Theta;
		TH2I* dHist_Thrownt_Recot[2];
		TH2I* dHist_ThrownDeltaThetaShower_DeltaThetaMissing[2];
		TH2I* dHist_PiPlusDeltaPOverP_Theta[2];

		TF1* dNeutronDeltaBeta;
		double dDeltaPhiCut, dDeltaThetaMinCut, dDeltaThetaMaxCut, dNeutronShowerThetaCut;
		double dBeamELowCut, dBeamEHighCut;

	ClassDef(DSelector_npi_missing, 0);
};

void DSelector_npi_missing::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
}

#endif // DSelector_npi_missing_h
