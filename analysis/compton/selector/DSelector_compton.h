#ifndef DSelector_compton_h
#define DSelector_compton_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_compton : public DSelector
{
	public:

		DSelector_compton(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_compton(){}

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
		DNeutralParticleHypothesis* dPhotonWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		// DEFINE YOUR HISTOGRAMS HERE

		// index is locRFBin=2 and locBeamEBin = 5
		TH1I* dHist_BeamEnergy[2];
		TH2I* dHist_PhotonP_Theta_Init[2][5], *dHist_ProtonP_Theta_Init[2][5], *dHist_PhotonTheta_ThetaDet_Init[2][5];
		TH1I* dHist_VertexZ[2][5], *dHist_VertexR[2][5];
		TH1I* dHist_MissingMassSquared[2][5], *dHist_MissingEnergy[2][5];
		TH2I* dHist_UnusedEnergyBCAL_t[2][5], *dHist_UnusedEnergyFCAL_t[2][5], *dHist_UnusedEnergyTotal_t[2][5];
		TH2I* dHist_PhotonThetaMeasure_PhotonThetaMissing[2][5], *dHist_DeltaE_DeltaTheta_BCAL[2][5];
		TH2I* dHist_DeltaPhi_DeltaTheta_BCAL[2][5], *dHist_DeltaPhi_DeltaTheta_FCAL[2][5];
		TH2I* dHist_DeltaPhi_t[2][5];
		
		// indices are locRFBin=2 and locBeamEBin = 5 (old locDeltaPhiBin=2 Signal and Sideband)
		TH2I* dHist_ProtonP_Theta_Final[2][5], *dHist_PhotonP_Theta_Final[2][5], *dHist_ThetaCM_ProtonTheta[2][5];
		TH2I* dHist_PhotonTheta_ThetaDet_Final[2][5];
		TH2I* dHist_ThetaCM_PhotonTheta[2][5], *dHist_ThetaCM_t[2][5];
		TH2I* dHist_u_t[2][5];

		TH2I* dHist_BCALSigTrans_SigLong[2][5], *dHist_BCALSigTrans_SigTheta[2][5];

		// indices are locRFBin=2 and locCutCounter=5
		TH2I* dHist_ProtonPhi_t[2][5], *dHist_ProtonPhi_ThetaCM[2][5];

	ClassDef(DSelector_compton, 0);
};

void DSelector_compton::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPhotonWrapper = static_cast<DNeutralParticleHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_compton_h
