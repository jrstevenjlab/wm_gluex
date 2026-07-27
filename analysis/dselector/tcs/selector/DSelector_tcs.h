#ifndef DSelector_tcs_h
#define DSelector_tcs_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_tcs : public DSelector
{
	public:

		DSelector_tcs(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_tcs(){}

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
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DChargedTrackHypothesis* dElectronWrapper;
		DChargedTrackHypothesis* dPositronWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_Mee, *dHist_Mee_CDCdEdx, *dHist_Mee_ShowerCut;
		TH1I* dHist_Mee_KinFit1, *dHist_Mee_KinFit2, *dHist_Mee_KinFit3, *dHist_Mee_KinFit4;
		TH2I* dHist_EleEOverP_PosEOverP_Jpsi, *dHist_EleEOverP_PosEOverP_Bkgd;
		TH2I* dHist_EOverP_Mee[2], *dHist_P_Mee[2], *dHist_Theta_Mee[2], *dHist_P_Mee_CDCdEdx[2], *dHist_Theta_Mee_CDCdEdx[2];
		TH2I* dHist_EOverP_Preshower_Jpsi[2], *dHist_EOverP_Preshower_Bkgd[2];

		TH2I *dHistEOverP_P_BCAL_Jpsi[2], *dHistEOverP_P_FCAL_Jpsi[2], *dHistEOverP_Theta_BCAL_Jpsi[2],*dHistEOverP_Theta_FCAL_Jpsi[2] , *dHistEOverP_Preshower_BCAL_Jpsi[2], *dHistFDCdEdx_P_Jpsi[2], *dHistCDCdEdx_P_Jpsi[2], *dHistCDCdEdx_Theta_Jpsi[2], *dHistCDCP_Theta_Jpsi[2], *dHistBCALDeltaT_P_Jpsi[2], *dHistTOFDeltaT_P_Jpsi[2];
		TH2I *dHistEOverP_P_BCAL_Bkgd[2], *dHistEOverP_P_FCAL_Bkgd[2], *dHistEOverP_Theta_BCAL_Bkgd[2],*dHistEOverP_Theta_FCAL_Bkgd[2] , *dHistEOverP_Preshower_BCAL_Bkgd[2], *dHistFDCdEdx_P_Bkgd[2], *dHistCDCdEdx_P_Bkgd[2], *dHistCDCdEdx_Theta_Bkgd[2], *dHistCDCP_Theta_Bkgd[2], *dHistBCALDeltaT_P_Bkgd[2], *dHistTOFDeltaT_P_Bkgd[2];

		TH2I* dHist_SigTrans_Preshower_Jpsi[2], *dHist_SigTrans_Preshower_Bkgd[2];
		TH2I *dHist_SigTheta_Preshower_Jpsi[2], *dHist_SigTheta_Preshower_Bkgd[2];
		TH2I *dHist_SigTrans_SigTheta_Jpsi[2], *dHist_SigTrans_SigTheta_Bkgd[2];
		TH2I *dHist_SigTrans_SigLong_Jpsi[2], *dHist_SigTrans_SigLong_Bkgd[2];
		TH2I *dHist_SigLong_SigTheta_Jpsi[2], *dHist_SigLong_SigTheta_Bkgd[2];

	ClassDef(DSelector_tcs, 0);
};

void DSelector_tcs::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dElectronWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPositronWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_tcs_h
