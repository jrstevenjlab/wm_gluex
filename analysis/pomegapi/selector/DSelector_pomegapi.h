#ifndef DSelector_pomegapi_h
#define DSelector_pomegapi_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_pomegapi : public DSelector
{
	public:

		DSelector_pomegapi(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_pomegapi(){}

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
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DKinematicData* dDecayingPi01Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		//Step 2
		DParticleComboStep* dStep2Wrapper;
		DKinematicData* dDecayingPi02Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;
		DNeutralParticleHypothesis* dPhoton4Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;
		TH1I* dHist_3PiMass_Measured;
		TH1I* dHist_MM2_Weighted;
		TH1I* dHist_lambda_peak;
		TH1I* dHist_lambda_wings;
		TH1I* dHist_lambda_uncut;
		TH1I* dHist_4PiMass;
		TH1I* dHist_OmegaPiMass;
		TH2I* dHist_3vs4;
		TH1I* dHist_Man_t;
		TH1I* dHist_costheta;
		TH1I* dHist_phi;
		TH1I* dHist_costhetaH;
		TH1I* dHist_phiH;
		TH2I* dHist_CosThetaVsMass;
		TH2I* dHist_PhiVsMass;
		TH2I* dHist_CosThetaHVsMass;
		TH2I* dHist_PhiHVsMass;
		TH1I* dHist_CosTheta_t1;
		TH1I* dHist_Phi_t1;
		TH1I* dHist_CosThetaH_t1;
		TH1I* dHist_PhiH_t1;
		TH1I* dHist_CosTheta_t2;
		TH1I* dHist_Phi_t2;
		TH1I* dHist_CosThetaH_t2;
		TH1I* dHist_PhiH_t2;
		TH2I* dHist_CosThetaVsMass_t1;
		TH2I* dHist_PhiVsMass_t1;
		TH2I* dHist_CosThetaHVsMass_t1;
		TH2I* dHist_PhiHVsMass_t1;
		TH2I* dHist_CosThetaVsMass_t2;
		TH2I* dHist_PhiVsMass_t2;
		TH2I* dHist_CosThetaHVsMass_t2;
		TH2I* dHist_PhiHVsMass_t2;

	ClassDef(DSelector_pomegapi, 0);
};

void DSelector_pomegapi::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(4));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dDecayingPi01Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
	dDecayingPi02Wrapper = dStep2Wrapper->Get_InitialParticle();
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_pomegapi_h
