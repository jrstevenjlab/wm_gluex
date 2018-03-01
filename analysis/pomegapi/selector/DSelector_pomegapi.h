#ifndef DSelector_pomegapi_h
#define DSelector_pomegapi_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"

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
		TH1F* dHist_MissingMassSquared;
		TH1F* dHist_BeamEnergy;
		TH1F* dHist_3PiMass_Measured;
		TH1F* dHist_MM2_Weighted;
		TH1F* dHist_lambda_peak;
		TH1F* dHist_lambda_wings;
		TH1F* dHist_lambda_uncut;
		TH1F* dHist_4PiMass;
		TH1F* dHist_OmegaPiMass;
		TH2F* dHist_3vs4;
		TH1F* dHist_Man_t;
		TH1F* dHist_costheta;
		TH1F* dHist_phi;
		TH1F* dHist_costhetaH;
		TH1F* dHist_phiH;
		TH2F* dHist_CosThetaVsMass;
		TH2F* dHist_PhiVsMass;
		TH2F* dHist_CosThetaHVsMass;
		TH2F* dHist_PhiHVsMass;
		TH1F* dHist_CosTheta_t1;
		TH1F* dHist_Phi_t1;
		TH1F* dHist_CosThetaH_t1;
		TH1F* dHist_PhiH_t1;
		TH1F* dHist_CosTheta_t2;
		TH1F* dHist_Phi_t2;
		TH1F* dHist_CosThetaH_t2;
		TH1F* dHist_PhiH_t2;
		TH2F* dHist_CosThetaVsMass_t1;
		TH2F* dHist_PhiVsMass_t1;
		TH2F* dHist_CosThetaHVsMass_t1;
		TH2F* dHist_PhiHVsMass_t1;
		TH2F* dHist_CosThetaVsMass_t2;
		TH2F* dHist_PhiVsMass_t2;
		TH2F* dHist_CosThetaHVsMass_t2;
		TH2F* dHist_PhiHVsMass_t2;
		//Moment sums for all t
		TH1F* dHist_H0000;
		TH1F* dHist_H0020;
		TH1F* dHist_H0021;
		TH1F* dHist_H0022;
		TH1F* dHist_H2000;
		TH1F* dHist_H2020;
		TH1F* dHist_H2021;
		TH1F* dHist_H2022;
		TH1F* dHist_H2120;
		TH1F* dHist_H2121_plus;
		TH1F* dHist_H2122_plus;
		TH1F* dHist_H2220;
		TH1F* dHist_H2221_plus;
		TH1F* dHist_H2222_plus;
		TH1F* dHist_H2111_plus;
		TH1F* dHist_H0010;
		TH1F* dHist_H0011;
		TH1F* dHist_H2110;
		TH1F* dHist_H2111_minus;
		TH1F* dHist_H2121_minus;
		TH1F* dHist_H2122_minus;
		TH1F* dHist_H2221_minus;
		TH1F* dHist_H2222_minus;
		TH1F* dHist_H2010;
		TH1F* dHist_H2011;
		//Moment sums in t bins
		TH1F* dHist_HlmLM_t[25][5];



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
