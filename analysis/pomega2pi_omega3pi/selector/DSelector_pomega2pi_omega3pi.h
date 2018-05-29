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

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPiPlus1Wrapper;
		DChargedTrackHypothesis* dPiMinus1Wrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DChargedTrackHypothesis* dPiPlus2Wrapper;
		DChargedTrackHypothesis* dPiMinus2Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;

		//Step 2
		DParticleComboStep* dStep2Wrapper;
		DKinematicData* dDecayingPi0Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;

	ClassDef(DSelector_pomega2pi_omega3pi, 0);
};

void DSelector_pomega2pi_omega3pi::Get_ComboWrappers(void)
{
	//Step 0
        dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
        dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
        dPiMinus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
        dPiMinus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
        dPiPlus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
        dPiPlus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(4));
        dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(5));

        //Step 1
        dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
        dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
        dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
        dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

}

#endif // DSelector_pomega2pi_omega3pi_h
