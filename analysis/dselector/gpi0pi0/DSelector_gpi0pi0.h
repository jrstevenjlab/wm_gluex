#ifndef DSelector_gpi0pi0_h
#define DSelector_gpi0pi0_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_gpi0pi0 : public DSelector
{
	public:

		DSelector_gpi0pi0(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_gpi0pi0(){}

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
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;

		//Step 2
		DParticleComboStep* dStep2Wrapper;
		DNeutralParticleHypothesis* dPhoton4Wrapper;
		DNeutralParticleHypothesis* dPhoton5Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;

	ClassDef(DSelector_gpi0pi0, 0);
};

void DSelector_gpi0pi0::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPhoton5Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_gpi0pi0_h
