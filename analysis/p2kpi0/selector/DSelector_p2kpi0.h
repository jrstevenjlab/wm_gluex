#ifndef DSelector_p2kpi0_h
#define DSelector_p2kpi0_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_p2kpi0 : public DSelector
{
	public:

		DSelector_p2kpi0(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_p2kpi0(){}

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
		DChargedTrackHypothesis* dKMinusWrapper;
		DChargedTrackHypothesis* dKPlusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DKinematicData* dDecayingPi0Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH2I *dHist_2kMass_t;
                TH2I *dHist_kppiMass_kmpiMass;

                TH2I *dHist_2kMass_ProtonPiMass_noKstar, *dHist_2kMass_2kpiMass_noKstar;

                TH2I *dHist_2kpiMass_t, *dHist_2kpiMass_t_SB;
                TH2I *dHist_2kpiMass_Egamma, *dHist_2kpiMass_Egamma_SB;

                TH1I* dHist_BeamEnergy;

	ClassDef(DSelector_p2kpi0, 0);
};

void DSelector_p2kpi0::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_p2kpi0_h
