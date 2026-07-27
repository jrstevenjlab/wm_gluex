#ifndef DSelector_p2k2pi_h
#define DSelector_p2k2pi_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_p2k2pi : public DSelector
{
	public:

		DSelector_p2k2pi(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_p2k2pi(){}

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
		DChargedTrackHypothesis* dKPlusWrapper;
		DChargedTrackHypothesis* dKMinusWrapper;
		DChargedTrackHypothesis* dPiPlusWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH2I *dHist_2piMass_deltaR, *dHist_2piMass_deltaZ, *dHist_deltaR_deltaZ;
		TH2I *dHist_2kMass_2piMass, *dHist_2kMass_t;
		TH2I *dHist_kppimMass_kmpipMass;

		TH2I *dHist_2kMass_2piMass_vertex, *dHist_2kMass_t_vertex;
		TH2I *dHist_kppimMass_kmpipMass_vertex;

		TH2I *dHist_2kMass_2piMass_noKstar, *dHist_2kMass_t_noKstar;

		TH2I *dHist_2kMass_ProtonPimMass_noKstar,  *dHist_2kMass_ProtonPipMass_noKstar, *dHist_2kMass_2kpimMass_noKstar;

		TH2I *dHist_2kMass_2piMass_noDelta, *dHist_phi2piMass_2piMass_noDelta, *dHist_phi2piMass_2piMass_noDelta_SB;

		TH2I *dHist_2k2piMass_t, *dHist_2k2piMass_t_SB;
		TH2I *dHist_2k2piMass_Egamma, *dHist_2k2piMass_Egamma_SB;

		TH1I* dHist_BeamEnergy;

		ClassDef(DSelector_p2k2pi, 0);
};

void DSelector_p2k2pi::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
        dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
        dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
        dKMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
        dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
        dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(4));
}

#endif // DSelector_p2k2pi_h
