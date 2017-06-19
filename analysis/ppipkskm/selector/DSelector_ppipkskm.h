#ifndef DSelector_ppipkskm_h
#define DSelector_ppipkskm_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_ppipkskm : public DSelector
{
	public:

		DSelector_ppipkskm(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_ppipkskm(){}

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
		DChargedTrackHypothesis* dPiPlus1Wrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DKinematicData* dDecayingKShortWrapper;
		DChargedTrackHypothesis* dPiPlus2Wrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH2I *dHist_2piMass_deltaR, *dHist_2piMass_deltaZ, *dHist_deltaR_deltaZ;
		TH2I *dHist_2piMass_kppimMass, *dHist_2kpiMass_kppimMass, *dHist_2kpiMass_kppimMass_SB1, *dHist_2kpiMass_kppimMass_SB2;
		TH2I *dHist_2piMass_kppimMass_vertex, *dHist_2kpiMass_kppimMass_vertex, *dHist_2kpiMass_kppimMass_vertex_SB1, *dHist_2kpiMass_kppimMass_vertex_SB2;

	ClassDef(DSelector_ppipkskm, 0);
};

void DSelector_ppipkskm::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiPlus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dKMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dDecayingKShortWrapper = dStep1Wrapper->Get_InitialParticle();
	dPiPlus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_ppipkskm_h
