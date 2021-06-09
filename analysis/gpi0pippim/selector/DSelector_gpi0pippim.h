#ifndef DSelector_gpi0pippim_h
#define DSelector_gpi0pippim_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_gpi0pippim : public DSelector
{
	public:

		DSelector_gpi0pippim(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_gpi0pippim(){}

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
		DChargedTrackHypothesis* dPiPlusWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;
		TH1I* dHist_Pi0Mass, *dHist_Pi0GammaMass, *dHist_OmegaPiMinusMass;
		TH2I* dHist_Pi0GammaVsPi0Mass;
		TH1I* dHist_3PiMass, *dHist_3Pi2Mass, *dHist_3Pi3Mass;
		TH1I* dHist_Pi02Mass, *dHist_Pi03Mass;
		TH1I* dHist_ProtonPiPlusMass, *dHist_Pi0PiMinusMass;
		TH2I* dHist_Pi0PiMinusGammaVsPi0Gamma, *dHist_Pi0PiMinusGammaDiffVsPi0Gamma;
		TH2I* dHist_Pi0PiMinusGammaVsPi0PiMinus, *dHist_KinFitChiSqVsPi0GammaMass;
		TH2I* dHist_Pi0PiMinusVsPi0Gamma, *dHist_Pi0PiMinusVsPi0GammaDiff;

	ClassDef(DSelector_gpi0pippim, 0);
};

void DSelector_gpi0pippim::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(4));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_gpi0pippim_h
