#ifndef DSelector_p2gamma_h
#define DSelector_p2gamma_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_p2gamma : public DSelector
{
	public:

		DSelector_p2gamma(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_p2gamma(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		// SEPARATE HISTOGRAM ACTIONS
		DHistogramAction_InvariantMass* dFinalMassAction;
		DHistogramAction_ParticleComboKinematics *dKinematicsActionu, *dKinematicsActiont;

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH2F* dHist_InvariantMassVsBeamE;
		TH2F* dHist_InvariantMassVsTAGMcolumn, *dHist_InvariantMassVsTAGHcounter;
		TH2F* dHist_InvariantMassVsBeamE_KinFit[6];
		TH2F* dHist_InvariantMassVst, *dHist_InvariantMassVsu;
		TH2F* dHist_InvariantMassKinFitVst, *dHist_InvariantMassKinFitVsu;
		TH2F* dHist_InvariantMassVst_BeamE[13], *dHist_InvariantMassVsu_BeamE[13];
		TH2F* dHist_PhiVst, *dHist_PhiVsu, *dHist_PhiVsEgamma;
		TH2F* dHist_Deltau_u, *dHist_Deltau_t, *dHist_ProtonPVsTheta, *dHist_PhotonPVsTheta, *dHist_ProtonNHitVsTheta, *dHist_ProtondEdxVsP;
		TH1F *dHist_EventPassCut;
		TH2F* dHist_PhotonDeltaPhiVsDeltaTheta;
		
		vector<double> dEventPassCut;

	ClassDef(DSelector_p2gamma, 0);
};

void DSelector_p2gamma::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
}

#endif // DSelector_p2gamma_h
