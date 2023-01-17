#ifndef DSelector_gpi0pippim_h
#define DSelector_gpi0pippim_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1F.h"
#include "TH2F.h"

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

		bool dIsMC; // will be false for data
		Bool_t isVarStudy = false; // set to false if not studying variable effects on OmegaPiMinus spectrum
		bool inTRange, inKinFitChiSqRange,  inNumUnusedShowersRange, inRejectWrongPi0Range, inReject3PiRange, inDeltaPlusPlusRange, inPi0Range, inPi0GammaRange, inSidebandRange; // for manipulating when cuts applied

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

		// Mass / other histograms:		
		TH1F* dHist_MM2;
		TH1F* dHist_MM2NoPhoton;

		TH1F* dHist_CosTheta;
		TH1F* dHist_Phi;
		TH1F* dHist_CosTheta_H;
		TH1F* dHist_Phi_H;
		TH2F* dHist_CosThetaVsPhi;
		TH2F* dHist_CosTheta_HVsPhi_H;
		TH2F* dHist_TVsCosTheta;
		TH2F* dHist_TVsPhi;
		TH2F* dHist_TVsCosTheta_H;
		TH2F* dHist_TVsPhi_H;

		TH1F* dHist_KinFitChiSq;
		TH1F* dHist_T;
		TH1F* dHist_ShowerQuality;
		TH1F* dHist_EnergyUnusedShowers;
		TH1F* dHist_Pi0Mass, *dHist_Pi02Mass, *dHist_Pi03Mass;
		TH1F* dHist_3PiMass, *dHist_3Pi2Mass, *dHist_3Pi3Mass;
		TH1F* dHist_ProtonPiPlusMass;
		TH2F* dHist_ProtonPiPlusMassVsPi0PiMinusMass, *dHist_ProtonPiPlusMassVsPi0PiMinusMass_PK, *dHist_ProtonPiPlusMassVsPi0PiMinusMass_SB;
		TH1F* dHist_Pi0PiMinusMass;
		TH1F* dHist_Pi0GammaMass;
		TH1F* dHist_OmegaPiMinusMass;
		
		TH2F* dHist_Pi0GammaMassVsMM2;
		TH2F* dHist_Pi0GammaMassVsMM2NoPhoton;
		TH2F* dHist_Pi0GammaMassVsCosTheta;
		TH2F* dHist_Pi0GammaMassVsPhi;
		TH2F* dHist_Pi0GammaMassVsCosTheta_H;
		TH2F* dHist_Pi0GammaMassVsPhi_H;
		TH2F* dHist_Pi0GammaMassVsKinFitChiSq;
		TH2F* dHist_Pi0GammaMassVsT;
		TH2F* dHist_Pi0GammaMassVsShowerQuality;
		TH2F* dHist_Pi0GammaMassVsEnergyUnusedShowers;
		TH2F* dHist_Pi0GammaMassVsPi0Mass;
		TH2F* dHist_Pi0GammaMassVs3PiMass;
		TH2F* dHist_Pi0GammaMassVsProtonPiPlusMass;
		TH2F* dHist_Pi0GammaMassVsPi0PiMinusMass;
		TH2F* dHist_Pi0GammaMassVsOmegaPiMinusMass;

		TH2F* dHist_OmegaPiMinusMassVsPi0PiMinusMass;
		TH2F* dHist_OmegaPiMinusMassVsCosTheta; 
		TH2F* dHist_OmegaPiMinusMassVsPhi;
		TH2F* dHist_OmegaPiMinusMassVsCosTheta_H; 
		TH2F* dHist_OmegaPiMinusMassVsPhi_H;
		
		/*
		  | *********************
		  | * MONTE CARLO PLOTS *
		  | *********************
		 */

		// Thrown Topologies Percentage Plots:
		TH1F* dHistTop_PercentStart;
		TH1F* dHistTop_PercentAngular;
		TH1F* dHistTop_PercentChi2;
		TH1F* dHistTop_PercentT;
		TH1F* dHistTop_PercentPhoton;
		TH1F* dHistTop_PercentPi0;
		TH1F* dHistTop_PercentPi0Combo;
		TH1F* dHistTop_Percent3PiCombo;
		TH1F* dHistTop_PercentProtonPiPlus;
		TH1F* dHistTop_PercentPi0PiMinus;
		TH1F* dHistTop_PercentPi0Gamma;
		TH1F* dHistTop_PercentOmega;
		
		// Mass-topology plots
		TH1F* dHistThrown_BeamEnergy; // no topology, but only fills for MC	
		TH1F* dHistThrown_BeamEnergySig;	

		map<TString, TH1F*> dHistThrown_Pi0GammaMass;
		
		map<TString, TH1F*> dHistTop_MM2;
		map<TString, TH1F*> dHistTop_MM2NoPhoton;
 
		map<TString, TH1F*> dHistTop_CosTheta; 
		map<TString, TH1F*> dHistTop_Phi;
		map<TString, TH1F*> dHistTop_CosTheta_H; 
		map<TString, TH1F*> dHistTop_Phi_H;
		map<TString, TH1F*> dHistThrown_CosTheta; 
		map<TString, TH1F*> dHistThrown_Phi;
		map<TString, TH1F*> dHistThrown_CosTheta_H; 
		map<TString, TH1F*> dHistThrown_Phi_H;

		map<TString, TH2F*> dHistTop_CosThetaVsPhi;
		map<TString, TH2F*> dHistTop_CosTheta_HVsPhi_H;
		map<TString, TH2F*> dHistThrown_CosThetaVsPhi;
		map<TString, TH2F*> dHistThrown_CosTheta_HVsPhi_H;

		map<TString, TH2F*> dHistTop_TVsCosTheta;
		map<TString, TH2F*> dHistTop_TVsPhi;
		map<TString, TH2F*> dHistTop_TVsCosTheta_H;
		map<TString, TH2F*> dHistTop_TVsPhi_H;

		map<TString, TH2F*> dHistThrown_TVsCosTheta;
		map<TString, TH2F*> dHistThrown_TVsPhi;
		map<TString, TH2F*> dHistThrown_TVsCosTheta_H;
		map<TString, TH2F*> dHistThrown_TVsPhi_H;
		
		map<TString, TH1F*> dHistTop_KinFitChiSq;
		map<TString, TH1F*> dHistTop_T;
		map<TString, TH1F*> dHistTop_EnergyUnusedShowers;
		map<TString, TH1F*> dHistTop_ShowerQuality;

		map<TString, TH1F*> dHistTop_Pi0Mass;  
		map<TString, TH1F*> dHistTop_Pi02Mass;  
		map<TString, TH1F*> dHistTop_Pi03Mass;

		map<TString, TH1F*> dHistTop_3PiMass;  
		map<TString, TH1F*> dHistTop_3Pi2Mass;  
		map<TString, TH1F*> dHistTop_3Pi3Mass;

		map<TString, TH1F*> dHistTop_ProtonPiPlusMass;
		map<TString, TH2F*> dHistTop_ProtonPiPlusMassVsPi0PiMinusMass;
		map<TString, TH2F*> dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_PK;
		map<TString, TH2F*> dHistTop_ProtonPiPlusMassVsPi0PiMinusMass_SB;
		
		map<TString, TH1F*> dHistTop_Pi0PiMinusMass;
		map<TString, TH1F*> dHistTop_Pi0GammaMass;
		map<TString, TH1F*> dHistTop_OmegaPiMinusMass;
		
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsMM2;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsMM2NoPhoton;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsCosTheta;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsPhi;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsCosTheta_H;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsPhi_H;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsKinFitChiSq;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsT;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsShowerQuality;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsEnergyUnusedShowers;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsPi0Mass;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVs3PiMass;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsProtonPiPlusMass;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsPi0PiMinusMass;
		map<TString, TH2F*> dHistTop_Pi0GammaMassVsOmegaPiMinusMass;

		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsPi0PiMinusMass;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsCosTheta;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsPhi;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsCosTheta_H;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsPhi_H;
		map<TString, TH2F*> dHistThrown_OmegaPiMinusMassVsCosTheta;
		map<TString, TH2F*> dHistThrown_OmegaPiMinusMassVsPhi;
		map<TString, TH2F*> dHistThrown_OmegaPiMinusMassVsCosTheta_H;
		map<TString, TH2F*> dHistThrown_OmegaPiMinusMassVsPhi_H;		

		// VarStudy 2D hists
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsPi0Gamma;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsKinFitChiSq;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsEnergyUnusedShowers;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsPi0Width;
		map<TString, TH2F*> dHistTop_OmegaPiMinusMassVsShowerQuality;		

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
