#ifndef DSelector_pi0kpkm_h
#define DSelector_pi0kpkm_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"

#include "TFolder.h"
#include "TDirectory.h"

class DSelector_pi0kpkm : public DSelector
{
	public:

		DSelector_pi0kpkm(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_pi0kpkm(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		float pi_const;
		static const int num_bars=22;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dKPlusWrapper;
		DChargedTrackHypothesis* dKMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DKinematicData* dDecayingPi0Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;

		TH1I* h2k_invariantmass;
		TH1I* h2kpi0_invariantmass;

		TH1I* h2kpi0_invariantmass_below_phi;
//		TFolder* f1;
//		TDirectory* g1;

		TH1I* h2kpi0_invariantmass_a;
		TH1I* h2kpi0_invariantmass_b;
		TH1I* h2kpi0_invariantmass_c;

		TH1I* hkppi0_invariantmass;
		TH1I* hkmpi0_invariantmass;

		TH1I* hpi0_mass;
		TH1I* hpi0_mass_measured;

		TH1I* dHist_2kpi_t;
		TH1I* dHist_2k_t;
		TH2I* dHist_2kpiMass_t;
		TH2I* dHist_t2k_2kpiMass;
		TH2I* dHist_t_t2k;

//		TH2I* dalitz_pikp_pikm_low;
//		TH2I* dalitz_pikp_pikm_mid; 
//		TH2I* dalitz_pikp_pikm_high;
//		TH2I* dalitz_pikp_pikm_ultra;

//		TH2I* dalitz_kp_km;
//		TH2I* dalitz_kp_pi;
//		TH2I* dalitz_km_pi;

		TH2I* dalitz_kppi0_kmpi0;
		TH2I* dalitz_kppi0_kmpi0_sq;

		TH2I* dalitz_kpkmpi0_kpkm;
		TH2I* dalitz_kpkmpi0_kpkm_sq;

		TH2I* dalitz_kpkm_kppi0;
		TH2I* dalitz_kpkm_kppi0_sq;

		TH2I* dalitz_kpkm_kmpi0;
		TH2I* dalitz_kpkm_kmpi0_sq;

		TH2I* dalitz_kppi0_kpkmpi0;
		TH2I* dalitz_kppi0_kpkmpi0_sq;

		TH2I* dalitz_kmpi0_kpkmpi0;
		TH2I* dalitz_kmpi0_kpkmpi0_sq;

		TH3I* dalitz_t2kpi_t2k_kpkmpi0;
		TH3I* dalitz_t2kpi_t2k_kpkm;

		TH3I* dalitz_kpkmpi0_kpkm_t2k;
		TH3I* dalitz_kpkmpi0_kpkm_t2kpi;

		TH3I* dalitz_kpkmpi0_kpkm_t2k_1;
		TH3I* dalitz_kpkmpi0_kpkm_t2k_2;
		TH3I* dalitz_kpkmpi0_kpkm_t2k_3;
		TH3I* dalitz_kpkmpi0_kpkm_t2k_4;

		TH3I* dalitz_kpkmpi0_kpkm_t2kpi_1;
		TH3I* dalitz_kpkmpi0_kpkm_t2kpi_2;
		TH3I* dalitz_kpkmpi0_kpkm_t2kpi_3;
		TH3I* dalitz_kpkmpi0_kpkm_t2kpi_4;

		Double_t boundary_1, boundary_2, boundary_3; 

		//-----------------------------------------------
		// DIRC specific histograms

		TH1I* hDIRC_l_photon;
		TH1I* hDIRC_l_pi;
		TH1I* hDIRC_l_K;

		TH2I* hPID_extrpl_track_X_Y;

		//-----------------------------------------------

		TH2I* hPID_l_pi_momentum;
		TH2I* hPID_l_K_momentum;
		TH2I* hPID_l_diff_momentum;
		
		TH2I* hPID_l_diff_momentum_rho;
		TH2I* hPID_l_diff_momentum_omega;

		TH2I* hPID_l_diff_momentum_phi;
		TH2I* hPID_l_diff_momentum_below_phi;
		TH2I* hPID_l_diff_momentum_above_phi;

		TH2I* hPID_l_diff_momentum_km;
		TH2I* hPID_l_diff_momentum_km_phi;
		TH2I* hPID_l_diff_momentum_km_below_phi;
		TH2I* hPID_l_diff_momentum_km_above_phi;

		TH1I* hPID_momentum;
		TH1I* hPID_momentum_km;

		TH2I* hPID_l_diff_momentum_kpcut_km;

		TH1I* h2kpi0_invariantmass_below_phi_dirc;
		TH1I* h2k_invariantmass_dirc;
		TH1I* h2kpi0_invariantmass_dirc;

		TH1I* h2kpi0_invariantmass_below_phi_pi0_cut;

		TH1I* h2kpi0_invariantmass_below_phi_pi0_sideband;
		TH1I* h2kpi0_invariantmass_below_phi_pi0_lo_sideband;
		TH1I* h2kpi0_invariantmass_below_phi_pi0_hi_sideband;

		TH2I* hPID_l_momentum_Kp_Km;

		TH2I* h2gamma_2kpi0_mass;

		TH2I* hist_deltaLL_momemtum_bar[num_bars];



		TH1I* hist_extrapolatedY_check;


		//-------------------------------
		// Angle

		TH2D* hAngle_GJ_Phi_pi0;
		TH2D* hAngle_GJ_a0_pi0;
		TH2D* hAngle_GJ_Kstar_Kplus;
		TH2D* hAngle_GJ_Kstar_Kminus;

		TH2D* hAngle_TY_Phi_pi0;

		ClassDef(DSelector_pi0kpkm, 0);

};

void DSelector_pi0kpkm::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dKMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_pi0kpkm_h
