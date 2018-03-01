#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include "clebschGordan.h"
#include "clebschGordan.cc"
#include "barrierFactor.h"
#include "barrierFactor.cc"
#include "breakupMomentum.h"
#include "breakupMomentum.cc"

#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"
 
TH1F* h[25]; 

//Define breakup momentum (x here is the measured mass of the omegapi)
double q1(double x) {
  double pi0mass = 0.1349766; //mass of the pi0 in GeV
  double omegamass = 0.78265; //mass of the omega in GeV
  if (x < pi0mass + omegamass)
    return 0.0;
  return 0.5 * TMath::Sqrt((TMath::Power(x, 2.0) - (pi0mass + omegamass)*(pi0mass + omegamass))*(TMath::Power(x, 2.0) - (omegamass - pi0mass)*(omegamass - pi0mass))) / x;
}

//Define barrier factor ratio
double barrierratio(double x, double resonancemass, int l) {
  return barrierFactor(q1(x), l) / barrierFactor(q1(resonancemass), l);
}

double Gamma_alpha(double x, double resonancewidth, double resonancemass, int l) {
  return resonancewidth * q1(x) / q1(resonancemass) * TMath::Power(barrierratio(x, resonancemass, l), 2);
}

std::complex<double> D_alpha(double x, double resonancemass, double resonancewidth, int l) {
  std::complex<double> denom = std::complex<double>(TMath::Power(x, 2) - TMath::Power(resonancemass, 2),-1.0*resonancemass * Gamma_alpha(x, resonancewidth, resonancemass, l));
  return resonancewidth * resonancemass / denom;
}

//Define J (spin) for a given ialpha
int J_spin(int ialpha) { //ialpha{0,1,2} -> JP{1+,1-,0-}
  if (ialpha == 2)
    return 0;
  return 1;
}

//Define a parity function
int eta_parity(int ialpha) {
  if (ialpha == 0)
    return 1;
  return -1;
}

//Define sum over l with Clebsch-Gordan coefficients and barrier factors
double lsum(double x, double resonancemass, int ialpha, int lambda, double DoverS) {
  double c_alpha; //partial wave amplitudes
  double lsum = 0;
  for (int l = 0; l < 3; l++) { 
    if (l == 1 && (ialpha == 1 || ialpha == 2))
      c_alpha = 1;
    else if (l == 2 && ialpha == 0)
      c_alpha = DoverS;
    else if (l == 0 && ialpha == 0)
      c_alpha = 1;
    else 
      c_alpha = 0; 
    lsum += TMath::Sqrt((2.0*l + 1.0)/(2.0*J_spin(ialpha) + 1.0)) * clebschGordan(l, 1, 0, lambda, J_spin(ialpha), lambda) * c_alpha * barrierratio(x, resonancemass, l);
  }
  return lsum;
}

std::complex<double> F_alpha_lambda(double x, double resonancemass, double resonancewidth, double G_alpha, int ialpha, int lambda, double DoverS, int l) {
  return D_alpha(x, resonancemass, resonancewidth, l) * G_alpha * lsum(x, resonancemass, ialpha, lambda, DoverS);
}

//Define sum over omega spin states
std::complex<double> f_Llm(double x, double resonancemass, double resonancewidth, double G_alpha, int ialpha, double resonancemass2, double resonancewidth2, double G_beta, int ibeta, double DoverS, int l, int m, int L) {
  std::complex<double> fsum = std::complex<double>(0, 0);
  for (int lambda = -1; lambda < 2; lambda++) {
    for (int lambdaprime = -1; lambdaprime < 2; lambdaprime++) {
      if (ibeta == 2 && lambdaprime != 0)
	continue;
      if (ialpha == 2 && lambda != 0)
	continue;
      fsum += F_alpha_lambda(x, resonancemass, resonancewidth, G_alpha, ialpha, lambda, DoverS, l) * std::conj(F_alpha_lambda(x, resonancemass2, resonancewidth2, G_beta, ibeta, lambdaprime, DoverS, l)) * clebschGordan(J_spin(ibeta), L, lambdaprime, m, J_spin(ialpha), lambda) * clebschGordan(1, l, lambdaprime, m, 1, lambda);
    }
  }
  return fsum;
}

const std::complex<double> ic(0, 1);

//Define complex "amplitudes"
std::complex<double> f_0(double phi0, double theta) { 
  return TMath::Sqrt(0.5) * std::exp(ic * phi0) * TMath::Cos(theta);
}
std::complex<double> f_plus(double phiplus, double theta, double psi) {
  return TMath::Sqrt(0.5) * std::exp(ic * phiplus) * TMath::Sin(theta) * TMath::Cos(psi);
}
std::complex<double> f_minus(double phiminus, double theta, double psi) {
  return TMath::Sqrt(0.5) * std::exp(ic * phiminus) * TMath::Sin(theta) * TMath::Sin(psi);
}

//Define a helicity function
int Lambda_H(int iH) { // {0,1,2}->{0,+1,-1}
  if (iH == 0)
    return 0;
  else if (iH == 1) 
    return 1;
  else
    return -1;
}


//Define production density matrix
std::complex<double> rho(double phi0, double theta, double phiplus, double phiminus, double psi, double phi02, double theta2, double phiplus2, double phiminus2, double psi2, int iH, int iHPrime, int ialpha, int ibeta) {
  std::complex<double> f_alpha = std::complex<double>(0,0);
  std::complex<double> f_beta = std::complex<double>(0,0); 
  std::complex<double> f_alpha_neg = std::complex<double>(0,0);
  std::complex<double> f_beta_neg = std::complex<double>(0,0); 
  if (iH == 0) {
    f_alpha = f_0(phi0, theta);
    f_alpha_neg = f_0(phi0, theta);
  }
  else if (iH == 1) {
    f_alpha = f_plus(phiplus, theta, psi);
    f_alpha_neg = f_minus(phiminus, theta, psi);
  }
  else { 
    f_alpha = f_minus(phiminus, theta, psi);
    f_alpha_neg = f_plus(phiplus, theta, psi);
  }
  if (iHPrime == 0) {
    f_beta = f_0(phi02, theta2);
    f_beta_neg = f_0(phi02, theta2);
  }
  else if (iHPrime == 1) {
    f_beta = f_plus(phiplus2, theta2, psi2);
    f_beta_neg = f_minus(phiminus2, theta2, psi2);
  }
  else { 
    f_beta = f_minus(phiminus2, theta2, psi2);
    f_beta_neg = f_plus(phiplus2, theta2, psi2);
  }
  return f_alpha * std::conj(f_beta) + eta_parity(ialpha) * eta_parity(ibeta) * TMath::Power(-1.0, J_spin(ialpha) - J_spin(ibeta)) * TMath::Power(-1.0, Lambda_H(iH) - Lambda_H(iHPrime)) * f_alpha_neg * std::conj(f_beta_neg);
}


//Define sum over helicity states
std::complex<double> HelicitySum(double phi0, double theta, double phiplus, double phiminus, double psi, double phi02, double theta2, double phiplus2, double phiminus2, double psi2, int ialpha, int ibeta, int L, int M) {
  std::complex<double> sumH(0, 0);
  for (int iH = 0; iH < 3; iH++) {
    for (int iHPrime = 0; iHPrime < 3; iHPrime++) {
      if (ialpha == 2 && iH > 0)
	continue;
      if (ibeta == 2 && iHPrime > 0)
	continue;
      sumH += rho(phi0, theta, phiplus, phiminus, psi, phi02, theta2, phiplus2, phiminus2, psi2, iH, iHPrime, ialpha, ibeta) * clebschGordan(J_spin(ibeta), L, Lambda_H(iHPrime), M, J_spin(ialpha), Lambda_H(iH));
    }
  }
  return sumH;
}

//Define t_star (sum of production density matrix)
std::complex<double> t_star_LM(double phi0, double theta, double phiplus, double phiminus, double psi, double phi02, double theta2, double phiplus2, double phiminus2, double psi2, int ialpha, int ibeta, int L, int M) {
  return TMath::Sqrt((2.0*J_spin(ibeta) + 1.0)/(2.0*J_spin(ialpha) + 1.0)) * HelicitySum(phi0, theta, phiplus, phiminus, psi, phi02, theta2, phiplus2, phiminus2, psi2, ialpha, ibeta, L, M);
}

//Create array of lmLM:
int lmLM[25][4] = {{0,0,0,0}, {0,0,2,0}, {0,0,2,1}, {0,0,2,2}, {2,0,0,0}, {2,0,2,0}, {2,0,2,1}, {2,0,2,2}, {2,1,2,0}, {2,1,2,1}, {2,1,2,2}, {2,2,2,0}, {2,2,2,1}, {2,2,2,2}, {2,1,1,1}, {0,0,1,0}, {0,0,1,1}, {2,1,1,0}, {2,1,1,1}, {2,1,2,1}, {2,1,2,2}, {2,2,2,1}, {2,2,2,2}, {2,0,1,0}, {2,0,1,1}};

double SingleIntensity(double *x, double *par, int l, int m, int L, int M, int ialpha, int ibeta) {
    double phiplus_alpha = 0;
    double phiminus_alpha = 0;
    double psi_alpha = 0;
    double phiplus_beta = 0;
    double phiminus_beta = 0;
    double psi_beta = 0;
    if (ialpha < 2) {
      phiplus_alpha = par[5*ialpha + 12];
      phiminus_alpha = par[5*ialpha + 13];
      psi_alpha = par[5*ialpha + 14];
    }
    if (ibeta < 2) {
      phiplus_beta = par[5*ibeta + 12];
      phiminus_beta = par[5*ibeta + 13];
      psi_beta = par[5*ibeta + 14];
    }
  return std::real(t_star_LM(par[5*ialpha + 10], par[5*ialpha + 11], phiplus_alpha, phiminus_alpha, psi_alpha, par[5*ibeta + 10], par[5*ibeta + 11], phiplus_beta, phiminus_beta, psi_beta, ialpha, ibeta, L, M) * f_Llm(x[0], par[3*ialpha + 0], par[3*ialpha + 1], par[3*ialpha + 2], ialpha, par[3*ibeta + 0], par[3*ibeta + 1], par[3*ibeta + 2], ibeta, par[9], l, m, L) * clebschGordan(1, l, 0, 0, 1, 0));
}

double IntensityPlus(double *x, double *par, int l, int m, int L, int M) {
  double IntensityPlus = 0;
  for (int ialpha = 0; ialpha < 3; ialpha++) {
    int ibeta = ialpha;
    IntensityPlus += SingleIntensity(x, par, l, m, L, M, ialpha, ibeta);
  } 
  return IntensityPlus;
}

double OneIntensity_1plus(double *x, double *par, int i) {
  int ialpha = 0;
  int ibeta = 0;
  return SingleIntensity(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
}

double OneIntensityPlot_1plus(double *x, double *par) {
 double par22[22];
  for (int i = 0; i < 22; i++) { //loop over 22 fit parameters
    par22[i] = par[i];
  }
  int index = (int)par[22];
  return OneIntensity_1plus(x, par22, index);
}

double OneIntensity_1minus(double *x, double *par, int i) {
  int ialpha = 1;
  int ibeta = 1;
  return SingleIntensity(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
}

double OneIntensityPlot_1minus(double *x, double *par) {
 double par22[22];
  for (int i = 0; i < 22; i++) { //loop over 22 fit parameters
    par22[i] = par[i];
  }
  int index = (int)par[22];
  return OneIntensity_1minus(x, par22, index);
}

double OneIntensity_0minus(double *x, double *par, int i) {
  int ialpha = 2;
  int ibeta = 2;
  return SingleIntensity(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
}

double OneIntensityPlot_0minus(double *x, double *par) {
 double par22[22];
  for (int i = 0; i < 22; i++) { //loop over 22 fit parameters
    par22[i] = par[i];
  }
  int index = (int)par[22];
  return OneIntensity_0minus(x, par22, index);
}

double IntensityMinus(double *x, double *par, int i) {
  double IntensityMinus = 0;
  for (int ialpha = 0; ialpha < 3; ialpha++) { //double sum over combination of JP states
    for (int ibeta = 0; ibeta < 3; ibeta++) {
      if (ialpha < ibeta) //avoid doublecounting
	continue;
      if (ialpha == ibeta) //Only want interference moments
	continue;
      if (eta_parity(ialpha) == eta_parity(ibeta)) //We only consider states with opposite parity
	continue;
      IntensityMinus += 2.0 * SingleIntensity(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
      //cout << "Interference of state " << J_spin(ialpha) << eta_parity(ialpha) << " with " << J_spin(ibeta) << eta_parity(ibeta) << endl;
    }
  }
  return IntensityMinus; 
}


double IntensityCombo_11(double *x, double *par, int i) {//ialpha = 0, ibeta = 1
  int ialpha = 1;
  int ibeta = 0;
  return 2.0 * SingleIntensity(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
}

double IntensityCombo_11plot(double *x, double *par) {
 double par22[22];
  for (int i = 0; i < 22; i++) { //loop over 22 fit parameters
    par22[i] = par[i];
  }
  int index = (int)par[22];
  return IntensityCombo_11(x, par22, index);
}

double IntensityCombo_10(double *x, double *par, int i) {//ialpha = 0, ibeta = 2
  int ialpha = 2;
  int ibeta = 0;
  return 2.0 * SingleIntensity(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3], ialpha, ibeta);
}

double IntensityCombo_10plot(double *x, double *par) {
 double par22[22];
  for (int i = 0; i < 22; i++) { //loop over 22 fit parameters
    par22[i] = par[i];
  }
  int index = (int)par[22];
  return IntensityCombo_10(x, par22, index);
}

double Intensity(double *x, double *par, int i) {
  if (i < 15) {
    return IntensityPlus(x, par, lmLM[i][0], lmLM[i][1], lmLM[i][2], lmLM[i][3]);
  }
  else {
    return IntensityMinus(x, par, i);
  }
}

double IntensityPlot(double *x, double *par) {
  double par22[22];
  for (int i = 0; i < 22; i++) { //loop over 22 fit parameters
    par22[i] = par[i];
  }
  int index = (int)par[22];
  return Intensity(x, par22, index); 
}

void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/, Double_t &fval, Double_t *par, Int_t /*iflag*/) {
  double chi2 = 0;
  double tmp,x;
  //loop over 25 moments
  for (int ih = 0; ih < 25; ih++) {
    int n = h[ih]->GetNbinsX();
    //n = n*3/5; //uncomment this to fit over 0.9-1.5 GeV only
    //loop over n bins
    for (int in = 1; in < n; in++) {
      x = h[ih]->GetXaxis()->GetBinCenter(in);
      tmp = (h[ih]->GetBinContent(in) - Intensity(&x, par, ih))/(h[ih]->GetBinError(in));
      chi2 += tmp*tmp;
    }  
  }
  fval = chi2;
}



void BreitWigner(bool performfit = true) { 
  TFile *fpomegapi = TFile::Open("/sciclone/data10/amschertz/pomegapi/hist_30ktotal.root");
  TString hname[25] = {"H0000", "H0020", "H0021", "H0022", "H2000", "H2020", "H2021", "H2022", "H2120", "H2121_plus", "H2122_plus", "H2220", "H2221_plus", "H2222_plus", "H2111_plus", "H0010", "H0011", "H2110", "H2111_minus", "H2121_minus", "H2122_minus", "H2221_minus", "H2222_minus", "H2010", "H2011"};
  TString parname[22] = {"1+ resonance mass", "1+ resonance width", "1+ normalization", "1- resonance mass", "1- resonance width", "1- normalization", "0- resonance mass", "0- resonance width", "0- normalization", "D/S ratio", "1+ phi0", "1+ theta", "1+ phi+", "1+ phi-", "1+ psi", "1- phi0", "1- theta", "1- phi+", "1- phi-", "1- psi", "0- phi0", "0- theta"};
  double parinitial[22] = {1.217, 0.231, 1.0, 1.2295, 0.142, 1.0, 1.2295, 0.142, 1.0, 0.26, 3.14, 0.78, 3.14, 3.14, 0.78, 3.14, 0.78, 3.14, 3.14, 0.78, 3.14, 0.78};
  double parmax[22] = {0, 0.5, 0, 0, 10000, 0, 0, 10000, 0, 0, 6.28, 1.57, 6.28, 6.28, 1.57, 6.28, 1.57, 6.28, 6.28, 1.57, 6.28, 1.57};
  TF1 *fit[25];

  for (int i = 0; i < 25; i++) { //loop over 25 moments
    h[i] = (TH1F*)fpomegapi->Get(hname[i]);
    TString ts = hname[i];
    ts += "_fit";
    fit[i] = new TF1(ts, IntensityPlot, 0.9, 1.9, 23); 
  }
  
  for (int i = 0; i < 25; i++) { //loop over 25 moments
    for (int j = 0; j < 22; j++) { //loop over 22 parameters
      fit[i]->SetParameter(j, parinitial[j]);
    }
    fit[i]->FixParameter(22, i); //index which histogram and fit to use
    
    //if (performfit)
      //h[i]->Fit(fit[i], "", "N");
  }

  
  TF1 *plot_[3][15]; //3->J^P, 15->lmLM
  TString jpc[3] = {"1plus", "1minus", "0minus"};

  for (int i = 0; i < 15; i++) { //loop over the 15 H+(lmLM) moments
    TString ts0 = "plot_";
    ts0 += jpc[0];
    ts0 += hname[i];
    plot_[0][i] = new TF1(ts0, OneIntensityPlot_1plus, 0.9, 1.9, 23);
    TString ts1 = "plot_";
    ts1 += jpc[1];
    ts1 += hname[i];
    plot_[1][i] = new TF1(ts1, OneIntensityPlot_1minus, 0.9, 1.9, 23);
    TString ts2 = "plot_";
    ts2 += jpc[2];
    ts2 += hname[i];
    plot_[2][i] = new TF1(ts2, OneIntensityPlot_0minus, 0.9, 1.9, 23);
    for (int j = 0; j < 3; j++) { //loop over 3 J^P states
      for (int p = 0; p < 22; p++) { //loop over 22 parameters
	plot_[j][i]->SetParameter(p, parinitial[p]);
      }
      plot_[j][i]->FixParameter(22, i); //index which histogram and fit to use
    }
  }
  
  TF1 *comboplot11[10];
  for (int i = 0; i < 10; i++) {
    TString ts11 = "1p1m_";
    ts11 += hname[i + 15];
    cout << "Comboplot name = " << ts11 << endl;
    comboplot11[i] = new TF1(ts11, IntensityCombo_11plot, 0.9, 1.9, 23);
    for (int p = 0; p < 22; p++) {
      comboplot11[i]->SetParameter(p, parinitial[p]);
    }
    comboplot11[i]->FixParameter(22, i+15);
  }

  TF1 *comboplot10[10];
  for (int i = 0; i < 10; i++) {
    TString ts10 = "1p0m_";
    ts10 += hname[i + 15];
    comboplot10[i] = new TF1(ts10, IntensityCombo_10plot, 0.9, 1.9, 23);
    for (int p = 0; p < 22; p++) {
      comboplot10[i]->SetParameter(p, parinitial[p]);
    }
    comboplot10[i]->FixParameter(22, i+15);
  }
  
  
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, 22);

  for (int j = 0; j < 22; j++) { //loop over 22 parameters
    minuit->SetParameter(j, parname[j], parinitial[j], 0.01, 0, parmax[j]);
  }
  minuit->FixParameter(0);  //fix 1+ resonance mass at b1 mass
  minuit->FixParameter(1);  //fix 1+ resonance width at b1 width
  //minuit->FixParameter(3);  //fix 1- resonance mass at rho0 mass
  //minuit->FixParameter(5);  //fix 1- normalization to zero to eliminate all but b1 state
  //minuit->FixParameter(8);  //fix 0- normalization to zero to eliminate 0- state
  //minuit->FixParameter(9);  //fix D/S ratio at 0.26
  minuit->FixParameter(10); //fix 1+ phi0 at some arbitrary value

  minuit->SetFCN(myFcn);
  
  if (!performfit)
    return;

  double arglist[100];
  arglist[0] = 0;
  //set print level
  minuit->ExecuteCommand("SET PRINT",arglist,2);

  //minimize
  arglist[0] = 5000000; // number of function calls
  arglist[1] = 0.0001; // tolerance
  minuit->ExecuteCommand("MIGRAD",arglist,2);

  //get result
  double minParams[22];
  double parErrors[22];
  for (int i = 0; i < 22; ++i) { //loop over 22 parameters
    minParams[i] = minuit->GetParameter(i);
    parErrors[i] = minuit->GetParError(i);
  }
  double chi2, edm, errdef;
  int nvpar, nparx;
  minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
 
  
  for (int i = 0; i < 25; i++) { //loop over 25 moments
    for (int p = 0; p < 22; p++) { //loop over 22 fit parameters
      fit[i]->SetParameter(p, minParams[p]);
    }
    fit[i]->SetChisquare(chi2);
    fit[i]->SetLineColor(kRed);
  } 

 
  for (int i = 0; i < 15; i++) {
    for (int j = 0; j < 3; j++) {
      for (int p = 0; p < 22; p++) { //loop over 22 fit parameters
	plot_[j][i]->SetParameter(p, minParams[p]);
      }
    }
    plot_[0][i]->SetLineColor(kBlue);
    plot_[1][i]->SetLineColor(kGreen);
    plot_[2][i]->SetLineColor(kBlack);
  }

  
  for (int i = 0; i < 10; i++) { //loop over the 10 interference moments
    for (int p; p < 22; p++) { //loop over the 22 fit parameters
      comboplot11[i]->SetParameter(p, minParams[p]);
      comboplot10[i]->SetParameter(p, minParams[p]);
    }
    comboplot11[i]->SetLineColor(kMagenta);
    comboplot10[i]->SetLineColor(kCyan);
  }

  //Set variables for fit parameters
  double resmass[3] = {minuit->GetParameter(0), minuit->GetParameter(3), minuit->GetParameter(6)};
  double reswidth[3] = {minuit->GetParameter(1), minuit->GetParameter(4), minuit->GetParameter(7)};
  double resnorm[3] = {minuit->GetParameter(2), minuit->GetParameter(5), minuit->GetParameter(8)};
  double DSratio = minuit->GetParameter(9);
  double phi0[3] = {minuit->GetParameter(10), minuit->GetParameter(15), minuit->GetParameter(20)};
  double theta[3] = {minuit->GetParameter(11), minuit->GetParameter(16), minuit->GetParameter(21)};
  double phiplus[3] = {minuit->GetParameter(12), minuit->GetParameter(17), 0};
  double phiminus[3] = {minuit->GetParameter(13), minuit->GetParameter(18), 0}; 
  double psi[3] = {minuit->GetParameter(14), minuit->GetParameter(19), 0}; 

 

  //Output values of the production density matrix to check vs Atkinson Table 4
  //rho_00(1+)
  cout << "rho_00(1+) = " << rho(phi0[0], theta[0], phiplus[0], phiminus[0], psi[0], phi0[0], theta[0], phiplus[0], phiminus[0], psi[0], 0, 0, 0, 0) << endl;
  //rho_+0(1+)
  cout << "rho_+0(1+) = " << rho(phi0[0], theta[0], phiplus[0], phiminus[0], psi[0], phi0[0], theta[0], phiplus[0], phiminus[0], psi[0], 1, 0, 0, 0) << endl;
  //rho_+-(1+)
  cout << "rho_+-(1+) = " << rho(phi0[0], theta[0], phiplus[0], phiminus[0], psi[0], phi0[0], theta[0], phiplus[0], phiminus[0], psi[0], 1, 2, 0, 0) << endl;
  //rho_00(1-)
  cout << "rho_00(1-) = " << rho(phi0[1], theta[1], phiplus[1], phiminus[1], psi[1], phi0[1], theta[1], phiplus[1], phiminus[1], psi[1], 0, 0, 1, 1) << endl;
  //rho_+0(1-)
  cout << "rho_+0(1-) = " << rho(phi0[1], theta[1], phiplus[1], phiminus[1], psi[1], phi0[1], theta[1], phiplus[1], phiminus[1], psi[1], 1, 0, 1, 1) << endl;
  //rho_+-(1-)
  cout << "rho_+-(1-) = " << rho(phi0[1], theta[1], phiplus[1], phiminus[1], psi[1], phi0[1], theta[1], phiplus[1], phiminus[1], psi[1], 1, 2, 1, 1)  << endl;
  
  
  TCanvas *c15 = new TCanvas("c15", "Nine Fits", 900, 900);
  c15->Divide(3, 3);

  /*
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(h[0], "GlueX Data");
  legend->AddEntry(fit[0], "Fit (3 states)");
  legend->AddEntry(plot_[0][0], "1+ Contribution");
  legend->AddEntry(plot_[1][0], "1- Contribution");
  legend->AddEntry(plot_[2][0], "0- Contribution");
  */

  c15->cd(1);
  h[0]->Draw();
  fit[0]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][0]->Draw("same");
  }
  //legend->Draw();

  c15->cd(2);
  h[1]->Draw();
  fit[1]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][1]->Draw("same");
  }

  c15->cd(3);
  h[4]->Draw();
  fit[4]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][4]->Draw("same");
  }

  c15->cd(4);
  h[5]->Draw();
  fit[5]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][5]->Draw("same");
  }

  c15->cd(5);
  h[8]->Draw();
  fit[8]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][8]->Draw("same");
  }

  c15->cd(6);
  h[11]->Draw();
  fit[11]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][11]->Draw("same");
  }

  c15->cd(7);
  h[3]->Draw();
  fit[3]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][3]->Draw("same");
  }

  c15->cd(8);
  h[7]->Draw();
  fit[7]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][7]->Draw("same");
  }

  c15->cd(9);
  h[13]->Draw();
  fit[13]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][13]->Draw("same");
  }

  TCanvas *c16 = new TCanvas("c16", "Nine More Fits", 900, 900);
  c16->Divide(3, 3);

  c16->cd(1);
  h[10]->Draw();
  fit[10]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][10]->Draw("same");
  }

  c16->cd(2);
  h[2]->Draw();
  fit[2]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][2]->Draw("same");
  }

  c16->cd(3);
  h[6]->Draw();
  fit[6]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][6]->Draw("same");
  }

  c16->cd(4);
  h[12]->Draw();
  fit[12]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][12]->Draw("same");
  }

  c16->cd(5);
  h[9]->Draw();
  fit[9]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][9]->Draw("same");
  }

  c16->cd(6);
  h[14]->Draw();
  fit[14]->Draw("same");
  for (int j = 0; j < 3; j++){
    plot_[j][14]->Draw("same");
  }
  
  c16->cd(7);
  h[15]->Draw();
  fit[15]->Draw("same");
  comboplot11[0]->Draw("same");
  comboplot10[0]->Draw("same");

  c16->cd(8);
  h[16]->Draw();
  fit[16]->Draw("same");
  comboplot11[1]->Draw("same");
  comboplot10[1]->Draw("same");

  c16->cd(9);
  h[17]->Draw();
  fit[17]->Draw("same");
  comboplot11[2]->Draw("same");
  comboplot10[2]->Draw("same");

  TCanvas *c17 = new TCanvas("c17", "Last Seven Fits", 900, 900);
  c17->Divide(3, 3);

  for (int i = 0; i < 7; i++) {
    c17->cd(i + 1);
    h[i + 18]->Draw();
    fit[i + 18]->Draw("same");
    comboplot11[i + 3]->Draw("same");
    comboplot10[i + 3]->Draw("same");
  }

  return;
}
