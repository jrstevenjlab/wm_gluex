#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"

TH2F* locHist_MMOP_phi2D[7];
TH2F* locHist_MMOP_theta2D[7];
TH2F* locHist_MMOP_p2D[7];
TH2F* locHist_MMOP_phi2D_0track;
TH2F* locHist_MMOP_theta2D_0track;
TH2F* locHist_MMOP_p2D_0track;

TH1F* locHist_MMOP_phi;
TH1F* locHist_MMOP_theta;
TH1F* locHist_MMOP_p;
TH1F* locHist_MMOP_phi_0track;
TH1F* locHist_MMOP_theta_0track;
TH1F* locHist_MMOP_p_0track;

double fitfn(double *x, double *par) { // two Gaussians + cubic background
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5], 2.0)) + par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0];
}

double fit1gaus(double *x, double *par) { // single Gaussian + cubic background
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
}

double fit2gaus(double *x, double *par) { //two gaussians with the same centroid, cubic polynomial background
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[4], 2.0)) + par[5] + par[6]*x[0] + par[7]*x[0]*x[0] + par[8]*x[0]*x[0]*x[0];
}

double omega2gaus(double *x, double *par) { //two gaussians with the same centroid, no background
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[4], 2.0));
}

double sgf1(double *x, double *par) {  // Skewed Gaussian + cubic background
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2)) + par[4] + par[5]*x[0] + par[6]*x[0]*x[0] + par[7]*x[0]*x[0]*x[0];
}

double sgf1plot(double *x, double *par) {  // Skewed Gaussian only
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2));
}

double sgf1quad(double *x, double *par) {  // Skewed Gaussian + quadratic background
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2)) + par[4] + par[5]*x[0] + par[6]*x[0]*x[0];
}

double sgf2cubeOld(double *x, double *par) { // Two skewed gaussians (omega and phi peaks) + cubic background
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2)) + par[4]*exp(-0.5*TMath::Power((x[0]-par[5])/(par[6]+(x[0]<par[5])*par[7]*(x[0]-par[5])),2)) + par[8] + par[9]*x[0] + par[10]*x[0]*x[0] + par[11]*x[0]*x[0]*x[0];
}

double sgf2cube(double *x, double *par) {  // Two skewed gaussians w/ same width and skew + cubic background
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2)) + par[4]*exp(-0.5*TMath::Power((x[0]-par[5])/(par[2]+(x[0]<par[5])*par[3]*(x[0]-par[5])),2)) + par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0];
}

double sgf2quad(double *x, double *par) {  // Two skewed gaussians w/ same width and skew + quadratic background
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2)) + par[4]*exp(-0.5*TMath::Power((x[0]-par[5])/(par[2]+(x[0]<par[5])*par[3]*(x[0]-par[5])),2)) + par[6] + par[7]*x[0] + par[8]*x[0]*x[0];
}

////// Assorted limits and such that are used in a lot of places //////

// Fit limits: limits used when fitting the histograms and/or plotting the fits
Double_t fitMin = 0.3; // GeV
Double_t fitMaxIfBG = 1.15; // GeV; use when fitting/plotting background or total fit
Double_t fitMaxFull = 1.3; // GeV; use when fitting/plotting gaussian peak(s) only (Full = full range)

// Parameter limits: limits on fit parameters
Double_t omegaMeanMin = 0.7; // lower limit on omega dist. mean
Double_t omegaMeanMax = 0.87; // upper limit on omega dist. mean
Double_t omegaConstMax = 1000000; // upper limit on omega dist. amplitude
Double_t phiMeanMin = 0.95; // lower limit on phi (meson) dist. mean
Double_t phiMeanMax = 1.05; // upper limit on phi dist. mean
Double_t phiConstMax = 10000; // upper limit on phi dist. amplitude
Double_t sigmaMin = 0.01; // lower limit on gaussian widths
Double_t sigmaMax = 0.07; // upper limit on gaussian widths
Double_t skewMin = -0.25; // lower limit on gaussian skews
Double_t skewInit = -0.15; // initial skewness value (when setting parameters; this isn't really a limit)


void fit(TString sample = "bggen_2017_ver03_misspim", bool save = true, const char *whichcut = "pv") {
  TFile* fhist;
  TFile* outROOT;

  fhist = TFile::Open(Form("/sciclone/gluex10/jrstevens01/omega_misspi/hist_%s.root", sample.Data());
  if(save == true)    
    outROOT = new TFile(Form("mmop_%s_params_%scut.root", sample.Data(), whichcut), "recreate");

  TString plotDir = sample += "/fits_paramTest";
  TString plotData = sample;

  for(int im = 0; im < 7; im++){ //loop over cuts on (P(chi^2) or theta, or 3-momentum)
    locHist_MMOP_phi2D[im] = (TH2F*)fhist->Get(Form("OmegaMassVsPhi_1track_%scut%d", whichcut, im+1));
    locHist_MMOP_phi2D[im]->Sumw2();
    locHist_MMOP_phi2D[im]->RebinX(10);
    locHist_MMOP_phi2D[im]->RebinY(30); // turns 600 bins into 20 bins

    locHist_MMOP_theta2D[im] = (TH2F*)fhist->Get(Form("OmegaMassVsTheta_1track_%scut%d", whichcut, im+1));
    locHist_MMOP_theta2D[im]->Sumw2();
    locHist_MMOP_theta2D[im]->RebinX(10);
    locHist_MMOP_theta2D[im]->RebinY(30);

    locHist_MMOP_p2D[im] = (TH2F*)fhist->Get(Form("OmegaMassVsP_1track_%scut%d", whichcut, im+1));
    locHist_MMOP_p2D[im]->Sumw2();
    locHist_MMOP_p2D[im]->RebinX(10);
    locHist_MMOP_p2D[im]->RebinY(30);
  }

  locHist_MMOP_phi2D_0track = (TH2F*)fhist->Get("OmegaMassVsPhi_0tracks");
  locHist_MMOP_phi2D_0track->Sumw2();
  locHist_MMOP_phi2D_0track->RebinX(10);
  locHist_MMOP_phi2D_0track->RebinY(30);

  locHist_MMOP_theta2D_0track = (TH2F*)fhist->Get("OmegaMassVsTheta_0tracks");
  locHist_MMOP_theta2D_0track->Sumw2();
  locHist_MMOP_theta2D_0track->RebinX(10);
  locHist_MMOP_theta2D_0track->RebinY(30);

  locHist_MMOP_p2D_0track = (TH2F*)fhist->Get("OmegaMassVsP_0tracks");
  locHist_MMOP_p2D_0track->Sumw2();
  locHist_MMOP_p2D_0track->RebinX(10);
  locHist_MMOP_p2D_0track->RebinY(30);

  double binsize = 1./60.;

  TCanvas *cphi[7];
  TCanvas *ctheta[7];
  TCanvas *cp[7];

  for(int im = 0; im < 7; im++){
    cphi[im] = new TCanvas(Form("cphi_%d", im+1), "Phi fits (1 candidate track)", 1200, 1200);
    cphi[im]->Divide(5,4);
    ctheta[im] = new TCanvas(Form("ctheta_%d", im+1), "Theta fits (1 candidate track)", 1200, 1200);
    ctheta[im]->Divide(5,4);
    cp[im] = new TCanvas(Form("cp_%d", im+1), "P fits (1 candidate track)", 1200, 1200);
    cp[im]->Divide(5,4);
  }

  TCanvas *cphi_0track = new TCanvas("cphi_0track", "Phi fits (0 candidate tracks)", 1200, 1200);
  cphi_0track->Divide(5,4);

  TCanvas *ctheta_0track = new TCanvas("ctheta_0track", "Theta fits (0 candidate tracks)", 1200, 1200);
  ctheta_0track->Divide(5,4);

  TCanvas *cp_0track = new TCanvas("cp_0track", "P fits (0 candidate tracks)", 1200, 1200);
  cp_0track->Divide(5,4);

  gStyle->SetOptStat(000000);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(1);

  TH1F* locHist_phi_reco[7];
  TH1F* locHist_theta_reco[7];
  TH1F* locHist_p_reco[7];
  //TH1F* locHist_phi_denom[7];
  //TH1F* locHist_theta_denom[7];
  //TH1F* locHist_p_denom[7];
  TH1F* locHist_phi_missing[7];
  TH1F* locHist_theta_missing[7];
  TH1F* locHist_p_missing[7];

  TH1F* locHist_phi_params[7][11]; // save all fit params, 11th one is chi squared
  TH1F* locHist_theta_params[7][11];
  TH1F* locHist_p_params[7][11];
  TH1F* locHist_phi_miss_params[7][11];
  TH1F* locHist_theta_miss_params[7][11];
  TH1F* locHist_p_miss_params[7][11];

  TString paramNames[11] = {"omega_const", "omega_mean", "sigma", "skew", "phi_const", "phi_mean", "constant", "linear", "quadratic", "cubic", "chiSquared"};

  for(int im = 0; im < 7; im++){
    locHist_phi_reco[im] = new TH1F(Form("phi_mmop_reco_%d", im+1), "#omega_{mmop} (1 #pi);#phi_{mmop} (deg);#omega yield / 18#circ", 20, -180., 180.);
    locHist_phi_reco[im]->Sumw2();
    locHist_theta_reco[im] = new TH1F(Form("theta_mmop_reco_%d", im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
    locHist_theta_reco[im]->Sumw2();
    locHist_p_reco[im] = new TH1F(Form("p_mmop_reco_%d", im+1), ";p_{mmop} (GeV);#omega yield / 0.5 GeV", 12, 0., 6.);
    locHist_p_reco[im]->Sumw2();
    //locHist_phi_denom[im] = new TH1F(Form("phi_mmop_denom_%d", im+1), "#omega_{mmop} (1 or 0 #pi);#phi_{mmop} (deg);#omega yield / 18#circ", 20, -180., 180.);
    //locHist_phi_denom[im]->Sumw2();
    //locHist_theta_denom[im] = new TH1F(Form("theta_mmop_denom_%d", im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
    //locHist_theta_denom[im]->Sumw2();
    //locHist_p_denom[im] = new TH1F(Form("p_mmop_denom_%d", im+1), ";p_{mmop} (GeV);#omega yield / 0.5 GeV", 12, 0., 6.);
    //locHist_p_denom[im]->Sumw2();
    locHist_phi_missing[im] = new TH1F(Form("phi_mmop_missing_%d", im+1), "#omega_{mmop} (0 #pi);#phi_{mmop} (deg);#omega yield / 18#circ", 20, -180., 180.);
    locHist_phi_missing[im]->Sumw2();
    locHist_theta_missing[im] = new TH1F(Form("theta_mmop_missing_%d", im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
    locHist_theta_missing[im]->Sumw2();
    locHist_p_missing[im] = new TH1F(Form("p_mmop_missing_%d", im+1), ";p_{mmop} (GeV);#omega yield / 0.5 GeV", 12, 0., 6.);
    locHist_p_missing[im]->Sumw2();

    for(int i = 0; i < 11; i++){
	locHist_phi_params[im][i] = new TH1F(Form("phi_mmop_param_%d_pv%d", i+1, im+1), Form("fit params (mmop, 1 #pi);#phi_{mmop} (deg);%s", paramNames[i].Data()), 20, -180., 180.);
	locHist_theta_params[im][i] = new TH1F(Form("theta_mmop_param_%d_pv%d", i+1, im+1), Form(";#theta_{mmop} (deg);%s", paramNames[i].Data()), 20, 0., 30.);
	locHist_p_params[im][i] = new TH1F(Form("p_mmop_param_%d_pv%d", i+1, im+1), Form(";#p_{mmop} (GeV);%s", paramNames[i].Data()), 12, 0., 6.);
	locHist_phi_miss_params[im][i] = new TH1F(Form("phi_mmop_miss_param_%d_pv%d", i+1, im+1), Form("fit params (mmop, 0 #pi);#phi_{mmop} (deg);%s", paramNames[i].Data()), 20, -180., 180.);
	locHist_theta_miss_params[im][i] = new TH1F(Form("theta_mmop_miss_param_%d_pv%d", i+1, im+1), Form(";#theta_{mmop} (deg);%s", paramNames[i].Data()), 20, 0., 30.);
	locHist_p_miss_params[im][i] = new TH1F(Form("p_mmop_miss_param_%d_pv%d", i+1, im+1), Form(";#p_{mmop} (GeV);%s", paramNames[i].Data()), 12, 0., 6.);
    }
  }

  double phi_errorD[20];
  double theta_errorD[20];
  double p_errorD[12];
  double phi_errorU[20];
  double theta_errorU[20];
  double p_errorU[12];

  TCanvas *cphi_1;
  cphi_1 = new TCanvas("cphi_single", "Phi fit (1 candidate track)", 600, 600);

  TCanvas *ctheta_1;
  ctheta_1 = new TCanvas("ctheta_single", "Theta fit (1 candidate track)", 600, 600);

  cout << "Phi fits" << endl;
  for(int im = 0; im < 7; im++){//loop over mm^2 cuts
    for(int i = 0; i < 20; i++){
      cphi[im]->cd(i+1);
      //if(im == 0 && i == 10){
	//cphi_1->cd();
      //}
      locHist_MMOP_phi = (TH1F*)locHist_MMOP_phi2D[im]->ProjectionX(Form("phi_%d", i), i+1, i+1);

      double phimin = -180. + 18.*i;
      double phimax = -162. + 18.*i;

      TString ts;
      ts += phimin;
      ts += "#circ < #phi < ";
      ts += phimax;
      ts += "#circ;Missing Mass off the Proton (GeV)";

      locHist_MMOP_phi->SetTitle(ts);

      // Preliminary fit function to phi (meson) peak
      TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
      gausphi->SetParameters(100., 1., 0.1);
      gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);
    
      TFitResultPtr r1 = locHist_MMOP_phi->Fit("gaus", "SQN", "", 0.7, 0.9); // rough omega fit
      TFitResultPtr rphi = locHist_MMOP_phi->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax); // rough phi fit
      //TF1 *sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
      TF1 *sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);

      //sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
      sgf->SetParameter(0, r1->Parameter(0));
      sgf->SetParameter(1, r1->Parameter(1));
      sgf->SetParameter(2, r1->Parameter(2));
      sgf->SetParameter(3, skewInit);
      sgf->SetParameter(4, rphi->Parameter(0));
      sgf->SetParameter(5, rphi->Parameter(1));
      //sgf->SetParameter(6, rphi->Parameter(2));
      //sgf->SetParameter(7, -0.15);
      sgf->SetParameter(6, 1.);
      sgf->SetParameter(7, 1.);
      sgf->SetParameter(8, 1.);
      sgf->SetParameter(9, 1.);

      sgf->SetParLimits(0, 0, omegaConstMax); // omega "height"
      sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean
      sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
      sgf->SetParLimits(3, skewMin, 0.); // skew
      sgf->SetParLimits(4, 0, phiConstMax);  // phi "height"
      sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
      //sgf->SetParLimits(6, 0.01, 0.08); // sigma
      //sgf->SetParLimits(7, -0.5, 0.); // skew
      
      sgf->SetParName(0, "#omega constant");
      sgf->SetParName(1, "#omega mean");
      sgf->SetParName(2, "#sigma");
      sgf->SetParName(3, "Skew");
      sgf->SetParName(4, "#phi constant");
      sgf->SetParName(5, "#phi mean");
      //sgf->SetParName(6, "#sigma_{#phi}");
      //sgf->SetParName(7, "#phi skew");
      sgf->SetParName(6, "Constant");
      sgf->SetParName(7, "Linear");
      sgf->SetParName(8, "Quadratic");
      sgf->SetParName(9, "Cubic");

      //TFitResultPtr r = locHist_MMOP_phi->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print parameters
      TFitResultPtr r = locHist_MMOP_phi->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram, don't print parameters
      Double_t chiSquare;
      chiSquare = locHist_MMOP_phi->Chisquare(sgf, "R");

      TMatrixD covMatrix(4,4);
 
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<4; i++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};
      
      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(fitMin, fitMaxFull);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
      error = error/binsize;

      TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
      sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));

      TF1 *polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
      polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
    
      sgfplot->SetLineColor(kBlue);
      sgfPhi->SetLineColor(kGreen);
      polBG->SetLineColor(kViolet);

      // Change line width for smaller plots
      sgfplot->SetLineWidth(2);
      sgfPhi->SetLineWidth(2);
      polBG->SetLineWidth(2);
      //locHist_MMOP_phi->GetFunction("sgf")->SetLineWidth(2);
      sgf->SetLineWidth(2);
    
      locHist_MMOP_phi->Draw();
      polBG->Draw("same");
      sgfPhi->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");

      locHist_phi_reco[im]->SetBinContent(i+1, intomega);
      locHist_phi_reco[im]->SetBinError(i+1, error);
      //locHist_phi_denom[im]->SetBinContent(i+1, intomega);

      for(int ipar = 0; ipar < 10; ipar++){
	locHist_phi_params[im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
      }
      locHist_phi_params[im][10]->SetBinContent(i+1, chiSquare);

      //if (im == 0 && i == 10){
	//cphi_1->Print("plots_S2018_bggen_new/fits_bggen/MMOP_bggen_PhiFit_1_eg_phiPeak1_S2018_new.pdf");
      //}
    }
    if(im == 0){
      cphi[im]->Print(Form("plots_%s/MMOP_%s_PhiFits_%d.pdf", plotDir.Data(), plotData.Data(), im+1));
    }
  }

  //fit phi (missing)
  for(int i = 0; i < 20; i++){
    cphi_0track->cd(i+1);
    locHist_MMOP_phi_0track = (TH1F*)locHist_MMOP_phi2D_0track->ProjectionX(Form("phi_0track_%d", i), i+1, i+1);
    
    double phimin = -180. + 18.*i;
    double phimax = -162. + 18.*i;
    
    TString ts;
    ts += phimin;
    ts += "#circ < #phi < ";
    ts += phimax;
    ts += "#circ;Missing Mass off the Proton (GeV)";
     
    locHist_MMOP_phi_0track->SetTitle(ts);

    // Preliminary fit function to phi (meson) peak
    TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
    gausphi->SetParameters(100., 1., 0.1);
    gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);
    
    TFitResultPtr r1 = locHist_MMOP_phi_0track->Fit("gaus", "SQN", "", 0.7, 0.9); // rough omega fit
    TFitResultPtr rphi = locHist_MMOP_phi_0track->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax); // rough phi fit
    //TF1 *sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
    TF1 *sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);

    sgf->SetParameter(0, r1->Parameter(0));
    sgf->SetParameter(1, r1->Parameter(1));
    sgf->SetParameter(2, r1->Parameter(2));
    sgf->SetParameter(3, skewInit);
    sgf->SetParameter(4, rphi->Parameter(0));
    sgf->SetParameter(5, rphi->Parameter(1));
    sgf->SetParameter(6, 1.);
    sgf->SetParameter(7, 1.);
    sgf->SetParameter(8, 1.);
    sgf->SetParameter(9, 1.);

    //sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), skewInit, 1., 1., 1., 1.);
    sgf->SetParLimits(0, 0, omegaConstMax); // omega "height"
    sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean
    sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
    sgf->SetParLimits(3, skewMin, 0.); // skew
    sgf->SetParLimits(4, 0, phiConstMax);  // phi "height"
    sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
      
    sgf->SetParName(0, "#omega constant");
    sgf->SetParName(1, "#omega mean");
    sgf->SetParName(2, "#sigma");
    sgf->SetParName(3, "Skew");
    sgf->SetParName(4, "#phi constant");
    sgf->SetParName(5, "#phi mean");
    sgf->SetParName(6, "Constant");
    sgf->SetParName(7, "Linear");
    sgf->SetParName(8, "Quadratic");
    sgf->SetParName(9, "Cubic");

    TFitResultPtr r = locHist_MMOP_phi_0track->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram or print stat box
    //TFitResultPtr r = locHist_MMOP_phi_0track->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print stat box
    Double_t chiSquare;
    chiSquare = locHist_MMOP_phi_0track->Chisquare(sgf, "R");

    TMatrixD covMatrix(4,4);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<4; i++){
      for (Int_t j=0; j<4; j++){
	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

    TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
    sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
    double intomega = sgfplot->Integral(fitMin, fitMaxFull);
    intomega = intomega/binsize;
    double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
    sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));

    TF1 *polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
    polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
    
    sgfplot->SetLineColor(kBlue);
    sgfPhi->SetLineColor(kGreen);
    polBG->SetLineColor(kViolet);

    // Change line width for smaller plots
    sgfplot->SetLineWidth(2);
    sgfPhi->SetLineWidth(2);
    polBG->SetLineWidth(2);
    //locHist_MMOP_phi_0track->GetFunction("sgf")->SetLineWidth(2);
    sgf->SetLineWidth(2);
    
    locHist_MMOP_phi_0track->Draw();
    polBG->Draw("same");
    sgfPhi->Draw("same");
    sgfplot->Draw("same");
    sgf->Draw("same");
    
    cout << "Integral under omega peak (from fit) = " << intomega << endl;
    cout << "Error = " << error << endl;

    phi_errorU[i] = error;

    for(int im = 0; im < 7; im++){
      //locHist_phi_denom[im]->AddBinContent(i+1, intomega);
      //locHist_phi_denom[im]->SetBinError(i+1, error);
      locHist_phi_missing[im]->SetBinContent(i+1, intomega);
      locHist_phi_missing[im]->SetBinError(i+1, error);

      for(int ipar = 0; ipar < 10; ipar++){
	locHist_phi_miss_params[im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
      }
      locHist_phi_miss_params[im][10]->SetBinContent(i+1, chiSquare);
    }
  }
  cphi_0track->Print(Form("plots_%s/MMOP_%s_PhiFits_missing.pdf", plotDir.Data(), plotData.Data()));


  cout << "Theta fits:" << endl;
  //fit theta (reco)
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 20; i++){
      ctheta[im]->cd(i+1);
      //if(im == 0 && i == 0){
	//ctheta_1->cd();
      //}

      locHist_MMOP_theta = (TH1F*)locHist_MMOP_theta2D[im]->ProjectionX(Form("theta_%d", i), i+1, i+1);
      double hsum = locHist_MMOP_theta->Integral();
      if(hsum < 75.)
	continue;

      double thetamin = 0. + 1.5*i;
      double thetamax = 1.5 + 1.5*i;

      TString ts;
      ts += Form("%.1f", thetamin);
      ts += "#circ < #theta < ";
      ts += Form("%.1f", thetamax);
      ts += "#circ;Missing Mass off the Proton (GeV)";

      locHist_MMOP_theta->SetTitle(ts);

      // Preliminary fit to phi (meson) peak
      TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
      gausphi->SetParameters(100., 1., 0.1);
      gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);

      TFitResultPtr r1 = locHist_MMOP_theta->Fit("gaus", "SQN", "", 0.7, 0.9);
      TFitResultPtr rphi = locHist_MMOP_theta->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax);
      TF1 *sgf;
      //if(i < 10){
      	//sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
        sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);
        //sgf = new TF1("sgf", sgf2quad, fitMin, fitMaxIfBG, 9);

        //sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);

      	sgf->SetParameter(0, r1->Parameter(0));
      	sgf->SetParameter(1, r1->Parameter(1));
      	sgf->SetParameter(2, r1->Parameter(2));
      	sgf->SetParameter(3, skewInit);
      	sgf->SetParameter(4, rphi->Parameter(0));
      	sgf->SetParameter(5, rphi->Parameter(1));
      	//sgf->SetParameter(6, rphi->Parameter(2));
      	//sgf->SetParameter(7, -0.15);
      	sgf->SetParameter(6, 1.);
      	sgf->SetParameter(7, 1.);
      	sgf->SetParameter(8, 1.);
      	sgf->SetParameter(9, 1.);

      	sgf->SetParLimits(0, 0, omegaConstMax); // omega "height"
      	sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean
      	sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
      	sgf->SetParLimits(3, skewMin, 0.); // skew
        sgf->SetParLimits(4, 0, phiConstMax);  // phi "height"
        sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
        //sgf->SetParLimits(6, 0.01, 0.08);
        //sgf->SetParLimits(7, -0.5, 0.);
      
        sgf->SetParName(0, "#omega constant");
        sgf->SetParName(1, "#omega mean");
        sgf->SetParName(2, "#sigma");
        sgf->SetParName(3, "Skew");
        sgf->SetParName(4, "#phi constant");
      	sgf->SetParName(5, "#phi mean");
      	//sgf->SetParName(6, "#sigma_{#phi}");
      	//sgf->SetParName(7, "#phi skew");
      	sgf->SetParName(6, "Constant");
      	sgf->SetParName(7, "Linear");
      	sgf->SetParName(8, "Quadratic");
      	sgf->SetParName(9, "Cubic");

        //sgf->SetParName(0, "Constant");
      	//sgf->SetParName(1, "Mean");
      	//sgf->SetParName(2, "Sigma");
      	//sgf->SetParName(3, "Skew");
      	//sgf->SetParName(4, "Constant");
      	//sgf->SetParName(5, "Linear");
      	//sgf->SetParName(6, "Quadratic");
      	//sgf->SetParName(7, "Cubic");
      //}
      /*else{
	sgf = new TF1("sgf", sgf1quad, 0.3, fitMaxIfBG, 7);
      	sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1.);
      	sgf->SetParLimits(0, 0, 1000000);
      	sgf->SetParLimits(1, 0.7, 0.9);
      	sgf->SetParLimits(2, 0.01, 0.1);
      	sgf->SetParLimits(3, -0.25, 0.);
      
      	sgf->SetParName(0, "Constant");
      	sgf->SetParName(1, "Mean");
      	sgf->SetParName(2, "Sigma");
      	sgf->SetParName(3, "Skew");
      	sgf->SetParName(4, "Constant");
      	sgf->SetParName(5, "Linear");
      	sgf->SetParName(6, "Quadratic");
      }*/

      TFitResultPtr r = locHist_MMOP_theta->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram, don't print parameters
      //TFitResultPtr r = locHist_MMOP_theta->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print parameters
      Double_t chiSquare;
      chiSquare = locHist_MMOP_theta->Chisquare(sgf, "R");
      
      TMatrixD covMatrix(4,4);
      
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<4; i++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(fitMin, fitMaxFull);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
      error = error/binsize;

      TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
      sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));

      TF1* polBG;
      //if (i < 10){
      	polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
      	polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
      //}
      //else{
	//polBG = new TF1("polBG", "pol2", fitMin, fitMaxIfBG);
	//polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8));
      //}
    
      sgfplot->SetLineColor(kBlue);
      polBG->SetLineColor(kViolet);
      sgfPhi->SetLineColor(kGreen);

      // Change line width for smaller plots
      sgfplot->SetLineWidth(2);
      polBG->SetLineWidth(2);
      sgfPhi->SetLineWidth(2);
      //locHist_MMOP_theta->GetFunction("sgf")->SetLineWidth(2);
      sgf->SetLineWidth(2);
    
      locHist_MMOP_theta->Draw();
      polBG->Draw("same");
      sgfPhi->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");

      locHist_theta_reco[im]->SetBinContent(i+1, intomega);
      locHist_theta_reco[im]->SetBinError(i+1, error);
      //locHist_theta_denom[im]->SetBinContent(i+1, intomega);

      for(int ipar = 0; ipar < 10; ipar++){
	locHist_theta_params[im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
      }
      locHist_theta_params[im][10]->SetBinContent(i+1, chiSquare);
      //if (im == 0 && i == 0){
	//ctheta_1->Print("plots_S2018_bggen_new/fits_bggen/MMOP_bggen_ThetaFit_1_eg_phiPeak1_S2018_new.pdf");
      //}
    }
    if(im == 0){
      ctheta[im]->Print(Form("plots_%s/MMOP_%s_ThetaFits_%d.pdf", plotDir.Data(), plotData.Data(), im+1));
    }
  }

  cout << "Theta (miss) fits" << endl;
  //fit theta (missing)
  for(int i = 0; i < 20; i++){
    ctheta_0track->cd(i+1);
    locHist_MMOP_theta_0track = (TH1F*)locHist_MMOP_theta2D_0track->ProjectionX(Form("theta_0track_%d", i), i+1, i+1);
    double hsum = locHist_MMOP_theta_0track->Integral();
    if(hsum < 75.)
      continue;

    double thetamin = 0. + 1.5*i;
    double thetamax = 1.5 + 1.5*i;

    TString ts;
    ts += Form("%.1f", thetamin);
    ts += "#circ < #theta < ";
    ts += Form("%.1f", thetamax);
    ts += "#circ;Missing Mass off the Proton (GeV)";

    locHist_MMOP_theta_0track->SetTitle(ts);

    // Preliminary fit function to phi (meson) peak
    TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
    gausphi->SetParameters(100., 1., 0.1);
    gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);

    TFitResultPtr r1 = locHist_MMOP_theta_0track->Fit("gaus", "SQN", "", 0.7, 0.9); // rough omega fit
    TFitResultPtr rphi = locHist_MMOP_theta_0track->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax); // rough phi fit
    //TF1 *sgf;
    //sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
    TF1 *sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);

    sgf->SetParameter(0, r1->Parameter(0));
    sgf->SetParameter(1, r1->Parameter(1));
    sgf->SetParameter(2, r1->Parameter(2));
    sgf->SetParameter(3, skewInit);
    sgf->SetParameter(4, rphi->Parameter(0));
    sgf->SetParameter(5, rphi->Parameter(1));
    sgf->SetParameter(6, 1.);
    sgf->SetParameter(7, 1.);
    sgf->SetParameter(8, 1.);
    sgf->SetParameter(9, 1.);

    //sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), skewInit, 1., 1., 1., 1.);
    sgf->SetParLimits(0, 0, omegaConstMax); // omega "height"
    sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean
    sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
    sgf->SetParLimits(3, skewMin, 0.); // skew
    sgf->SetParLimits(4, 0, phiConstMax);  // phi "height"
    sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
      
    sgf->SetParName(0, "#omega constant");
    sgf->SetParName(1, "#omega mean");
    sgf->SetParName(2, "#sigma");
    sgf->SetParName(3, "Skew");
    sgf->SetParName(4, "#phi constant");
    sgf->SetParName(5, "#phi mean");
    sgf->SetParName(6, "Constant");
    sgf->SetParName(7, "Linear");
    sgf->SetParName(8, "Quadratic");
    sgf->SetParName(9, "Cubic");

    TFitResultPtr r = locHist_MMOP_theta_0track->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram or print stat box
    //TFitResultPtr r = locHist_MMOP_theta_0track->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print stat box
    Double_t chiSquare;
    chiSquare = locHist_MMOP_theta_0track->Chisquare(sgf, "R");

    TMatrixD covMatrix(4,4);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<4; i++){
      for (Int_t j=0; j<4; j++){
  	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

    TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
    sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
    double intomega = sgfplot->Integral(fitMin, fitMaxFull);
    intomega = intomega/binsize;
    double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
    sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));

    TF1* polBG;
    //if (i < 10){
      	polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
      	polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
    /*}
    else{
	polBG = new TF1("polBG", "pol2", fitMin, fitMaxIfBG);
	polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6));
    }*/
    
    sgfplot->SetLineColor(kBlue);
    sgfPhi->SetLineColor(kGreen);
    polBG->SetLineColor(kViolet);

    // Change line width for smaller plots
    sgfplot->SetLineWidth(2);
    sgfPhi->SetLineWidth(2);
    polBG->SetLineWidth(2);
    //locHist_MMOP_theta_0track->GetFunction("sgf")->SetLineWidth(2);
    sgf->SetLineWidth(2);
    
    locHist_MMOP_theta_0track->Draw();
    polBG->Draw("same");
    sgfPhi->Draw("same");
    sgfplot->Draw("same");
    sgf->Draw("same");
    
    //Time to check errors:
    double intHist = locHist_MMOP_theta_0track->Integral();

    cout << "Integral under omega peak (from fit) = " << intomega << endl;
    cout << "Error (from fit) = " << error << endl;
    cout << "Integral under omega peak (from hist) = " << intHist << endl;
    cout << "Error (from sqrt(N)) = " << TMath::Sqrt(intHist) << endl;

    for(int im = 0; im < 7; im++){
      //locHist_theta_denom[im]->AddBinContent(i+1, intomega);
      //locHist_theta_denom[im]->SetBinError(i+1, error);
      locHist_theta_missing[im]->SetBinContent(i+1, intomega);
      locHist_theta_missing[im]->SetBinError(i+1, error);

      for(int ipar = 0; ipar < 10; ipar++){
	locHist_theta_miss_params[im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
      }
      locHist_theta_miss_params[im][10]->SetBinContent(i+1, chiSquare);
    }
  }
  ctheta_0track->Print(Form("plots_%s/MMOP_%s_ThetaFits_missing.pdf", plotDir.Data(), plotData.Data()));


  cout << "p fits:" << endl;
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 12; i++){
      cp[im]->cd(i+1);
      locHist_MMOP_p = (TH1F*)locHist_MMOP_p2D[im]->ProjectionX(Form("p_%d", i), i+1, i+1);

      double pmin = 0. + 0.5*i;
      double pmax = 0.5 + 0.5*i;

      TString ts;
      ts += Form("%.1f", pmin);
      ts += " < p < ";
      ts += Form("%.1f", pmax);
      ts += " GeV;Missing Mass off the Proton (GeV)";

      locHist_MMOP_p->SetTitle(ts);

      // Preliminary fit to phi (meson) peak
      TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
      gausphi->SetParameters(100., 1., 0.1);
      gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);
      
      TFitResultPtr r1 = locHist_MMOP_p->Fit("gaus", "SQN", "", 0.7, 0.9);
      TFitResultPtr rphi = locHist_MMOP_p->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax);
      //TF1 *sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
      TF1 *sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);

      //sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);

      sgf->SetParameter(0, r1->Parameter(0));
      sgf->SetParameter(1, r1->Parameter(1));
      sgf->SetParameter(2, r1->Parameter(2));
      sgf->SetParameter(3, skewInit);
      sgf->SetParameter(4, rphi->Parameter(0));
      sgf->SetParameter(5, rphi->Parameter(1));
      //sgf->SetParameter(6, rphi->Parameter(2));
      //sgf->SetParameter(7, -0.15);
      sgf->SetParameter(6, 1.);
      sgf->SetParameter(7, 1.);
      sgf->SetParameter(8, 1.);
      sgf->SetParameter(9, 1.);

      sgf->SetParLimits(0, 0, omegaConstMax); // omega constant
      sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean 
      sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
      sgf->SetParLimits(3, skewMin, 0.); // skew
      sgf->SetParLimits(4, 0, phiConstMax);  // phi constant
      sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
      //sgf->SetParLimits(6, 0.01, 0.08);
      //sgf->SetParLimits(7, -0.5, 0.);
      
      sgf->SetParName(0, "#omega constant");
      sgf->SetParName(1, "#omega mean");
      sgf->SetParName(2, "#sigma");
      sgf->SetParName(3, "Skew");
      sgf->SetParName(4, "#phi constant");
      sgf->SetParName(5, "#phi mean");
      //sgf->SetParName(6, "#sigma_{#phi}");
      //sgf->SetParName(7, "#phi skew");
      sgf->SetParName(6, "Constant");
      sgf->SetParName(7, "Linear");
      sgf->SetParName(8, "Quadratic");
      sgf->SetParName(9, "Cubic");


      //sgf->SetParName(0, "Constant");
      //sgf->SetParName(1, "Mean");
      //sgf->SetParName(2, "Sigma");
      //sgf->SetParName(3, "Skew");
      //sgf->SetParName(4, "Constant");
      //sgf->SetParName(5, "Linear");
      //sgf->SetParName(6, "Quadratic");
      //sgf->SetParName(7, "Cubic");

      TFitResultPtr r = locHist_MMOP_p->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram, don't print parameters
      //TFitResultPtr r = locHist_MMOP_p->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print parameters
      Double_t chiSquare;
      chiSquare = locHist_MMOP_p->Chisquare(sgf, "R");

      TMatrixD covMatrix(4,4);
      
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<4; i++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(fitMin, fitMaxFull);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
      error = error/binsize;

      TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
      sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));

      TF1 *polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
      polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
      
      sgfplot->SetLineColor(kBlue);
      polBG->SetLineColor(kViolet);
      sgfPhi->SetLineColor(kGreen);

      // Change line width for smaller plots
      sgfplot->SetLineWidth(2);
      polBG->SetLineWidth(2);
      sgfPhi->SetLineWidth(2);
      sgf->SetLineWidth(2);
      //locHist_MMOP_p->GetFunction("sgf")->SetLineWidth(2);

    
      locHist_MMOP_p->Draw();
      polBG->Draw("same");
      sgfPhi->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");
    
      cout << "Integral under omega peak = " << intomega << endl;
      cout << "Error = " << error << endl;

      locHist_p_reco[im]->SetBinContent(i+1, intomega);
      locHist_p_reco[im]->SetBinError(i+1, error);
      //locHist_p_denom[im]->SetBinContent(i+1, intomega);

      for(int ipar = 0; ipar < 10; ipar++){
	locHist_p_params[im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
      }
      locHist_p_params[im][10]->SetBinContent(i+1, chiSquare);
    }
    if(im == 0){
      cp[im]->Print(Form("plots_%s/MMOP_%s_p_fits_%d.pdf", plotDir.Data(), plotData.Data(), im+1));
    }
  }

  //fit p (missing)
  for(int i = 0; i < 12; i++){
    cp_0track->cd(i+1);
    locHist_MMOP_p_0track = (TH1F*)locHist_MMOP_p2D_0track->ProjectionX(Form("p_0track_%d", i), i+1, i+1);

    double pmin = 0. + 0.5*i;
    double pmax = 0.5 + 0.5*i;

    TString ts;
    ts += Form("%.1f", pmin);
    ts += " < p < ";
    ts += Form("%.1f", pmax);
    ts += " GeV;Missing Mass off the Proton (GeV)";

    locHist_MMOP_p_0track->SetTitle(ts);

    // Preliminary fit function to phi (meson) peak
    TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
    gausphi->SetParameters(100., 1., 0.1);
    gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);
    
    TFitResultPtr r1 = locHist_MMOP_p_0track->Fit("gaus", "SQN", "", 0.7, 0.9); // rough omega fit
    TFitResultPtr rphi = locHist_MMOP_p_0track->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax); // rough phi fit
    //TF1 *sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
    TF1 *sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);

    sgf->SetParameter(0, r1->Parameter(0));
    sgf->SetParameter(1, r1->Parameter(1));
    sgf->SetParameter(2, r1->Parameter(2));
    sgf->SetParameter(3, skewInit);
    sgf->SetParameter(4, rphi->Parameter(0));
    sgf->SetParameter(5, rphi->Parameter(1));
    sgf->SetParameter(6, 1.);
    sgf->SetParameter(7, 1.);
    sgf->SetParameter(8, 1.);
    sgf->SetParameter(9, 1.);

    //sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), skewInit, 1., 1., 1., 1.);
    sgf->SetParLimits(0, 0, omegaConstMax); // omega "height"
    sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean
    sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
    sgf->SetParLimits(3, skewMin, 0.); // skew
    sgf->SetParLimits(4, 0, phiConstMax);  // phi "height"
    sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
      
    sgf->SetParName(0, "#omega constant");
    sgf->SetParName(1, "#omega mean");
    sgf->SetParName(2, "#sigma");
    sgf->SetParName(3, "Skew");
    sgf->SetParName(4, "#phi constant");
    sgf->SetParName(5, "#phi mean");
    sgf->SetParName(6, "Constant");
    sgf->SetParName(7, "Linear");
    sgf->SetParName(8, "Quadratic");
    sgf->SetParName(9, "Cubic");

    TFitResultPtr r = locHist_MMOP_p_0track->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram or print stat box
    //TFitResultPtr r = locHist_MMOP_p_0track->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print stat box
    Double_t chiSquare;
    chiSquare = locHist_MMOP_p_0track->Chisquare(sgf, "R");

    TMatrixD covMatrix(4,4);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<4; i++){
      for (Int_t j=0; j<4; j++){
  	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

    TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
    sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
    double intomega = sgfplot->Integral(fitMin, fitMaxFull);
    intomega = intomega/binsize;
    double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
    sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));

    TF1 *polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
    polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
    
    sgfplot->SetLineColor(kBlue);
    sgfPhi->SetLineColor(kGreen);
    polBG->SetLineColor(kViolet);

    // Change line width for smaller plots
    sgfplot->SetLineWidth(2);
    sgfPhi->SetLineWidth(2);
    polBG->SetLineWidth(2);
    //locHist_MMOP_p_0track->GetFunction("sgf")->SetLineWidth(2);
    sgf->SetLineWidth(2);
    
    locHist_MMOP_p_0track->Draw();
    polBG->Draw("same");
    sgfPhi->Draw("same");
    sgfplot->Draw("same");
    sgf->Draw("same");
    
    cout << "Integral under omega peak (from fit) = " << intomega << endl;
    cout << "Error = " << error << endl;

    for(int im = 0; im < 7; im++){
      //locHist_p_denom[im]->AddBinContent(i+1, intomega);
      //locHist_p_denom[im]->SetBinError(i+1, error);
      locHist_p_missing[im]->SetBinContent(i+1, intomega);
      locHist_p_missing[im]->SetBinError(i+1, error);

      for(int ipar = 0; ipar < 10; ipar++){
	locHist_p_miss_params[im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
      }
      locHist_p_miss_params[im][10]->SetBinContent(i+1, chiSquare);
    }
  }
  cp_0track->Print(Form("plots_%s/MMOP_%s_p_fits_missing.pdf", plotDir.Data(), plotData.Data()));

  /************************************ FIT MASS VS THETA IN p BINS (with P-value cuts) ***************************************/
  TH2F* locHist_theta2D_1track_p[9][7];
  TH2F* locHist_theta2D_0track_p[9];
  for(int i = 0; i < 9; i++){
    for(int im = 0; im < 7; im++){
      locHist_theta2D_1track_p[i][im] = (TH2F*)fhist->Get(Form("OmegaMassVsTheta_%scut%d_p%d", whichcut, im+1, i+1));
      locHist_theta2D_1track_p[i][im]->Sumw2();
      locHist_theta2D_1track_p[i][im]->RebinX(10);
      locHist_theta2D_1track_p[i][im]->RebinY(30);
    }
    locHist_theta2D_0track_p[i] = (TH2F*)fhist->Get(Form("OmegaMassVsTheta_0track_p%d", i+1));
    locHist_theta2D_0track_p[i]->Sumw2();
    locHist_theta2D_0track_p[i]->RebinX(10);
    locHist_theta2D_0track_p[i]->RebinY(30);
  }

  TCanvas* cthetap1[9][7];
  TCanvas* cthetap0[9];
  for(int i = 0; i < 9; i++){
    for(int im = 0; im < 7; im++){
      cthetap1[i][im] = new TCanvas(Form("cthetap1%d_m%d", i+1, im+1), Form("Omega Mass vs Theta in %dth p Bin (1 track)(%dth cut)", i+1, im), 1);
      cthetap1[i][im]->Divide(5,4);
    }
    cthetap0[i] = new TCanvas(Form("cthetap0%d", i+1), Form("Omega Mass vs Theta in %dth p Bin (0 tracks)", i+1), 1);
    cthetap0[i]->Divide(5,4);
  }

  TH1F* locHist_yield_theta_p[9][7];
  //TH1F* locHist_yield_theta_denom_p[9][7];
  TH1F* locHist_yield_theta_missing_p[9][7];
  //TH1F* locHist_error_1track[9];
  //TH1F* locHist_error_0track[9];

  TH1F* locHist_thetap_params[9][7][11];
  TH1F* locHist_thetap_miss_params[9][7][11];
  for(int i = 0; i < 9; i++){
    for(int im = 0; im < 7; im++){
      locHist_yield_theta_p[i][im] = new TH1F(Form("theta_yield_%d_m%d", i+1, im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
      locHist_yield_theta_p[i][im]->Sumw2();
      //locHist_yield_theta_denom_p[i][im] = new TH1F(Form("theta_yield_denom_%d_m%d", i+1, im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
      //locHist_yield_theta_denom_p[i][im]->Sumw2();
      locHist_yield_theta_missing_p[i][im] = new TH1F(Form("theta_yield_missing_%d_m%d", i+1, im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
      locHist_yield_theta_missing_p[i][im]->Sumw2();
      for(int ipar = 0; ipar < 11; ipar++){
	locHist_thetap_params[i][im][ipar] = new TH1F(Form("thetap_mmop_param_%d_p%d_pv%d", ipar+1, i+1, im+1), Form(";#theta_{mmop} (deg);%s", paramNames[ipar].Data()), 20, 0., 30.);
	locHist_thetap_miss_params[i][im][ipar] = new TH1F(Form("thetap_mmop_miss_param_%d_p%d_pv%d", ipar+1, i+1, im+1), Form(";#theta_{mmop} (deg);%s", paramNames[ipar].Data()), 20, 0., 30.);
      }
    }
    //locHist_error_1track[i] = new TH1F(Form("error_1track_%d", i+1), ";#theta_{mmop} (deg); error", 20, 0., 30.);
    //locHist_error_0track[i] = new TH1F(Form("error_0track_%d", i+1), ";#theta_{mmop} (deg); error", 20, 0., 30.);
  }

  //Fit theta (1 track):
  for(int ip = 0; ip < 9; ip++){
    for(int im = 0; im < 7; im++){
      for(int i = 0; i < 20; i++){
	cthetap1[ip][im]->cd(i+1);
	TH1F* locHist_Mass = (TH1F*)locHist_theta2D_1track_p[ip][im]->ProjectionX(Form("theta1_p%d_%d_%d", ip, im, i), i+1, i+1);
	double hsum = locHist_Mass->Integral();
	if(hsum < 75.)
	  continue;

	double thetamin = 0. + 1.5*i;
	double thetamax = 1.5 + 1.5*i;

	TString ts;
	ts += Form("%.1f", thetamin);
	ts += "#circ < #theta < ";
	ts += Form("%.1f", thetamax);
	ts += "#circ;Missing Mass off the Proton (GeV)";

	locHist_Mass->SetTitle(ts);

      	// Preliminary fit to phi (meson) peak
      	TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
      	gausphi->SetParameters(100., 1., 0.1);
      	gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);
	
	TFitResultPtr r1 = locHist_Mass->Fit("gaus", "SQN", "", 0.7, 0.9); // rough omega fit
    	TFitResultPtr rphi = locHist_Mass->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax); // rough phi fit
	//TF1* sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
    	TF1 *sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);

    	sgf->SetParameter(0, r1->Parameter(0));
    	sgf->SetParameter(1, r1->Parameter(1));
    	sgf->SetParameter(2, r1->Parameter(2));
    	sgf->SetParameter(3, skewInit);
    	sgf->SetParameter(4, rphi->Parameter(0));
    	sgf->SetParameter(5, rphi->Parameter(1));
    	sgf->SetParameter(6, 1.);
    	sgf->SetParameter(7, 1.);
    	sgf->SetParameter(8, 1.);
    	sgf->SetParameter(9, 1.);
	//sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), skewInit, 1., 1., 1., 1.);

    	sgf->SetParLimits(0, 0, omegaConstMax); // omega "height"
    	sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean
    	sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
    	sgf->SetParLimits(3, skewMin, 0.); // skew
    	sgf->SetParLimits(4, 0, phiConstMax);  // phi "height"
    	sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
      
    	sgf->SetParName(0, "#omega constant");
    	sgf->SetParName(1, "#omega mean");
    	sgf->SetParName(2, "#sigma");
    	sgf->SetParName(3, "Skew");
    	sgf->SetParName(4, "#phi constant");
    	sgf->SetParName(5, "#phi mean");
    	sgf->SetParName(6, "Constant");
    	sgf->SetParName(7, "Linear");
    	sgf->SetParName(8, "Quadratic");
    	sgf->SetParName(9, "Cubic");


	TFitResultPtr r = locHist_Mass->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram or print stat box
	//TFitResultPtr r = locHist_Mass->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print stat box 
      	Double_t chiSquare;
      	chiSquare = locHist_Mass->Chisquare(sgf, "R");
	
	TMatrixD covMatrix(4,4);
 
	TMatrixDSym fCovar = r->GetCovarianceMatrix();
	for (Int_t k=0; k<4; k++){
	  for (Int_t j=0; j<4; j++){
	    covMatrix[k][j] = fCovar[k][j];
	  }
	}

	Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};
      
	TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
	sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
	double intomega = sgfplot->Integral(fitMin, fitMaxFull);
	intomega = intomega/binsize;
	double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
	if(error == 0)
	  cout << "WARNING: IntegralError failed for p-bin = " << ip << ", hist bin = " << i << ". (1-track fit)" << endl;
	error = error/binsize;

    	TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
    	sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));
      
	TF1 *polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
	polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
      
	sgfplot->SetLineColor(kBlue);
    	sgfPhi->SetLineColor(kGreen);
	polBG->SetLineColor(kViolet);

    	// Change line width for smaller plots
    	sgfplot->SetLineWidth(2);
    	sgfPhi->SetLineWidth(2);
    	polBG->SetLineWidth(2);
    	//locHist_Mass->GetFunction("sgf")->SetLineWidth(2);
    	sgf->SetLineWidth(2);

	locHist_Mass->Draw();
	polBG->Draw("same");
    	sgfPhi->Draw("same");
	sgfplot->Draw("same");
	sgf->Draw("same");

	locHist_yield_theta_p[ip][im]->SetBinContent(i+1, intomega);
	locHist_yield_theta_p[ip][im]->SetBinError(i+1, error);
	//locHist_yield_theta_denom_p[ip][im]->SetBinContent(i+1, intomega);
	//locHist_error_1track[ip]->SetBinContent(i+1, error/intomega);

	for(int ipar = 0; ipar < 10; ipar++){
	  locHist_thetap_params[ip][im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
	}
	locHist_thetap_params[ip][im][10]->SetBinContent(i+1, chiSquare);
      }
      if(im == 0){
      	cthetap1[ip][im]->Print(Form("plots_%s/MMOP_p_bins/MMOP_%s_ThetaFits1_p%d.pdf", plotDir.Data(), plotData.Data(), ip+1));
      }
    }

    for(int i = 0; i < 20; i++){
      cthetap0[ip]->cd(i+1);
      TH1F* locHist_Mass = (TH1F*)locHist_theta2D_0track_p[ip]->ProjectionX(Form("theta0_p%d_%d", ip, i), i+1, i+1);
      double hsum = locHist_Mass->Integral();
      if(hsum < 75.)
	continue;

      double thetamin = 0. + 1.5*i;
      double thetamax = 1.5 + 1.5*i;

      TString ts;
      ts += Form("%.1f", thetamin);
      ts += "#circ < #theta < ";
      ts += Form("%.1f", thetamax);
      ts += "#circ;Missing Mass off the Proton (GeV)";

      locHist_Mass->SetTitle(ts);

      // Preliminary fit to phi (meson) peak
      TF1 *gausphi = new TF1("gausphi", "gaus", fitMin, fitMaxFull);
      gausphi->SetParameters(100., 1., 0.1);
      gausphi->SetParLimits(1, phiMeanMin, phiMeanMax);

      TFitResultPtr r1 = locHist_Mass->Fit("gaus", "SQN", "", 0.7, 0.9); // rough omega fit
      TFitResultPtr rphi = locHist_Mass->Fit(gausphi, "SQN", "", phiMeanMin, phiMeanMax); // rough phi fit
      //TF1* sgf = new TF1("sgf", sgf1, fitMin, fitMaxIfBG, 8);
      TF1 *sgf = new TF1("sgf", sgf2cube, fitMin, fitMaxIfBG, 10);

      sgf->SetParameter(0, r1->Parameter(0));
      sgf->SetParameter(1, r1->Parameter(1));
      sgf->SetParameter(2, r1->Parameter(2));
      sgf->SetParameter(3, skewInit);
      sgf->SetParameter(4, rphi->Parameter(0));
      sgf->SetParameter(5, rphi->Parameter(1));
      sgf->SetParameter(6, 1.);
      sgf->SetParameter(7, 1.);
      sgf->SetParameter(8, 1.);
      sgf->SetParameter(9, 1.);

      //sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), skewInit, 1., 1., 1., 1.);

      sgf->SetParLimits(0, 0, omegaConstMax); // omega dist. amplitude
      sgf->SetParLimits(1, omegaMeanMin, omegaMeanMax); // omega mean
      sgf->SetParLimits(2, sigmaMin, sigmaMax); // sigma
      sgf->SetParLimits(3, skewMin, 0.); // skew
      sgf->SetParLimits(4, 0, phiConstMax);  // phi "height"
      sgf->SetParLimits(5, phiMeanMin, phiMeanMax); // phi mean
      
      sgf->SetParName(0, "#omega constant");
      sgf->SetParName(1, "#omega mean");
      sgf->SetParName(2, "#sigma");
      sgf->SetParName(3, "Skew");
      sgf->SetParName(4, "#phi constant");
      sgf->SetParName(5, "#phi mean");
      sgf->SetParName(6, "Constant");
      sgf->SetParName(7, "Linear");
      sgf->SetParName(8, "Quadratic");
      sgf->SetParName(9, "Cubic");

      TFitResultPtr r = locHist_Mass->Fit(sgf, "SQN", "", fitMin, fitMaxIfBG); // Don't attach fit to histogram or print stat box
      //TFitResultPtr r = locHist_Mass->Fit(sgf, "SQ", "", fitMin, fitMaxIfBG); // Attach fit to histogram and print stat box
      Double_t chiSquare;
      chiSquare = locHist_Mass->Chisquare(sgf, "R");

      TMatrixD covMatrix(4,4);
 
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t k=0; k<4; k++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[k][j] = fCovar[k][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};
      
      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, fitMin, fitMaxFull, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(fitMin, fitMaxFull);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(fitMin, fitMaxFull, params, covMatrix.GetMatrixArray(), .0001);
      if(error == 0)
	cout << "WARNING: IntegralError failed for p-bin = " << ip << ", hist bin = " << i << ". (0-track fit)" << endl;
      error = error/binsize;

      TF1 *sgfPhi = new TF1("sgfPhi", sgf1plot, fitMin, fitMaxFull, 4);
      sgfPhi->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(2), r->Parameter(3));
      
      TF1 *polBG = new TF1("polBG", "pol3", fitMin, fitMaxIfBG);
      polBG->SetParameters(r->Parameter(6), r->Parameter(7), r->Parameter(8), r->Parameter(9));
      
      sgfplot->SetLineColor(kBlue);
      sgfPhi->SetLineColor(kGreen);
      polBG->SetLineColor(kViolet);

      // Change line width for smaller plots
      sgfplot->SetLineWidth(2);
      sgfPhi->SetLineWidth(2);
      polBG->SetLineWidth(2);
      //locHist_Mass->GetFunction("sgf")->SetLineWidth(2);
      sgf->SetLineWidth(2);

      locHist_Mass->Draw();
      polBG->Draw("same");
      sgfPhi->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");

      for(int im = 0; im < 7; im++){
	//locHist_yield_theta_denom_p[ip][im]->AddBinContent(i+1, intomega);
	//locHist_yield_theta_denom_p[ip][im]->SetBinError(i+1, error);
	locHist_yield_theta_missing_p[ip][im]->SetBinContent(i+1, intomega);
	locHist_yield_theta_missing_p[ip][im]->SetBinError(i+1, error);
	//locHist_error_0track[ip]->SetBinContent(i+1, error/intomega);

	for(int ipar = 0; ipar < 10; ipar++){
	  locHist_thetap_miss_params[ip][im][ipar]->SetBinContent(i+1, r->Parameter(ipar));
	}
	locHist_thetap_miss_params[ip][im][10]->SetBinContent(i+1, chiSquare);
      }
    }
    cthetap0[ip]->Print(Form("plots_%s/MMOP_p_bins/MMOP_%s_ThetaFitsMissing_p%d.pdf", plotDir.Data(), plotData.Data(), ip+1));
  }



  // TCanvas* c1 = new TCanvas("c1", "omega yields", 900, 300);
  // c1->Divide(3,1);
  // c1->cd(1);
  // locHist_phi_denom->SetMinimum(0);
  // locHist_phi_reco->SetLineColor(kGreen+3);
  // locHist_phi_denom->SetLineColor(kMagenta+3);
  // locHist_phi_reco->SetMarkerColor(kGreen+3);
  // locHist_phi_denom->SetMarkerColor(kMagenta+3);
  // locHist_phi_reco->SetMarkerStyle(kOpenCircle);
  // locHist_phi_denom->SetMarkerStyle(kFullCircle);
  // locHist_phi_denom->Draw();
  // locHist_phi_reco->Draw("same");
  // c1->cd(2);
  // locHist_theta_denom->SetMinimum(0);
  // locHist_theta_reco->SetLineColor(kGreen+3);
  // locHist_theta_denom->SetLineColor(kMagenta+3);
  // locHist_theta_reco->SetMarkerColor(kGreen+3);
  // locHist_theta_denom->SetMarkerColor(kMagenta+3);
  // locHist_theta_reco->SetMarkerStyle(kOpenCircle);
  // locHist_theta_denom->SetMarkerStyle(kFullCircle);
  // locHist_theta_denom->Draw();
  // locHist_theta_reco->Draw("same");
  // c1->cd(3);
  // locHist_p_denom->SetMinimum(0);
  // locHist_p_reco->SetLineColor(kGreen+3);
  // locHist_p_denom->SetLineColor(kMagenta+3);
  // locHist_p_reco->SetMarkerColor(kGreen+3);
  // locHist_p_denom->SetMarkerColor(kMagenta+3);
  // locHist_p_reco->SetMarkerStyle(kOpenCircle);
  // locHist_p_denom->SetMarkerStyle(kFullCircle);
  // locHist_p_denom->Draw();
  // locHist_p_reco->Draw("same");

  // TCanvas* c2 = new TCanvas("c2", "Omega yield vs theta, p", 900, 900);
  // c2->Divide(3, 3);
  // for(int i = 0; i < 9; i++){
  //   c2->cd(i+1);
  //   locHist_yield_theta_denom_p[i]->SetMinimum(0);
  //   locHist_yield_theta_denom_p[i]->SetMarkerStyle(kFullCircle);
  //   locHist_yield_theta_denom_p[i]->SetMarkerColor(kMagenta+3);
  //   locHist_yield_theta_denom_p[i]->SetLineColor(kMagenta+3);
  //   locHist_yield_theta_p[i]->SetMarkerStyle(kOpenCircle);
  //   locHist_yield_theta_p[i]->SetMarkerColor(kGreen+3);
  //   locHist_yield_theta_p[i]->SetLineColor(kGreen+3);
  //   locHist_yield_theta_denom_p[i]->Draw();
  //   locHist_yield_theta_p[i]->Draw("same");
  // }

  if(save == true){
    for(int im = 0; im < 7; im++){
      locHist_phi_reco[im]->Write();
      locHist_theta_reco[im]->Write();
      locHist_p_reco[im]->Write();
      locHist_phi_missing[im]->Write();
      locHist_theta_missing[im]->Write();
      locHist_p_missing[im]->Write();
      for(int ipar = 0; ipar < 11; ipar++){
	locHist_phi_params[im][ipar]->Write();
	locHist_phi_miss_params[im][ipar]->Write();
	locHist_theta_params[im][ipar]->Write();
	locHist_theta_miss_params[im][ipar]->Write();
	locHist_p_params[im][ipar]->Write();
	locHist_p_miss_params[im][ipar]->Write();
      }
      for(int i = 0; i < 9; i++){
	locHist_yield_theta_p[i][im]->Write();
	locHist_yield_theta_missing_p[i][im]->Write();
	for(int ipar = 0; ipar < 11; ipar++){
	  locHist_thetap_params[i][im][ipar]->Write();
	  locHist_thetap_miss_params[i][im][ipar]->Write();
	}
      }
    }
  }
 
  return;
}
