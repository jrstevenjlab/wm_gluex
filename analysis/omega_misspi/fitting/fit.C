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

double fitfn(double *x, double *par) {
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5], 2.0)) + par[6] + par[7]*x[0] + par[8]*x[0]*x[0] + par[9]*x[0]*x[0]*x[0];
}

double fit1gaus(double *x, double *par) {
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
}

double fit2gaus(double *x, double *par) { //two gaussians with the same centroid, cubic polynomial background
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[4], 2.0)) + par[5] + par[6]*x[0] + par[7]*x[0]*x[0] + par[8]*x[0]*x[0]*x[0];
}

double omega2gaus(double *x, double *par) { //two gaussians with the same centroid, no background
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[4], 2.0));
}

double sgf1(double *x, double *par) {
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2)) + par[4] + par[5]*x[0] + par[6]*x[0]*x[0] + par[7]*x[0]*x[0]*x[0];
}

double sgf1plot(double *x, double *par) {
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2));
}

void fit(bool mc = true, bool save = true, const char *whichcut = "pv", TString datadir = "/sciclone/gluex10/amschertz/omega_misspim/test7_full/", TString mcdir = "/sciclone/gluex10/amschertz/omega_misspim/geant4/test2/", TString name = "30274_31057") {
  TFile* fhist;
  TFile* outROOT;

  if(mc == false) {
    fhist = TFile::Open(datadir + Form("hist_sum_%s.root", name.Data()));
    if(save == true)    
      outROOT = new TFile(Form("mmop_data_%scut_%s.root", whichcut, name.Data()), "recreate");
  }
  else {
    fhist = TFile::Open(mcdir + Form("hist_sum_%s.root", name.Data()));
    if(save == true)
      outROOT = new TFile(Form("mmop_gen_%scut_%s.root", whichcut, name.Data()), "recreate");
  }

  for(int im = 0; im < 7; im++){ //loop over cuts on (P(chi^2) or theta, or 3-momentum)
    locHist_MMOP_phi2D[im] = (TH2F*)fhist->Get(Form("OmegaMassVsPhi_1track_%scut%d", whichcut, im+1));
    locHist_MMOP_phi2D[im]->Sumw2();
    locHist_MMOP_phi2D[im]->RebinX(10);
    locHist_MMOP_phi2D[im]->RebinY(30);

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
    cphi[im] = new TCanvas(Form("cphi_%d", im+1), "Phi fits (1 candidate track)", 600, 600);
    cphi[im]->Divide(5,4);
    ctheta[im] = new TCanvas(Form("ctheta_%d", im+1), "Theta fits (1 candidate track)", 600, 600);
    ctheta[im]->Divide(5,4);
    cp[im] = new TCanvas(Form("cp_%d", im+1), "P fits (1 candidate track)", 600, 600);
    cp[im]->Divide(5,4);
  }

  TCanvas *cphi_0track = new TCanvas("cphi_0track", "Phi fits (0 candidate tracks)", 600, 600);
  cphi_0track->Divide(5,4);

  TCanvas *ctheta_0track = new TCanvas("ctheta_0track", "Theta fits (0 candidate tracks)", 600, 600);
  ctheta_0track->Divide(5,4);

  TCanvas *cp_0track = new TCanvas("cp_0track", "P fits (0 candidate tracks)", 600, 600);
  cp_0track->Divide(5,4);

  gStyle->SetOptStat(000000);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(1);

  TH1F* locHist_phi_reco[7];
  TH1F* locHist_theta_reco[7];
  TH1F* locHist_p_reco[7];
  TH1F* locHist_phi_denom[7];
  TH1F* locHist_theta_denom[7];
  TH1F* locHist_p_denom[7];

  for(int im = 0; im < 7; im++){
    locHist_phi_reco[im] = new TH1F(Form("phi_mmop_reco_%d", im+1), "#omega_{mmop} (1 #pi^{+});#phi_{mmop} (deg);#omega yield / 18#circ", 20, -180., 180.);
    locHist_phi_reco[im]->Sumw2();
    locHist_theta_reco[im] = new TH1F(Form("theta_mmop_reco_%d", im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
    locHist_theta_reco[im]->Sumw2();
    locHist_p_reco[im] = new TH1F(Form("p_mmop_reco_%d", im+1), ";p_{mmop} (GeV);#omega yield / 0.5 GeV", 12, 0., 6.);
    locHist_p_reco[im]->Sumw2();
    locHist_phi_denom[im] = new TH1F(Form("phi_mmop_denom_%d", im+1), "#omega_{mmop} (1 or 0 #pi^{+});#phi_{mmop} (deg);#omega yield / 18#circ", 20, -180., 180.);
    locHist_phi_denom[im]->Sumw2();
    locHist_theta_denom[im] = new TH1F(Form("theta_mmop_denom_%d", im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
    locHist_theta_denom[im]->Sumw2();
    locHist_p_denom[im] = new TH1F(Form("p_mmop_denom_%d", im+1), ";p_{mmop} (GeV);#omega yield / 0.5 GeV", 12, 0., 6.);
    locHist_p_denom[im]->Sumw2();
  }

  double phi_errorD[20];
  double theta_errorD[20];
  double p_errorD[12];
  double phi_errorU[20];
  double theta_errorU[20];
  double p_errorU[12];

  cout << "Phi fits" << endl;
  for(int im = 0; im < 7; im++){//loop over mm^2 cuts
    for(int i = 0; i < 20; i++){
      cphi[im]->cd(i+1);
      locHist_MMOP_phi = (TH1F*)locHist_MMOP_phi2D[im]->ProjectionX(Form("phi_%d", i), i+1, i+1);

      double phimin = -180. + 18.*i;
      double phimax = -162. + 18.*i;

      TString ts;
      ts += phimin;
      ts += "#circ < #phi < ";
      ts += phimax;
      ts += "#circ;Missing Mass off the Proton (GeV)";

      locHist_MMOP_phi->SetTitle(ts);
    
      TFitResultPtr r1 = locHist_MMOP_phi->Fit("gaus", "SQ", "", 0.7, 0.9);
      TF1 *sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);

      sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
      sgf->SetParName(7, "Cubic");

      TFitResultPtr r = locHist_MMOP_phi->Fit(sgf, "SQ");

      TMatrixD covMatrix(4,4);
 
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<4; i++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};
      
      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(0.3, 1.3);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
      error = error/binsize;

      TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
      polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
    
      sgfplot->SetLineColor(kBlue);
      polBG->SetLineColor(kViolet);
    
      locHist_MMOP_phi->Draw();
      polBG->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");

      locHist_phi_reco[im]->SetBinContent(i+1, intomega);
      locHist_phi_reco[im]->SetBinError(i+1, error);
      locHist_phi_denom[im]->SetBinContent(i+1, intomega);
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
    
    TFitResultPtr r1 = locHist_MMOP_phi_0track->Fit("gaus", "SQ", "", 0.7, 0.9);
    TF1 *sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);

    sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
    sgf->SetParName(7, "Cubic");

    TFitResultPtr r = locHist_MMOP_phi_0track->Fit(sgf, "SQ");

    TMatrixD covMatrix(4,4);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<4; i++){
      for (Int_t j=0; j<4; j++){
	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

    TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
    sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
    double intomega = sgfplot->Integral(0.3, 1.3);
    intomega = intomega/binsize;
    double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
    polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
    
    sgfplot->SetLineColor(kBlue);
    polBG->SetLineColor(kViolet);
    
    locHist_MMOP_phi_0track->Draw();
    polBG->Draw("same");
    sgfplot->Draw("same");
    sgf->Draw("same");
    
    cout << "Integral under omega peak (from fit) = " << intomega << endl;
    cout << "Error = " << error << endl;

    phi_errorU[i] = error;

    for(int im = 0; im < 7; im++){
      locHist_phi_denom[im]->AddBinContent(i+1, intomega);
      locHist_phi_denom[im]->SetBinError(i+1, error);
    }
  }


  cout << "Theta fits:" << endl;
  //fit theta (reco)
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 20; i++){
      ctheta[im]->cd(i+1);
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

      TFitResultPtr r1 = locHist_MMOP_theta->Fit("gaus", "SQ", "", 0.7, 0.9);
      TF1 *sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);
      sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
      sgf->SetParName(7, "Cubic");

      TFitResultPtr r = locHist_MMOP_theta->Fit(sgf, "SQ");
      
      TMatrixD covMatrix(4,4);
      
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<4; i++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(0.3, 1.3);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
      error = error/binsize;

      TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
      polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
    
      sgfplot->SetLineColor(kBlue);
      polBG->SetLineColor(kViolet);
    
      locHist_MMOP_theta->Draw();
      polBG->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");

      locHist_theta_reco[im]->SetBinContent(i+1, intomega);
      locHist_theta_reco[im]->SetBinError(i+1, error);
      locHist_theta_denom[im]->SetBinContent(i+1, intomega);
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
    TFitResultPtr r1 = locHist_MMOP_theta_0track->Fit("gaus", "SQ", "", 0.7, 0.9);
    TF1 *sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);
    sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
    sgf->SetParName(7, "Cubic");

    TFitResultPtr r = locHist_MMOP_theta_0track->Fit(sgf, "SQ");

    TMatrixD covMatrix(4,4);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<4; i++){
      for (Int_t j=0; j<4; j++){
  	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

    TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
    sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
    double intomega = sgfplot->Integral(0.3, 1.3);
    intomega = intomega/binsize;
    double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
    polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
    
    sgfplot->SetLineColor(kBlue);
    polBG->SetLineColor(kViolet);
    
    locHist_MMOP_theta_0track->Draw();
    polBG->Draw("same");
    sgfplot->Draw("same");
    sgf->Draw("same");
    
    //Time to check errors:
    double intHist = locHist_MMOP_theta_0track->Integral();

    cout << "Integral under omega peak (from fit) = " << intomega << endl;
    cout << "Error (from fit) = " << error << endl;
    cout << "Integral under omega peak (from hist) = " << intHist << endl;
    cout << "Error (from sqrt(N)) = " << TMath::Sqrt(intHist) << endl;

    for(int im = 0; im < 7; im++){
      locHist_theta_denom[im]->AddBinContent(i+1, intomega);
      locHist_theta_denom[im]->SetBinError(i+1, error);
    }
  }


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
      
      TFitResultPtr r1 = locHist_MMOP_p->Fit("gaus", "SQ", "", 0.7, 0.9);
      TF1 *sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);
      sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
      sgf->SetParName(7, "Cubic");

      TFitResultPtr r = locHist_MMOP_p->Fit(sgf, "SQ");

      TMatrixD covMatrix(4,4);
      
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<4; i++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(0.3, 1.3);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
      error = error/binsize;

      TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
      polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
      
      sgfplot->SetLineColor(kBlue);
      polBG->SetLineColor(kViolet);
    
      locHist_MMOP_p->Draw();
      polBG->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");
    
      cout << "Integral under omega peak = " << intomega << endl;
      cout << "Error = " << error << endl;

      locHist_p_reco[im]->SetBinContent(i+1, intomega);
      locHist_p_reco[im]->SetBinError(i+1, error);
      locHist_p_denom[im]->SetBinContent(i+1, intomega);
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
    
    TFitResultPtr r1 = locHist_MMOP_p_0track->Fit("gaus", "SQ", "", 0.7, 0.9);
    TF1 *sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);
    sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
    sgf->SetParName(7, "Cubic");

    TFitResultPtr r = locHist_MMOP_p_0track->Fit(sgf, "SQ");

    TMatrixD covMatrix(4,4);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<4; i++){
      for (Int_t j=0; j<4; j++){
  	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};

    TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
    sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
    double intomega = sgfplot->Integral(0.3, 1.3);
    intomega = intomega/binsize;
    double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
    polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
    
    sgfplot->SetLineColor(kBlue);
    polBG->SetLineColor(kViolet);
    
    locHist_MMOP_p_0track->Draw();
    polBG->Draw("same");
    sgfplot->Draw("same");
    sgf->Draw("same");
    
    cout << "Integral under omega peak (from fit) = " << intomega << endl;
    cout << "Error = " << error << endl;

    for(int im = 0; im < 7; im++){
      locHist_p_denom[im]->AddBinContent(i+1, intomega);
      locHist_p_denom[im]->SetBinError(i+1, error);
    }
  }

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
  TH1F* locHist_yield_theta_denom_p[9][7];
  //TH1F* locHist_error_1track[9];
  //TH1F* locHist_error_0track[9];
  for(int i = 0; i < 9; i++){
    for(int im = 0; im < 7; im++){
      locHist_yield_theta_p[i][im] = new TH1F(Form("theta_yield_%d_m%d", i+1, im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
      locHist_yield_theta_p[i][im]->Sumw2();
      locHist_yield_theta_denom_p[i][im] = new TH1F(Form("theta_yield_denom_%d_m%d", i+1, im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
      locHist_yield_theta_denom_p[i][im]->Sumw2();
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
	
	TFitResultPtr r1 = locHist_Mass->Fit("gaus", "SQ0", "", 0.7, 0.9);
	TF1* sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);
	sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
	sgf->SetParName(7, "Cubic");

	TFitResultPtr r = locHist_Mass->Fit(sgf, "SQ0");
	
	TMatrixD covMatrix(4,4);
 
	TMatrixDSym fCovar = r->GetCovarianceMatrix();
	for (Int_t k=0; k<4; k++){
	  for (Int_t j=0; j<4; j++){
	    covMatrix[k][j] = fCovar[k][j];
	  }
	}

	Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};
      
	TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
	sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
	double intomega = sgfplot->Integral(0.3, 1.3);
	intomega = intomega/binsize;
	double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
	if(error == 0)
	  cout << "WARNING: IntegralError failed for p-bin = " << ip << ", hist bin = " << i << ". (1-track fit)" << endl;
	error = error/binsize;
      
	TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
	polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
      
	sgfplot->SetLineColor(kBlue);
	polBG->SetLineColor(kViolet);

	locHist_Mass->Draw();
	polBG->Draw("same");
	sgfplot->Draw("same");
	sgf->Draw("same");

	locHist_yield_theta_p[ip][im]->SetBinContent(i+1, intomega);
	locHist_yield_theta_p[ip][im]->SetBinError(i+1, error);
	locHist_yield_theta_denom_p[ip][im]->SetBinContent(i+1, intomega);
	//locHist_error_1track[ip]->SetBinContent(i+1, error/intomega);
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

      TFitResultPtr r1 = locHist_Mass->Fit("gaus", "SQ0", "", 0.7, 0.9);
      TF1* sgf = new TF1("sgf", sgf1, 0.3, 1.3, 8);
      sgf->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), -0.5, 1., 1., 1., 1.);
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
      sgf->SetParName(7, "Cubic");

      TFitResultPtr r = locHist_Mass->Fit(sgf, "SQ0");

      TMatrixD covMatrix(4,4);
 
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t k=0; k<4; k++){
	for (Int_t j=0; j<4; j++){
	  covMatrix[k][j] = fCovar[k][j];
	}
      }

      Double_t params[4] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3)};
      
      TF1 *sgfplot = new TF1("sgfplot", sgf1plot, 0.3, 1.3, 4);
      sgfplot->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3));
      double intomega = sgfplot->Integral(0.3, 1.3);
      intomega = intomega/binsize;
      double error = sgfplot->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
      if(error == 0)
	cout << "WARNING: IntegralError failed for p-bin = " << ip << ", hist bin = " << i << ". (0-track fit)" << endl;
      error = error/binsize;
      
      TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
      polBG->SetParameters(r->Parameter(4), r->Parameter(5), r->Parameter(6), r->Parameter(7));
      
      sgfplot->SetLineColor(kBlue);
      polBG->SetLineColor(kViolet);

      locHist_Mass->Draw();
      polBG->Draw("same");
      sgfplot->Draw("same");
      sgf->Draw("same");

      for(int im = 0; im < 7; im++){
	locHist_yield_theta_denom_p[ip][im]->AddBinContent(i+1, intomega);
	locHist_yield_theta_denom_p[ip][im]->SetBinError(i+1, error);
	//locHist_error_0track[ip]->SetBinContent(i+1, error/intomega);
      }
    }
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
      locHist_phi_denom[im]->Write();
      locHist_theta_denom[im]->Write();
      locHist_p_denom[im]->Write();
      for(int i = 0; i < 9; i++){
	locHist_yield_theta_p[i][im]->Write();
	locHist_yield_theta_denom_p[i][im]->Write();
      }
    }
  }
 
  return;
}
