#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"

double fit1gaus(double *x, double *par) {
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
}

double fit4gaus(double *x, double *par) {
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[4], 2.0)) + par[5]*TMath::Exp(-0.5*TMath::Power((x[0]-par[6])/par[7], 2.0)) + par[8]*TMath::Exp(-0.5*TMath::Power((x[0]-par[9])/par[10], 2.0)) + par[11] + par[12]*x[0] + par[13]*x[0]*x[0] + par[14]*x[0]*x[0]*x[0];
}

double omega2gaus(double *x, double *par) {
  return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.0)) + par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[4], 2.0));
}

double sgf1(double *x, double *par) {
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2)) + par[4] + par[5]*x[0] + par[6]*x[0]*x[0] + par[7]*x[0]*x[0]*x[0];
}

double sgf1plot(double *x, double *par) {
  return par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2));
}

void threepifit(bool mc = true, bool save = true, const char *whichcut = "pv", TString datadir = "/sciclone/gluex10/amschertz/omega_misspim/test7_full/", TString mcdir = "/sciclone/gluex10/amschertz/omega_misspim/geant4/test2/",  TString name = "30274_31057") {
  TFile* fhist;
  TFile* outROOT;
  if(mc == false) {
    fhist = TFile::Open(datadir + Form("hist_sum_%s.root", name.Data()));
    if(save == true)
      outROOT = new TFile(Form("method2_data_%scut_%s.root", whichcut, name.Data()), "recreate");
  }
  else {
    fhist = TFile::Open(mcdir + Form("hist_sum_%s.root", name.Data()));
    if(save == true)
      outROOT = new TFile(Form("method2_gen_%scut_%s.root", whichcut, name.Data()), "recreate");
  }

  TH2F* locHist_3Pi_phi2D[7];
  TH2F* locHist_3Pi_theta2D[7];
  TH2F* locHist_3Pi_p2D[7];

  for(int im = 0; im < 7; im++){
    locHist_3Pi_phi2D[im] = (TH2F*)fhist->Get(Form("3PiMassVsPhi_%scut%d", whichcut, im+1));
    locHist_3Pi_phi2D[im]->Sumw2();
    locHist_3Pi_phi2D[im]->RebinX(10);
    locHist_3Pi_phi2D[im]->RebinY(30);
    locHist_3Pi_theta2D[im] = (TH2F*)fhist->Get(Form("3PiMassVsTheta_%scut%d", whichcut, im+1));
    locHist_3Pi_theta2D[im]->Sumw2();
    locHist_3Pi_theta2D[im]->RebinX(10);
    locHist_3Pi_theta2D[im]->RebinY(30);
    locHist_3Pi_p2D[im] = (TH2F*)fhist->Get(Form("3PiMassVsP_%scut%d", whichcut, im+1));
    locHist_3Pi_p2D[im]->Sumw2();
    locHist_3Pi_p2D[im]->RebinX(10);
    locHist_3Pi_p2D[im]->RebinY(30);
  }

  TH2F* locHist_diff_phi2D = (TH2F*)fhist->Get("OmegaMassVsPhi_0tracks");
  locHist_diff_phi2D->Sumw2();
  locHist_diff_phi2D->RebinX(10);
  locHist_diff_phi2D->RebinY(30);
  TH2F* locHist_diff_theta2D = (TH2F*)fhist->Get("OmegaMassVsTheta_0tracks");
  locHist_diff_theta2D->Sumw2();
  locHist_diff_theta2D->RebinX(10);
  locHist_diff_theta2D->RebinY(30);
  TH2F* locHist_diff_p2D = (TH2F*)fhist->Get("OmegaMassVsP_0tracks");
  locHist_diff_p2D->Sumw2();
  locHist_diff_p2D->RebinX(10);
  locHist_diff_p2D->RebinY(30);

  double binsize = 1./60.;

  TCanvas* cphi[7];
  TCanvas* ctheta[7];
  TCanvas* cp[7];

  for(int im = 0; im < 7; im++){
    cphi[im] = new TCanvas(Form("cphi_%d", im+1), "Phi fits", 600, 600);
    cphi[im]->Divide(5,4);
    ctheta[im] = new TCanvas(Form("ctheta_%d", im+1), "Theta fits", 600, 600);
    ctheta[im]->Divide(5,4);
    cp[im] = new TCanvas(Form("cp_%d", im+1), "3-Momentum fits", 600, 600);
    cp[im]->Divide(5,4);
  }

  TCanvas *cphi_diff = new TCanvas("cphi_diff", "Phi fits", 600, 600);
  cphi_diff->Divide(5,4);
  TCanvas *ctheta_diff = new TCanvas("ctheta_diff", "Theta fits", 600, 600);
  ctheta_diff->Divide(5,4);
  TCanvas *cp_diff = new TCanvas("cp_diff", "3-Momentum fits", 600, 600);
  cp_diff->Divide(5,4);

  TH1F* locHist_yield_phi[7];
  TH1F* locHist_yield_theta[7];
  TH1F* locHist_yield_p[7];
  TH1F* locHist_yield_phi_denom[7];
  TH1F* locHist_yield_theta_denom[7];
  TH1F* locHist_yield_p_denom[7];

  for(int im = 0; im < 7; im++){
    locHist_yield_phi[im] = new TH1F(Form("phi_yield_%d", im+1), "#omega_{inv};#phi_{mmop} (deg);#omega yield / 18#circ", 20, -180., 180.);
    locHist_yield_theta[im] = new TH1F(Form("theta_yield_%d", im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
    locHist_yield_p[im] = new TH1F(Form("p_yield_%d", im+1), ";p_{mmop} (GeV);#omega yield / 0.5 GeV", 12, 0., 6.);
    locHist_yield_phi_denom[im] = new TH1F(Form("phi_yield_denom_%d", im+1), "#omega_{inv} + #omega_{mmop}(0 #pi^{+});#phi_{mmop} (deg);#omega yield / 18#circ", 20, -180., 180.);
    locHist_yield_theta_denom[im] = new TH1F(Form("theta_yield_denom_%d", im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
    locHist_yield_p_denom[im] = new TH1F(Form("p_yield_denom_%d", im+1), ";p_{mmop} (GeV);#omega yield / 0.5 GeV", 12, 0., 6.);
    locHist_yield_phi[im]->Sumw2();
    locHist_yield_theta[im]->Sumw2();
    locHist_yield_p[im]->Sumw2();
    locHist_yield_phi_denom[im]->Sumw2();
    locHist_yield_theta_denom[im]->Sumw2();
    locHist_yield_p_denom[im]->Sumw2();
  }

  gStyle->SetOptStat(000000);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(1);

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(7000);

  double phi_error[20];
  double theta_error[20];
  double p_error[20];


  
  //loop over phi bins
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 20; i++) {
      cphi[im]->cd(i+1);
      //cout << i+1 << "th phi bin:" << endl;
      TH1F* locHist_3Pi = (TH1F*)locHist_3Pi_phi2D[im]->ProjectionX(Form("phi_%d", i), i+1, i+1);

      double phimin = -180. + 18.*i;
      double phimax = -162. + 18.*i;

      TString ts;
      ts += phimin;
      ts += "#circ < #phi < ";
      ts += phimax;
      ts += "#circ;M_{#pi^{+}#pi^{-}#pi^{0}} (GeV)";

      locHist_3Pi->SetTitle(ts);
     
      TF1 *gauseta = new TF1("gauseta", "gaus", 0.4, 1.3);
      gauseta->SetParameters(100., 0.55, 0.1);
      gauseta->SetParLimits(1, 0.5, 0.6);
      
      TF1 *gausphi = new TF1("gausphi", "gaus", 0.4, 1.3);
      gausphi->SetParameters(100., 1., 0.1);
      gausphi->SetParLimits(1, 0.95, 1.05);

      TF1 *gausomegafit = new TF1("gausomegafit", "gaus", 0.4, 1.3);
      gausomegafit->SetParameters(100., 0.8, 0.1);

      TFitResultPtr reta = locHist_3Pi->Fit(gauseta, "SQ", "", 0.5, 0.6);
      TFitResultPtr rphi = locHist_3Pi->Fit(gausphi, "SQ", "", 0.95, 1.05);
      TFitResultPtr romega = locHist_3Pi->Fit(gausomegafit, "SQ", "", 0.7, 0.9);
     
      TF1* gaus4pol3 = new TF1("gaus4pol3", fit4gaus, 0.4, 1.3, 15);
      
      gaus4pol3->SetParameter(0, romega->Parameter(0));
      gaus4pol3->SetParameter(1, romega->Parameter(1));
      gaus4pol3->SetParameter(2, romega->Parameter(2));
      gaus4pol3->SetParameter(3, romega->Parameter(0));
      gaus4pol3->SetParameter(4, 2.0 * romega->Parameter(2));
      gaus4pol3->SetParameter(5, reta->Parameter(0));
      gaus4pol3->FixParameter(6, reta->Parameter(1));
      gaus4pol3->SetParameter(7, reta->Parameter(2));
      gaus4pol3->SetParameter(8, rphi->Parameter(0));
      gaus4pol3->FixParameter(9, rphi->Parameter(1));
      gaus4pol3->SetParameter(10, rphi->Parameter(2));
      gaus4pol3->SetParameter(11, 1.);
      gaus4pol3->SetParameter(12, 1.);
      gaus4pol3->SetParameter(13, 1.);
      gaus4pol3->SetParameter(14, 1.);
     
      gaus4pol3->SetParLimits(0, 0, 1000000);
      gaus4pol3->SetParLimits(1, 0.7, 0.9);
      gaus4pol3->SetParLimits(2, 0.01, 0.1);
      gaus4pol3->SetParLimits(3, 0, 1000000);
      gaus4pol3->SetParLimits(4, 0.01, 0.1);
      gaus4pol3->SetParLimits(5, 0, 10000);
      gaus4pol3->SetParLimits(7, 0.01, 0.1);
      gaus4pol3->SetParLimits(8, 0, 10000);
      gaus4pol3->SetParLimits(10, 0.01, 0.03);
     
      gaus4pol3->SetParName(0, "#omega_{1}");
      gaus4pol3->SetParName(1, "#bar{x}_{#omega}");
      gaus4pol3->SetParName(2, "#sigma_{#omega_{1}}");
      gaus4pol3->SetParName(3, "#omega_{2}");
      gaus4pol3->SetParName(4, "#sigma_{#omega_{2}}");
      gaus4pol3->SetParName(5, "#eta");
      gaus4pol3->SetParName(6, "#eta mean");
      gaus4pol3->SetParName(7, "#sigma_{#eta}");
      gaus4pol3->SetParName(8, "#phi");
      gaus4pol3->SetParName(9, "#phi mean");
      gaus4pol3->SetParName(10, "#sigma_{#phi}");
      gaus4pol3->SetParName(11, "const");
      gaus4pol3->SetParName(12, "lin");
      gaus4pol3->SetParName(13, "quad");
      gaus4pol3->SetParName(14, "cube");
     
      TFitResultPtr r;
      if(i == 11)
	r = locHist_3Pi->Fit(gaus4pol3, "SQ");
      else 
	r = locHist_3Pi->Fit(gaus4pol3, "SQ");


      TMatrixD covMatrix(5,5);
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<5; i++){
	for (Int_t j=0; j<5; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }
     
      Double_t params[5] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4)};

      TF1 *gausomega = new TF1("gausomega", omega2gaus, 0.4, 1.3, 5);
      gausomega->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4));
      double intomega = gausomega->Integral(0.4, 1.3);
      intomega = intomega/binsize;
      double error = gausomega->IntegralError(0.4, 1.3, params, covMatrix.GetMatrixArray(), 0.0001);
      error = error/binsize;
     

      TF1 *gaus1 = new TF1("gaus1", "gaus", 0.4, 1.3);
      gaus1->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2));
     
      TF1 *gaus2 = new TF1("gaus2", "gaus", 0.4, 1.3);
      gaus2->SetParameters(r->Parameter(3), r->Parameter(1), r->Parameter(4));

      TF1 *polBG = new TF1("polBG", "pol3", 0.4, 1.3);
      polBG->SetParameters(r->Parameter(11), r->Parameter(12), r->Parameter(13), r->Parameter(14));
 
      gauseta->SetParameters(r->Parameter(5), r->Parameter(6), r->Parameter(7));
      
      gausphi->SetParameters(r->Parameter(8), r->Parameter(9), r->Parameter(10));

      gaus1->SetLineColor(kBlue);
      gaus2->SetLineColor(kAzure+10);
      polBG->SetLineColor(kViolet);
      gauseta->SetLineColor(kGreen);
      gausphi->SetLineColor(kGreen);
    
      locHist_3Pi->SetMarkerStyle(kFullCircle);
      locHist_3Pi->Draw();
      polBG->Draw("same");
      gauseta->Draw("same");
      gausphi->Draw("same");
      gaus1->Draw("same");
      gaus2->Draw("same");
      gaus4pol3->Draw("same");
      
    if(i == 10) {
      TCanvas* c10 = new TCanvas("c10", "c10", 600, 600);
      c10->cd(1);

      locHist_3Pi->Draw();
      polBG->Draw("same");
      gauseta->Draw("same");
      gausphi->Draw("same");
      gaus1->Draw("same");
      gaus2->Draw("same");
      gaus4pol3->Draw("same");
    }

     // cout << "Integral under omega peak (from fit) = " << intomega << endl;
     // cout << "Error = " << error << endl;

      locHist_yield_phi[im]->SetBinContent(i+1, intomega);
      locHist_yield_phi[im]->SetBinError(i+1, error);
      locHist_yield_phi_denom[im]->SetBinContent(i+1, intomega);
     //phi_error[i] = error;
    }
  }
  
  cout << "Theta fits (inv):" << endl;

  //loop over theta bins
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 20; i++) {
      ctheta[im]->cd(i+1);
      //cout << i+1 << "th theta bin:" << endl;
      TH1F* locHist_3Pi = (TH1F*)locHist_3Pi_theta2D[im]->ProjectionX(Form("theta_%d", i), i+1, i+1);
      double hsum = locHist_3Pi->Integral();
      if(hsum < 75.)
	continue;

      double thetamin = 0. + 1.5*i;
      double thetamax = 1.5 + 1.5*i;

      TString ts;
      ts += Form("%.1f", thetamin);
      ts += "#circ < #theta < ";
      ts += Form("%.1f", thetamax);
      ts += "#circ;M_{#pi^{+}#pi^{-}#pi^{0}} (GeV)";

      locHist_3Pi->SetTitle(ts);
     
      TF1 *gauseta = new TF1("gauseta", "gaus", 0.4, 1.3);
      gauseta->SetParameters(100., 0.55, 0.1);

      TF1 *gausphi = new TF1("gausphi", "gaus", 0.4, 1.3);
      gausphi->SetParameters(100., 1., 0.1);

      TF1 *gausomegafit = new TF1("gausomegafit", "gaus", 0.4, 1.3);
      gausomegafit->SetParameters(100., 0.8, 0.1);

      TFitResultPtr reta = locHist_3Pi->Fit(gauseta, "SQ", "", 0.5, 0.6);
      TFitResultPtr rphi = locHist_3Pi->Fit(gausphi, "SQ", "", 0.95, 1.05);
      TFitResultPtr romega = locHist_3Pi->Fit(gausomegafit, "SQ", "", 0.7, 0.9);
     
      TF1* gaus4pol3 = new TF1("gaus4pol3", fit4gaus, 0.4, 1.3, 15);

      gaus4pol3->SetParameter(0, romega->Parameter(0));
      gaus4pol3->SetParameter(1, romega->Parameter(1));
      gaus4pol3->SetParameter(2, romega->Parameter(2));
      gaus4pol3->SetParameter(3, romega->Parameter(0));
      gaus4pol3->SetParameter(4, 2.0 * romega->Parameter(2));
      if(reta.Get() == nullptr){
        gaus4pol3->SetParameter(5, 100);
        gaus4pol3->FixParameter(6, .55);
        gaus4pol3->SetParameter(7, 0.05);
      }
      else{
        gaus4pol3->SetParameter(5, reta->Parameter(0));
        gaus4pol3->FixParameter(6, reta->Parameter(1));
        gaus4pol3->SetParameter(7, reta->Parameter(2));
      }
      gaus4pol3->SetParameter(8, rphi->Parameter(0));
      gaus4pol3->FixParameter(9, rphi->Parameter(1));
      gaus4pol3->SetParameter(10, rphi->Parameter(2));
      gaus4pol3->SetParameter(11, 1.);
      gaus4pol3->SetParameter(12, 1.);
      gaus4pol3->SetParameter(13, 1.);
      gaus4pol3->SetParameter(14, 1.);
     
      gaus4pol3->SetParLimits(0, 0, 1000000);
      gaus4pol3->SetParLimits(1, 0.7, 0.9);
      gaus4pol3->SetParLimits(2, 0.01, 0.1);
      gaus4pol3->SetParLimits(3, 0, 1000000);
      gaus4pol3->SetParLimits(4, 0.01, 0.1);
      gaus4pol3->SetParLimits(5, 0, 10000);
      gaus4pol3->SetParLimits(7, 0.01, 0.1);
      gaus4pol3->SetParLimits(8, 0, 10000);
      gaus4pol3->SetParLimits(10, 0.01, 0.03);
     
      gaus4pol3->SetParName(0, "#omega_{1}");
      gaus4pol3->SetParName(1, "#bar{x}_{#omega}");
      gaus4pol3->SetParName(2, "#sigma_{#omega_{1}}");
      gaus4pol3->SetParName(3, "#omega_{2}");
      gaus4pol3->SetParName(4, "#sigma_{#omega_{2}}");
      gaus4pol3->SetParName(5, "#eta");
      gaus4pol3->SetParName(6, "#eta mean");
      gaus4pol3->SetParName(7, "#sigma_{#eta}");
      gaus4pol3->SetParName(8, "#phi");
      gaus4pol3->SetParName(9, "#phi mean");
      gaus4pol3->SetParName(10, "#sigma_{#phi}");
      gaus4pol3->SetParName(11, "const");
      gaus4pol3->SetParName(12, "lin");
      gaus4pol3->SetParName(13, "quad");
      gaus4pol3->SetParName(14, "cube");
     
      TFitResultPtr r;
      //if(i == 7)
      r = locHist_3Pi->Fit(gaus4pol3, "SQ");
      //else 
      //r = locHist_3Pi->Fit(gaus4pol3, "SQ");


      TMatrixD covMatrix(5,5);
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<5; i++){
	for (Int_t j=0; j<5; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }
     
      Double_t params[5] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4)};

      TF1 *gausomega = new TF1("gausomega", omega2gaus, 0.4, 1.3, 5);
      gausomega->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4));
      double intomega = gausomega->Integral(0.4, 1.3);
      intomega = intomega/binsize;
      double error = gausomega->IntegralError(0.4, 1.3, params, covMatrix.GetMatrixArray(), 0.0001);
      error = error/binsize;
     

      TF1 *gaus1 = new TF1("gaus1", "gaus", 0.4, 1.3);
      gaus1->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2));
     
      TF1 *gaus2 = new TF1("gaus2", "gaus", 0.4, 1.3);
      gaus2->SetParameters(r->Parameter(3), r->Parameter(1), r->Parameter(4));

      TF1 *polBG = new TF1("polBG", "pol3", 0.4, 1.3);
      polBG->SetParameters(r->Parameter(11), r->Parameter(12), r->Parameter(13), r->Parameter(14));

 
      gauseta->SetParameters(r->Parameter(5), r->Parameter(6), r->Parameter(7));
      
      gausphi->SetParameters(r->Parameter(8), r->Parameter(9), r->Parameter(10));

    
      gaus1->SetLineColor(kBlue);
      gaus2->SetLineColor(kAzure+10);
      polBG->SetLineColor(kViolet);
      gauseta->SetLineColor(kGreen);
      gausphi->SetLineColor(kGreen);
    
      locHist_3Pi->SetMarkerStyle(kFullCircle);
      locHist_3Pi->Draw();
      polBG->Draw("same");
      gauseta->Draw("same");
      gausphi->Draw("same");
      gaus1->Draw("same");
      gaus2->Draw("same");
      gaus4pol3->Draw("same");
    
      double intHist = locHist_3Pi->Integral();

      //cout << i+1 << "th theta (inv) bin:" << endl;
      //cout << "Integral under omega peak (from fit) = " << intomega << endl;
      //cout << "Error = " << error << endl;
      //cout << "Integral under omega peak (from hist) = " << intHist << endl;
      //cout << "Error (from sqrt(N)) = " << TMath::Sqrt(intHist) << endl;



      locHist_yield_theta[im]->SetBinContent(i+1, intomega);
      locHist_yield_theta[im]->SetBinError(i+1, error);
      locHist_yield_theta_denom[im]->SetBinContent(i+1, intomega);
      //theta_error[i] = error;
    }
   }
  
  //loop over p bins
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 12; i++) {
      cp[im]->cd(i+1);
      //cout << i+1 << "th p bin:" << endl;
      TH1F* locHist_3Pi = (TH1F*)locHist_3Pi_p2D[im]->ProjectionX(Form("p_%d", i), i+1, i+1);

      double pmin = 0. + 0.5*i;
      double pmax = 0.5 + 0.5*i;

      TString ts;
      ts += Form("%.1f", pmin);
      ts += " < p < ";
      ts += Form("%.1f", pmax);
      ts += " GeV;M_{#pi^{+}#pi^{-}#pi^{0}} (GeV)";

      locHist_3Pi->SetTitle(ts);
     
      TF1 *gauseta = new TF1("gauseta", "gaus", 0.4, 1.3);
      gauseta->SetParameters(100., 0.55, 0.1);

      TF1 *gausphi = new TF1("gausphi", "gaus", 0.4, 1.3);
      gausphi->SetParameters(100., 1., 0.1);

      TF1 *gausomegafit = new TF1("gausomegafit", "gaus", 0.4, 1.3);
      gausomegafit->SetParameters(100., 0.8, 0.1);

      TFitResultPtr reta = locHist_3Pi->Fit(gauseta, "SQ", "", 0.5, 0.6);
      TFitResultPtr rphi = locHist_3Pi->Fit(gausphi, "SQ", "", 0.95, 1.05);
      TFitResultPtr romega = locHist_3Pi->Fit(gausomegafit, "SQ", "", 0.7, 0.9);
     
      TF1* gaus4pol3 = new TF1("gaus4pol3", fit4gaus, 0.4, 1.3, 15);

      gaus4pol3->SetParameter(0, romega->Parameter(0));
      gaus4pol3->SetParameter(1, romega->Parameter(1));
      gaus4pol3->SetParameter(2, romega->Parameter(2));
      gaus4pol3->SetParameter(3, romega->Parameter(0));
      gaus4pol3->SetParameter(4, 2.0 * romega->Parameter(2));
      gaus4pol3->SetParameter(5, reta->Parameter(0));
      gaus4pol3->FixParameter(6, reta->Parameter(1));
      gaus4pol3->SetParameter(7, reta->Parameter(2));
      gaus4pol3->SetParameter(8, rphi->Parameter(0));
      gaus4pol3->FixParameter(9, rphi->Parameter(1));
      gaus4pol3->SetParameter(10, rphi->Parameter(2));
      gaus4pol3->SetParameter(11, 1.);
      gaus4pol3->SetParameter(12, 1.);
      gaus4pol3->SetParameter(13, 1.);
      gaus4pol3->SetParameter(14, 1.);
      
      gaus4pol3->SetParLimits(0, 0, 1000000);
      gaus4pol3->SetParLimits(1, 0.7, 0.9);
      gaus4pol3->SetParLimits(2, 0.01, 0.1);
      gaus4pol3->SetParLimits(3, 0, 1000000);
      gaus4pol3->SetParLimits(4, 0.01, 0.1);
      gaus4pol3->SetParLimits(5, 0, 10000);
      gaus4pol3->SetParLimits(7, 0.01, 0.1);
      gaus4pol3->SetParLimits(8, 0, 10000);
      gaus4pol3->SetParLimits(10, 0.01, 0.03);
      
      gaus4pol3->SetParName(0, "#omega_{1}");
      gaus4pol3->SetParName(1, "#bar{x}_{#omega}");
      gaus4pol3->SetParName(2, "#sigma_{#omega_{1}}");
      gaus4pol3->SetParName(3, "#omega_{2}");
      gaus4pol3->SetParName(4, "#sigma_{#omega_{2}}");
      gaus4pol3->SetParName(5, "#eta");
      gaus4pol3->SetParName(6, "#eta mean");
      gaus4pol3->SetParName(7, "#sigma_{#eta}");
      gaus4pol3->SetParName(8, "#phi");
      gaus4pol3->SetParName(9, "#phi mean");
      gaus4pol3->SetParName(10, "#sigma_{#phi}");
      gaus4pol3->SetParName(11, "const");
      gaus4pol3->SetParName(12, "lin");
      gaus4pol3->SetParName(13, "quad");
      gaus4pol3->SetParName(14, "cube");
      
      TFitResultPtr r;
      //if(i == 7)
      r = locHist_3Pi->Fit(gaus4pol3, "SQ");
      //else 
      //r = locHist_3Pi->Fit(gaus4pol3, "SQ");


      TMatrixD covMatrix(5,5);
      TMatrixDSym fCovar = r->GetCovarianceMatrix();
      for (Int_t i=0; i<5; i++){
	for (Int_t j=0; j<5; j++){
	  covMatrix[i][j] = fCovar[i][j];
	}
      }
      
      Double_t params[5] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4)};

      TF1 *gausomega = new TF1("gausomega", omega2gaus, 0.4, 1.3, 5);
      gausomega->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4));
      double intomega = gausomega->Integral(0.4, 1.3);
      intomega = intomega/binsize;
      double error = gausomega->IntegralError(0.4, 1.3, params, covMatrix.GetMatrixArray(), 0.0001);
      error = error/binsize;
     

      TF1 *gaus1 = new TF1("gaus1", "gaus", 0.4, 1.3);
      gaus1->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2));
     
      TF1 *gaus2 = new TF1("gaus2", "gaus", 0.4, 1.3);
      gaus2->SetParameters(r->Parameter(3), r->Parameter(1), r->Parameter(4));

      TF1 *polBG = new TF1("polBG", "pol3", 0.4, 1.3);
      polBG->SetParameters(r->Parameter(11), r->Parameter(12), r->Parameter(13), r->Parameter(14));

 
      gauseta->SetParameters(r->Parameter(5), r->Parameter(6), r->Parameter(7));
      
      gausphi->SetParameters(r->Parameter(8), r->Parameter(9), r->Parameter(10));

    
      gaus1->SetLineColor(kBlue);
      gaus2->SetLineColor(kAzure+10);
      polBG->SetLineColor(kViolet);
      gauseta->SetLineColor(kGreen);
      gausphi->SetLineColor(kGreen);
      
      locHist_3Pi->SetMarkerStyle(kFullCircle);
      locHist_3Pi->Draw();
      polBG->Draw("same");
      gauseta->Draw("same");
      gausphi->Draw("same");
      gaus1->Draw("same");
      gaus2->Draw("same");
      gaus4pol3->Draw("same");
    

      //cout << "Integral under omega peak (from fit) = " << intomega << endl;
      //cout << "Error = " << error << endl;

      locHist_yield_p[im]->SetBinContent(i+1, intomega);
      locHist_yield_p[im]->SetBinError(i+1, error);
      locHist_yield_p_denom[im]->SetBinContent(i+1, intomega);
      //p_error[i] = error;
    }
  }
  
  //loop over phi bins
  for(int i = 0; i < 20; i++) {
    cphi_diff->cd(i+1);
    TH1F* locHist_diff = (TH1F*)locHist_diff_phi2D->ProjectionX(Form("phi_diff_%d", i), i+1, i+1);

    double phimin = -180. + 18.*i;
    double phimax = -162. + 18.*i;

     TString ts;
     ts += phimin;
     ts += "#circ < #phi < ";
     ts += phimax;
     ts += "#circ;Missing Mass off the Proton (GeV)";

     locHist_diff->SetTitle(ts);

    TFitResultPtr r1 = locHist_diff->Fit("gaus", "SQ", "", 0.7, 0.9);
    TF1* gaus1pol3 = new TF1("gaus1pol3", fit1gaus, 0.3, 1.3, 7);
    gaus1pol3->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), 1., 1., 1., 1.);
    gaus1pol3->SetParLimits(0, 0, 1000000);
    gaus1pol3->SetParLimits(1, 0.7, 0.9);
    gaus1pol3->SetParLimits(2, 0, 0.1);

    TFitResultPtr r = locHist_diff->Fit(gaus1pol3, "SQ");
    TMatrixD covMatrix(3,3);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<3; i++){
      for (Int_t j=0; j<3; j++){
	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[3] = {r->Parameter(0), r->Parameter(1), r->Parameter(2)};

    TF1 *gaus1 = new TF1("gaus1", "gaus", 0.3, 1.3);
    gaus1->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2));
    double intgaus1 = gaus1->Integral(0.3, 1.3);
    intgaus1 = intgaus1/binsize;
    double error = gaus1->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
    polBG->SetParameters(r->Parameter(3), r->Parameter(4), r->Parameter(5), r->Parameter(6));
    
    gaus1->SetLineColor(kBlue);
    polBG->SetLineColor(kViolet);




    locHist_diff->Draw();
    polBG->Draw("same");
    gaus1->Draw("same");
    gaus1pol3->Draw("same");

    //cout << "Integral under omega peak (from fit) = " << intgaus1 << endl;
    //cout << "Error = " << error << endl;

    double prop_err = TMath::Sqrt(phi_error[i]*phi_error[i] + error*error);

    for(int im = 0; im < 7; im++){
      locHist_yield_phi_denom[im]->AddBinContent(i+1, intgaus1);
      locHist_yield_phi_denom[im]->SetBinError(i+1, error);
    }
  }

  cout << "Theta fits (0 candidates):" << endl;
  //loop over theta bins
  for(int i = 0; i < 20; i++) {
    ctheta_diff->cd(i+1);
    TH1F* locHist_diff = (TH1F*)locHist_diff_theta2D->ProjectionX(Form("theta_diff_%d", i), i+1, i+1);
      double hsum = locHist_diff->Integral();
      if(hsum < 75.)
	continue;

     double thetamin = 0. + 1.5*i;
     double thetamax = 1.5 + 1.5*i;

     TString ts;
     ts += Form("%.1f", thetamin);
     ts += "#circ < #theta < ";
     ts += Form("%.1f", thetamax);
     ts += "#circ;Missing Mass off the Proton (GeV)";

     locHist_diff->SetTitle(ts);

    TFitResultPtr r1 = locHist_diff->Fit("gaus", "SQ", "", 0.7, 0.9);
    TF1* gaus1pol3 = new TF1("gaus1pol3", fit1gaus, 0.3, 1.3, 7);
    gaus1pol3->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), 1., 1., 1., 1.);
    gaus1pol3->SetParLimits(0, 0, 1000000);
    gaus1pol3->SetParLimits(1, 0.7, 0.9);
    gaus1pol3->SetParLimits(2, 0, 0.1);

    TFitResultPtr r = locHist_diff->Fit(gaus1pol3, "SQ");
    TMatrixD covMatrix(3,3);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<3; i++){
      for (Int_t j=0; j<3; j++){
  	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[3] = {r->Parameter(0), r->Parameter(1), r->Parameter(2)};

    TF1 *gaus1 = new TF1("gaus1", "gaus", 0.3, 1.3);
    gaus1->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2));
    double intgaus1 = gaus1->Integral(0.3, 1.3);
    intgaus1 = intgaus1/binsize;
    double error = gaus1->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
    polBG->SetParameters(r->Parameter(3), r->Parameter(4), r->Parameter(5), r->Parameter(6));
    
    gaus1->SetLineColor(kBlue);
    polBG->SetLineColor(kViolet);

    locHist_diff->Draw();
    polBG->Draw("same");
    gaus1->Draw("same");
    gaus1pol3->Draw("same");

    double intHist = locHist_diff->Integral();

    cout << "Integral under omega peak (from fit) = " << intgaus1 << endl;
    cout << "Error = " << error << endl;
    //cout << "Integral under omega peak (from hist) = " << intHist << endl;
    //cout << "Error (from sqrt(N)) = " << TMath::Sqrt(intHist) << endl;

    double prop_err = TMath::Sqrt(theta_error[i]*theta_error[i] + error*error);

    for(int im = 0; im < 7; im++){
      locHist_yield_theta_denom[im]->AddBinContent(i+1, intgaus1);
      locHist_yield_theta_denom[im]->SetBinError(i+1, error);
    }
  }
  
  //loop over p bins
  for(int i = 0; i < 12; i++) {
    cp_diff->cd(i+1);
    TH1F* locHist_diff = (TH1F*)locHist_diff_p2D->ProjectionX(Form("p_diff_%d", i), i+1, i+1);

     double pmin = 0. + 0.5*i;
     double pmax = 0.5 + 0.5*i;

     TString ts;
     ts += Form("%.1f", pmin);
     ts += " < p < ";
     ts += Form("%.1f", pmax);
     ts += " GeV;Missing Mass off the Proton (GeV)";

     locHist_diff->SetTitle(ts);

    TFitResultPtr r1 = locHist_diff->Fit("gaus", "SQ", "", 0.7, 0.9);
    TF1* gaus1pol3 = new TF1("gaus1pol3", fit1gaus, 0.3, 1.3, 7);
    gaus1pol3->SetParameters(r1->Parameter(0), r1->Parameter(1), r1->Parameter(2), 1., 1., 1., 1.);
    gaus1pol3->SetParLimits(0, 0, 1000000);
    gaus1pol3->SetParLimits(1, 0.7, 0.9);
    gaus1pol3->SetParLimits(2, 0, 0.1);

    TFitResultPtr r = locHist_diff->Fit(gaus1pol3, "SQ");
    TMatrixD covMatrix(3,3);
 
    TMatrixDSym fCovar = r->GetCovarianceMatrix();
    for (Int_t i=0; i<3; i++){
      for (Int_t j=0; j<3; j++){
  	covMatrix[i][j] = fCovar[i][j];
      }
    }

    Double_t params[3] = {r->Parameter(0), r->Parameter(1), r->Parameter(2)};

    TF1 *gaus1 = new TF1("gaus1", "gaus", 0.3, 1.3);
    gaus1->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2));
    double intgaus1 = gaus1->Integral(0.3, 1.3);
    intgaus1 = intgaus1/binsize;
    double error = gaus1->IntegralError(0.3, 1.3, params, covMatrix.GetMatrixArray(), .0001);
    error = error/binsize;

    TF1 *polBG = new TF1("polBG", "pol3", 0.3, 1.3);
    polBG->SetParameters(r->Parameter(3), r->Parameter(4), r->Parameter(5), r->Parameter(6));
    
    gaus1->SetLineColor(kBlue);
    polBG->SetLineColor(kViolet);

    locHist_diff->Draw();
    polBG->Draw("same");
    gaus1->Draw("same");
    gaus1pol3->Draw("same");

    //cout << "Integral under omega peak (from fit) = " << intgaus1 << endl;
    //cout << "Error = " << error << endl;

    for(int im = 0; im < 7; im++){
      locHist_yield_p_denom[im]->AddBinContent(i+1, intgaus1);
      locHist_yield_p_denom[im]->SetBinError(i+1, error);
    }
  }
  

  /************************************************** FIT MASSES BINNED IN THETA, P **********************************/
  
  TH2F* locHist_theta2D_3pi_p[9][7];
  TH2F* locHist_theta2D_0track_p[9];
  for(int i = 0; i < 9; i++){
    for(int im = 0; im < 7; im++){
      locHist_theta2D_3pi_p[i][im] = (TH2F*)fhist->Get(Form("3PiMassVsTheta_%scut%d_p%d", whichcut, im+1, i+1));
      locHist_theta2D_3pi_p[i][im]->Sumw2();
      locHist_theta2D_3pi_p[i][im]->RebinX(10);
      locHist_theta2D_3pi_p[i][im]->RebinY(30);
    }
    locHist_theta2D_0track_p[i] = (TH2F*)fhist->Get(Form("OmegaMassVsTheta_0track_p%d", i+1));
    locHist_theta2D_0track_p[i]->Sumw2();
    locHist_theta2D_0track_p[i]->RebinX(10);
    locHist_theta2D_0track_p[i]->RebinY(30);
  }

  TCanvas* cthetap3pi[9][7];
  TCanvas* cthetap0[9];
  for(int i = 0; i < 9; i++){
    for(int im = 0; im < 7; im++){
      cthetap3pi[i][im] = new TCanvas(Form("cthetap3pi%d_m%d", i+1, im+1), Form("3 Pion Mass vs Theta in %dth p Bin", i+1), 1);
      cthetap3pi[i][im]->Divide(5,4);
    }
    cthetap0[i] = new TCanvas(Form("cthetap0%d", i+1), Form("Omega Mass vs Theta in %dth p Bin (0 tracks)", i+1), 1);
    cthetap0[i]->Divide(5,4);
  }

  TH1F* locHist_yield_theta_p[9][7];
  TH1F* locHist_yield_theta_denom_p[9][7];
  for(int i = 0; i < 9; i++){
    for(int im = 0; im < 7; im++){
      locHist_yield_theta_p[i][im] = new TH1F(Form("theta_yield_%d_m%d", i+1, im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
      locHist_yield_theta_denom_p[i][im] = new TH1F(Form("theta_yield_denom_%d_m%d", i+1, im+1), ";#theta_{mmop} (deg);#omega yield / 1.5#circ", 20, 0., 30.);
      locHist_yield_theta_p[i][im]->Sumw2();
      locHist_yield_theta_denom_p[i][im]->Sumw2();
    }
  }

  //Fit theta (3 pions):
  for(int ip = 0; ip < 9; ip++){
    for(int im = 0; im < 7; im++){
      for(int i = 0; i < 20; i++){
	cthetap3pi[ip][im]->cd(i+1);
	TH1F* locHist_3Pi = (TH1F*)locHist_theta2D_3pi_p[ip][im]->ProjectionX(Form("theta3pi_p%d_%d", ip, i), i+1, i+1);
	double hsum = locHist_3Pi->Integral();
	if(hsum < 75.)
	  continue;

	double thetamin = 0. + 1.5*i;
	double thetamax = 1.5 + 1.5*i;

	TString ts;
	ts += Form("%.1f", thetamin);
	ts += "#circ < #theta < ";
	ts += Form("%.1f", thetamax);
	ts += "#circ;M_{#pi^{+}#pi^{-}#pi^{0}} (GeV)";
	
	locHist_3Pi->SetTitle(ts);

	TF1 *gauseta = new TF1("gauseta", "gaus", 0.4, 1.3);
	gauseta->SetParameters(100., 0.55, 0.1);
	
	TF1 *gausphi = new TF1("gausphi", "gaus", 0.4, 1.3);
	gausphi->SetParameters(100., 1., 0.1);
      
	TF1 *gausomegafit = new TF1("gausomegafit", "gaus", 0.4, 1.3);
	gausomegafit->SetParameters(100., 0.8, 0.1);

	TFitResultPtr reta = locHist_3Pi->Fit(gauseta, "SQ", "", 0.5, 0.6);
	TFitResultPtr rphi = locHist_3Pi->Fit(gausphi, "SQ", "", 0.95, 1.05);
	TFitResultPtr romega = locHist_3Pi->Fit(gausomegafit, "SQ", "", 0.7, 0.9);
	
	TF1* gaus4pol3 = new TF1("gaus4pol3", fit4gaus, 0.4, 1.3, 15);

	gaus4pol3->SetParameter(0, romega->Parameter(0));
	gaus4pol3->SetParameter(1, romega->Parameter(1));
	gaus4pol3->SetParameter(2, romega->Parameter(2));
	gaus4pol3->SetParameter(3, romega->Parameter(0));
	gaus4pol3->SetParameter(4, 2.0 * romega->Parameter(2));
	if(reta.Get() == nullptr){
	  gaus4pol3->SetParameter(5, 100);
	  gaus4pol3->FixParameter(6, .55);
	  gaus4pol3->SetParameter(7, 0.05);
	}
	else{
	  gaus4pol3->SetParameter(5, reta->Parameter(0));
	  gaus4pol3->FixParameter(6, reta->Parameter(1));	 
	  gaus4pol3->SetParameter(7, reta->Parameter(2));
	}
	if(rphi.Get() == nullptr){
  	  gaus4pol3->SetParameter(8, 100);
	  gaus4pol3->FixParameter(9, 1.02);
	  gaus4pol3->SetParameter(10, 0.05);
	}
	else{
	  gaus4pol3->SetParameter(8, rphi->Parameter(0));
          gaus4pol3->FixParameter(9, rphi->Parameter(1));
          gaus4pol3->SetParameter(10, rphi->Parameter(2));
	}
	gaus4pol3->SetParameter(11, 1.);
	gaus4pol3->SetParameter(12, 1.);
	gaus4pol3->SetParameter(13, 1.);
	gaus4pol3->SetParameter(14, 1.);
     
	gaus4pol3->SetParLimits(0, 0, 1000000);
	gaus4pol3->SetParLimits(1, 0.7, 0.9);
	gaus4pol3->SetParLimits(2, 0.01, 0.1);
	gaus4pol3->SetParLimits(3, 0, 1000000);
	gaus4pol3->SetParLimits(4, 0.01, 0.1);
	gaus4pol3->SetParLimits(5, 0, 10000);
	gaus4pol3->SetParLimits(7, 0.01, 0.1);
	gaus4pol3->SetParLimits(8, 0, 10000);
	gaus4pol3->SetParLimits(10, 0.01, 0.03);
	
	gaus4pol3->SetParName(0, "#omega_{1}");
	gaus4pol3->SetParName(1, "#bar{x}_{#omega}");
	gaus4pol3->SetParName(2, "#sigma_{#omega_{1}}");
	gaus4pol3->SetParName(3, "#omega_{2}");
	gaus4pol3->SetParName(4, "#sigma_{#omega_{2}}");
	gaus4pol3->SetParName(5, "#eta");
	gaus4pol3->SetParName(6, "#eta mean");
	gaus4pol3->SetParName(7, "#sigma_{#eta}");
	gaus4pol3->SetParName(8, "#phi");
	gaus4pol3->SetParName(9, "#phi mean");
	gaus4pol3->SetParName(10, "#sigma_{#phi}");
	gaus4pol3->SetParName(11, "const");
	gaus4pol3->SetParName(12, "lin");
	gaus4pol3->SetParName(13, "quad");
	gaus4pol3->SetParName(14, "cube");
      
	TFitResultPtr r;
	r = locHist_3Pi->Fit(gaus4pol3, "SQ");

	TMatrixD covMatrix(5,5);
	TMatrixDSym fCovar = r->GetCovarianceMatrix();
	for (Int_t k=0; k<5; k++){
	  for (Int_t j=0; j<5; j++){
	    covMatrix[k][j] = fCovar[k][j];
	  }
	}
     
	Double_t params[5] = {r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4)};
	
	TF1 *gausomega = new TF1("gausomega", omega2gaus, 0.4, 1.3, 5);
	gausomega->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2), r->Parameter(3), r->Parameter(4));
	double intomega = gausomega->Integral(0.4, 1.3);
	intomega = intomega/binsize;
	double error = gausomega->IntegralError(0.4, 1.3, params, covMatrix.GetMatrixArray(), 0.0001);
	if(error == 0)
	  cout << "WARNING: IntegralError failed for p-bin = " << ip << ", hist bin = " << i << ". (3Pi fit)" << endl;
	error = error/binsize;

	TF1 *gaus1 = new TF1("gaus1", "gaus", 0.4, 1.3);
	gaus1->SetParameters(r->Parameter(0), r->Parameter(1), r->Parameter(2));
      
	TF1 *gaus2 = new TF1("gaus2", "gaus", 0.4, 1.3);
	gaus2->SetParameters(r->Parameter(3), r->Parameter(1), r->Parameter(4));
	
	TF1 *polBG = new TF1("polBG", "pol3", 0.4, 1.3);
	polBG->SetParameters(r->Parameter(11), r->Parameter(12), r->Parameter(13), r->Parameter(14));
	
	gauseta->SetParameters(r->Parameter(5), r->Parameter(6), r->Parameter(7));
      
	gausphi->SetParameters(r->Parameter(8), r->Parameter(9), r->Parameter(10));

	gaus1->SetLineColor(kBlue);
	gaus2->SetLineColor(kAzure+10);
	polBG->SetLineColor(kViolet);
	gauseta->SetLineColor(kGreen);
	gausphi->SetLineColor(kGreen);
    
	locHist_3Pi->SetMarkerStyle(kFullCircle);
	locHist_3Pi->Draw();
	polBG->Draw("same");
	gauseta->Draw("same");
	gausphi->Draw("same");
	gaus1->Draw("same");
	gaus2->Draw("same");
	gaus4pol3->Draw("same");
	
	locHist_yield_theta_p[ip][im]->SetBinContent(i+1, intomega);
	locHist_yield_theta_p[ip][im]->SetBinError(i+1, error);
	locHist_yield_theta_denom_p[ip][im]->SetBinContent(i+1, intomega);
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
      }
    }
  }
  

  
  TCanvas* c1[7];

  for(int im = 0; im < 7; im++){
    c1[im] = new TCanvas(Form("c1_%d", im+1), "omega yield", 900, 300);
    c1[im]->Divide(3,1);
    c1[im]->cd(1);
    locHist_yield_phi_denom[im]->SetMinimum(0);
    locHist_yield_phi_denom[im]->SetMarkerStyle(kFullCircle);
    locHist_yield_phi_denom[im]->SetMarkerColor(kAzure);
    locHist_yield_phi_denom[im]->SetLineColor(kAzure);
    locHist_yield_phi[im]->SetMarkerStyle(kOpenCircle);
    locHist_yield_phi[im]->SetMarkerColor(kOrange+7);
    locHist_yield_phi[im]->SetLineColor(kOrange+7);
    locHist_yield_phi_denom[im]->Draw();
    locHist_yield_phi[im]->Draw("same");
    c1[im]->cd(2);
    locHist_yield_theta_denom[im]->SetMinimum(0);
    locHist_yield_theta_denom[im]->SetMarkerStyle(kFullCircle);
    locHist_yield_theta_denom[im]->SetMarkerColor(kAzure);
    locHist_yield_theta_denom[im]->SetLineColor(kAzure);
    locHist_yield_theta[im]->SetMarkerStyle(kOpenCircle);
    locHist_yield_theta[im]->SetMarkerColor(kOrange+7);
    locHist_yield_theta[im]->SetLineColor(kOrange+7);
    locHist_yield_theta_denom[im]->Draw();
    locHist_yield_theta[im]->Draw("same");
    c1[im]->cd(3);
    locHist_yield_p_denom[im]->SetMinimum(0);
    locHist_yield_p_denom[im]->SetMarkerStyle(kFullCircle);
    locHist_yield_p_denom[im]->SetMarkerColor(kAzure);
    locHist_yield_p_denom[im]->SetLineColor(kAzure);
    locHist_yield_p[im]->SetMarkerStyle(kOpenCircle);
    locHist_yield_p[im]->SetMarkerColor(kOrange+7);
    locHist_yield_p[im]->SetLineColor(kOrange+7);
    locHist_yield_p_denom[im]->Draw();
    locHist_yield_p[im]->Draw("same");
  }

  TCanvas* c2[7];

  for(int im = 0; im < 7; im++){
    c2[im] = new TCanvas(Form("c2_%d", im+1), "Omega yield vs theta, p", 900, 900);
    c2[im]->Divide(3, 3);
    for(int i = 0; i < 9; i++){
      c2[im]->cd(i+1);
      locHist_yield_theta_denom_p[i][im]->SetMinimum(0);
      locHist_yield_theta_denom_p[i][im]->SetMarkerStyle(kFullCircle);
      locHist_yield_theta_denom_p[i][im]->SetMarkerColor(kAzure);
      locHist_yield_theta_denom_p[i][im]->SetLineColor(kAzure);
      locHist_yield_theta_p[i][im]->SetMarkerStyle(kOpenCircle);
      locHist_yield_theta_p[i][im]->SetMarkerColor(kOrange+7);
      locHist_yield_theta_p[i][im]->SetLineColor(kOrange+7);
      locHist_yield_theta_denom_p[i][im]->Draw();
      locHist_yield_theta_p[i][im]->Draw("same");
    }
  }

  


  if(save == true) {
    for(int im = 0; im < 7; im++){
      locHist_yield_phi[im]->Write();
      locHist_yield_theta[im]->Write();
      locHist_yield_p[im]->Write();
      locHist_yield_phi_denom[im]->Write();
      locHist_yield_theta_denom[im]->Write();
      locHist_yield_p_denom[im]->Write();
      for(int i = 0; i < 9; i++){
	locHist_yield_theta_p[i][im]->Write();
	locHist_yield_theta_denom_p[i][im]->Write();
      }
    }
  }

   return;
}
