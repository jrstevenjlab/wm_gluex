void efficiency(const char *whichcut = "pv", bool print = false, bool save = false){
  TFile *fdata = TFile::Open(Form("mmop_data_%scut.root", whichcut));
  TFile *fgen = TFile::Open(Form("mmop_g4_%scut.root", whichcut));

  TFile *fdata_m2 = TFile::Open(Form("method2_data_%scut.root", whichcut));
  TFile *fgen_m2 = TFile::Open(Form("method2_g4_%scut.root", whichcut));

  if(save == true)
    TFile *outROOT = new TFile(Form("piplus_efficiency_%scut.root", whichcut), "recreate");

  //New histograms of sigma_D and sigma_U
  //TH1F* locHist_sigmaD_theta_gen_m2 = new TH1F("sigmaD_theta_gen_m2", ";#theta (deg)", 20, 0., 30.);
  //TH1F* locHist_sigmaU_theta_gen_m2 = new TH1F("sigmaU_theta_gen_m2", ";#theta (deg)", 20, 0., 30.);

  //Method 1
  TH1F* locHist_phi_data[7];
  TH1F* locHist_phi_denom_data[7];
  double phi_errorD_data[20][7];
  double phi_data[20][7];
  double phi_errorU_data[20][7];
  double phi_denom_data[20][7];
  for(int im = 0; im < 7; im++){
    locHist_phi_data[im] = (TH1F*)fdata->Get(Form("phi_mmop_reco_%d", im+1));
    locHist_phi_denom_data[im] = (TH1F*)fdata->Get(Form("phi_mmop_denom_%d", im+1));
    for(int i = 0; i < 20; i++) {
      phi_errorD_data[i][im] = locHist_phi_data[im]->GetBinError(i+1);
      phi_data[i][im] = locHist_phi_data[im]->GetBinContent(i+1);
      phi_errorU_data[i][im] = locHist_phi_denom_data[im]->GetBinError(i+1);
      phi_denom_data[i][im] = locHist_phi_denom_data[im]->GetBinContent(i+1);
    }
    locHist_phi_data[im]->Sumw2();
    locHist_phi_denom_data[im]->Sumw2();
    locHist_phi_data[im]->Divide(locHist_phi_denom_data[im]);
  }

  
  TH1F* locHist_theta_num_data[7];
  TH1F* locHist_theta_data[7];
  double theta_errorD_data[20][7];
  double theta_data[20][7];
  TH1F* locHist_theta_denom_data[7];
  double theta_errorU_data[20][7];
  double theta_denom_data[20][7];
  for(int im = 0; im < 7; im++){
    locHist_theta_num_data[im] = (TH1F*)fdata->Get(Form("theta_mmop_reco_%d", im+1));
    locHist_theta_denom_data[im] = (TH1F*)fdata->Get(Form("theta_mmop_denom_%d", im+1));
    locHist_theta_data[im] = new TH1F(Form("theta_data_efficiency_%d", im+1), ";#theta (deg)", 20, 0., 30.);
    for(int i = 0; i < 20; i++) {
      if(locHist_theta_num_data[im]->GetBinContent(i+1) == locHist_theta_denom_data[im]->GetBinContent(i+1))
	continue;
      if(locHist_theta_num_data[im]->GetBinContent(i+1) == 0)
	continue;
      theta_errorD_data[i][im] = locHist_theta_num_data[im]->GetBinError(i+1);
      theta_data[i][im] = locHist_theta_num_data[im]->GetBinContent(i+1);
      theta_errorU_data[i][im] = locHist_theta_denom_data[im]->GetBinError(i+1);
      theta_denom_data[i][im] = locHist_theta_denom_data[im]->GetBinContent(i+1);
      double content = theta_data[i][im] / theta_denom_data[i][im];
      double error = TMath::Sqrt((TMath::Power((theta_denom_data[i][im] - theta_data[i][im])/(theta_denom_data[i][im]*theta_denom_data[i][im]), 2.)*theta_errorD_data[i][im]*theta_errorD_data[i][im] + TMath::Power(theta_data[i][im]/(theta_denom_data[i][im]*theta_denom_data[i][im]), 2.)*theta_errorU_data[i][im]*theta_errorU_data[i][im]));
      locHist_theta_data[im]->SetBinContent(i+1, content);
      locHist_theta_data[im]->SetBinError(i+1, error);
    }
    locHist_theta_data[im]->Sumw2();
  }

  TH1F* locHist_p_data[7];
  TH1F* locHist_p_denom_data[7];
  double p_errorD_data[20][7];
  double p_data[20][7];
  double p_errorU_data[20][7];
  double p_denom_data[20][7];
  for(int im = 0; im < 7; im++){
    locHist_p_data[im] = (TH1F*)fdata->Get(Form("p_mmop_reco_%d", im+1));
    locHist_p_denom_data[im] = (TH1F*)fdata->Get(Form("p_mmop_denom_%d", im+1));
    for(int i = 0; i < 20; i++) {
      p_errorD_data[i][im] = locHist_p_data[im]->GetBinError(i+1);
      p_data[i][im] = locHist_p_data[im]->GetBinContent(i+1);
      p_errorU_data[i][im] = locHist_p_denom_data[im]->GetBinError(i+1);
      p_denom_data[i][im] = locHist_p_denom_data[im]->GetBinContent(i+1);
    }
    locHist_p_data[im]->Sumw2();
    locHist_p_denom_data[im]->Sumw2();
    locHist_p_data[im]->Divide(locHist_p_denom_data[im]);
  }
  
  TH1F* locHist_phi_gen[7];
  TH1F* locHist_phi_denom_gen[7];
  double phi_errorD_gen[20][7];
  double phi_gen[20][7];
  double phi_errorU_gen[20][7];
  double phi_denom_gen[20][7];
  for(int im = 0; im < 7; im++){
    locHist_phi_gen[im] = (TH1F*)fgen->Get(Form("phi_mmop_reco_%d", im+1));
    locHist_phi_denom_gen[im] = (TH1F*)fgen->Get(Form("phi_mmop_denom_%d", im+1));
    for(int i = 0; i < 20; i++) {
      phi_errorD_gen[i][im] = locHist_phi_gen[im]->GetBinError(i+1);
      phi_gen[i][im] = locHist_phi_gen[im]->GetBinContent(i+1);
      phi_errorU_gen[i][im] = locHist_phi_denom_gen[im]->GetBinError(i+1);
      phi_denom_gen[i][im] = locHist_phi_denom_gen[im]->GetBinContent(i+1);
    }
    locHist_phi_gen[im]->Sumw2();
    locHist_phi_denom_gen[im]->Sumw2();
    locHist_phi_gen[im]->Divide(locHist_phi_denom_gen[im]);
  }
  
  TH1F* locHist_theta_num_gen[7];
  TH1F* locHist_theta_gen[7];
  double theta_errorD_gen[20][7];
  double theta_gen[20][7];
  TH1F* locHist_theta_denom_gen[7];
  double theta_errorU_gen[20][7];
  double theta_denom_gen[20][7];
  for(int im = 0; im < 7; im++){
    locHist_theta_num_gen[im] = (TH1F*)fgen->Get(Form("theta_mmop_reco_%d", im+1));
    locHist_theta_denom_gen[im] = (TH1F*)fgen->Get(Form("theta_mmop_denom_%d", im+1));
    locHist_theta_gen[im] = new TH1F(Form("theta_gen_efficiency_%d", im+1), ";#theta (deg)", 20, 0., 30.);
    for(int i = 0; i < 20; i++) {
      if(locHist_theta_num_gen[im]->GetBinContent(i+1) == locHist_theta_denom_gen[im]->GetBinContent(i+1))
	continue;
      if(locHist_theta_num_gen[im]->GetBinContent(i+1) == 0)
	continue;
      theta_errorD_gen[i][im] = locHist_theta_num_gen[im]->GetBinError(i+1);
      theta_gen[i][im] = locHist_theta_num_gen[im]->GetBinContent(i+1);
      theta_errorU_gen[i][im] = locHist_theta_denom_gen[im]->GetBinError(i+1);
      theta_denom_gen[i][im] = locHist_theta_denom_gen[im]->GetBinContent(i+1);
      double content = theta_gen[i][im] / theta_denom_gen[i][im];
      double error = TMath::Sqrt((TMath::Power((theta_denom_gen[i][im] - theta_gen[i][im])/(theta_denom_gen[i][im]*theta_denom_gen[i][im]), 2.)*theta_errorD_gen[i][im]*theta_errorD_gen[i][im] + TMath::Power(theta_gen[i][im]/(theta_denom_gen[i][im]*theta_denom_gen[i][im]), 2.)*theta_errorU_gen[i][im]*theta_errorU_gen[i][im]));
      locHist_theta_gen[im]->SetBinContent(i+1, content);
      locHist_theta_gen[im]->SetBinError(i+1, error);
    }
    locHist_theta_gen[im]->Sumw2();
  }

  TH1F* locHist_p_gen[7];
  TH1F* locHist_p_denom_gen[7];
  double p_errorD_gen[12][7];
  double p_gen[12][7];
  double p_errorU_gen[12][7];
  double p_denom_gen[12][7];
  for(int im = 0; im < 7; im++){
    locHist_p_gen[im] = (TH1F*)fgen->Get(Form("p_mmop_reco_%d", im+1));
    locHist_p_denom_gen[im] = (TH1F*)fgen->Get(Form("p_mmop_denom_%d", im+1));
    for(int i = 0; i < 12; i++) {
      p_errorD_gen[i][im] = locHist_p_gen[im]->GetBinError(i+1);
      p_gen[i][im] = locHist_p_gen[im]->GetBinContent(i+1);
      p_errorU_gen[i][im] = locHist_p_denom_gen[im]->GetBinError(i+1);
      p_denom_gen[i][im] = locHist_p_denom_gen[im]->GetBinContent(i+1);
    }
    locHist_p_gen[im]->Sumw2();
    locHist_p_denom_gen[im]->Sumw2();
    locHist_p_gen[im]->Divide(locHist_p_denom_gen[im]);
  }
  
  //Propagate errors:
  double phi_err_data[7];
  double theta_err_data[7];
  double p_err_data[7];
  double phi_err_gen[7];
  double theta_err_gen[7];
  double p_err_gen[7];
  
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 20; i++){
      
      phi_err_data[im] = TMath::Sqrt((TMath::Power((phi_denom_data[i][im] - phi_data[i][im])/(phi_denom_data[i][im]*phi_denom_data[i][im]), 2.)*phi_errorD_data[i][im]*phi_errorD_data[i][im] + TMath::Power(phi_data[i][im]/(phi_denom_data[i][im]*phi_denom_data[i][im]), 2.)*phi_errorU_data[i][im]*phi_errorU_data[i][im]));
      locHist_phi_data[im]->SetBinError(i+1, phi_err_data[im]);
      
      p_err_data[im] = TMath::Sqrt((TMath::Power((p_denom_data[i][im] - p_data[i][im])/(p_denom_data[i][im]*p_denom_data[i][im]), 2.)*p_errorD_data[i][im]*p_errorD_data[i][im] + TMath::Power(p_data[i][im]/(p_denom_data[i][im]*p_denom_data[i][im]), 2.)*p_errorU_data[i][im]*p_errorU_data[i][im]));
      locHist_p_data[im]->SetBinError(i+1, p_err_data[im]);

      phi_err_gen[im] = TMath::Sqrt((TMath::Power((phi_denom_gen[i][im] - phi_gen[i][im])/(phi_denom_gen[i][im]*phi_denom_gen[i][im]), 2.)*phi_errorD_gen[i][im]*phi_errorD_gen[i][im] + TMath::Power(phi_gen[i][im]/(phi_denom_gen[i][im]*phi_denom_gen[i][im]), 2.)*phi_errorU_gen[i][im]*phi_errorU_gen[i][im]));
      locHist_phi_gen[im]->SetBinError(i+1, phi_err_gen[im]);
   
      p_err_gen[im] = TMath::Sqrt((TMath::Power((p_denom_gen[i][im] - p_gen[i][im])/(p_denom_gen[i][im]*p_denom_gen[i][im]), 2.)*p_errorD_gen[i][im]*p_errorD_gen[i][im] + TMath::Power(p_gen[i][im]/(p_denom_gen[i][im]*p_denom_gen[i][im]), 2.)*p_errorU_gen[i][im]*p_errorU_gen[i][im]));
      locHist_p_gen[im]->SetBinError(i+1, p_err_gen[im]);
    }
  }
  
  //Find Data/MC ratio (Method 1):
  TH1F* locHist_phi_ratio[7];
  TH1F* locHist_theta_ratio[7];
  TH1F* locHist_p_ratio[7];
  for(int im = 0; im < 7; im++){
    locHist_phi_ratio[im] = (TH1F*)locHist_phi_data[im]->Clone(Form("phi_ratio_%d", im+1));
    locHist_phi_ratio[im]->Sumw2();
    locHist_phi_ratio[im]->Divide(locHist_phi_gen[im]);
    locHist_theta_ratio[im] = (TH1F*)locHist_theta_data[im]->Clone(Form("theta_ratio_%d", im+1));
    locHist_theta_ratio[im]->Sumw2();
    locHist_theta_ratio[im]->Divide(locHist_theta_gen[im]);
    locHist_p_ratio[im] = (TH1F*)locHist_p_data[im]->Clone(Form("p_ratio_%d", im+1));
    locHist_p_ratio[im]->Sumw2();
    locHist_p_ratio[im]->Divide(locHist_p_gen[im]);
  }
  

  //Method 2
  TH1F* locHist_phi_data_m2[7];
  TH1F* locHist_phi_denom_data_m2[7];
  double phi_errorD_data_m2[20][7];
  double phi_data_m2[20][7];
  double phi_errorU_data_m2[20][7];
  double phi_denom_data_m2[20][7];
  for(int im = 0; im < 7; im++){
    locHist_phi_data_m2[im] = (TH1F*)fdata_m2->Get(Form("phi_yield_%d", im+1));
    locHist_phi_denom_data_m2[im] = (TH1F*)fdata_m2->Get(Form("phi_yield_denom_%d", im+1));
    for(int i = 0; i < 20; i++) {
      phi_errorD_data_m2[i][im] = locHist_phi_data_m2[im]->GetBinError(i+1);
      phi_data_m2[i][im] = locHist_phi_data_m2[im]->GetBinContent(i+1);
      phi_errorU_data_m2[i][im] = locHist_phi_denom_data_m2[im]->GetBinError(i+1);
      phi_denom_data_m2[i][im] = locHist_phi_denom_data_m2[im]->GetBinContent(i+1);
    }
    locHist_phi_data_m2[im]->Sumw2();
    locHist_phi_denom_data_m2[im]->Sumw2();
    locHist_phi_data_m2[im]->Divide(locHist_phi_denom_data_m2[im]);
  }

  TH1F* locHist_theta_num_data_m2[7];
  TH1F* locHist_theta_data_m2[7];
  double theta_errorD_data_m2[20][7];
  double theta_data_m2[20][7];
  TH1F* locHist_theta_denom_data_m2[7];
  double theta_errorU_data_m2[20][7];
  double theta_denom_data_m2[20][7];
  for(int im = 0; im < 7; im++){
    locHist_theta_num_data_m2[im] = (TH1F*)fdata_m2->Get(Form("theta_yield_%d", im+1));
    locHist_theta_denom_data_m2[im] = (TH1F*)fdata_m2->Get(Form("theta_yield_denom_%d", im+1));
    locHist_theta_data_m2[im] = new TH1F(Form("theta_data_m2_efficiency_%d", im+1), ";#theta (deg)", 20, 0., 30.);
    for(int i = 0; i < 20; i++) {
      if(locHist_theta_num_data_m2[im]->GetBinContent(i+1) == locHist_theta_denom_data_m2[im]->GetBinContent(i+1))
	continue;
      if(locHist_theta_num_data_m2[im]->GetBinContent(i+1) == 0)
	continue;
      theta_errorD_data_m2[i][im] = locHist_theta_num_data_m2[im]->GetBinError(i+1);
      theta_data_m2[i][im] = locHist_theta_num_data_m2[im]->GetBinContent(i+1);
      theta_errorU_data_m2[i][im] = locHist_theta_denom_data_m2[im]->GetBinError(i+1);
      theta_denom_data_m2[i][im] = locHist_theta_denom_data_m2[im]->GetBinContent(i+1);
      double content = theta_data_m2[i][im] / theta_denom_data_m2[i][im];
      double error = TMath::Sqrt((TMath::Power((theta_denom_data_m2[i][im] - theta_data_m2[i][im])/(theta_denom_data_m2[i][im]*theta_denom_data_m2[i][im]), 2.)*theta_errorD_data_m2[i][im]*theta_errorD_data_m2[i][im] + TMath::Power(theta_data_m2[i][im]/(theta_denom_data_m2[i][im]*theta_denom_data_m2[i][im]), 2.)*theta_errorU_data_m2[i][im]*theta_errorU_data_m2[i][im]));
      locHist_theta_data_m2[im]->SetBinContent(i+1, content);
      locHist_theta_data_m2[im]->SetBinError(i+1, error);
    }
    locHist_theta_data_m2[im]->Sumw2();
  }


  TH1F* locHist_p_data_m2[7];
  TH1F* locHist_p_denom_data_m2[7];
  double p_errorD_data_m2[20][7];
  double p_data_m2[20][7];
  double p_errorU_data_m2[20][7];
  double p_denom_data_m2[20][7];
  for(int im = 0; im < 7; im++){
    locHist_p_data_m2[im] = (TH1F*)fdata_m2->Get(Form("p_yield_%d", im+1));
    locHist_p_denom_data_m2[im] = (TH1F*)fdata_m2->Get(Form("p_yield_denom_%d", im+1));
    for(int i = 0; i < 20; i++) {
      p_errorD_data_m2[i][im] = locHist_p_data_m2[im]->GetBinError(i+1);
      p_data_m2[i][im] = locHist_p_data_m2[im]->GetBinContent(i+1);
      p_errorU_data_m2[i][im] = locHist_p_denom_data_m2[im]->GetBinError(i+1);
      p_denom_data_m2[i][im] = locHist_p_denom_data_m2[im]->GetBinContent(i+1);
    }
    locHist_p_data_m2[im]->Sumw2();
    locHist_p_denom_data_m2[im]->Sumw2();
    locHist_p_data_m2[im]->Divide(locHist_p_denom_data_m2[im]);
  }

  TH1F* locHist_phi_gen_m2[7];
  TH1F* locHist_phi_denom_gen_m2[7];
  double phi_errorD_gen_m2[20][7];
  double phi_gen_m2[20][7];
  double phi_errorU_gen_m2[20][7];
  double phi_denom_gen_m2[20][7];
  for(int im = 0; im < 7; im++){
    locHist_phi_gen_m2[im] = (TH1F*)fgen_m2->Get(Form("phi_yield_%d", im+1));
    locHist_phi_denom_gen_m2[im] = (TH1F*)fgen_m2->Get(Form("phi_yield_denom_%d", im+1));
    for(int i = 0; i < 20; i++) {
      phi_errorD_gen_m2[i][im] = locHist_phi_gen_m2[im]->GetBinError(i+1);
      phi_gen_m2[i][im] = locHist_phi_gen_m2[im]->GetBinContent(i+1);
      phi_errorU_gen_m2[i][im] = locHist_phi_denom_gen_m2[im]->GetBinError(i+1);
      phi_denom_gen_m2[i][im] = locHist_phi_denom_gen_m2[im]->GetBinContent(i+1);
    }
    locHist_phi_gen_m2[im]->Sumw2();
    locHist_phi_denom_gen_m2[im]->Sumw2();
    locHist_phi_gen_m2[im]->Divide(locHist_phi_denom_gen_m2[im]);
  }

  TH1F* locHist_theta_num_gen_m2[7];
  TH1F* locHist_theta_gen_m2[7];
  double theta_errorD_gen_m2[20][7];
  double theta_gen_m2[20][7];
  TH1F* locHist_theta_denom_gen_m2[7];
  double theta_errorU_gen_m2[20][7];
  double theta_denom_gen_m2[20][7];
  for(int im = 0; im < 7; im++){
    locHist_theta_num_gen_m2[im] = (TH1F*)fgen_m2->Get(Form("theta_yield_%d", im+1));
    locHist_theta_denom_gen_m2[im] = (TH1F*)fgen_m2->Get(Form("theta_yield_denom_%d", im+1));
    locHist_theta_gen_m2[im] = new TH1F(Form("theta_gen_m2_efficiency_%d", im+1), ";#theta (deg)", 20, 0., 30.);
    for(int i = 0; i < 20; i++) {
      if(locHist_theta_num_gen_m2[im]->GetBinContent(i+1) == locHist_theta_denom_gen_m2[im]->GetBinContent(i+1))
	continue;
      if(locHist_theta_num_gen_m2[im]->GetBinContent(i+1) == 0)
	continue;
      theta_errorD_gen_m2[i][im] = locHist_theta_num_gen_m2[im]->GetBinError(i+1);
      theta_gen_m2[i][im] = locHist_theta_num_gen_m2[im]->GetBinContent(i+1);
      theta_errorU_gen_m2[i][im] = locHist_theta_denom_gen_m2[im]->GetBinError(i+1);
      theta_denom_gen_m2[i][im] = locHist_theta_denom_gen_m2[im]->GetBinContent(i+1);
      double content = theta_gen_m2[i][im] / theta_denom_gen_m2[i][im];
      double error = TMath::Sqrt((TMath::Power((theta_denom_gen_m2[i][im] - theta_gen_m2[i][im])/(theta_denom_gen_m2[i][im]*theta_denom_gen_m2[i][im]), 2.)*theta_errorD_gen_m2[i][im]*theta_errorD_gen_m2[i][im] + TMath::Power(theta_gen_m2[i][im]/(theta_denom_gen_m2[i][im]*theta_denom_gen_m2[i][im]), 2.)*theta_errorU_gen_m2[i][im]*theta_errorU_gen_m2[i][im]));
      locHist_theta_gen_m2[im]->SetBinContent(i+1, content);
      locHist_theta_gen_m2[im]->SetBinError(i+1, error);
    }
    locHist_theta_gen_m2[im]->Sumw2();
  }
  
  TH1F* locHist_p_gen_m2[7];
  TH1F* locHist_p_denom_gen_m2[7];
  double p_errorD_gen_m2[20][7];
  double p_gen_m2[20][7];
  double p_errorU_gen_m2[20][7];
  double p_denom_gen_m2[20][7];
  for(int im = 0; im < 7; im++){
    locHist_p_gen_m2[im] = (TH1F*)fgen_m2->Get(Form("p_yield_%d", im+1));
    locHist_p_denom_gen_m2[im] = (TH1F*)fgen_m2->Get(Form("p_yield_denom_%d", im+1));
    for(int i = 0; i < 20; i++) {
      p_errorD_gen_m2[i][im] = locHist_p_gen_m2[im]->GetBinError(i+1);
      p_gen_m2[i][im] = locHist_p_gen_m2[im]->GetBinContent(i+1);
      p_errorU_gen_m2[i][im] = locHist_p_denom_gen_m2[im]->GetBinError(i+1);
      p_denom_gen_m2[i][im] = locHist_p_denom_gen_m2[im]->GetBinContent(i+1);
    }
    locHist_p_gen_m2[im]->Sumw2();
    locHist_p_denom_gen_m2[im]->Sumw2();
    locHist_p_gen_m2[im]->Divide(locHist_p_denom_gen_m2[im]); 
  }
  
  //Propagate errors:
  double phi_err_data_m2[7];
  double theta_err_data_m2[7];
  double p_err_data_m2[7];
  double phi_err_gen_m2[7];
  double theta_err_gen_m2[7];
  double p_err_gen_m2[7];
  for(int im = 0; im < 7; im++){
    for(int i = 0; i < 20; i++) {
      phi_err_data_m2[im] = TMath::Sqrt((TMath::Power((phi_denom_data_m2[i][im] - phi_data_m2[i][im])/(phi_denom_data_m2[i][im]*phi_denom_data_m2[i][im]), 2.)*phi_errorD_data_m2[i][im]*phi_errorD_data_m2[i][im] + TMath::Power(phi_data_m2[i][im]/(phi_denom_data_m2[i][im]*phi_denom_data_m2[i][im]), 2.)*phi_errorU_data_m2[i][im]*phi_errorU_data_m2[i][im]));
      locHist_phi_data_m2[im]->SetBinError(i+1, phi_err_data_m2[im]);

      p_err_data_m2[im] = TMath::Sqrt((TMath::Power((p_denom_data_m2[i][im] - p_data_m2[i][im])/(p_denom_data_m2[i][im]*p_denom_data_m2[i][im]), 2.)*p_errorD_data_m2[i][im]*p_errorD_data_m2[i][im] + TMath::Power(p_data_m2[i][im]/(p_denom_data_m2[i][im]*p_denom_data_m2[i][im]), 2.)*p_errorU_data_m2[i][im]*p_errorU_data_m2[i][im]));
      locHist_p_data_m2[im]->SetBinError(i+1, p_err_data_m2[im]);

      phi_err_gen_m2[im] = TMath::Sqrt((TMath::Power((phi_denom_gen_m2[i][im] - phi_gen_m2[i][im])/(phi_denom_gen_m2[i][im]*phi_denom_gen_m2[i][im]), 2.)*phi_errorD_gen_m2[i][im]*phi_errorD_gen_m2[i][im] + TMath::Power(phi_gen_m2[i][im]/(phi_denom_gen_m2[i][im]*phi_denom_gen_m2[i][im]), 2.)*phi_errorU_gen_m2[i][im]*phi_errorU_gen_m2[i][im]));
      locHist_phi_gen_m2[im]->SetBinError(i+1, phi_err_gen_m2[im]);

      p_err_gen_m2[im] = TMath::Sqrt((TMath::Power((p_denom_gen_m2[i][im] - p_gen_m2[i][im])/(p_denom_gen_m2[i][im]*p_denom_gen_m2[i][im]), 2.)*p_errorD_gen_m2[i][im]*p_errorD_gen_m2[i][im] + TMath::Power(p_gen_m2[i][im]/(p_denom_gen_m2[i][im]*p_denom_gen_m2[i][im]), 2.)*p_errorU_gen_m2[i][im]*p_errorU_gen_m2[i][im]));
      locHist_p_gen_m2[im]->SetBinError(i+1, p_err_gen_m2[im]);
    }
  }
  
  //Repeat method 1 for theta vs p:
  TH1F* locHist_thetap_data[9][7];
  TH1F* locHist_thetap_denom_data[9][7];
  TH1F* locHist_thetap_gen[9][7];
  TH1F* locHist_thetap_denom_gen[9][7];
  TH1F* locHist_thetap_eff_data[9][7];
  TH1F* locHist_thetap_eff_gen[9][7];

  //these are for error calculations
  double thetap_errorD_data[9][7][20];
  double thetap_data[9][7][20];
  double thetap_errorU_data[9][7][20];
  double thetap_denom_data[9][7][20];
  double thetap_errorD_gen[9][7][20];
  double thetap_gen[9][7][20];
  double thetap_errorU_gen[9][7][20];
  double thetap_denom_gen[9][7][20];
  for(int im = 0; im < 7; im++){
    for(int ip = 0; ip < 9; ip++){
      locHist_thetap_data[ip][im] = (TH1F*)fdata->Get(Form("theta_yield_%d_m%d", ip+1, im+1));
      locHist_thetap_denom_data[ip][im] = (TH1F*)fdata->Get(Form("theta_yield_denom_%d_m%d", ip+1, im+1));
      locHist_thetap_gen[ip][im] = (TH1F*)fgen->Get(Form("theta_yield_%d_m%d", ip+1, im+1));
      locHist_thetap_denom_gen[ip][im] = (TH1F*)fgen->Get(Form("theta_yield_denom_%d_m%d", ip+1, im+1));
      locHist_thetap_eff_data[ip][im] = new TH1F(Form("thetap_efficiency_data_%d_m%d", ip+1, im+1), ";#theta (deg);#pi^{+} efficiency / 1.5#circ", 20, 0., 30.);
      locHist_thetap_eff_gen[ip][im] = new TH1F(Form("thetap_efficiency_gen_%d_m%d", ip+1, im+1), ";#theta (deg);#pi^{+} efficiency / 1.5#circ", 20, 0., 30.);
      for(int i = 0; i < 20; i++){
	if(locHist_thetap_data[ip][im]->GetBinContent(i+1) == locHist_thetap_denom_data[ip][im]->GetBinContent(i+1))
	  continue;
	thetap_data[ip][im][i] = locHist_thetap_data[ip][im]->GetBinContent(i+1);
	thetap_errorD_data[ip][im][i] = locHist_thetap_data[ip][im]->GetBinError(i+1);
	thetap_denom_data[ip][im][i] = locHist_thetap_denom_data[ip][im]->GetBinContent(i+1);
	thetap_errorU_data[ip][im][i] = locHist_thetap_denom_data[ip][im]->GetBinError(i+1); 
	double content = thetap_data[ip][im][i] / thetap_denom_data[ip][im][i];
	double error = TMath::Sqrt((TMath::Power((thetap_denom_data[ip][im][i] - thetap_data[ip][im][i])/(thetap_denom_data[ip][im][i]*thetap_denom_data[ip][im][i]), 2.)*thetap_errorD_data[ip][im][i]*thetap_errorD_data[ip][im][i] + TMath::Power(thetap_data[ip][im][i]/(thetap_denom_data[ip][im][i]*thetap_denom_data[ip][im][i]), 2.)*thetap_errorU_data[ip][im][i]*thetap_errorU_data[ip][im][i]));
	locHist_thetap_eff_data[ip][im]->SetBinContent(i+1, content);
	locHist_thetap_eff_data[ip][im]->SetBinError(i+1, error);   
      }
      for(int i = 0; i < 20; i++){
	if(locHist_thetap_gen[ip][im]->GetBinContent(i+1) == locHist_thetap_denom_gen[ip][im]->GetBinContent(i+1))
	  continue;
	thetap_gen[ip][im][i] = locHist_thetap_gen[ip][im]->GetBinContent(i+1);
	thetap_errorD_gen[ip][im][i] = locHist_thetap_gen[ip][im]->GetBinError(i+1);
	thetap_denom_gen[ip][im][i] = locHist_thetap_denom_gen[ip][im]->GetBinContent(i+1);
	thetap_errorU_gen[ip][im][i] = locHist_thetap_denom_gen[ip][im]->GetBinError(i+1);
	double content = thetap_gen[ip][im][i] / thetap_denom_gen[ip][im][i];
	double error = TMath::Sqrt((TMath::Power((thetap_denom_gen[ip][im][i] - thetap_gen[ip][im][i])/(thetap_denom_gen[ip][im][i]*thetap_denom_gen[ip][im][i]), 2.)*thetap_errorD_gen[ip][im][i]*thetap_errorD_gen[ip][im][i] + TMath::Power(thetap_gen[ip][im][i]/(thetap_denom_gen[ip][im][i]*thetap_denom_gen[ip][im][i]), 2.)*thetap_errorU_gen[ip][im][i]*thetap_errorU_gen[ip][im][i]));
	locHist_thetap_eff_gen[ip][im]->SetBinContent(i+1, content);
	locHist_thetap_eff_gen[ip][im]->SetBinError(i+1, error);
      }
    }
  }

  //Put theta vs p in a TH2F (method 1):
  const Int_t NBINS = 9;
  Double_t edges[NBINS + 1] = {0.5, 0.75, 1., 1.25, 1.5, 2., 3., 4., 5., 6.};
  TH2F* locHist_thetap_2D_data[7];
  TH2F* locHist_thetap_2D_gen[7];

  for(int im = 0; im < 7; im++){//loop over cuts
    locHist_thetap_2D_data[im] = new TH2F(Form("thetap_2D_data_%d", im+1), ";#theta (deg);p (GeV)", 20, 0., 30., NBINS, edges);
    locHist_thetap_2D_gen[im] = new TH2F(Form("thetap_2D_gen_%d", im+1), ";#theta (deg);p (GeV)", 20, 0., 30., NBINS, edges);

    for(int iy = 0; iy < 9; iy++){//loop over p (y-axis)
      for(int ix = 0; ix < 20; ix++){//loop over theta (x-axis)
	if(!locHist_thetap_eff_data[iy][im]->GetBinContent(ix+1))
	  continue;
	if(locHist_thetap_eff_data[iy][im]->GetBinError(ix+1) > 0.05*locHist_thetap_eff_data[iy][im]->GetBinContent(ix+1))
	  continue;
	locHist_thetap_2D_data[im]->SetBinContent(ix+1, iy+1, locHist_thetap_eff_data[iy][im]->GetBinContent(ix+1));
	locHist_thetap_2D_data[im]->SetBinError(ix+1, iy+1, locHist_thetap_eff_data[iy][im]->GetBinError(ix+1));
      }
    }
    for(int iy = 0; iy < 9; iy++){//loop over p (y-axis)
      for(int ix = 0; ix < 20; ix++){//loop over theta (x-axis)
	if(!locHist_thetap_eff_gen[iy][im]->GetBinContent(ix+1))
	  continue;
	if(locHist_thetap_eff_gen[iy][im]->GetBinError(ix+1) > 0.05*locHist_thetap_eff_gen[iy][im]->GetBinContent(ix+1))
	  continue;
	locHist_thetap_2D_gen[im]->SetBinContent(ix+1, iy+1, locHist_thetap_eff_gen[iy][im]->GetBinContent(ix+1));
	locHist_thetap_2D_gen[im]->SetBinError(ix+1, iy+1, locHist_thetap_eff_gen[iy][im]->GetBinError(ix+1));
      }
    }
  }
  
  //Repeat method 2 for theta vs p:
  TH1F* locHist_thetap_data_m2[9][7];
  TH1F* locHist_thetap_denom_data_m2[9][7];
  TH1F* locHist_thetap_gen_m2[9][7];
  TH1F* locHist_thetap_denom_gen_m2[9][7];
  TH1F* locHist_thetap_eff_data_m2[9][7];
  TH1F* locHist_thetap_eff_gen_m2[9][7];
  double thetap_errorD_data_m2[9][7][20];
  double thetap_data_m2[9][7][20];
  double thetap_errorU_data_m2[9][7][20];
  double thetap_denom_data_m2[9][7][20];
  double thetap_errorD_gen_m2[9][7][20];
  double thetap_gen_m2[9][7][20];
  double thetap_errorU_gen_m2[9][7][20];
  double thetap_denom_gen_m2[9][7][20];
  for(int im = 0; im < 7; im++){
    for(int ip = 0; ip < 9; ip++){
      locHist_thetap_data_m2[ip][im] = (TH1F*)fdata_m2->Get(Form("theta_yield_%d_m%d", ip+1, im+1));
      locHist_thetap_denom_data_m2[ip][im] = (TH1F*)fdata_m2->Get(Form("theta_yield_denom_%d_m%d", ip+1, im+1));
      locHist_thetap_gen_m2[ip][im] = (TH1F*)fgen_m2->Get(Form("theta_yield_%d_m%d", ip+1, im+1));
      locHist_thetap_denom_gen_m2[ip][im] = (TH1F*)fgen_m2->Get(Form("theta_yield_denom_%d_m%d", ip+1, im+1));
      locHist_thetap_eff_data_m2[ip][im] = new TH1F(Form("thetap_efficiency_data_m2_%d_m%d", ip+1, im+1), ";#theta (deg);#pi^{+} efficiency / 1.5#circ", 20, 0., 30.);
      locHist_thetap_eff_gen_m2[ip][im] = new TH1F(Form("thetap_efficiency_gen_m2_%d_m%d", ip+1, im+1), ";#theta (deg);#pi^{+} efficiency / 1.5#circ", 20, 0., 30.);
    for(int i = 0; i < 20; i++){
      if(locHist_thetap_data_m2[ip][im]->GetBinContent(i+1) == locHist_thetap_denom_data_m2[ip][im]->GetBinContent(i+1))
	continue;
      thetap_data_m2[ip][im][i] = locHist_thetap_data_m2[ip][im]->GetBinContent(i+1);
      thetap_errorD_data_m2[ip][im][i] = locHist_thetap_data_m2[ip][im]->GetBinError(i+1);
      thetap_denom_data_m2[ip][im][i] = locHist_thetap_denom_data_m2[ip][im]->GetBinContent(i+1);
      thetap_errorU_data_m2[ip][im][i] = locHist_thetap_denom_data_m2[ip][im]->GetBinError(i+1); 
      double content = thetap_data_m2[ip][im][i] / thetap_denom_data_m2[ip][im][i];
      double error = TMath::Sqrt((TMath::Power((thetap_denom_data_m2[ip][im][i] - thetap_data_m2[ip][im][i])/(thetap_denom_data_m2[ip][im][i]*thetap_denom_data_m2[ip][im][i]), 2.)*thetap_errorD_data_m2[ip][im][i]*thetap_errorD_data_m2[ip][im][i] + TMath::Power(thetap_data_m2[ip][im][i]/(thetap_denom_data_m2[ip][im][i]*thetap_denom_data_m2[ip][im][i]), 2.)*thetap_errorU_data_m2[ip][im][i]*thetap_errorU_data_m2[ip][im][i]));
      locHist_thetap_eff_data_m2[ip][im]->SetBinContent(i+1, content);
      locHist_thetap_eff_data_m2[ip][im]->SetBinError(i+1, error);   
    }//
    for(int i = 0; i < 20; i++){
      if(locHist_thetap_gen_m2[ip][im]->GetBinContent(i+1) == locHist_thetap_denom_gen_m2[ip][im]->GetBinContent(i+1))
	continue;
      thetap_gen_m2[ip][im][i] = locHist_thetap_gen_m2[ip][im]->GetBinContent(i+1);
      thetap_errorD_gen_m2[ip][im][i] = locHist_thetap_gen_m2[ip][im]->GetBinError(i+1);
      thetap_denom_gen_m2[ip][im][i] = locHist_thetap_denom_gen_m2[ip][im]->GetBinContent(i+1);
      thetap_errorU_gen_m2[ip][im][i] = locHist_thetap_denom_gen_m2[ip][im]->GetBinError(i+1);
      double content = thetap_gen_m2[ip][im][i] / thetap_denom_gen_m2[ip][im][i];
      double error = TMath::Sqrt((TMath::Power((thetap_denom_gen_m2[ip][im][i] - thetap_gen_m2[ip][im][i])/(thetap_denom_gen_m2[ip][im][i]*thetap_denom_gen_m2[ip][im][i]), 2.)*thetap_errorD_gen_m2[ip][im][i]*thetap_errorD_gen_m2[ip][im][i] + TMath::Power(thetap_gen_m2[ip][im][i]/(thetap_denom_gen_m2[ip][im][i]*thetap_denom_gen_m2[ip][im][i]), 2.)*thetap_errorU_gen_m2[ip][im][i]*thetap_errorU_gen_m2[ip][im][i]));
      locHist_thetap_eff_gen_m2[ip][im]->SetBinContent(i+1, content);
      locHist_thetap_eff_gen_m2[ip][im]->SetBinError(i+1, error);
    }
    }
  }

  //Put theta vs p in a TH2F (method 2)
  TH2F* locHist_thetap_2D_data_m2[7];
  TH2F* locHist_thetap_2D_gen_m2[7];

  for(int im = 0; im < 7; im++){//loop over cuts
    locHist_thetap_2D_data_m2[im] = new TH2F(Form("thetap_2D_data_m2_%d", im+1), ";#theta (deg);p (GeV)", 20, 0., 30., NBINS, edges);
    locHist_thetap_2D_gen_m2[im] = new TH2F(Form("thetap_2D_gen_m2_%d", im+1), ";#theta (deg);p (GeV)", 20, 0., 30., NBINS, edges);

    for(int iy = 0; iy < 9; iy++){//loop over p (y-axis)
      for(int ix = 0; ix < 20; ix++){//loop over theta (x-axis)
	if(!locHist_thetap_eff_data_m2[iy][im]->GetBinContent(ix+1))
	  continue;
	if(locHist_thetap_eff_data_m2[iy][im]->GetBinError(ix+1) > 0.05*locHist_thetap_eff_data_m2[iy][im]->GetBinContent(ix+1))
	  continue;
	locHist_thetap_2D_data_m2[im]->SetBinContent(ix+1, iy+1, locHist_thetap_eff_data_m2[iy][im]->GetBinContent(ix+1));
	locHist_thetap_2D_data_m2[im]->SetBinError(ix+1, iy+1, locHist_thetap_eff_data_m2[iy][im]->GetBinError(ix+1));
      }
    }
    for(int iy = 0; iy < 9; iy++){//loop over p (y-axis)
      for(int ix = 0; ix < 20; ix++){//loop over theta (x-axis)
	if(!locHist_thetap_eff_gen_m2[iy][im]->GetBinContent(ix+1))
	  continue;
	if(locHist_thetap_eff_gen_m2[iy][im]->GetBinError(ix+1) > 0.05*locHist_thetap_eff_gen_m2[iy][im]->GetBinContent(ix+1))
	  continue;
	locHist_thetap_2D_gen_m2[im]->SetBinContent(ix+1, iy+1, locHist_thetap_eff_gen_m2[iy][im]->GetBinContent(ix+1));
	locHist_thetap_2D_gen_m2[im]->SetBinError(ix+1, iy+1, locHist_thetap_eff_gen_m2[iy][im]->GetBinError(ix+1));
      }
    }
  }

  //Calculate ratio TH2Fs:
  TH2F* locHist_thetap_ratio_2D_data[7];
  TH2F* locHist_thetap_ratio_2D_gen[7];
  TH2F* locHist_thetap_ratio_2D_m1[7];
  TH2F* locHist_thetap_ratio_2D_m2[7];
  for(int im = 0; im < 7; im++){
    locHist_thetap_ratio_2D_data[im] = (TH2F*)locHist_thetap_2D_data[im]->Clone(Form("thetap_ratio_data_%d", im+1));
    locHist_thetap_ratio_2D_data[im]->Sumw2();
    locHist_thetap_ratio_2D_data[im]->Divide(locHist_thetap_2D_data_m2[im]);

    locHist_thetap_ratio_2D_gen[im] = (TH2F*)locHist_thetap_2D_gen[im]->Clone(Form("thetap_ratio_gen_%d", im+1));
    locHist_thetap_ratio_2D_gen[im]->Sumw2();
    locHist_thetap_ratio_2D_gen[im]->Divide(locHist_thetap_2D_gen_m2[im]);

    locHist_thetap_ratio_2D_m1[im] = (TH2F*)locHist_thetap_2D_data[im]->Clone(Form("thetap_ratio_m1_%d", im+1));
    locHist_thetap_ratio_2D_m1[im]->Sumw2();
    locHist_thetap_ratio_2D_m1[im]->Divide(locHist_thetap_2D_gen[im]);
    locHist_thetap_ratio_2D_m2[im] = (TH2F*)locHist_thetap_2D_data_m2[im]->Clone(Form("thetap_ratio_m2_%d", im+1));
    locHist_thetap_ratio_2D_m2[im]->Sumw2();
    locHist_thetap_ratio_2D_m2[im]->Divide(locHist_thetap_2D_gen_m2[im]);
  }
  
  //Find Data/MC ratio (method 2):
  TH1F* locHist_phi_ratio_m2[7];
  TH1F* locHist_theta_ratio_m2[7];
  TH1F* locHist_p_ratio_m2[7];
  for(int im = 0; im < 7; im++){
    locHist_phi_ratio_m2[im] = (TH1F*)locHist_phi_data_m2[im]->Clone(Form("phi_ratio_m2_%d", im+1));
    locHist_phi_ratio_m2[im]->Sumw2();
    locHist_phi_ratio_m2[im]->Divide(locHist_phi_gen_m2[im]);
    locHist_theta_ratio_m2[im] = (TH1F*)locHist_theta_data_m2[im]->Clone(Form("theta_ratio_m2_%d", im+1));
    locHist_theta_ratio_m2[im]->Sumw2();
    locHist_theta_ratio_m2[im]->Divide(locHist_theta_gen_m2[im]);
    locHist_p_ratio_m2[im] = (TH1F*)locHist_p_data_m2[im]->Clone(Form("p_ratio_m2_%d", im+1));
    locHist_p_ratio_m2[im]->Sumw2();
    locHist_p_ratio_m2[im]->Divide(locHist_p_gen_m2[im]);
  }

  
  //Find method 1/2 ratio (data)
  TH1F* locHist_phi_12ratio_data[7];
  TH1F* locHist_theta_12ratio_data[7];
  TH1F* locHist_p_12ratio_data[7];
  TH1F* locHist_phi_12ratio_gen[7];
  TH1F* locHist_theta_12ratio_gen[7];
  TH1F* locHist_p_12ratio_gen[7];
  for(int im = 0; im < 7; im++){
    locHist_phi_12ratio_data[im] = (TH1F*)locHist_phi_data[im]->Clone();
    locHist_phi_12ratio_data[im]->Sumw2();
    locHist_phi_12ratio_data[im]->Divide(locHist_phi_data_m2[im]);
    locHist_theta_12ratio_data[im] = (TH1F*)locHist_theta_data[im]->Clone();
    locHist_theta_12ratio_data[im]->Sumw2();
    locHist_theta_12ratio_data[im]->Divide(locHist_theta_data_m2[im]);  
    locHist_p_12ratio_data[im] = (TH1F*)locHist_p_data[im]->Clone();
    locHist_p_12ratio_data[im]->Sumw2();
    locHist_p_12ratio_data[im]->Divide(locHist_p_data_m2[im]);
    
    locHist_phi_12ratio_gen[im] = (TH1F*)locHist_phi_gen[im]->Clone();
    locHist_phi_12ratio_gen[im]->Sumw2();
    locHist_phi_12ratio_gen[im]->Divide(locHist_phi_gen_m2[im]);
    locHist_theta_12ratio_gen[im] = (TH1F*)locHist_theta_gen[im]->Clone();
    locHist_theta_12ratio_gen[im]->Sumw2();
    locHist_theta_12ratio_gen[im]->Divide(locHist_theta_gen_m2[im]);  
    locHist_p_12ratio_gen[im] = (TH1F*)locHist_p_gen[im]->Clone();
    locHist_p_12ratio_gen[im]->Sumw2();
    locHist_p_12ratio_gen[im]->Divide(locHist_p_gen_m2[im]);
  }

  //Find cut/uncut ratios
  TH1F* locHist_theta_cutratio_data[7];
  TH1F* locHist_theta_cutratio_gen[7];
  TH1F* locHist_theta_cutratio_data_m2[7];
  TH1F* locHist_theta_cutratio_gen_m2[7];
  TH1F* locHist_p_cutratio_data[7];
  TH1F* locHist_p_cutratio_gen[7];
  TH1F* locHist_p_cutratio_data_m2[7];
  TH1F* locHist_p_cutratio_gen_m2[7];
  for(int im = 0; im < 7; im++){
    locHist_theta_cutratio_data[im] = (TH1F*)locHist_theta_data[im]->Clone();
    locHist_theta_cutratio_data[im]->Sumw2();
    locHist_theta_cutratio_data[im]->Divide(locHist_theta_data[0]);
    locHist_theta_cutratio_gen[im] = (TH1F*)locHist_theta_gen[im]->Clone();
    locHist_theta_cutratio_gen[im]->Sumw2();
    locHist_theta_cutratio_gen[im]->Divide(locHist_theta_gen[0]);
    locHist_theta_cutratio_data_m2[im] = (TH1F*)locHist_theta_data_m2[im]->Clone();
    locHist_theta_cutratio_data_m2[im]->Sumw2();
    locHist_theta_cutratio_data_m2[im]->Divide(locHist_theta_data_m2[0]);
    locHist_theta_cutratio_gen_m2[im] = (TH1F*)locHist_theta_gen_m2[im]->Clone();
    locHist_theta_cutratio_gen_m2[im]->Sumw2();
    locHist_theta_cutratio_gen_m2[im]->Divide(locHist_theta_gen_m2[0]);
    locHist_p_cutratio_data[im] = (TH1F*)locHist_p_data[im]->Clone();
    locHist_p_cutratio_data[im]->Sumw2();
    locHist_p_cutratio_data[im]->Divide(locHist_p_data[0]);
    locHist_p_cutratio_gen[im] = (TH1F*)locHist_p_gen[im]->Clone();
    locHist_p_cutratio_gen[im]->Sumw2();
    locHist_p_cutratio_gen[im]->Divide(locHist_p_gen[0]);
    locHist_p_cutratio_data_m2[im] = (TH1F*)locHist_p_data_m2[im]->Clone();
    locHist_p_cutratio_data_m2[im]->Sumw2();
    locHist_p_cutratio_data_m2[im]->Divide(locHist_p_data_m2[0]);
    locHist_p_cutratio_gen_m2[im] = (TH1F*)locHist_p_gen_m2[im]->Clone();
    locHist_p_cutratio_gen_m2[im]->Sumw2();
    locHist_p_cutratio_gen_m2[im]->Divide(locHist_p_gen_m2[0]);
  }

  
  Double_t shiftphi = 0.1*locHist_phi_ratio[0]->GetBinWidth(1);
  Double_t shifttheta = 0.1*locHist_theta_ratio[0]->GetBinWidth(1);
  Double_t shiftp = 0.1*locHist_p_ratio[0]->GetBinWidth(1);
  

  TString masscut[7] = {"0.", "0.1", "0.05", "0.02", "0.01", "0.0075", "0.005"};
  TString pcut[7] = {"0", "0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1"};

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat("");


  
  TCanvas* c1[7];
  for(int im = 0; im < 7; im++){
    TString cutname = "Efficiencies (Method 1)";

    c1[im] = new TCanvas(Form("c1_%d", im+1), cutname, 900, 300);
    c1[im]->Divide(3,1);

    c1[im]->cd(1);
    locHist_phi_data[im]->SetTitle("Data (Method 1);#phi (deg);#pi^{+} efficiency (method 1)");
    locHist_phi_gen[im]->SetTitle("MC (Method 1)");
    locHist_phi_data[im]->SetMaximum(1.05);
    locHist_phi_data[im]->SetMinimum(0.0);
    locHist_phi_data[im]->SetMarkerStyle(kFullCircle);
    locHist_phi_gen[im]->SetMarkerStyle(kOpenSquare);
    locHist_phi_data[im]->SetMarkerColor(kBlue);
    locHist_phi_gen[im]->SetMarkerColor(kRed);
    locHist_phi_data[im]->SetLineColor(kBlue);
    locHist_phi_gen[im]->SetLineColor(kRed);
    locHist_phi_data[im]->Draw();
    locHist_phi_gen[im]->Draw("same");
    c1[im]->cd(2);
    locHist_theta_data[im]->SetTitle(";#theta (deg);;");
    locHist_theta_gen[im]->SetTitle(";#theta (deg);;");
    locHist_theta_data[im]->SetMaximum(1.05);
    locHist_theta_data[im]->SetMinimum(0.0);
    locHist_theta_data[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_gen[im]->SetMarkerStyle(kOpenSquare);
    locHist_theta_data[im]->SetMarkerColor(kBlue);
    locHist_theta_gen[im]->SetMarkerColor(kRed);
    locHist_theta_data[im]->SetLineColor(kBlue);
    locHist_theta_gen[im]->SetLineColor(kRed);
    locHist_theta_data[im]->Draw();
    locHist_theta_gen[im]->Draw("same");
    c1[im]->cd(3);
    locHist_p_data[im]->SetTitle(";p (GeV);;");
    locHist_p_data[im]->SetMaximum(1.05);
    locHist_p_data[im]->SetMinimum(0.0);
    locHist_p_data[im]->SetMarkerStyle(kFullCircle);
    locHist_p_gen[im]->SetMarkerStyle(kOpenSquare);
    locHist_p_data[im]->SetMarkerColor(kBlue);
    locHist_p_gen[im]->SetMarkerColor(kRed);
    locHist_p_data[im]->SetLineColor(kBlue);
    locHist_p_gen[im]->SetLineColor(kRed);
    locHist_p_data[im]->Draw();
    locHist_p_gen[im]->Draw("same");
  }



  
  TCanvas* c2[7];
  for(int im = 0; im < 7; im++){
    TString cutname = "Efficiencies (Method 2) (";
    //if(im == 0)
    //cutname += "No MM^2 cut)";
    //else{
      cutname += "P(Chi2) > ";
      cutname += pcut[im];
      cutname += ")";
      //}

    c2[im] = new TCanvas(Form("c2_%d", im+1), cutname, 900, 300);
    c2[im]->Divide(3,1);

    c2[im]->cd(1);
    locHist_phi_data_m2[im]->SetTitle("Data (Method 2);#phi (deg);#pi^{+} efficiency (method 2)");
    locHist_phi_gen_m2[im]->SetTitle("MC (Method 2)");
    locHist_phi_data_m2[im]->SetMaximum(1.05);
    locHist_phi_data_m2[im]->SetMinimum(0.0);
    locHist_phi_data_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_phi_gen_m2[im]->SetMarkerStyle(kOpenSquare);
    locHist_phi_data_m2[im]->SetMarkerColor(kAzure+10);
    locHist_phi_gen_m2[im]->SetMarkerColor(kPink+1);
    locHist_phi_data_m2[im]->SetLineColor(kAzure+10);
    locHist_phi_gen_m2[im]->SetLineColor(kPink+1);
    locHist_phi_data_m2[im]->Draw();
    locHist_phi_gen_m2[im]->Draw("same");
    c2[im]->cd(2);
    locHist_theta_data_m2[im]->SetTitle(";#theta (deg);;");
    locHist_theta_gen_m2[im]->SetTitle(";#theta (deg);;");
    locHist_theta_data_m2[im]->SetMaximum(1.05);
    locHist_theta_data_m2[im]->SetMinimum(0.0);
    locHist_theta_data_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_gen_m2[im]->SetMarkerStyle(kOpenSquare);
    locHist_theta_data_m2[im]->SetMarkerColor(kAzure+10);
    locHist_theta_gen_m2[im]->SetMarkerColor(kPink+1);
    locHist_theta_data_m2[im]->SetLineColor(kAzure+10);
    locHist_theta_gen_m2[im]->SetLineColor(kPink+1);
    locHist_theta_data_m2[im]->Draw();
    locHist_theta_gen_m2[im]->Draw("same");
    c2[im]->cd(3);
    locHist_p_data_m2[im]->SetTitle(";p (GeV);;");
    locHist_p_data_m2[im]->SetMaximum(1.05);
    locHist_p_data_m2[im]->SetMinimum(0.0);
    locHist_p_data_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_p_gen_m2[im]->SetMarkerStyle(kOpenSquare);
    locHist_p_data_m2[im]->SetMarkerColor(kAzure+10);
    locHist_p_gen_m2[im]->SetMarkerColor(kPink+1);
    locHist_p_data_m2[im]->SetLineColor(kAzure+10);
    locHist_p_gen_m2[im]->SetLineColor(kPink+1);
    locHist_p_data_m2[im]->Draw();
    locHist_p_gen_m2[im]->Draw("same");
  }
 


  
  Int_t color[7];
  for(int i = 0; i < 7; i++) {
    color[i] = 2 + i;
  }
  TString pvcut[7] = {"No P(#chi^{2}) cuts", "P(#chi^{2}) > 0.000001", "P(#chi^{2}) > 0.00001", "P(#chi^{2}) > 0.0001", "P(#chi^{2}) > 0.001", "P(#chi^{2}) > 0.01", "P(#chi^{2}) > 0.1"};

  /* 
  TCanvas* ccuts = new TCanvas("ccuts", "#theta Efficiencies with P(Chi2) cuts", 800, 800);
  ccuts->Divide(2,2);
  for(int im = 0; im < 7; im++){
    ccuts->cd(1);
    locHist_theta_data[im]->SetMinimum(0);
    locHist_theta_data[im]->SetMaximum(1.05);
    locHist_theta_data[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_data[im]->SetTitle("Data (Method 1);#theta (deg);#pi^{+} efficiency");
    locHist_theta_data[im]->SetMarkerColor(color[im]);
    locHist_theta_data[im]->SetLineColor(color[im]);
    locHist_theta_data[im]->Draw("same");
    ccuts->cd(2);
    locHist_theta_gen[im]->SetMinimum(0);
    locHist_theta_gen[im]->SetMaximum(1.05);
    locHist_theta_gen[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_gen[im]->SetTitle("MC (Method 1);#theta (deg);#pi^{+} efficiency");
    locHist_theta_gen[im]->SetMarkerColor(color[im]);
    locHist_theta_gen[im]->SetLineColor(color[im]);
    locHist_theta_gen[im]->Draw("same");
    ccuts->cd(3);
    locHist_theta_data_m2[im]->SetMinimum(0);
    locHist_theta_data_m2[im]->SetMaximum(1.05);
    locHist_theta_data_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_data_m2[im]->SetTitle("Data (Method 2);#theta (deg);#pi^{+} efficiency");
    locHist_theta_data_m2[im]->SetMarkerColor(color[im]);
    locHist_theta_data_m2[im]->SetLineColor(color[im]);
    locHist_theta_data_m2[im]->Draw("same");
    ccuts->cd(4);
    locHist_theta_gen_m2[im]->SetMinimum(0);
    locHist_theta_gen_m2[im]->SetMaximum(1.05);
    locHist_theta_gen_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_gen_m2[im]->SetTitle("MC (Method 2);#theta (deg);#pi^{+} efficiency");
    locHist_theta_gen_m2[im]->SetMarkerColor(color[im]);
    locHist_theta_gen_m2[im]->SetLineColor(color[im]);
    locHist_theta_gen_m2[im]->Draw("same");
  }
  auto lcuts = new TLegend(0.4, 0.2, 0.9, 0.5);
  for(int im = 0; im < 7; im++){
    lcuts->AddEntry(locHist_theta_gen_m2[im], pvcut[im]);
  }
  ccuts->cd(4);
  lcuts->Draw();

  ccuts->Print("plots/ThetaEff_PValue_cuts.pdf");
  */
  /*
  TCanvas* ccuts_ratio = new TCanvas("ccuts_ratio", "#theta efficiency ratios with cuts", 800, 800);
  ccuts_ratio->Divide(2,2);
  for(int im = 0; im < 7; im++){
    ccuts_ratio->cd(1);
    locHist_theta_cutratio_data[im]->SetMinimum(0);
    locHist_theta_cutratio_data[im]->SetMaximum(1.05);
    locHist_theta_cutratio_data[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_cutratio_data[im]->SetTitle("Data (Method 1);#theta (deg);#pi^{+} efficiency (cut/uncut)");
    locHist_theta_cutratio_data[im]->SetMarkerColor(color[im]);
    locHist_theta_cutratio_data[im]->SetLineColor(color[im]);
    locHist_theta_cutratio_data[im]->Draw("same");
    ccuts_ratio->cd(2);
    locHist_theta_cutratio_gen[im]->SetMinimum(0);
    locHist_theta_cutratio_gen[im]->SetMaximum(1.05);
    locHist_theta_cutratio_gen[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_cutratio_gen[im]->SetTitle("MC (Method 1);#theta (deg);#pi^{+} efficiency (cut/uncut)");
    locHist_theta_cutratio_gen[im]->SetMarkerColor(color[im]);
    locHist_theta_cutratio_gen[im]->SetLineColor(color[im]);
    locHist_theta_cutratio_gen[im]->Draw("same");
    ccuts_ratio->cd(3);
    locHist_theta_cutratio_data_m2[im]->SetMinimum(0);
    locHist_theta_cutratio_data_m2[im]->SetMaximum(1.05);
    locHist_theta_cutratio_data_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_cutratio_data_m2[im]->SetTitle("Data (Method 2);#theta (deg);#pi^{+} efficiency (cut/uncut)");
    locHist_theta_cutratio_data_m2[im]->SetMarkerColor(color[im]);
    locHist_theta_cutratio_data_m2[im]->SetLineColor(color[im]);
    locHist_theta_cutratio_data_m2[im]->Draw("same");
    ccuts_ratio->cd(4);
    locHist_theta_cutratio_gen_m2[im]->SetMinimum(0);
    locHist_theta_cutratio_gen_m2[im]->SetMaximum(1.05);
    locHist_theta_cutratio_gen_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_cutratio_gen_m2[im]->SetTitle("MC (Method 2);#theta (deg);#pi^{+} efficiency (cut/uncut)");
    locHist_theta_cutratio_gen_m2[im]->SetMarkerColor(color[im]);
    locHist_theta_cutratio_gen_m2[im]->SetLineColor(color[im]);
    locHist_theta_cutratio_gen_m2[im]->Draw("same");
  }
  */
  /*
  TCanvas* ccuts_p = new TCanvas("ccuts_p", "p Efficiencies with P(Chi2) cuts", 800, 800);
  ccuts_p->Divide(2,2);
  for(int im = 0; im < 7; im++){
    ccuts_p->cd(1);
    locHist_p_data[im]->SetMinimum(0);
    locHist_p_data[im]->SetMaximum(1.05);
    locHist_p_data[im]->SetMarkerStyle(kFullCircle);
    locHist_p_data[im]->SetTitle("Data (Method 1);p (GeV);#pi^{+} efficiency");
    locHist_p_data[im]->SetMarkerColor(color[im]);
    locHist_p_data[im]->SetLineColor(color[im]);
    locHist_p_data[im]->Draw("same");
    ccuts_p->cd(2);
    locHist_p_gen[im]->SetMinimum(0);
    locHist_p_gen[im]->SetMaximum(1.05);
    locHist_p_gen[im]->SetMarkerStyle(kFullCircle);
    locHist_p_gen[im]->SetTitle("MC (Method 1);p (GeV);#pi^{+} efficiency");
    locHist_p_gen[im]->SetMarkerColor(color[im]);
    locHist_p_gen[im]->SetLineColor(color[im]);
    locHist_p_gen[im]->Draw("same");
    ccuts_p->cd(3);
    locHist_p_data_m2[im]->SetMinimum(0);
    locHist_p_data_m2[im]->SetMaximum(1.05);
    locHist_p_data_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_p_data_m2[im]->SetTitle("Data (Method 2);p (GeV);#pi^{+} efficiency");
    locHist_p_data_m2[im]->SetMarkerColor(color[im]);
    locHist_p_data_m2[im]->SetLineColor(color[im]);
    locHist_p_data_m2[im]->Draw("same");
    ccuts_p->cd(4);
    locHist_p_gen_m2[im]->SetMinimum(0);
    locHist_p_gen_m2[im]->SetMaximum(1.05);
    locHist_p_gen_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_p_gen_m2[im]->SetTitle("MC (Method 2);p (GeV);#pi^{+} efficiency");
    locHist_p_gen_m2[im]->SetMarkerColor(color[im]);
    locHist_p_gen_m2[im]->SetLineColor(color[im]);
    locHist_p_gen_m2[im]->Draw("same");
  }
  ccuts_p->cd(4);
  lcuts->Draw();

  ccuts_p->Print("plots/pEff_PValue_cuts.pdf");
  */

  /*
  TCanvas* ccuts_ratio_p = new TCanvas("ccuts_ratio_p", "p efficiency ratios with cuts", 800, 800);
  ccuts_ratio_p->Divide(2,2);
  for(int im = 0; im < 7; im++){
    ccuts_ratio_p->cd(1);
    locHist_p_cutratio_data[im]->SetMinimum(0);
    locHist_p_cutratio_data[im]->SetMaximum(1.05);
    locHist_p_cutratio_data[im]->SetMarkerStyle(kFullCircle);
    locHist_p_cutratio_data[im]->SetTitle("Data (Method 1);p (GeV);#pi^{+} efficiency (cut/uncut)");
    locHist_p_cutratio_data[im]->SetMarkerColor(color[im]);
    locHist_p_cutratio_data[im]->SetLineColor(color[im]);
    locHist_p_cutratio_data[im]->Draw("same");
    ccuts_ratio_p->cd(2);
    locHist_p_cutratio_gen[im]->SetMinimum(0);
    locHist_p_cutratio_gen[im]->SetMaximum(1.05);
    locHist_p_cutratio_gen[im]->SetMarkerStyle(kFullCircle);
    locHist_p_cutratio_gen[im]->SetTitle("MC (Method 1);p (GeV);#pi^{+} efficiency (cut/uncut)");
    locHist_p_cutratio_gen[im]->SetMarkerColor(color[im]);
    locHist_p_cutratio_gen[im]->SetLineColor(color[im]);
    locHist_p_cutratio_gen[im]->Draw("same");
    ccuts_ratio_p->cd(3);
    locHist_p_cutratio_data_m2[im]->SetMinimum(0);
    locHist_p_cutratio_data_m2[im]->SetMaximum(1.05);
    locHist_p_cutratio_data_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_p_cutratio_data_m2[im]->SetTitle("Data (Method 2);p (GeV);#pi^{+} efficiency (cut/uncut)");
    locHist_p_cutratio_data_m2[im]->SetMarkerColor(color[im]);
    locHist_p_cutratio_data_m2[im]->SetLineColor(color[im]);
    locHist_p_cutratio_data_m2[im]->Draw("same");
    ccuts_ratio_p->cd(4);
    locHist_p_cutratio_gen_m2[im]->SetMinimum(0);
    locHist_p_cutratio_gen_m2[im]->SetMaximum(1.05);
    locHist_p_cutratio_gen_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_p_cutratio_gen_m2[im]->SetTitle("MC (Method 2);p (GeV);#pi^{+} efficiency (cut/uncut)");
    locHist_p_cutratio_gen_m2[im]->SetMarkerColor(color[im]);
    locHist_p_cutratio_gen_m2[im]->SetLineColor(color[im]);
    locHist_p_cutratio_gen_m2[im]->Draw("same");
  }
  */
  
  /*
  // TCanvas* cl = new TCanvas("cl", "legend", 600, 600);
  // cl->cd(1);
  // locHist_phi_data_m2->Draw();
  // locHist_phi_gen_m2->Draw("same");
  // cl->BuildLegend();
  */

  TCanvas* c3 = new TCanvas("c3", "Efficiencies (Both methods)", 900, 300);
  c3->Divide(3,1);

  c3->cd(1);
  locHist_phi_data[0]->SetTitle("Data (Method 1);#phi (deg);#pi^{+} efficiency");
  locHist_phi_data[0]->GetXaxis()->SetLimits(-180-2*shiftphi, 180-2*shiftphi);
  locHist_phi_data_m2[0]->GetXaxis()->SetLimits(-180-shiftphi, 180-shiftphi);
  locHist_phi_gen[0]->GetXaxis()->SetLimits(-180+2*shiftphi, 180+2*shiftphi);
  locHist_phi_gen_m2[0]->GetXaxis()->SetLimits(-180+shiftphi, 180+shiftphi);
  locHist_phi_data[0]->Draw();
  locHist_phi_data_m2[0]->Draw("same");
  locHist_phi_gen[0]->Draw("same");
  locHist_phi_gen_m2[0]->Draw("same");
  c3->cd(2);
  locHist_theta_data[0]->GetXaxis()->SetLimits(0-2*shifttheta, 30-2*shifttheta);
  locHist_theta_data_m2[0]->GetXaxis()->SetLimits(0-shifttheta, 30-shifttheta);
  locHist_theta_gen[0]->GetXaxis()->SetLimits(0+2*shifttheta, 30+2*shifttheta);
  locHist_theta_gen_m2[0]->GetXaxis()->SetLimits(0+shifttheta, 30+shifttheta);
  locHist_theta_data[0]->Draw();
  locHist_theta_data_m2[0]->Draw("same");
  locHist_theta_gen[0]->Draw("same");
  locHist_theta_gen_m2[0]->Draw("same");
  c3->cd(3);
  locHist_p_data[0]->GetXaxis()->SetLimits(0-2*shiftp, 6-2*shiftp);
  locHist_p_data_m2[0]->GetXaxis()->SetLimits(0-shiftp, 6-shiftp);
  locHist_p_gen[0]->GetXaxis()->SetLimits(0+2*shiftp, 6+2*shiftp);
  locHist_p_gen_m2[0]->GetXaxis()->SetLimits(0+shiftp, 6+shiftp);
  locHist_p_data[0]->Draw();
  locHist_p_data_m2[0]->Draw("same");
  locHist_p_gen[0]->Draw("same");
  locHist_p_gen_m2[0]->Draw("same");
  auto l3 = new TLegend(0.4, 0.2, 0.9, 0.5);
  l3->AddEntry(locHist_p_data[0], "Data (Method 1)");
  l3->AddEntry(locHist_p_data_m2[0], "Data (Method 2)");
  l3->AddEntry(locHist_p_gen[0], "MC (Method 1)");
  l3->AddEntry(locHist_p_gen_m2[0], "MC (Method 2)");
  l3->Draw();
  c3->Print("plots/PiPlus_Efficiency.C");
  /*
  TCanvas* c5 = new TCanvas("c5", "c5", 900, 300);
  c5->Divide(3,1);
  c5->cd(1);
  locHist_phi_ratio[0]->SetTitle(";#phi (deg);Data/MC Ratio");
  locHist_phi_ratio[0]->SetMinimum(0.9);
  locHist_phi_ratio[0]->SetMaximum(1.1);
  locHist_phi_ratio[0]->SetMarkerStyle(kFullCircle);
  locHist_phi_ratio[0]->SetMarkerColor(kViolet+10);
  locHist_phi_ratio[0]->SetLineColor(kViolet+10);
  locHist_phi_ratio_m2[0]->SetMarkerStyle(kFullCircle);
  locHist_phi_ratio_m2[0]->SetMarkerColor(kViolet+6);
  locHist_phi_ratio_m2[0]->SetLineColor(kViolet+6);
  locHist_phi_ratio[0]->Draw();
  locHist_phi_ratio_m2[0]->Draw("same");
  c5->cd(2);
  locHist_theta_ratio[0]->SetTitle(";#theta (deg);;");
  locHist_theta_ratio[0]->SetMinimum(0.9);
  locHist_theta_ratio[0]->SetMaximum(1.1);
  locHist_theta_ratio[0]->SetMarkerStyle(kFullCircle);
  locHist_theta_ratio[0]->SetMarkerColor(kViolet+10);
  locHist_theta_ratio[0]->SetLineColor(kViolet+10);
  locHist_theta_ratio_m2[0]->SetMarkerStyle(kFullCircle);
  locHist_theta_ratio_m2[0]->SetMarkerColor(kViolet+6);
  locHist_theta_ratio_m2[0]->SetLineColor(kViolet+6);
  locHist_theta_ratio[0]->Draw();
  locHist_theta_ratio_m2[0]->Draw("same");
  c5->cd(3);
  locHist_p_ratio[0]->SetTitle(";p (GeV);;");
  locHist_p_ratio[0]->SetMinimum(0.9);
  locHist_p_ratio[0]->SetMaximum(1.1);
  locHist_p_ratio[0]->SetMarkerStyle(kFullCircle);
  locHist_p_ratio[0]->SetMarkerColor(kViolet+10);
  locHist_p_ratio[0]->SetLineColor(kViolet+10);
  locHist_p_ratio_m2[0]->SetMarkerStyle(kFullCircle);
  locHist_p_ratio_m2[0]->SetMarkerColor(kViolet+6);
  locHist_p_ratio_m2[0]->SetLineColor(kViolet+6);
  locHist_p_ratio[0]->Draw();
  locHist_p_ratio_m2[0]->Draw("same");
  auto* l5 = new TLegend(0.5, 0.2, 0.9, 0.4);
  l5->AddEntry(locHist_p_ratio[0], "Method 1");
  l5->AddEntry(locHist_p_ratio_m2[0], "Method 2");
  l5->Draw();

  //c5->Print("plots/PiPlus_Ratio_Int.pdf");
  */
  /*
  TCanvas* call = new TCanvas("call", "efficiencies with cuts", 1);
  call->Divide(6, 4);

  for(int i = 0; i < 7; i++){
    call->cd(1 + 3*i);
    locHist_phi_data[i]->SetTitle("Data (Method 1);#phi (deg);#pi^{+} efficiency");
    locHist_phi_data[i]->GetXaxis()->SetLimits(-180-2*shiftphi, 180-2*shiftphi);
    locHist_phi_data_m2[i]->GetXaxis()->SetLimits(-180-shiftphi, 180-shiftphi);
    locHist_phi_gen[i]->GetXaxis()->SetLimits(-180+2*shiftphi, 180+2*shiftphi);
    locHist_phi_gen_m2[i]->GetXaxis()->SetLimits(-180+shiftphi, 180+shiftphi);
    locHist_phi_data[i]->Draw();
    locHist_phi_data_m2[i]->Draw("same");
    locHist_phi_gen[i]->Draw("same");
    locHist_phi_gen_m2[i]->Draw("same");
    call->cd(2 + 3*i);
    locHist_theta_data[i]->GetXaxis()->SetLimits(0-2*shifttheta, 30-2*shifttheta);
    locHist_theta_data_m2[i]->GetXaxis()->SetLimits(0-shifttheta, 30-shifttheta);
    locHist_theta_gen[i]->GetXaxis()->SetLimits(0+2*shifttheta, 30+2*shifttheta);
    locHist_theta_gen_m2[i]->GetXaxis()->SetLimits(0+shifttheta, 30+shifttheta);
    locHist_theta_data[i]->Draw();
    locHist_theta_data_m2[i]->Draw("same");
    locHist_theta_gen[i]->Draw("same");
    locHist_theta_gen_m2[i]->Draw("same");
    call->cd(3 + 3*i);
    locHist_p_data[i]->GetXaxis()->SetLimits(0-2*shiftp, 6-2*shiftp);
    locHist_p_data_m2[i]->GetXaxis()->SetLimits(0-shiftp, 6-shiftp);
    locHist_p_gen[i]->GetXaxis()->SetLimits(0+2*shiftp, 6+2*shiftp);
    locHist_p_gen_m2[i]->GetXaxis()->SetLimits(0+shiftp, 6+shiftp);
    locHist_p_data[i]->Draw();
    locHist_p_data_m2[i]->Draw("same");
    locHist_p_gen[i]->Draw("same");
    locHist_p_gen_m2[i]->Draw("same");
  }
  */


  locHist_thetap_eff_data[0][0]->SetTitle("0.5 < p < 0.75 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[1][0]->SetTitle("0.75 < p < 1.0 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[2][0]->SetTitle("1.0 < p < 1.25 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[3][0]->SetTitle("1.25 < p < 1.5 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[4][0]->SetTitle("1.5 < p < 2.0 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[5][0]->SetTitle("2.0 < p < 3.0 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[6][0]->SetTitle("3.0 < p < 4.0 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[7][0]->SetTitle("4.0 < p < 5.0 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data[8][0]->SetTitle("5.0 < p < 6.0 GeV;#theta (deg);#pi^{+} efficiency");



  /*
  TCanvas* cgrid = new TCanvas("cgrid", "pi+ efficiencies vs theta, p", 900, 900);
  cgrid->Divide(3,3);

  for(int i = 0; i < 9; i++){
    cgrid->cd(i+1);
    //locHist_thetap_eff_data[i][0]->SetTitle(";#theta (deg);#pi^{+} efficiency");
    locHist_thetap_eff_data[i][0]->SetMaximum(1.05);
    locHist_thetap_eff_data[i][0]->SetMinimum(0);
    locHist_thetap_eff_data[i][0]->SetMarkerColor(kBlue);
    locHist_thetap_eff_data[i][0]->SetLineColor(kBlue);
    locHist_thetap_eff_data[i][0]->SetMarkerStyle(kFullCircle);
    locHist_thetap_eff_data[i][0]->GetXaxis()->SetLimits(0-2*shifttheta, 30-2*shifttheta);
    locHist_thetap_eff_gen[i][0]->SetMarkerColor(kRed);
    locHist_thetap_eff_gen[i][0]->SetLineColor(kRed);
    locHist_thetap_eff_gen[i][0]->SetMarkerStyle(kOpenSquare);
    locHist_thetap_eff_gen[i][0]->GetXaxis()->SetLimits(0+2*shifttheta, 30+2*shifttheta);
    locHist_thetap_eff_data_m2[i][0]->SetMarkerColor(kAzure+10);
    locHist_thetap_eff_data_m2[i][0]->SetLineColor(kAzure+10);
    locHist_thetap_eff_data_m2[i][0]->SetMarkerStyle(kFullCircle);
    locHist_thetap_eff_data_m2[i][0]->GetXaxis()->SetLimits(0-shifttheta, 30-shifttheta);
    locHist_thetap_eff_gen_m2[i][0]->SetMarkerColor(kPink+1);
    locHist_thetap_eff_gen_m2[i][0]->SetLineColor(kPink+1);
    locHist_thetap_eff_gen_m2[i][0]->SetMarkerStyle(kOpenSquare);
    locHist_thetap_eff_gen_m2[i][0]->GetXaxis()->SetLimits(0+shifttheta, 30+shifttheta);
    locHist_thetap_eff_data[i][0]->Draw();
    locHist_thetap_eff_gen[i][0]->Draw("same");
    locHist_thetap_eff_data_m2[i][0]->Draw("same");
    locHist_thetap_eff_gen_m2[i][0]->Draw("same");
  }
  cgrid->cd(9);
  auto lgrid = new TLegend(0.4, 0.2, 0.9, 0.5);
  lgrid->AddEntry(locHist_thetap_eff_data[8][0], "Data (Method 1)");
  lgrid->AddEntry(locHist_thetap_eff_data_m2[8][0], "Data (Method 2)");
  lgrid->AddEntry(locHist_thetap_eff_gen[8][0], "MC (Method 1)");
  lgrid->AddEntry(locHist_thetap_eff_gen_m2[8][0], "MC (Method 2)");
  lgrid->Draw();

  //cgrid->Print("plots/PiPlus_PostageStamp.pdf");
  */
  /*
  locHist_thetap_eff_data_m2[0][0]->SetTitle("0.5 < p < 0.75 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data_m2[1][0]->SetTitle("0.75 < p < 1.0 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data_m2[2][0]->SetTitle("1.0 < p < 1.25 GeV;#theta (deg);#pi^{+} efficiency");
  locHist_thetap_eff_data_m2[3][0]->SetTitle("1.25 < p < 1.5 GeV;#theta (deg);#pi^{+} efficiency");

  const char* plabel[4] = {"0.5 < p < 0.75 GeV", "0.75 < p < 1.0 GeV", "1.0 < p < 1.25 GeV", "1.25 < p < 1.5 GeV"};

  TCanvas* c4 = new TCanvas("c4", "c4", 800, 800);
  c4->Divide(2,2);
  for(int i = 0; i < 4; i++){
    c4->cd(i+1);
    //locHist_thetap_eff_data[i][0]->Draw();
    //locHist_thetap_eff_gen[i][0]->Draw("same");
    locHist_thetap_eff_data_m2[i][0]->SetMaximum(1.05);
    locHist_thetap_eff_data_m2[i][0]->SetMinimum(0);
    locHist_thetap_eff_data_m2[i][0]->SetMarkerColor(kBlack);
    locHist_thetap_eff_data_m2[i][0]->SetLineColor(kBlack);
    locHist_thetap_eff_data_m2[i][0]->SetMarkerStyle(kOpenCircle);
    locHist_thetap_eff_gen_m2[i][0]->SetMarkerColor(kRed);
    locHist_thetap_eff_gen_m2[i][0]->SetLineColor(kRed);
    locHist_thetap_eff_data_m2[i][0]->Draw();
    locHist_thetap_eff_gen_m2[i][0]->Draw("same");
    TLatex *tx = new TLatex(11, 0.6, plabel[i]);
    tx->SetTextSize(0.055);
    tx->Draw("same");

  }
  c4->cd(4);
  auto* legend = new TLegend(0.5, 0.2, 0.9, 0.4);
  legend->AddEntry(locHist_thetap_eff_data_m2[3][0], "Data");
  legend->AddEntry(locHist_thetap_eff_gen_m2[3][0], "MC");
  legend->Draw();

  if(print == true)
    c4->Print("plots/PiPlusEfficiency.pdf");
  */
  /*
  TCanvas* c4[7];
  for(int im = 0; im < 1; im++){
    TString cutname = "Efficiency ratios (Data/MC) (";
    //if(im == 0)
    //cutname += "No MM^2 cut)";
    //else{
      cutname += "P(Chi2) > ";
      cutname += pcut[im];
      cutname += ")";
      //}
    c4[im] = new TCanvas(Form("c4_%d", im+1), cutname, 900, 300);
    c4[im]->Divide(3,1);

    c4[im]->cd(1);
    locHist_phi_ratio[im]->SetTitle("Data/MC;#phi (deg);#pi^{+} efficiency ratio (data/MC)");
    locHist_phi_ratio_m2[im]->SetTitle("Data/MC (Method 2)");
    locHist_phi_ratio[im]->SetMarkerStyle(kFullCircle);
    locHist_phi_ratio_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_phi_ratio[im]->SetMarkerColor(kViolet+5);
    locHist_phi_ratio[im]->SetLineColor(kViolet+5);
    locHist_phi_ratio[im]->SetMaximum(1.2);
    locHist_phi_ratio[im]->SetMinimum(0.8);
    locHist_phi_ratio_m2[im]->SetMarkerColor(kMagenta-8);
    locHist_phi_ratio_m2[im]->SetLineColor(kMagenta-8);
    //locHist_phi_ratio[im]->GetXaxis()->SetLimits(-180-shiftphi,180-shiftphi);
    //locHist_phi_ratio_m2[im]->GetXaxis()->SetLimits(-180+shiftphi,180+shiftphi);
    locHist_phi_ratio[im]->Draw();
    locHist_phi_ratio_m2[im]->Draw("same");
    //c4->BuildLegend();
    c4[im]->cd(2);
    locHist_theta_ratio[im]->SetTitle(";#theta (deg);;");
    locHist_theta_ratio[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_ratio_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_ratio[im]->SetMarkerColor(kViolet+5);
    locHist_theta_ratio[im]->SetLineColor(kViolet+5);
    locHist_theta_ratio_m2[im]->SetMarkerColor(kMagenta-8);
    locHist_theta_ratio_m2[im]->SetLineColor(kMagenta-8);
    locHist_theta_ratio[im]->SetMaximum(1.2);
    locHist_theta_ratio[im]->SetMinimum(0.8);
    //locHist_theta_ratio[im]->GetXaxis()->SetLimits(0-shifttheta,30-shifttheta);
    //locHist_theta_ratio_m2[im]->GetXaxis()->SetLimits(0+shifttheta,30+shifttheta);
    locHist_theta_ratio[im]->Draw();
    locHist_theta_ratio_m2[im]->Draw("same");
    c4[im]->cd(3);
    locHist_p_ratio[im]->SetTitle(";p (GeV);;");
    locHist_p_ratio[im]->SetMarkerStyle(kFullCircle);
    locHist_p_ratio_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_p_ratio[im]->SetMarkerColor(kViolet+5);
    locHist_p_ratio[im]->SetLineColor(kViolet+5);
    locHist_p_ratio_m2[im]->SetMarkerColor(kMagenta-8);
    locHist_p_ratio_m2[im]->SetLineColor(kMagenta-8);
    locHist_p_ratio[im]->SetMaximum(1.2);
    locHist_p_ratio[im]->SetMinimum(0.8);
    //locHist_p_ratio[im]->GetXaxis()->SetLimits(0-shiftp,6-shiftp);
    //locHist_p_ratio_m2[im]->GetXaxis()->SetLimits(0+shiftp,6+shiftp);
    locHist_p_ratio[im]->Draw();
    locHist_p_ratio_m2[im]->Draw("same");
  }

  
  TCanvas* c5[7];

  for(int im = 0; im < 1; im++){
    TString cutname = "Efficiency ratios (Method 1/2) (";
    //if(im == 0)
    //cutname += "No MM^2 cut)";
    //else{
      cutname += "P(Chi2) > ";
      cutname += pcut[im];
      cutname += ")";
      //}
      c5[im] = new TCanvas(Form("c5_%d", im+1), cutname, 900, 300);
      c5[im]->Divide(3,1);
      
      c5[im]->cd(1);
      locHist_phi_12ratio_data[im]->SetTitle("Method 1/Method 2;#phi (deg);#pi^{+} efficiency ratio (m1/m2)");
      locHist_phi_12ratio_gen[im]->SetTitle("Method 1/Method 2 (MC)");
      //locHist_phi_12ratio_data[im]->GetXaxis()->SetLimits(-180-shiftphi, 180-shiftphi);
      //locHist_phi_12ratio_gen[im]->GetXaxis()->SetLimits(-180+shiftphi, 180+shiftphi);
      locHist_phi_12ratio_data[im]->SetMarkerColor(kGreen-3);
      locHist_phi_12ratio_data[im]->SetMarkerStyle(kFullCircle);
      locHist_phi_12ratio_data[im]->SetLineColor(kGreen-3);
      locHist_phi_12ratio_data[im]->SetMinimum(0.8);
      locHist_phi_12ratio_data[im]->SetMaximum(1.2);
      locHist_phi_12ratio_data[im]->Draw();
      locHist_phi_12ratio_gen[im]->SetMarkerColor(kGreen-9);
      locHist_phi_12ratio_gen[im]->SetMarkerStyle(kFullCircle);
      locHist_phi_12ratio_gen[im]->SetLineColor(kGreen-9);
      locHist_phi_12ratio_gen[im]->Draw("same");
      //c5->BuildLegend();
      c5[im]->cd(2);
      locHist_theta_12ratio_data[im]->SetTitle(";#theta (deg);;");
      //locHist_theta_12ratio_data[im]->GetXaxis()->SetLimits(0-shifttheta, 30-shifttheta);
      //locHist_theta_12ratio_gen[im]->GetXaxis()->SetLimits(0+shifttheta, 30+shifttheta);
      locHist_theta_12ratio_data[im]->SetMarkerColor(kGreen-3);
      locHist_theta_12ratio_data[im]->SetMarkerStyle(kFullCircle);
      locHist_theta_12ratio_data[im]->SetLineColor(kGreen-3);
      locHist_theta_12ratio_data[im]->SetMinimum(0.8);
      locHist_theta_12ratio_data[im]->SetMaximum(1.2);
      locHist_theta_12ratio_data[im]->Draw();
      locHist_theta_12ratio_gen[im]->SetMarkerColor(kGreen-9);
      locHist_theta_12ratio_gen[im]->SetMarkerStyle(kFullCircle);
      locHist_theta_12ratio_gen[im]->SetLineColor(kGreen-9);
      locHist_theta_12ratio_gen[im]->Draw("same");
      c5[im]->cd(3);
      locHist_p_12ratio_data[im]->SetTitle(";p (GeV);;");
      //locHist_p_12ratio_data[im]->GetXaxis()->SetLimits(0-shiftp, 6-shiftp);
      //locHist_p_12ratio_gen[im]->GetXaxis()->SetLimits(0+shiftp, 6+shiftp);
      locHist_p_12ratio_data[im]->SetMarkerColor(kGreen-3);
      locHist_p_12ratio_data[im]->SetMarkerStyle(kFullCircle);
      locHist_p_12ratio_data[im]->SetLineColor(kGreen-3);
      locHist_p_12ratio_data[im]->SetMinimum(0.8);
      locHist_p_12ratio_data[im]->SetMaximum(1.2);
      locHist_p_12ratio_data[im]->Draw();
      locHist_p_12ratio_gen[im]->SetMarkerColor(kGreen-9);
      locHist_p_12ratio_gen[im]->SetMarkerStyle(kFullCircle);
      locHist_p_12ratio_gen[im]->SetLineColor(kGreen-9);
      locHist_p_12ratio_gen[im]->Draw("same");
  }
  */

  /*
  TCanvas* cratio_cuts = new TCanvas("cratio_cuts", "cratio_cuts", 800, 800);
  cratio_cuts->Divide(2,2);
  for(int im = 0; im < 7; im++){
    cratio_cuts->cd(1); //data/MC (method 1)
    locHist_theta_ratio[im]->SetMinimum(0.8);
    locHist_theta_ratio[im]->SetMaximum(1.2);
    locHist_theta_ratio[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_ratio[im]->SetMarkerColor(color[im]);
    locHist_theta_ratio[im]->SetLineColor(color[im]);
    locHist_theta_ratio[im]->SetTitle("Data/MC (Method 1);#theta (deg);;");
    locHist_theta_ratio[im]->Draw("same");
    cratio_cuts->cd(2); //data/MC (method 2)
    locHist_theta_ratio_m2[im]->SetMinimum(0.8);
    locHist_theta_ratio_m2[im]->SetMaximum(1.2);
    locHist_theta_ratio_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_ratio_m2[im]->SetMarkerColor(color[im]);
    locHist_theta_ratio_m2[im]->SetLineColor(color[im]);
    locHist_theta_ratio_m2[im]->SetTitle("Data/MC (Method 2);#theta (deg);;");
    locHist_theta_ratio_m2[im]->Draw("same");
    cratio_cuts->cd(3); //m1/m2 (data)
    locHist_theta_12ratio_data[im]->SetMinimum(0.8);
    locHist_theta_12ratio_data[im]->SetMaximum(1.2);
    locHist_theta_12ratio_data[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_12ratio_data[im]->SetMarkerColor(color[im]);
    locHist_theta_12ratio_data[im]->SetLineColor(color[im]);
    locHist_theta_12ratio_data[im]->SetTitle("Method 1/Method 2 (Data);#theta (deg);;");
    locHist_theta_12ratio_data[im]->Draw("same");
    cratio_cuts->cd(4); //m1/m2 (MC)
    locHist_theta_12ratio_gen[im]->SetMinimum(0.8);
    locHist_theta_12ratio_gen[im]->SetMaximum(1.2);
    locHist_theta_12ratio_gen[im]->SetMarkerStyle(kFullCircle);
    locHist_theta_12ratio_gen[im]->SetMarkerColor(color[im]);
    locHist_theta_12ratio_gen[im]->SetLineColor(color[im]);
    locHist_theta_12ratio_gen[im]->SetTitle("Method 1/Method 2 (MC);#theta (deg);;");
    locHist_theta_12ratio_gen[im]->Draw("same");
  }

  TCanvas* cratio_cuts_p = new TCanvas("cratio_cuts_p", "cratio_cuts_p", 800, 800);
  cratio_cuts_p->Divide(2,2);
  for(int im = 0; im < 7; im++){
    cratio_cuts_p->cd(1); //data/MC (method 1)
    locHist_p_ratio[im]->SetMinimum(0.8);
    locHist_p_ratio[im]->SetMaximum(1.2);
    locHist_p_ratio[im]->SetMarkerStyle(kFullCircle);
    locHist_p_ratio[im]->SetMarkerColor(color[im]);
    locHist_p_ratio[im]->SetLineColor(color[im]);
    locHist_p_ratio[im]->SetTitle("Data/MC (Method 1);p (GeV);;");
    locHist_p_ratio[im]->Draw("same");
    cratio_cuts_p->cd(2); //data/MC (method 2)
    locHist_p_ratio_m2[im]->SetMinimum(0.8);
    locHist_p_ratio_m2[im]->SetMaximum(1.2);
    locHist_p_ratio_m2[im]->SetMarkerStyle(kFullCircle);
    locHist_p_ratio_m2[im]->SetMarkerColor(color[im]);
    locHist_p_ratio_m2[im]->SetLineColor(color[im]);
    locHist_p_ratio_m2[im]->SetTitle("Data/MC (Method 2);p (GeV);;");
    locHist_p_ratio_m2[im]->Draw("same");
    cratio_cuts_p->cd(3); //m1/m2 (data)
    locHist_p_12ratio_data[im]->SetMinimum(0.8);
    locHist_p_12ratio_data[im]->SetMaximum(1.2);
    locHist_p_12ratio_data[im]->SetMarkerStyle(kFullCircle);
    locHist_p_12ratio_data[im]->SetMarkerColor(color[im]);
    locHist_p_12ratio_data[im]->SetLineColor(color[im]);
    locHist_p_12ratio_data[im]->SetTitle("Method 1/Method 2 (Data);p (GeV);;");
    locHist_p_12ratio_data[im]->Draw("same");
    cratio_cuts_p->cd(4); //m1/m2 (MC)
    locHist_p_12ratio_gen[im]->SetMinimum(0.8);
    locHist_p_12ratio_gen[im]->SetMaximum(1.2);
    locHist_p_12ratio_gen[im]->SetMarkerStyle(kFullCircle);
    locHist_p_12ratio_gen[im]->SetMarkerColor(color[im]);
    locHist_p_12ratio_gen[im]->SetLineColor(color[im]);
    locHist_p_12ratio_gen[im]->SetTitle("Method 1/Method 2 (MC);p (GeV);;");
    locHist_p_12ratio_gen[im]->Draw("same");
  }
  */
  
  double pmin[9] = {0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0};
  double pmax[9] = {0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0};
  /*
  TCanvas* cptheta_data_m1 = new TCanvas("cptheta_data_m1", "Efficiencies vs theta, p (Data) (Method 1)", 900, 900);
  cptheta_data_m1->Divide(3,3);
  for(int i = 0; i < 9; i++){ //loop over momentum bins
    cptheta_data_m1->cd(i+1);
    TString ts;
    ts += pmin[i];
    ts += " < p < ";
    ts += pmax[i];
    ts += " GeV;#theta (deg);#pi^{+} efficiency";
    for(int im = 0; im < 7; im++){ // loop over cuts
      locHist_thetap_eff_data[i][im]->SetTitle(ts);
      locHist_thetap_eff_data[i][im]->SetMaximum(1.05);
      locHist_thetap_eff_data[i][im]->SetMinimum(0);
      locHist_thetap_eff_data[i][im]->SetMarkerColor(color[im]);
      locHist_thetap_eff_data[i][im]->SetLineColor(color[im]);
      locHist_thetap_eff_data[i][im]->SetMarkerStyle(kFullCircle);
      locHist_thetap_eff_data[i][im]->Draw("same");
    }
  }
  
  TCanvas* cptheta_gen_m1 = new TCanvas("cptheta_gen_m1", "Efficiencies vs theta, p (MC) (Method 1)", 900, 900);
  cptheta_gen_m1->Divide(3,3);
  for(int i = 0; i < 9; i++){ //loop over momentum bins
    cptheta_gen_m1->cd(i+1);
    TString ts;
    ts += pmin[i];
    ts += " < p < ";
    ts += pmax[i];
    ts += " GeV;#theta (deg);#pi^{+} efficiency";
    for(int im = 0; im < 7; im++){ // loop over cuts
      locHist_thetap_eff_gen[i][im]->SetTitle(ts);
      locHist_thetap_eff_gen[i][im]->SetMaximum(1.05);
      locHist_thetap_eff_gen[i][im]->SetMinimum(0);
      locHist_thetap_eff_gen[i][im]->SetMarkerColor(color[im]);
      locHist_thetap_eff_gen[i][im]->SetLineColor(color[im]);
      locHist_thetap_eff_gen[i][im]->SetMarkerStyle(kFullCircle);
      locHist_thetap_eff_gen[i][im]->Draw("same");
    }
  }

  TCanvas* cptheta_data_m2 = new TCanvas("cptheta_data_m2", "Efficiencies vs theta, p (Data) (Method 2)", 900, 900);
  cptheta_data_m2->Divide(3,3);
  for(int i = 0; i < 9; i++){ //loop over momentum bins
    cptheta_data_m2->cd(i+1);
    TString ts;
    ts += pmin[i];
    ts += " < p < ";
    ts += pmax[i];
    ts += " GeV;#theta (deg);#pi^{+} efficiency";
    for(int im = 0; im < 7; im++){ // loop over cuts
      locHist_thetap_eff_data_m2[i][im]->SetTitle(ts);
      locHist_thetap_eff_data_m2[i][im]->SetMaximum(1.05);
      locHist_thetap_eff_data_m2[i][im]->SetMinimum(0);
      locHist_thetap_eff_data_m2[i][im]->SetMarkerColor(color[im]);
      locHist_thetap_eff_data_m2[i][im]->SetLineColor(color[im]);
      locHist_thetap_eff_data_m2[i][im]->SetMarkerStyle(kFullCircle);
      locHist_thetap_eff_data_m2[i][im]->Draw("same");
    }
  }

  TCanvas* cptheta_gen_m2 = new TCanvas("cptheta_gen_m2", "Efficiencies vs theta, p (MC) (Method 2)", 900, 900);
  cptheta_gen_m2->Divide(3,3);
  for(int i = 0; i < 9; i++){ //loop over momentum bins
    cptheta_gen_m2->cd(i+1);
    TString ts;
    ts += pmin[i];
    ts += " < p < ";
    ts += pmax[i];
    ts += " GeV;#theta (deg);#pi^{+} efficiency";
    for(int im = 0; im < 7; im++){ // loop over cuts
      locHist_thetap_eff_gen_m2[i][im]->SetTitle(ts);
      locHist_thetap_eff_gen_m2[i][im]->SetMaximum(1.05);
      locHist_thetap_eff_gen_m2[i][im]->SetMinimum(0);
      locHist_thetap_eff_gen_m2[i][im]->SetMarkerColor(color[im]);
      locHist_thetap_eff_gen_m2[i][im]->SetLineColor(color[im]);
      locHist_thetap_eff_gen_m2[i][im]->SetMarkerStyle(kFullCircle);
      locHist_thetap_eff_gen_m2[i][im]->Draw("same");
    }
  }

  
  TCanvas *c6[7];
  for(int im = 0; im < 7; im++){
    TString cutname = "Efficiencies vs theta, p (Method 2) (";
    if(im == 0)
      cutname += "No MM^2 cut)";
    else{
      cutname += "|MM^2| < ";
      cutname += masscut[im];
      cutname += " GeV^2)";
    }
    c6[im] = new TCanvas(Form("c6_%d", im+1), cutname, 900, 900);
    c6[im]->Divide(3,3);
    for(int i = 0; i < 9; i++){
      c6[im]->cd(i+1);
      TString ts;
      ts += pmin[i];
      ts += " < p < ";
      ts += pmax[i];
      ts += " GeV";
      locHist_thetap_eff_data_m2[i][im]->SetTitle(ts);
      locHist_thetap_eff_data_m2[i][im]->SetMaximum(1.2);
      locHist_thetap_eff_data_m2[i][im]->SetMinimum(0);
      locHist_thetap_eff_data_m2[i][im]->SetMarkerColor(kOrange);
      locHist_thetap_eff_data_m2[i][im]->SetLineColor(kOrange);
      locHist_thetap_eff_data_m2[i][im]->SetMarkerStyle(kFullCircle);
      locHist_thetap_eff_gen_m2[i][im]->SetMarkerColor(kBlue);
      locHist_thetap_eff_gen_m2[i][im]->SetLineColor(kBlue);
      locHist_thetap_eff_gen_m2[i][im]->SetMarkerStyle(kOpenCircle);
      locHist_thetap_eff_data_m2[i][im]->Draw();
      locHist_thetap_eff_gen_m2[i][im]->Draw("same");
    }
  }
  

  TCanvas *c7[7];
  for(int im = 0; im < 7; im++){
    TString cutname = "Efficiencies vs theta, p (Method 1) (";
    if(im == 0)
      cutname += "No MM^2 cut)";
    else{
      cutname += "|MM^2| < ";
      cutname += masscut[im];
      cutname += " GeV^2)";
    }
    c7[im] = new TCanvas(Form("c7_%d", im+1), cutname, 900, 900);
    c7[im]->Divide(3,3);
    for(int i = 0; i < 9; i++){
      c7[im]->cd(i+1);
      TString ts;
      ts += pmin[i];
      ts += " < p < ";
      ts += pmax[i];
      ts += " GeV";
      locHist_thetap_eff_data[i][im]->SetTitle(ts);
      locHist_thetap_eff_data[i][im]->SetMaximum(1.2);
      locHist_thetap_eff_data[i][im]->SetMinimum(0);
      locHist_thetap_eff_data[i][im]->SetMarkerColor(kOrange);
      locHist_thetap_eff_data[i][im]->SetLineColor(kOrange);
      locHist_thetap_eff_data[i][im]->SetMarkerStyle(kFullCircle);
      locHist_thetap_eff_gen[i][im]->SetMarkerColor(kBlue);
      locHist_thetap_eff_gen[i][im]->SetLineColor(kBlue);
      locHist_thetap_eff_gen[i][im]->SetMarkerStyle(kOpenCircle);
      locHist_thetap_eff_data[i][im]->Draw();
      locHist_thetap_eff_gen[i][im]->Draw("same");
    }
  }
  */
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat("4.2f");

  TCanvas *c2D = new TCanvas("c2D", "c2D", 1);
  c2D->Divide(2,2);

  c2D->cd(1);
  gPad->SetRightMargin(0.2);
  locHist_thetap_2D_data[0]->SetTitle("Data (Method 1)");
  locHist_thetap_2D_data[0]->SetMinimum(0.5);
  locHist_thetap_2D_data[0]->SetMaximum(1.0);
  locHist_thetap_2D_data[0]->Draw("colz");
  c2D->cd(2);
  gPad->SetRightMargin(0.2);
  locHist_thetap_2D_gen[0]->SetTitle("MC (Method 1)");
  locHist_thetap_2D_gen[0]->SetMinimum(0.5);
  locHist_thetap_2D_gen[0]->SetMaximum(1.0);
  locHist_thetap_2D_gen[0]->Draw("colz");
  c2D->cd(3);
  gPad->SetRightMargin(0.2);
  locHist_thetap_2D_data_m2[0]->SetTitle("Data (Method 2)");
  locHist_thetap_2D_data_m2[0]->SetMinimum(0.5);
  locHist_thetap_2D_data_m2[0]->SetMaximum(1.0);
  locHist_thetap_2D_data_m2[0]->Draw("colz");
  c2D->cd(4);
  gPad->SetRightMargin(0.2);
  locHist_thetap_2D_gen_m2[0]->SetTitle("MC (Method 2)");
  locHist_thetap_2D_gen_m2[0]->SetMinimum(0.5);
  locHist_thetap_2D_gen_m2[0]->SetMaximum(1.0);
  locHist_thetap_2D_gen_m2[0]->Draw("colz");

  TCanvas* c2D_ratio = new TCanvas("c2D_ratio", "c2D_ratio", 1);
  c2D_ratio->Divide(2,2);

  c2D_ratio->cd(1);
  gPad->SetRightMargin(0.2);
  locHist_thetap_ratio_2D_m1[0]->SetMinimum(0.9);
  locHist_thetap_ratio_2D_m1[0]->SetMaximum(1.1);
  locHist_thetap_ratio_2D_m1[0]->Draw("colz text");
  c2D_ratio->cd(2);
  gPad->SetRightMargin(0.2);
  locHist_thetap_ratio_2D_m2[0]->SetMinimum(0.9);
  locHist_thetap_ratio_2D_m2[0]->SetMaximum(1.1);
  locHist_thetap_ratio_2D_m2[0]->Draw("colz text");
  c2D_ratio->cd(3);
  gPad->SetRightMargin(0.2);
  locHist_thetap_ratio_2D_data[0]->SetMinimum(0.9);
  locHist_thetap_ratio_2D_data[0]->SetMaximum(1.1);
  locHist_thetap_ratio_2D_data[0]->Draw("colz text");
  c2D_ratio->cd(4);
  gPad->SetRightMargin(0.2);
  locHist_thetap_ratio_2D_gen[0]->SetMinimum(0.9);
  locHist_thetap_ratio_2D_gen[0]->SetMaximum(1.1);
  locHist_thetap_ratio_2D_gen[0]->Draw("colz text");

  TCanvas* cm1 = new TCanvas("cm1", "cm1", 1);
  cm1->cd(1);
  gPad->SetRightMargin(0.2);
  locHist_thetap_ratio_2D_m1[0]->Draw("colz text");

  TCanvas* cm2 = new TCanvas("cm2", "cm2", 1);
  cm2->cd(1);
  gPad->SetRightMargin(0.2);
  locHist_thetap_ratio_2D_m2[0]->Draw("colz text");

  // c2D->Print("plots/PiPlusEff_Absolute.pdf");
  // cm1->Print("plots/PiPlusEff_Ratio_m1.pdf");
  // cm2->Print("plots/PiPlusEff_Ratio_m2.pdf");

  // TCanvas *c2D_gen = new TCanvas("c2D_gen", "c2D_gen", 1);
  // c2D_gen->cd(1);
  // locHist_thetap_2D_gen[0]->Draw("colz");

  // TCanvas *c2D_m2 = new TCanvas("c2D_m2", "c2D_m2", 1);
  // c2D_m2->cd(1);
  // locHist_thetap_2D_data_m2[0]->Draw("colz");

  // TCanvas *c2D_gen_m2 = new TCanvas("c2D_gen_m2", "c2D_gen_m2", 1);
  // c2D_gen_m2->cd(1);
  // locHist_thetap_2D_gen_m2[0]->Draw("colz");

  if(save == true){
    for(int i = 0; i < 7; i++){
      locHist_phi_data[i]->SetName(Form("phi_data_efficiency_%d", i+1));
      locHist_phi_data_m2[i]->SetName(Form("phi_data_efficiency_m2_%d", i+1));
      locHist_phi_gen[i]->SetName(Form("phi_gen_efficiency_%d", i+1));
      locHist_phi_gen_m2[i]->SetName(Form("phi_gen_efficiency_m2_%d", i+1));
      locHist_phi_data[i]->Write();
      locHist_phi_data_m2[i]->Write();
      locHist_phi_gen[i]->Write();
      locHist_phi_gen_m2[i]->Write();
      
      locHist_theta_data[i]->Write();
      locHist_theta_data_m2[i]->Write();
      locHist_theta_gen[i]->Write();
      locHist_theta_gen_m2[i]->Write();
      
      locHist_p_data[i]->SetName(Form("p_data_efficiency_%d", i+1));
      locHist_p_data_m2[i]->SetName(Form("p_data_efficiency_m2_%d", i+1));
      locHist_p_gen[i]->SetName(Form("p_gen_efficiency_%d", i+1));
      locHist_p_gen_m2[i]->SetName(Form("p_gen_efficiency_m2_%d", i+1));
      locHist_p_data[i]->Write();
      locHist_p_data_m2[i]->Write();
      locHist_p_gen[i]->Write();
      locHist_p_gen_m2[i]->Write();
      
      locHist_thetap_2D_data[i]->Write();
      locHist_thetap_2D_gen[i]->Write();
      locHist_thetap_2D_data_m2[i]->Write();
      locHist_thetap_2D_gen_m2[i]->Write();
    }
  }

  return;
}
