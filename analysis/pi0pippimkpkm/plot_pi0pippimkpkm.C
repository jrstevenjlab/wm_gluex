
TString NT("ntFSGlueX_MODECODE");
TString DEFAULT_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("100_110111")->addCategory("pi0pippimkpkm");

  // FIXED CUTS
  FSCut::defineCut("unusedE","EnUnusedSh<0.1");
  FSCut::defineCut("unusedTracks","NumUnusedTracks<1");
  FSCut::defineCut("z","ProdVz>=51.2&&ProdVz<=78.8");
  FSCut::defineCut("MM2","abs(RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5,-6))<0.05");
  FSCut::defineCut("eBeam","EnPB>8.0");
  FSCut::defineCut("chi2","Chi2DOF<5");
  for(int i=6; i<=10; i++) FSCut::defineCut(Form("chi2_%d",i),Form("Chi2DOF<%d",i));
  FSCut::defineCut("chi2rank","Chi2Rank==1");
  
  FSCut::defineCut("t","abs(-1*MASS2([proton],-GLUEXTARGET))<1.0");
  
  FSCut::defineCut("Phi","abs(MASS([K+],[K-]) - 1.019) < 0.015");
  FSCut::defineCut("PhiSB","abs(MASS([K+],[K-]) - 1.019) > 0.025 && abs(MASS([K+],[K-]) - 1.019) < 0.055");
  FSCut::defineCut("Phi2DSB","abs(MASS([K+],[K-]) - 1.019) < 0.015","abs(MASS([K+],[K-]) - 1.019) > 0.025 && abs(MASS([K+],[K-]) - 1.019) < 0.055",0.5);

  FSCut::defineCut("Eta","abs(MASS([pi0],[pi+],[pi-]) - 0.548) < 0.020");

  FSCut::defineCut("KshortVeto","abs(MASS([pi+],[pi-]) - 0.4976) > 0.04");
  FSCut::defineCut("DeltaVeto","MASS([pi+],[proton]) > 1.35");

  FSCut::defineCut("GlueXI","Run<70000");
  
  // Select RF signal region
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");

  // SET SOME DEFAULT SKIM CUTS
  DEFAULT_CUTS = "eBeam,chi2,rf,z,MM2,t,chi2rank,Phi,KshortVeto,DeltaVeto,Eta";
  // ,unusedE,unusedTracks
  // ,Lambda1520,KpipiMass
}

void plot_pi0pippimkpkm(bool bggen=false){

  setup();

  // with mass constraints
  TString FND_DATA = "tree_pi0pippimkpkm__B4_BestChi2_SKIM.root";
  TString FND_BGGEN = "tree_pi0pippimkpkm__B4_BGGEN_BestChi2_SKIM.root";
  TString FND_MC = "tree_pi0pippimkpkm__B4_ETAMC_BestChi2_SKIM.root";
  TString FND_ETAMC = "tree_pi0pippimkpkm__B4_ETAMC_BestChi2_SKIM.root";
  TString FND_OMEGAMC = "tree_pi0pippimkpkm__B4_OMEGAMC_BestChi2_SKIM.root";
  TString CATEGORY = "pi0pippimkpkm";

  FSTree::addFriendTree("Chi2Rank");

  setup();
  system("rm -rf plots");  system("mkdir plots");

  TCanvas* cb = new TCanvas("cb","cb",1000,600);
  cb->Divide(3,2);

    TString CUTS;
    TCanvas* c1 = new TCanvas("c1","c1",1000,600);
    c1->Divide(3,2);
    c1->cd(1);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",unusedE","");
    TH1F* hEnUnusedSh = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","EnUnusedSh","(100,0.0,1.0)",Form("CUT(%s)",CUTS.Data()));
    hEnUnusedSh->SetXTitle("E_{unused}  [GeV/c^{2}]");
    hEnUnusedSh->SetYTitle("Events");
    hEnUnusedSh->Draw();
    TH1F* hEnUnusedShMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","EnUnusedSh","(100,0.0,1.0)",Form("CUT(%s)",CUTS.Data()));
    hEnUnusedShMC->Scale(hEnUnusedSh->GetMaximum()/hEnUnusedShMC->GetMaximum());
    hEnUnusedShMC->SetMarkerColor(kMagenta);
    hEnUnusedShMC->Draw("same");
    TLine* cutUnusedE = new TLine(0.1,0,0.1,hEnUnusedSh->GetMaximum());
    cutUnusedE->SetLineColor(kRed);
    cutUnusedE->Draw("same");

    c1->cd(2);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",z","");
    TH1F* hProdVz = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","ProdVz","(100,0.,100.0)",Form("CUT(%s)",CUTS.Data()));
    hProdVz->SetXTitle("ProdVz  [cm]");
    hProdVz->SetYTitle("Events");
    hProdVz->Draw();
    TH1F* hProdVzMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","ProdVz","(100,0.,100.0)",Form("CUT(%s)",CUTS.Data()));
    hProdVzMC->Scale(hProdVz->GetMaximum()/hProdVzMC->GetMaximum());
    hProdVzMC->SetMarkerColor(kMagenta);
    hProdVzMC->Draw("same");
    TLine* cutVz_low = new TLine(52,0,52,hProdVz->GetMaximum());
    cutVz_low->SetLineColor(kRed);
    cutVz_low->Draw("same");
    TLine* cutVz_hi = new TLine(78,0,78,hProdVz->GetMaximum());
    cutVz_hi->SetLineColor(kRed);
    cutVz_hi->Draw("same");

    c1->cd(3);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",t","");
    TH1F* htk = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","abs(-1*MASS2([proton],-GLUEXTARGET))","(100,0,5)",Form("CUT(%s)",CUTS.Data()));
    htk->SetXTitle("|-t| [GeV^{2}]");
    htk->SetYTitle("Entries");
    htk->Draw();
    htk->Fit("expo","","",0.5,1.0);
    TH1F* htkMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","abs(-1*MASS2([proton],-GLUEXTARGET))","(100,0,5)",Form("CUT(%s)",CUTS.Data()));
    htkMC->Scale(htk->GetMaximum()/htkMC->GetMaximum());
    htkMC->SetMarkerColor(kMagenta);
    htkMC->Draw("same");
    htkMC->Fit("expo","","",0.5,1.0);

    c1->cd(4);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll("eBeam,","");
    TH1F* hEnPB = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","EnPB","(125,5,12)",Form("CUT(%s)",CUTS.Data()));
    hEnPB->SetXTitle("E_{beam} [GeV]");
    hEnPB->SetYTitle("Entries");
    hEnPB->Draw();
    TH1F* hEnPBMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","EnPB","(125,5,12)",Form("CUT(%s)",CUTS.Data()));
    hEnPBMC->Scale(hEnPB->GetMaximum()/hEnPBMC->GetMaximum());
    hEnPBMC->SetMarkerColor(kMagenta);
    hEnPBMC->Draw("same");

    c1->cd(5);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",MM2","");
    TH1F* hMM2 = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5,-6)","(100,-0.1,0.1)",Form("CUT(%s)",CUTS.Data()));
    hMM2->SetXTitle("Missing Mass Squared [GeV^{2}]");
    hMM2->SetYTitle("Entries");
    hMM2->Draw();
    TH1F* hMM2MC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5,-6)","(100,-0.1,0.1)",Form("CUT(%s)",CUTS.Data()));
    hMM2MC->Scale(hMM2->GetMaximum()/hMM2MC->GetMaximum());
    hMM2MC->SetMarkerColor(kMagenta);
    hMM2MC->Draw("same");

    c1->cd(6);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",chi2,",",");
    TH1F* hChi2DOF = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    hChi2DOF->SetXTitle("#chi^{2}/dof");
    hChi2DOF->SetYTitle("Events");
    hChi2DOF->Draw();
    TH1F* hChi2DOFMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    hChi2DOFMC->Scale(hChi2DOF->GetMaximum()/hChi2DOFMC->GetMaximum());
    hChi2DOFMC->SetMarkerColor(kMagenta);
    hChi2DOFMC->Draw("same");

    if(bggen) {
	    cb->cd(1);    
	    TH1F* hChi2DOF_BGGEN = FSModeHistogram::getTH1F(FND_BGGEN,NT,"pi0pippimkpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
	    hChi2DOF_BGGEN->Draw();
	    FSModeHistogram::drawMCComponentsSame(FND_BGGEN,NT,"pi0pippimkpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    }
    
    // Some mass spectra...
    CUTS = DEFAULT_CUTS;
    TH1F* hMkpkmSignal = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([K+],[K-])","(200,0.95,1.2)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMkpkmMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","MASS([K+],[K-])","(200,0.95,1.2)",Form("CUT(%s)",CUTS.Data()));
    hMkpkmMC->Scale(hMkpkmSignal->GetMaximum()/hMkpkmMC->GetMaximum());
    hMkpkmMC->SetMarkerColor(kMagenta);

    CUTS.ReplaceAll(",Phi",",PhiSB");
    TH1F* hMkpkmSB = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([K+],[K-])","(200,0.95,1.2)",Form("CUT(%s)",CUTS.Data()));
    hMkpkmSB->SetMarkerColor(kBlue);

    CUTS.ReplaceAll(",PhiSB","");
    TH1F* hMkpkm = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([K+],[K-])","(200,0.95,1.2)",Form("CUT(%s)",CUTS.Data()));

    if(bggen) {
	    cb->cd(2);    
	    TH1F* hMkpkm_BGGEN = FSModeHistogram::getTH1F(FND_BGGEN,NT,"pi0pippimkpkm","MASS([K+],[K-])","(200,0.95,1.2)",Form("CUT(%s)",CUTS.Data()));
	    hMkpkm_BGGEN->Draw();
	    FSModeHistogram::drawMCComponentsSame(FND_BGGEN,NT,"pi0pippimkpkm","MASS([K+],[K-])","(200,0.95,1.2)",Form("CUT(%s)",CUTS.Data()));
    }
    
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",KshortVeto","");
    TH1F* hMpi0pippimSignal = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMpi0pippimMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));
    hMpi0pippimMC->Scale(hMpi0pippimSignal->GetMaximum()/hMpi0pippimMC->GetMaximum());
    hMpi0pippimMC->SetMarkerColor(kMagenta);
    hMpi0pippimSignal->SetMarkerColor(kRed);

    // background MC
    TH1F* hMpi0pippimEtaMC = FSModeHistogram::getTH1F(FND_ETAMC,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));
    hMpi0pippimEtaMC->SetMarkerColor(kMagenta);
    TH1F* hMpi0pippimOmegaMC = FSModeHistogram::getTH1F(FND_OMEGAMC,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));
    hMpi0pippimOmegaMC->SetMarkerColor(kCyan);
    
    CUTS.ReplaceAll(",Phi",",PhiSB");
    TH1F* hMpi0pippimSB = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));
    hMpi0pippimSB->SetMarkerColor(kBlue);

    CUTS.ReplaceAll(",PhiSB","");
    TH1F* hMpi0pippim = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));

    if(bggen) {
	    cb->cd(3);    
	    TH1F* hMpi0pippim_BGGEN = FSModeHistogram::getTH1F(FND_BGGEN,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));
	    hMpi0pippim_BGGEN->Draw();
	    FSModeHistogram::drawMCComponentsSame(FND_BGGEN,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-])","(50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));
    }
    
    CUTS = DEFAULT_CUTS;
    TH1F* hMkpkmpi0pippim = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([K+],[K-],[pi0],[pi+],[pi-])","(100,1.5,2.5)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMkpkmpi0pippimMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","MASS([K+],[K-],[pi0],[pi+],[pi-])","(100,1.5,2.5)",Form("CUT(%s)",CUTS.Data()));
    hMkpkmpi0pippimMC->Scale(hMkpkmpi0pippim->GetMaximum()/hMkpkmpi0pippimMC->GetMaximum());
    hMkpkmpi0pippimMC->SetMarkerColor(kMagenta);
    
    CUTS.ReplaceAll(",Phi",",PhiSB");
    TH1F* hMkpkmpi0pippimSB = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([K+],[K-],[pi0],[pi+],[pi-])","(100,1.5,2.5)",Form("CUT(%s)",CUTS.Data()));
    hMkpkmpi0pippimSB->SetMarkerColor(kBlue);
    
    ////////////////////////
    // baryon backgrounds //
    ////////////////////////
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",DeltaVeto","");
    TH1F* hMpipp = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([proton],[pi+])","(100,1.0,2.0)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMpippMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","MASS([proton],[pi+])","(100,1.0,2.0)",Form("CUT(%s)",CUTS.Data()));
    hMpippMC->Scale(hMpipp->GetMaximum()/hMpippMC->GetMaximum());
    hMpippMC->SetMarkerColor(kMagenta);

    CUTS = DEFAULT_CUTS;
    TH1F* hMpimp = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([proton],[pi-])","(100,1.0,2.0)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMpimpMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","MASS([proton],[pi-])","(100,1.0,2.0)",Form("CUT(%s)",CUTS.Data()));
    hMpimpMC->Scale(hMpimp->GetMaximum()/hMpimpMC->GetMaximum());
    hMpimpMC->SetMarkerColor(kMagenta);

    CUTS = DEFAULT_CUTS;
    TH1F* hMkmp = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimkpkm","MASS([proton],[K-])","(100,1.2,2.2)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMkmpMC = FSModeHistogram::getTH1F(FND_MC,NT,"pi0pippimkpkm","MASS([proton],[K-])","(100,1.2,2.2)",Form("CUT(%s)",CUTS.Data()));
    hMkmpMC->Scale(hMkmp->GetMaximum()/hMkmpMC->GetMaximum());
    hMkmpMC->SetMarkerColor(kMagenta);

    // correlations
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",Phi","");
    CUTS.ReplaceAll(",Eta","");
    CUTS.ReplaceAll(",KshortVeto","");
    TH2F* hMpi0pippimkpkm_kpkm = FSModeHistogram::getTH2F(FND_DATA,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-],[K+],[K-]):MASS([K+],[K-])","(100,0.98,1.1,40,1.4,3.2)",Form("CUT(%s)",CUTS.Data()));
    TH2F* hMpi0pippim_kpkm = FSModeHistogram::getTH2F(FND_DATA,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-]):MASS([K+],[K-])","(100,0.98,1.1,50,0.4,0.9)",Form("CUT(%s)",CUTS.Data()));

    CUTS += ",Eta";
    TH2F* hMetakpkm_kpkm = FSModeHistogram::getTH2F(FND_DATA,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-],[K+],[K-]):MASS([K+],[K-])","(100,0.98,1.1,80,1.4,3.2)",Form("CUT(%s)",CUTS.Data()));
    CUTS.ReplaceAll(",Eta",",EtaSB");
    TH2F* hMetaSBkpkm_kpkm = FSModeHistogram::getTH2F(FND_DATA,NT,"pi0pippimkpkm","MASS([pi0],[pi+],[pi-],[K+],[K-]):MASS([K+],[K-])","(100,0.98,1.1,80,1.4,3.2)",Form("CUT(%s)",CUTS.Data()));
    
    FSHistogram::dumpHistogramCache();
    
    TCanvas* c11 = new TCanvas("c11","c11",1000,600);
    c11->Divide(3,2);
    c11->cd(1);
    hMkpkm->Draw();
    hMkpkmSignal->Draw("same");
    hMkpkmMC->Draw("same");
    hMkpkmSB->Draw("same");
    
    c11->cd(2);
    hMpi0pippim->Draw();
    hMpi0pippimSignal->Draw("same");
    hMpi0pippimMC->Draw("same");
    hMpi0pippimSB->Draw("same");
    hMpi0pippimSignal->Draw("same");
    hMpi0pippimEtaMC->Draw("same");
    hMpi0pippimOmegaMC->Draw("same");
    
    c11->cd(3);
    hMkpkmpi0pippim->Draw();
    hMkpkmpi0pippimMC->Draw("same");
    hMkpkmpi0pippimSB->Draw("same");

    c11->cd(4);
    hMpipp->Draw();
    hMpippMC->Draw("same");

    c11->cd(5);
    hMpimp->Draw();
    hMpimpMC->Draw("same");

    c11->cd(6);
    hMkmp->Draw();
    hMkmpMC->Draw("same");

    // Save histograms for making plots in separate macro
    // Set axis labels before writing to file
    hMpi0pippimkpkm_kpkm->SetXTitle("Mass(K^{+}K^{-}) [GeV]");
    hMpi0pippimkpkm_kpkm->SetYTitle("Mass(#pi^{0}#pi^{+}#pi^{-}K^{+}K^{-}) [GeV]");
    hMpi0pippimkpkm_kpkm->SetZTitle("Counts");
    hMpi0pippim_kpkm->SetXTitle("Mass(K^{+}K^{-}) [GeV]");
    hMpi0pippim_kpkm->SetYTitle("Mass(#pi^{0}#pi^{+}#pi^{-}) [GeV]");
    hMpi0pippim_kpkm->SetZTitle("Counts");
    hMpi0pippimEtaMC->SetXTitle("Mass(#pi^{0}#pi^{+}#pi^{-}) [GeV]");
    hMpi0pippimEtaMC->SetYTitle("Counts");
    hMpi0pippimOmegaMC->SetXTitle("Mass(#pi^{0}#pi^{+}#pi^{-}) [GeV]");
    hMpi0pippimOmegaMC->SetYTitle("Counts");

    TFile *f = new TFile("hist_pi0pippimkpkm.root","recreate");
    hMpi0pippimkpkm_kpkm->SetName("hMpi0pippimkpkm_kpkm"); hMpi0pippimkpkm_kpkm->Write();
    hMpi0pippim_kpkm->SetName("hMpi0pippim_kpkm"); hMpi0pippim_kpkm->Write();

    hMetakpkm_kpkm->SetName("hMeta3pikpkm_kpkm"); hMetakpkm_kpkm->Write();
    hMetaSBkpkm_kpkm->SetName("hMetaSBkpkm_kpkm"); hMetaSBkpkm_kpkm->Write();

    hMpi0pippimEtaMC->SetName("hM3pi_EtaMC"); hMpi0pippimEtaMC->Write();
    hMpi0pippimOmegaMC->SetName("hM3pi_OmegaMC"); hMpi0pippimOmegaMC->Write();
/*
    hMpippim_kpkm_PhiPiPiMass1->SetName("hMpipi_kpkm_PhiPiPi1.60_1.80"); hMpippim_kpkm_PhiPiPiMass1->Write();
    hMpippim_kpkm_PhiPiPiMass2->SetName("hMpipi_kpkm_PhiPiPi1.80_2.00"); hMpippim_kpkm_PhiPiPiMass2->Write();
    hMpippim_kpkm_PhiPiPiMass3->SetName("hMpipi_kpkm_PhiPiPi2.00_2.20"); hMpippim_kpkm_PhiPiPiMass3->Write();
    hMpippim_kpkm_PhiPiPiMass4->SetName("hMpipi_kpkm_PhiPiPi2.20_2.40"); hMpippim_kpkm_PhiPiPiMass4->Write();
    hMpippim_kpkm_PhiPiPiMass5->SetName("hMpipi_kpkm_PhiPiPi2.40_2.60"); hMpippim_kpkm_PhiPiPiMass5->Write();
    hMpippim_kpkm_PhiPiPiMass6->SetName("hMpipi_kpkm_PhiPiPi2.60_2.80"); hMpippim_kpkm_PhiPiPiMass6->Write();

    hMpi0pippimkpkm_kpkm_PiPiMass1->SetName("hMpipikpkm_kpkm_PiPi0.45_0.55"); hMpi0pippimkpkm_kpkm_PiPiMass1->Write();
    hMpi0pippimkpkm_kpkm_PiPiMass2->SetName("hMpipikpkm_kpkm_PiPi0.55_0.65"); hMpi0pippimkpkm_kpkm_PiPiMass2->Write();
    hMpi0pippimkpkm_kpkm_PiPiMass3->SetName("hMpipikpkm_kpkm_PiPi0.65_0.90"); hMpi0pippimkpkm_kpkm_PiPiMass3->Write();
    hMpi0pippimkpkm_kpkm_PiPiMass4->SetName("hMpipikpkm_kpkm_PiPi0.90_1.00"); hMpi0pippimkpkm_kpkm_PiPiMass4->Write();
    hMpi0pippimkpkm_kpkm_PiPiMass5->SetName("hMpipikpkm_kpkm_PiPi1.00_1.10"); hMpi0pippimkpkm_kpkm_PiPiMass5->Write();
    hMpi0pippimkpkm_kpkm_PiPiMass6->SetName("hMpipikpkm_kpkm_PiPi1.10_1.20"); hMpi0pippimkpkm_kpkm_PiPiMass6->Write();
    hMpi0pippimkpkm_kpkm_PiPiMass7->SetName("hMpipikpkm_kpkm_PiPi1.20_1.30"); hMpi0pippimkpkm_kpkm_PiPiMass7->Write();
    hMpi0pippimkpkm_kpkm_PiPiMass8->SetName("hMpipikpkm_kpkm_PiPi1.30_1.40"); hMpi0pippimkpkm_kpkm_PiPiMass8->Write();

    hMpippimEtaMC->SetName("hMpipi_EtaMC"); hMpippimEtaMC->Write();

    for(int i=6; i<=10; i++) {
	    hMpi0pippimkpkm_kpkm_chi2[i-6]->SetName(Form("hMpipikpkm_kpkm_chi2_%d",i));
	    hMpi0pippimkpkm_kpkm_chi2[i-6]->Write();
	    hMpippim_kpkm_chi2[i-6]->SetName(Form("hMpipi_kpkm_chi2_%d",i));
	    hMpippim_kpkm_chi2[i-6]->Write();
    }
*/
    
#if 0
    hMkmp->SetName("hMkmp"); hMkmp->Write();
    hMkmpMC->SetName("hMkmpMC"); hMkmpMC->Write();
    hMkmpSB->SetName("hMkmpSB"); hMkmpSB->Write();
    
    hMpippim->SetName("hMpippim"); hMpippim->Write();
    hMpippimMC->SetName("hMpippimMC"); hMpippimMC->Write();
    hMpippimSB->SetName("hMpippimSB"); hMpippimSB->Write();
    
    hMkppim->SetName("hMkppim"); hMkppim->Write();
    hMkppimMC->SetName("hMkppimMC"); hMkppimMC->Write();
    hMkppimSB->SetName("hMkppimSB"); hMkppimSB->Write();
    
    hMkppippim->SetName("hMkppippim"); hMkppippim->Write();
    hMkppippimMC->SetName("hMkppippimMC"); hMkppippimMC->Write();
    hMkppippimSB->SetName("hMkppippimSB"); hMkppippimSB->Write();
#endif    

    return;
}
