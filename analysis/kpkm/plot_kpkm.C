
TString NT("ntFSGlueX_MODECODE");
TString DEFAULT_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("100_110000")->addCategory("kpkm");

  // FIXED CUTS
  FSCut::defineCut("unusedE","EnUnusedSh<0.1");
  FSCut::defineCut("unusedTracks","NumUnusedTracks<1");
  FSCut::defineCut("z","ProdVz>=51.2&&ProdVz<=78.8");
  FSCut::defineCut("MM2","abs(RMASS2(GLUEXTARGET,B,-1,-2,-3))<0.05");
  FSCut::defineCut("eBeam","EnPB>8.0");
  FSCut::defineCut("chi2","Chi2DOF<5");
  FSCut::defineCut("chi2rank","Chi2RankHybrid==1");
  
  FSCut::defineCut("t","abs(-1*MASS2([proton],-GLUEXTARGET))<1.0");
  
  FSCut::defineCut("GlueXI","Run<70000");
  
  // Select RF signal region
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");

  // SET SOME DEFAULT SKIM CUTS
  DEFAULT_CUTS = "eBeam,chi2,rf,unusedE,unusedTracks,z,MM2,t,chi2rank";
}

void plot_kpkm(bool bggen=false){

  setup();

  // with mass constraints
  TString FND_DATA = "tree_kpkm__B4_Hybrid_SKIM_*.root";
  TString FND_BGGEN = "tree_kpkm__B4_BGGEN_Hybrid_SKIM.root";
  TString FND_MC = "tree_kpkm__B4_Hybrid_SKIM_*.root"; // temporarily plot data twice until MC available "tree_kpkm__B4_PSMC_Hybrid_SKIM.root";
  TString CATEGORY = "kpkm";

  FSTree::addFriendTree("Chi2RankHybrid");

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
    TH1F* hEnUnusedSh = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","EnUnusedSh","(100,0.0,1.0)",Form("CUT(%s)",CUTS.Data()));
    hEnUnusedSh->SetXTitle("E_{unused}  [GeV/c^{2}]");
    hEnUnusedSh->SetYTitle("Events");
    hEnUnusedSh->Draw();
    TH1F* hEnUnusedShMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","EnUnusedSh","(100,0.0,1.0)",Form("CUT(%s)",CUTS.Data()));
    hEnUnusedShMC->Scale(hEnUnusedSh->GetMaximum()/hEnUnusedShMC->GetMaximum());
    hEnUnusedShMC->SetMarkerColor(kMagenta);
    hEnUnusedShMC->Draw("same");
    TLine* cutUnusedE = new TLine(0.1,0,0.1,hEnUnusedSh->GetMaximum());
    cutUnusedE->SetLineColor(kRed);
    cutUnusedE->Draw("same");

    c1->cd(2);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",z","");
    TH1F* hProdVz = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","ProdVz","(100,0.,100.0)",Form("CUT(%s)",CUTS.Data()));
    hProdVz->SetXTitle("ProdVz  [cm]");
    hProdVz->SetYTitle("Events");
    hProdVz->Draw();
    TH1F* hProdVzMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","ProdVz","(100,0.,100.0)",Form("CUT(%s)",CUTS.Data()));
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
    TH1F* htk = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","abs(-1*MASS2([proton],-GLUEXTARGET))","(100,0,5)",Form("CUT(%s)",CUTS.Data()));
    htk->SetXTitle("|-t| [GeV^{2}]");
    htk->SetYTitle("Entries");
    htk->Draw();
    htk->Fit("expo","","",0.5,1.0);
    TH1F* htkMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","abs(-1*MASS2([proton],-GLUEXTARGET))","(100,0,5)",Form("CUT(%s)",CUTS.Data()));
    htkMC->Scale(htk->GetMaximum()/htkMC->GetMaximum());
    htkMC->SetMarkerColor(kMagenta);
    htkMC->Draw("same");
    htkMC->Fit("expo","","",0.5,1.0);

    c1->cd(4);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll("eBeam,","");
    TH1F* hEnPB = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","EnPB","(125,5,12)",Form("CUT(%s)",CUTS.Data()));
    hEnPB->SetXTitle("E_{beam} [GeV]");
    hEnPB->SetYTitle("Entries");
    hEnPB->Draw();
    TH1F* hEnPBMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","EnPB","(125,5,12)",Form("CUT(%s)",CUTS.Data()));
    hEnPBMC->Scale(hEnPB->GetMaximum()/hEnPBMC->GetMaximum());
    hEnPBMC->SetMarkerColor(kMagenta);
    hEnPBMC->Draw("same");

    c1->cd(5);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",MM2","");
    TH1F* hMM2 = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","RMASS2(GLUEXTARGET,B,-1,-2,-3)","(100,-0.1,0.1)",Form("CUT(%s)",CUTS.Data()));
    hMM2->SetXTitle("Missing Mass Squared [GeV^{2}]");
    hMM2->SetYTitle("Entries");
    hMM2->Draw();
    TH1F* hMM2MC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","RMASS2(GLUEXTARGET,B,-1,-2,-3)","(100,-0.1,0.1)",Form("CUT(%s)",CUTS.Data()));
    hMM2MC->Scale(hMM2->GetMaximum()/hMM2MC->GetMaximum());
    hMM2MC->SetMarkerColor(kMagenta);
    hMM2MC->Draw("same");

    c1->cd(6);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",chi2,",",");
    TH1F* hChi2DOF = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    hChi2DOF->SetXTitle("#chi^{2}/dof");
    hChi2DOF->SetYTitle("Events");
    hChi2DOF->Draw();
    TH1F* hChi2DOFMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    hChi2DOFMC->Scale(hChi2DOF->GetMaximum()/hChi2DOFMC->GetMaximum());
    hChi2DOFMC->SetMarkerColor(kMagenta);
    hChi2DOFMC->Draw("same");

    if(bggen) {
	    cb->cd(1);    
	    TH1F* hChi2DOF_BGGEN = FSModeHistogram::getTH1F(FND_BGGEN,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
	    hChi2DOF_BGGEN->Draw();
	    FSModeHistogram::drawMCComponentsSame(FND_BGGEN,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    }
    
    // Some mass spectra...
    CUTS = DEFAULT_CUTS;
    TH1F* hMkpkm = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","MASS([K+],[K-])","(200,0.95,1.5)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMkpkmMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","MASS([K+],[K-])","(200,0.95,1.5)",Form("CUT(%s)",CUTS.Data()));
    hMkpkmMC->Scale(hMkpkm->GetMaximum()/hMkpkmMC->GetMaximum());
    hMkpkmMC->SetMarkerColor(kMagenta);

    if(bggen) {
	    cb->cd(2);    
	    TH1F* hMkpkm_BGGEN = FSModeHistogram::getTH1F(FND_BGGEN,NT,"kpkm","MASS([K+],[K-])","(200,0.95,1.5)",Form("CUT(%s)",CUTS.Data()));
	    hMkpkm_BGGEN->Draw();
	    FSModeHistogram::drawMCComponentsSame(FND_BGGEN,NT,"kpkm","MASS([K+],[K-])","(200,0.95,1.5)",Form("CUT(%s)",CUTS.Data()));
    }

    CUTS = DEFAULT_CUTS;

    TH1F* hMkmp = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","MASS([proton],[K-])","(100,1.2,2.2)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMkmpMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","MASS([proton],[K-])","(100,1.2,2.2)",Form("CUT(%s)",CUTS.Data()));
    hMkmpMC->Scale(hMkmp->GetMaximum()/hMkmpMC->GetMaximum());
    hMkmpMC->SetMarkerColor(kMagenta);
    
    FSHistogram::dumpHistogramCache();
    
    TCanvas* c11 = new TCanvas("c11","c11",1000,500);
    c11->Divide(2,1);
    c11->cd(1);
    hMkpkm->Draw();
    hMkpkmMC->Draw("same");

    c11->cd(2);
    hMkmp->Draw();
    hMkmpMC->Draw("same");

    // Save histograms for making plots in separate macro
    TFile *f = new TFile("out_kpkm.root","recreate");

    hMkpkm->SetName("hMkpkm"); hMkpkm->Write();
    hMkpkmMC->SetName("hMkpkmMC"); hMkpkmMC->Write();

    hMkmp->SetName("hMkmp"); hMkmp->Write();
    hMkmpMC->SetName("hMkmpMC"); hMkmpMC->Write();

    return;
}
