
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
  FSCut::defineCut("chi2rank","Chi2Rank==1");
  FSCut::defineCut("chi2rankglobal","Chi2RankGlobal==1");
  FSCut::defineCut("hybridchi2rank","HybridChi2Rank==1");
  FSCut::defineCut("hybridchi2rankglobal","HybridChi2RankGlobal==1");
  
  FSCut::defineCut("t","abs(-1*MASS2([proton],-GLUEXTARGET))<1.0");
  
  FSCut::defineCut("GlueXI","Run<70000");
  
  // Hybrid accidental subtraction cut
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004","abs(RFDeltaT) > 2.004 && abs(RFDeltaT) < 18.036",0.125);

  // SET SOME DEFAULT SKIM CUTS
  DEFAULT_CUTS = "eBeam,chi2,unusedE,unusedTracks,z,MM2,t";
}

void plot_hybrid_bestChi2(){

  setup();

  TString FND_DATA_BestChi2 = "tree_kpkm__B4_BestChi2_SKIM_*.root";
  TString FND_MC_BestChi2 = "tree_kpkm__B4_SIGMC_BestChi2_SKIM_*.root";

  TString FND_DATA_Hybrid = "tree_kpkm__B4_BestHybridChi2_SKIM_*.root";
  TString FND_MC_Hybrid = "tree_kpkm__B4_SIGMC_BestHybridChi2_SKIM_*.root";


  FSTree::addFriendTree("Chi2Rank");
  FSTree::addFriendTree("HybridChi2Rank");

  setup();
  system("rm -rf plots");  system("mkdir plots");

  TString CUTS;
  TCanvas* c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(3,2);
  c1->cd(1);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",chi2,",",");
  TH1F* hChi2DOF_BestChi2 = FSModeHistogram::getTH1F(FND_DATA_BestChi2,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
  hChi2DOF_BestChi2->SetXTitle("#chi^{2}/dof");
  hChi2DOF_BestChi2->SetYTitle("Events");
  hChi2DOF_BestChi2->Draw();
  TH1F* hChi2DOFMC_BestChi2 = FSModeHistogram::getTH1F(FND_MC_BestChi2,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
  double scale_bestChi2 = hChi2DOF_BestChi2->GetMaximum()/hChi2DOFMC_BestChi2->GetMaximum();
  hChi2DOFMC_BestChi2->Scale(scale_bestChi2);
  hChi2DOFMC_BestChi2->SetMarkerColor(kMagenta);
  hChi2DOFMC_BestChi2->Draw("same");

  // add global chi2 rank cut
  CUTS += ",chi2rankglobal";
  TH1F* hChi2DOF_BestChi2_Global = FSModeHistogram::getTH1F(FND_DATA_BestChi2,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
  hChi2DOF_BestChi2_Global->SetMarkerColor(kGreen);
  hChi2DOF_BestChi2_Global->Draw("same");
  TH1F* hChi2DOFMC_BestChi2_Global = FSModeHistogram::getTH1F(FND_MC_BestChi2,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
  hChi2DOFMC_BestChi2_Global->Scale(scale_bestChi2);
  hChi2DOFMC_BestChi2_Global->SetMarkerColor(kBlue);
  hChi2DOFMC_BestChi2_Global->Draw("same");
  
  c1->cd(2);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",chi2,",",");
  TH1F* hChi2DOF_Hybrid = FSModeHistogram::getTH1F(FND_DATA_Hybrid,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)*CUTWT(rf)",CUTS.Data()));
  hChi2DOF_Hybrid->SetXTitle("#chi^{2}/dof");
  hChi2DOF_Hybrid->SetYTitle("Events");
  hChi2DOF_Hybrid->Draw();
  TH1F* hChi2DOFMC_Hybrid = FSModeHistogram::getTH1F(FND_MC_Hybrid,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)*CUTWT(rf)",CUTS.Data()));
  double scale_hybrid = hChi2DOF_Hybrid->GetMaximum()/hChi2DOFMC_Hybrid->GetMaximum();
  hChi2DOFMC_Hybrid->Scale(scale_hybrid);
  hChi2DOFMC_Hybrid->SetMarkerColor(kMagenta);
  hChi2DOFMC_Hybrid->Draw("same");

  // add global chi2 rank cut
  CUTS += ",hybridchi2rankglobal";
  TH1F* hChi2DOF_Hybrid_Global = FSModeHistogram::getTH1F(FND_DATA_Hybrid,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)*CUTWT(rf)",CUTS.Data()));
  hChi2DOF_Hybrid_Global->SetMarkerColor(kGreen);
  hChi2DOF_Hybrid_Global->Draw("same");
  TH1F* hChi2DOFMC_Hybrid_Global = FSModeHistogram::getTH1F(FND_MC_Hybrid,NT,"kpkm","Chi2DOF","(40,0,20)",Form("CUT(%s)*CUTWT(rf)",CUTS.Data()));
  hChi2DOFMC_Hybrid_Global->Scale(scale_hybrid);
  hChi2DOFMC_Hybrid_Global->SetMarkerColor(kBlue);
  hChi2DOFMC_Hybrid_Global->Draw("same");

/*
  // Some mass spectra...
  CUTS = DEFAULT_CUTS;
  TH1F* hMkpkm = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","MASS([K+],[K-])","(200,0.95,1.5)",Form("CUT(%s)",CUTS.Data()));
  TH1F* hMkpkmMC = FSModeHistogram::getTH1F(FND_MC,NT,"kpkm","MASS([K+],[K-])","(200,0.95,1.5)",Form("CUT(%s)",CUTS.Data()));
  hMkpkmMC->Scale(hMkpkm->GetMaximum()/hMkpkmMC->GetMaximum());
  hMkpkmMC->SetMarkerColor(kMagenta);

  // Add cut on global chi2 rank to see how much background is removed by vetoing cases where pi+pi- has lower chi2 than K+K- (cross-hypothesis ranking)
  TH1F* hMkpkm_global = FSModeHistogram::getTH1F(FND_DATA,NT,"kpkm","MASS([K+],[K-])","(200,0.95,1.5)",Form("CUT(%s,chi2rankglobal)",CUTS.Data()));
  hMkpkm_global->SetMarkerColor(kBlue);
*/

    
  FSHistogram::dumpHistogramCache();
    
  return;
}