
TString NT("ntFSGlueX_MODECODE");
TString DEFAULT_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("101_111")->addCategory("pi0pippimeta");

  // DEFINITION OF CUTS:

  // FIXED CUTS
  FSCut::defineCut("unusedE","EnUnusedSh<0.1");
  FSCut::defineCut("unusedTracks","NumUnusedTracks<1");
  FSCut::defineCut("z","ProdVz>=51.2&&ProdVz<=78.8");
  FSCut::defineCut("MM2","abs(RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5))<0.05");
  FSCut::defineCut("eBeam","(EnPB>8.0)");
  FSCut::defineCut("chi2","Chi2DOF<5");
  FSCut::defineCut("chi2rank","Chi2Rank==1");
 
  FSCut::defineCut("t","abs(-1*MASS2([proton],-GLUEXTARGET))<0.5");
  FSCut::defineCut("Metagg","abs(MASS([eta]) - 0.547) < 0.025");
  FSCut::defineCut("Meta3pi","abs(MASS([pi0],[pi+],[pi-]) - 0.547) < 0.025");

  // wrong combo pi0 veto
  FSCut::defineCut("pi0veto","abs(MASS([eta]a,[pi0]a)-0.135)>0.02 && abs(MASS([eta]a,[pi0]b)-0.135)>0.02 && abs(MASS([eta]b,[pi0]a)-0.135)>0.02 && abs(MASS([eta]b,[pi0]b)-0.135)>0.02");
  
  // SET SOME DEFAULT CUTS
  DEFAULT_CUTS = "unusedTracks,unusedE,z,eBeam,chi2,MM2,chi2rank,t,Metagg,Meta3pi,pi0veto";
}

void plots(bool bggen=false){
  bool mc=true;
	
  // Basic plots for flattened and skimmed DATA trees:
  TString FND_DATA = "tree_pi0pippimeta__B4_M17_GENERAL_SKIM.root";

  FSTree::addFriendTree("Chi2Rank");

  setup();
  system("rm -rf plots");  system("mkdir plots");

  TString CUTS;
  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(3,2);
  c1->cd(1);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",unusedE","");
  TH1F* hEnUnusedSh = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","EnUnusedSh","(100,0.0,1.0)",Form("CUT(%s)",CUTS.Data()));
  hEnUnusedSh->SetXTitle("E_{unused}  [GeV/c^{2}]");
  hEnUnusedSh->SetYTitle("Events");
  hEnUnusedSh->Draw();

  TLine* cutUnusedE = new TLine(0.1,0,0.1,hEnUnusedSh->GetMaximum());
  cutUnusedE->SetLineColor(kRed);
  cutUnusedE->Draw("same");

  c1->cd(2);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",z","");
  TH1F* hProdVz = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","ProdVz","(100,0.,100.0)",Form("CUT(%s)",CUTS.Data()));
  hProdVz->SetXTitle("ProdVz  [GeV/c^{2}]");
  hProdVz->SetYTitle("Events");
  hProdVz->Draw();

  TLine* cutVz_low = new TLine(52,0,52,hProdVz->GetMaximum());
  cutVz_low->SetLineColor(kRed);
  cutVz_low->Draw("same");
  TLine* cutVz_hi = new TLine(78,0,78,hProdVz->GetMaximum());
  cutVz_hi->SetLineColor(kRed);
  cutVz_hi->Draw("same");

  c1->cd(3);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",t","");
  TH1F* ht = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","abs(-1*MASS2([proton],-GLUEXTARGET))","(150,0,15)",Form("CUT(%s)",CUTS.Data()));
  ht->SetXTitle("|t| [GeV^{2}]");
  ht->SetYTitle("Entries");
  ht->Draw();

  c1->cd(4);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",eBeam","");
  TH1F* hEnPB = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","EnPB","(125,5,12)",Form("CUT(%s)",CUTS.Data()));
  hEnPB->SetXTitle("E_{beam} [GeV]");
  hEnPB->SetYTitle("Entries");
  hEnPB->Draw();

  c1->cd(5);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",MM2","");
  TH1F* hMM2 = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","RMASS2(GLUEXTARGET,B,-1,-2,-3,-4,-5)","(100,-0.1,0.1)",Form("CUT(%s)",CUTS.Data()));
  hMM2->SetXTitle("Missing Mass Squared [GeV^{2}]");
  hMM2->SetYTitle("Entries");
  hMM2->Draw();

  c1->cd(6);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",chi2,",",");
  TH1F* hChi2DOF = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","Chi2DOF","(100,0,10)",Form("CUT(%s)",CUTS.Data()));
  hChi2DOF->SetTitle("#gamma p #rightarrow #gamma#eta p");
  hChi2DOF->SetXTitle("#chi^{2}/dof");
  hChi2DOF->SetYTitle("Events");
  hChi2DOF->Draw();

  // Some mass spectra...
  TCanvas* c11 = new TCanvas("c11","c11",1200,800);
  c11->Divide(3,2);

  c11->cd(1);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",Metagg","");
  TH1F* hMgg = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","MASS([eta]a,[eta]b)","(50,0.35,0.75)",Form("CUT(%s)",CUTS.Data()));
  hMgg->Draw();

  c11->cd(2);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",Meta3pi","");
  TH1F* hMpi0pippim = FSModeHistogram::getTH1F(FND_DATA,NT,"pi0pippimeta","MASS([pi0],[pi+],[pi-])","(60,0.35,0.85)",Form("CUT(%s)",CUTS.Data()));
  hMpi0pippim->Draw();

  c11->cd(3);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",Metagg","");
  CUTS.ReplaceAll(",Meta3pi","");
  TH2F* hMggVsMpi0pippim = FSModeHistogram::getTH2F(FND_DATA,NT,"pi0pippimeta","MASS([pi0],[pi+],[pi-]):MASS([eta]a,[eta]b)","(50,0.35,0.75,60,0.35,0.85)",Form("CUT(%s)",CUTS.Data()));
  hMggVsMpi0pippim->Draw("colz");
 
  c11->cd(4);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",pi0veto","");
  TH2F* hMgg_alt_etaa = FSModeHistogram::getTH2F(FND_DATA,NT,"pi0pippimeta","MASS([eta]a,[pi0]a):MASS([eta]a,[pi0]b)","(50,0.0,0.5,50,0,0.5)",Form("CUT(%s)",CUTS.Data()));
  hMgg_alt_etaa->Draw("colz");

  c11->cd(5);
  CUTS = DEFAULT_CUTS;
  CUTS.ReplaceAll(",pi0veto","");
  TH2F* hMgg_alt_etab = FSModeHistogram::getTH2F(FND_DATA,NT,"pi0pippimeta","MASS([eta]b,[pi0]a):MASS([eta]b,[pi0]b)","(50,0.0,0.5,50,0,0.5)",Form("CUT(%s)",CUTS.Data()));
  hMgg_alt_etab->Draw("colz");

  FSHistogram::dumpHistogramCache();

  return;
}
