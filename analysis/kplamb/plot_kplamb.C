
TString NT("ntFSGlueX_MODECODE");
TString DEFAULT_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("100000000_100000")->addCategory("kplamb");

  // FIXED CUTS
  FSCut::defineCut("unusedE","EnUnusedSh<0.1");
  FSCut::defineCut("unusedTracks","NumUnusedTracks<1");
  FSCut::defineCut("z","ProdVz>=51.2&&ProdVz<=78.8");
  FSCut::defineCut("MM2","abs(RMASS2(GLUEXTARGET,B,-1,-2))<0.05");
  FSCut::defineCut("eBeam","(EnPB>8.2&&EnPB<8.8)");
  FSCut::defineCut("chi2","Chi2DOF<5");
  FSCut::defineCut("chi2rank","Chi2Rank==1");
    
  FSCut::defineCut("t","abs(-1*MASS2(1,-GLUEXTARGET))<1.0");
  FSCut::defineCut("Lambda","abs(MASS(1) - 1.115) < 0.010");
  FSCut::defineCut("LambdaSB","abs(MASS(1) - 1.090) < 0.010 || abs(MASS(1) - 1.140) < 0.010");

  FSCut::defineCut("Lambda2DSB","abs(MASS(1) - 1.115) < 0.010","abs(MASS(1) - 1.115) > 0.015 && abs(MASS(1) - 1.115) < 0.035",0.5);
  
  // Select RF signal region
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");

  // SET SOME DEFAULT SKIM CUTS
  DEFAULT_CUTS = "eBeam,chi2,rf,unusedE,unusedTracks,z,MM2,t,Lambda,chi2rank";
}

void plot_kplamb(bool bggen=false){

	//gStyle->SetCanvasPreferGL();
  gStyle->SetPadBottomMargin(0.22);
  gStyle->SetPadLeftMargin(0.25);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.15);

  setup();

  // with mass constraints
  TString FND_DATA = "tree_kplamb__B4_M18_BestChi2_SKIM*.root";
  TString FND_BGGEN = "tree_kplamb__B4_M18_BestChi2_SKIM*.root";
  TString FND_MC = "tree_kplamb__B4_M18_PSMC_BestChi2_SKIM*.root"; 
  TString CATEGORY = "kplamb";

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
    TH1F* hEnUnusedSh = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","EnUnusedSh","(100,0.0,1.0)",Form("CUT(%s)",CUTS.Data()));
    hEnUnusedSh->SetXTitle("E_{unused}  [GeV/c^{2}]");
    hEnUnusedSh->SetYTitle("Events");
    hEnUnusedSh->Draw();
    TH1F* hEnUnusedShMC = FSModeHistogram::getTH1F(FND_MC,NT,"kplamb","EnUnusedSh","(100,0.0,1.0)",Form("CUT(%s)",CUTS.Data()));
    hEnUnusedShMC->Scale(hEnUnusedSh->GetMaximum()/hEnUnusedShMC->GetMaximum());
    hEnUnusedShMC->SetMarkerColor(kMagenta);
    hEnUnusedShMC->Draw("same");
    TLine* cutUnusedE = new TLine(0.1,0,0.1,hEnUnusedSh->GetMaximum());
    cutUnusedE->SetLineColor(kRed);
    cutUnusedE->Draw("same");
    
    c1->cd(2);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",z","");
    TH1F* hProdVz = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","ProdVz","(100,0.,100.0)",Form("CUT(%s)",CUTS.Data()));
    hProdVz->SetXTitle("ProdVz  [cm]");
    hProdVz->SetYTitle("Events");
    hProdVz->Draw();
    TH1F* hProdVzMC = FSModeHistogram::getTH1F(FND_MC,NT,"kplamb","ProdVz","(100,0.,100.0)",Form("CUT(%s)",CUTS.Data()));
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
    TH1F* htk = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","abs(-1*MASS2(1,-GLUEXTARGET))","(100,0,5)",Form("CUT(%s)",CUTS.Data()));
    htk->SetXTitle("|-t| [GeV^{2}]");
    htk->SetYTitle("Entries");
    htk->Draw();
    htk->Fit("expo","","",0.5,1.0);
    TH1F* htkMC = FSModeHistogram::getTH1F(FND_MC,NT,"kplamb","abs(-1*MASS2(1,-GLUEXTARGET))","(100,0,5)",Form("CUT(%s)",CUTS.Data()));
    htkMC->Scale(htk->GetMaximum()/htkMC->GetMaximum());
    htkMC->SetMarkerColor(kMagenta);
    htkMC->Draw("same");
    htkMC->Fit("expo","","",0.5,1.0);

    c1->cd(4);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll("eBeam,","");
    TH1F* hEnPB = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","EnPB","(125,5,12)",Form("CUT(%s)",CUTS.Data()));
    hEnPB->SetXTitle("E_{beam} [GeV]");
    hEnPB->SetYTitle("Entries");
    hEnPB->Draw();
    TH1F* hEnPBMC = FSModeHistogram::getTH1F(FND_MC,NT,"kplamb","EnPB","(125,5,12)",Form("CUT(%s)",CUTS.Data()));
    hEnPBMC->Scale(hEnPB->GetMaximum()/hEnPBMC->GetMaximum());
    hEnPBMC->SetMarkerColor(kMagenta);
    hEnPBMC->Draw("same");

    c1->cd(5);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",MM2","");
    TH1F* hMM2 = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","RMASS2(GLUEXTARGET,B,-1,-2)","(100,-0.1,0.1)",Form("CUT(%s)",CUTS.Data()));
    hMM2->SetXTitle("Missing Mass Squared [GeV^{2}]");
    hMM2->SetYTitle("Entries");
    hMM2->Draw();
    TH1F* hMM2MC = FSModeHistogram::getTH1F(FND_MC,NT,"kplamb","RMASS2(GLUEXTARGET,B,-1,-2)","(100,-0.1,0.1)",Form("CUT(%s)",CUTS.Data()));
    hMM2MC->Scale(hMM2->GetMaximum()/hMM2MC->GetMaximum());
    hMM2MC->SetMarkerColor(kMagenta);
    hMM2MC->Draw("same");

    c1->cd(6);
    CUTS = DEFAULT_CUTS;
    CUTS.ReplaceAll(",chi2,",",");
    TH1F* hChi2DOF = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    hChi2DOF->SetXTitle("#chi^{2}/dof");
    hChi2DOF->SetYTitle("Events");
    hChi2DOF->Draw();
    TH1F* hChi2DOFMC = FSModeHistogram::getTH1F(FND_MC,NT,"kplamb","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    hChi2DOFMC->Scale(hChi2DOF->GetMaximum()/hChi2DOFMC->GetMaximum());
    hChi2DOFMC->SetMarkerColor(kMagenta);
    hChi2DOFMC->Draw("same");
    
    CUTS.ReplaceAll(",Lambda",",LambdaSB");
    TH1F* hChi2DOFSB = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    hChi2DOFSB->SetMarkerColor(kBlue);
    hChi2DOFSB->Draw("same");

    TLegend *legCuts = new TLegend(0.5,0.5,0.85,0.88);
    legCuts->AddEntry(hChi2DOF,"GlueX Data","pe");
    legCuts->AddEntry(hChi2DOFSB,"Sideband Data","pe");
    legCuts->AddEntry(hChi2DOFMC,"Phasespace MC","pe");
    legCuts->Draw("same");

    c1->Print("plots/kplamb_cuts.pdf");

    if(bggen) {
	    cb->cd(1);   
	    CUTS = DEFAULT_CUTS;
	    CUTS.ReplaceAll(",chi2,",",");
	    TH1F* hChi2DOF_BGGEN = FSModeHistogram::getTH1F(FND_BGGEN,NT,"kplamb","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
	    hChi2DOF_BGGEN->SetXTitle("#chi^{2}/dof");
	    hChi2DOF_BGGEN->SetYTitle("Events");
	    hChi2DOF_BGGEN->Draw();
	    FSModeHistogram::drawMCComponentsSame(FND_BGGEN,NT,"kplamb","Chi2DOF","(40,0,20)",Form("CUT(%s)",CUTS.Data()));
    }

    // Momentum and angle distributions
    CUTS = "eBeam,rf,Lambda"; //DEFAULT_CUTS;
    FSCut::defineCut("t05","abs(-1*MASS2(1,-GLUEXTARGET))<0.5");
    CUTS += ",t05";

    TH2F* hLambdaMomentumVsPathLengthSignal = FSModeHistogram::getTH2F(FND_DATA,NT,"kplamb","RMOMENTUM(1):VeeLP1","(100,-50,50,100,0,1)",Form("CUT(%s)",CUTS.Data()));
    TH2F* hLambdaMomentumVsPathLengthMC = FSModeHistogram::getTH2F(FND_MC,NT,"kplamb","RMOMENTUM(1):VeeLP1","(100,-50,50,100,0,1)",Form("CUT(%s)",CUTS.Data()));

    TH2F* hProtonMomentumVsThetaSignal = FSModeHistogram::getTH2F(FND_DATA,NT,"kplamb","RMOMENTUM(1a):acos(RCOSINE(1a))*180/3.141","(180,0,90,100,0,1)",Form("CUT(%s)",CUTS.Data()));
    TH2F* hProtonMomentumVsThetaMC = FSModeHistogram::getTH2F(FND_MC,NT,"kplamb","RMOMENTUM(1a):acos(RCOSINE(1a))*180/3.141","(180,0,90,100,0,1)",Form("CUT(%s)",CUTS.Data()));
    TH2F* hPiMinusMomentumVsThetaSignal = FSModeHistogram::getTH2F(FND_DATA,NT,"kplamb","RMOMENTUM(1b):acos(RCOSINE(1b))*180/3.141","(180,0,90,100,0,0.5)",Form("CUT(%s)",CUTS.Data()));
    TH2F* hPiMinusMomentumVsThetaMC = FSModeHistogram::getTH2F(FND_MC,NT,"kplamb","RMOMENTUM(1b):acos(RCOSINE(1b))*180/3.141","(180,0,90,100,0,0.5)",Form("CUT(%s)",CUTS.Data()));

    TH2F* hProtonMomentumVsPiMinusMomentumSignal = FSModeHistogram::getTH2F(FND_DATA,NT,"kplamb","RMOMENTUM(1a):RMOMENTUM(1b)","(100,0,0.5,100,0,1)",Form("CUT(%s)",CUTS.Data()));
    TH2F* hProtonMomentumVsPiMinusMomentumMC = FSModeHistogram::getTH2F(FND_MC,NT,"kplamb","RMOMENTUM(1a):RMOMENTUM(1b)","(100,0,0.5,100,0,1)",Form("CUT(%s)",CUTS.Data()));

    TH2F* hProtonThetaVsPiMinusThetaSignal = FSModeHistogram::getTH2F(FND_DATA,NT,"kplamb","acos(COSINE(1a))*180/3.141:acos(COSINE(1b))*180/3.141","(100,0,90,100,0,90)",Form("CUT(%s)",CUTS.Data()));
    TH2F* hProtonThetaVsPiMinusThetaMC = FSModeHistogram::getTH2F(FND_MC,NT,"kplamb","acos(COSINE(1a))*180/3.141:acos(COSINE(1b))*180/3.141","(100,0,90,100,0,90)",Form("CUT(%s)",CUTS.Data()));

    //hProtonThetaMC->Scale(hProtonThetaSignal->GetMaximum()/hProtonThetaMC->GetMaximum());
    //hProtonThetaMC->SetMarkerColor(kMagenta);
    

    // Some mass spectra...
    
    // protonpi mass
    CUTS = DEFAULT_CUTS;
    TH1F* hMpimpSignal = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","MASS(1)","(100,1.08,1.155)",Form("CUT(%s)",CUTS.Data()));
    TH1F* hMpimpMC = FSModeHistogram::getTH1F(FND_MC,NT,"kplamb","MASS(1)","(100,1.08,1.155)",Form("CUT(%s)",CUTS.Data()));
    hMpimpMC->Scale(hMpimpSignal->GetMaximum()/hMpimpMC->GetMaximum());
    hMpimpMC->SetMarkerColor(kMagenta);
    
    CUTS.ReplaceAll(",Lambda",",LambdaSB");
    TH1F* hMpimpSB = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","MASS(1)","(100,1.08,1.155)",Form("CUT(%s)",CUTS.Data()));
    hMpimpSB->SetFillColor(kBlue);

    CUTS.ReplaceAll(",LambdaSB","");
    TH1F* hMpimp = FSModeHistogram::getTH1F(FND_DATA,NT,"kplamb","MASS(1)","(100,1.08,1.155)",Form("CUT(%s)",CUTS.Data()));

    if(bggen) {
	    cb->cd(2);    
	    TH1F* hMpimp_BGGEN = FSModeHistogram::getTH1F(FND_BGGEN,NT,"kplamb","MASS(1)","(100,1.08,1.155)",Form("CUT(%s)",CUTS.Data()));
	    hMpimp_BGGEN->SetXTitle("Mass(p#pi^{-}) [GeV]");
	    hMpimp_BGGEN->SetYTitle("Events");
	    hMpimp_BGGEN->Draw();
	    FSModeHistogram::drawMCComponentsSame(FND_BGGEN,NT,"kplamb","MASS(1)","(100,1.08,1.155)",Form("CUT(%s)",CUTS.Data()));
    }

    FSHistogram::dumpHistogramCache();

    if(bggen) cb->Print("plots/kplamb_bggen.pdf");
    
    hMpimp->SetTitle("; Mass(p#pi^{-}) [GeV]; Counts");
    hMpimp->GetXaxis()->SetNdivisions(5);
    hMpimp->GetYaxis()->SetNdivisions(5);
    hMpimp->GetXaxis()->CenterTitle(true);
    hMpimp->GetYaxis()->CenterTitle(true);
    hMpimp->SetLabelSize(0.07,"xyz");
    hMpimp->SetTitleSize(0.07,"xyz");

    TLegend *leg = new TLegend(0.6,0.6,0.85,0.88);
    leg->AddEntry(hMpimp,"Data","pe");
    leg->AddEntry(hMpimpSB,"Sideband","f");
    leg->AddEntry(hMpimpMC,"PS MC","pe");
    
    TCanvas *c22 = new TCanvas("c22","c22",1000,600);
    c22->Divide(3,2);
    c22->cd(1);
    hMpimp->Draw();
    hMpimpMC->Draw("same");
    hMpimpSignal->Draw("same");
    hMpimpSB->SetFillColorAlpha(kBlue, 0.2);
    hMpimpSB->Draw("h same");
    hMpimp->Draw("same");
    leg->Draw("same");

    c22->Print("plots/kplamb_lambda.pdf");
    
    return;
}
