
TString NT("ntFSGlueX_MODECODE");
TString SKIM_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("100000000_100000")->addCategory("kplamb");

  // FIXED CUTS
  FSCut::defineCut("eBeam","EnPB>8.0");
  FSCut::defineCut("chi2","Chi2DOF<20");
  
  // Select Delta++ region
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");

  // SET SOME DEFAULT SKIM CUTS
  SKIM_CUTS = "eBeam,chi2,rf";
}

void skim_period(int period=3){

  // with mass constraints
  TString FND_DATA = Form("/sciclone/gluex10/flattened/kpi_gx1_pwa/ver01/tree_kplamb__B4_M18/data/tree_kplamb__B4_M18_FSROOT_%02d*.root", period);
  TString FND_BGGEN = Form("/sciclone/gluex10/flattened/kpi_gx1_pwa/ver01/tree_kplamb__B4_M18/*bggen*/tree_kplamb__B4_M18_FSROOT_%02d*.root", period);
  TString FND_PSMC = Form("/sciclone/gluex10/flattened/kpi_gx1_pwa/ver01/tree_kplamb__B4_M18/*phasespace*/tree_kplamb__B4_M18_FSROOT_%02d*.root", period);
  TString CATEGORY = "kplamb";

  FSModeTree::skimTree(FND_DATA,NT,CATEGORY,Form("tree_kplamb__B4_M18_GENERAL_SKIM_%02d.root", period),"CUT("+SKIM_CUTS+")");
  //FSModeTree::skimTree(FND_BGGEN,NT,CATEGORY,Form("tree_kplamb__B4_M18_BGGEN_GENERAL_SKIM_%02d.root", period),"CUT("+SKIM_CUTS+")");
  FSModeTree::skimTree(FND_PSMC,NT,CATEGORY,Form("tree_kplamb__B4_M18_PSMC_GENERAL_SKIM_%02d.root",period),"CUT("+SKIM_CUTS+")");

  FSModeTree::createChi2RankingTree(Form("./tree_kplamb__B4_M18_GENERAL_SKIM_%02d.root",period),NT,CATEGORY,"CUT(rf)");
  //FSModeTree::createChi2RankingTree(Form("./tree_kplamb__B4_M18_BGGEN_GENERAL_SKIM_%02d.root",period),NT,CATEGORY,"CUT(rf)");
  FSModeTree::createChi2RankingTree(Form("./tree_kplamb__B4_M18_PSMC_GENERAL_SKIM_%02d.root",period),NT,CATEGORY,"CUT(rf)");

  // skim only best chi2 candidates
  FSCut::defineCut("chi2rank","Chi2Rank==1");
  FSTree::addFriendTree("Chi2Rank");
  FSModeTree::skimTree(Form("tree_kplamb__B4_M18_GENERAL_SKIM_%02d.root",period),NT,CATEGORY,Form("tree_kplamb__B4_M18_BestChi2_SKIM_%02d.root",period),"CUT(chi2rank)");
  //FSModeTree::skimTree(Form("tree_kplamb__B4_M18_BGGEN_GENERAL_SKIM_%02d.root",period),NT,CATEGORY,Form("tree_kplamb__B4_M18_BGGEN_BestChi2_SKIM_%02d.root",period),"CUT(chi2rank)");
  FSModeTree::skimTree(Form("tree_kplamb__B4_M18_PSMC_GENERAL_SKIM_%02d.root",period),NT,CATEGORY,Form("tree_kplamb__B4_M18_PSMC_BestChi2_SKIM_%02d.root",period),"CUT(chi2rank)");

  return;
}

void skim(){
	setup();

	skim_period(3);
	skim_period(4);
	skim_period(5);
	skim_period(7);	
}
