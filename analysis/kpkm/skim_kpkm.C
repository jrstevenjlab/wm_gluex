
TString NT("ntFSGlueX_MODECODE");
TString SKIM_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();

  // define modes to consider in Mode Collection
  FSModeCollection::addModeInfo("100_110000")->addCategory("kpkm");
  FSModeCollection::addModeInfo("100_110")->addCategory("pippim");
  FSModeCollection::display();

  // FIXED CUTS
  FSCut::defineCut("chi2","Chi2DOF<20");
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");

  // RANKING CUTS
  FSCut::defineCut("chi2rankhybrid","Chi2RankHybrid==1");
  FSCut::defineCut("chi2rank","Chi2Rank==1 && Chi2RankGlobal==1");
}

void skim_period(int period=3){

  TString FND_DATA_kpkm = "/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/data/tree_kpkm__B4_FSROOT_0506*.root"; //subset of data for testing
  TString FND_DATA_pippim = "/volatile/halld/home/jrsteven/flattened/tree_pippim__B4/data/tree_pippim__B4_FSROOT_0506*.root"; //subset of data for testing
  TString FND_DATA_both = "/volatile/halld/home/jrsteven/flattened/tree_*p*m__B4/data/tree_*__B4_FSROOT_0506*.root"; //kpkm + pippim data for cross-hypothesis ranking

  //TString FND_BGGEN_kpkm = Form("/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/*bggen*/tree_kpkm__B4_FSROOT_%02d*.root", period);
  //TString FND_BGGEN_pippim = Form("/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/*bggen*/tree_kpkm__B4_FSROOT_%02d*.root", period);
  //TString FND_SIGMC_kpkm = Form("/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/*phi*/tree_kpkm__B4_FSROOT_%02d*.root", period);
  //TString FND_SIGMC_pippim = Form("/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/*phi*/tree_kpkm__B4_FSROOT_%02d*.root", period);

  // Rank both hypotheses together so Chi2RankGlobal reflects the best
  // kpkm-vs-pippim choice for each Run/Event group
  FSModeTree::createRankingTree(FND_DATA_both,NT,"","Chi2Rank","Chi2DOF*1000","CUT(rf)");
  
  FSModeTree::skimTree(FND_DATA_kpkm,NT,"kpkm",Form("tree_kpkm__B4_BestChi2_SKIM_%02d.root",period),"CUT(chi2,rf)");

  return;
}

void skim_kpkm(){
  setup();

  //skim_period(3); // 2017-01
  //skim_period(4); // 2018-01
  skim_period(5); // 2018-08
  //skim_period(7); // 2019-11
}
