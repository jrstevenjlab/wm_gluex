
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
  FSCut::defineCut("chi2rank","Chi2Rank==1");
  FSCut::defineCut("hybridchi2Rank","HybridChi2Rank==1");

}

void skim_period(int period=3){

  TString FND_DATA_kpkm = "/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/data/tree_kpkm__B4_FSROOT_0506*.root"; //subset of data for testing
  TString FND_DATA_pippim = "/volatile/halld/home/jrsteven/flattened/tree_pippim__B4/data/tree_pippim__B4_FSROOT_0506*.root"; //subset of data for testing
  TString FND_DATA_both = "/volatile/halld/home/jrsteven/flattened/tree_*p*m__B4/data/tree_*__B4_FSROOT_0506*.root"; //kpkm + pippim data for cross-hypothesis ranking

  TString FND_SIGMC_kpkm = "/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/akovatsb_kpkmMC__B4_4890/tree_kpkm__B4_FSROOT_050*.root";
  TString FND_SIGMC_pippim = "/volatile/halld/home/jrsteven/flattened/tree_pippim__B4/akovatsb_kpkmMC__B4_4890/tree_pippim__B4_FSROOT_050*.root";
  TString FND_SIGMC_both = "/volatile/halld/home/jrsteven/flattened/tree_*p*m__B4/akovatsb_kpkmMC__B4_4890/tree_*__B4_FSROOT_050*.root";

  ///////////////////////
  // Best Chi2 method: //
  ///////////////////////

  // Rank both hypotheses together so Chi2RankGlobal reflects the best kpkm-vs-pippim choice for each Run/Event group
  FSModeTree::createRankingTree(FND_DATA_both,NT,"","Chi2Rank","Chi2DOF*1000","CUT(rf)","Run","Event");
  FSModeTree::createRankingTree(FND_SIGMC_both,NT,"","Chi2Rank","Chi2DOF*1000","CUT(rf)","Run","Event","MCPxP2*10000");

  // Make skim where kpkm hypothesis had a better Chi2 than pippim hypothesis (one combination for each event)
  FSTree::addFriendTree("Chi2Rank");
  FSModeTree::skimTree(FND_DATA_kpkm,NT,"kpkm",Form("tree_kpkm__B4_BestChi2_SKIM_%02d.root",period),"CUT(chi2,chi2rank,rf)");
  FSModeTree::skimTree(FND_SIGMC_kpkm,NT,"kpkm",Form("tree_kpkm__B4_SIGMC_BestChi2_SKIM_%02d.root",period),"CUT(chi2,chi2rank,rf)");

  ////////////////////
  // Hybrid method: //
  ////////////////////

  // Rank both hypotheses together so Chi2RankGlobal reflects the best kpkm-vs-pippim choice for each Run/Event/Photon beam group
  FSModeTree::createRankingTree(FND_DATA_both,NT,"","HybridChi2Rank","Chi2DOF*1000","1==1","Run","Event","EnPB*100000");
  FSModeTree::createRankingTree(FND_SIGMC_both,NT,"","HybridChi2Rank","Chi2DOF*1000","1==1","Run","Event","MCPxP2*10000+EnPB*100000");

  // Make skim where kpkm hypothesis had a better Chi2 than pippim hypothesis (one combination for each beam photon, requires RFDeltaT subtraction in analysis)
  FSTree::addFriendTree("HybridChi2Rank");
  FSModeTree::skimTree(FND_DATA_kpkm,NT,"kpkm",Form("tree_kpkm__B4_BestHybridChi2_SKIM_%02d.root",period),"CUT(chi2,hybridchi2Rank)");
  FSModeTree::skimTree(FND_SIGMC_kpkm,NT,"kpkm",Form("tree_kpkm__B4_SIGMC_BestHybridChi2_SKIM_%02d.root",period),"CUT(chi2,hybridchi2Rank)");
  
  return;
}

void skim_kpkm(){
  setup();

  //skim_period(3); // 2017-01
  //skim_period(4); // 2018-01
  skim_period(5); // 2018-08
  //skim_period(7); // 2019-11
}
