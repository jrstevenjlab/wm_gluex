
TString NT("ntFSGlueX_MODECODE");
TString SKIM_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("100_110111")->addCategory("pi0pippimkpkm");

  // FIXED CUTS
  FSCut::defineCut("eBeam","(EnPB>8.0&&EnPB<11.6)");
  FSCut::defineCut("chi2","Chi2DOF<20");
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");

  // RANKING CUTS
  FSCut::defineCut("chi2rank","Chi2Rank==1");

  FSCut::defineCut("Phi","MASS([K+],[K-]) < 1.1");
  FSCut::defineCut("Omega","MASS([pi0],[pi+],[pi-]) < 0.9");

  // SET SOME DEFAULT SKIM CUTS
  SKIM_CUTS = "eBeam,chi2,rf,Phi,Omega";
}

void skim(bool bggen=false){

  setup();

  // with mass constraints
  TString FND_DATA = "/volatile/halld/home/jrsteven/flattened/tree_pi0pippimkpkm__B4/data/tree_pi0pippimkpkm__B4_FSROOT_0*.root";
  TString FND_BGGEN = "/cache/halld/home/jrsteven/flattened/phi2pi_gx1/tree_pi0pippimkpkm__B4/*bggen*/tree_pi0pippimkpkm__B4_FSROOT_0*.root";
  TString FND_SIGMC = "/cache/halld/home/jrsteven/flattened/phi2pi_gx1/tree_pi0pippimkpkm__B4/phi2pi_python*/tree_pi0pippimkpkm__B4_FSROOT_0*.root";
  TString FND_ETAMC = "/volatile/halld/home/jrsteven/flattened/tree_pi0pippimkpkm__B4/kpkmeta_mc*/tree_pi0pippimkpkm__B4_FSROOT*.root";
  TString FND_OMEGAMC = "/volatile/halld/home/jrsteven/flattened/tree_pi0pippimkpkm__B4/phiomega_3pi*/tree_pi0pippimkpkm__B4_FSROOT*.root";
  TString CATEGORY = "pi0pippimkpkm";

  FSTree::addFriendTree("Chi2Rank");
  
  //FSModeTree::skimTree(FND_DATA,NT,CATEGORY,"tree_pi0pippimkpkm__B4_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  //FSModeTree::skimTree(FND_BGGEN,NT,CATEGORY,"tree_pi0pippimkpkm__B4_BGGEN_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  //FSModeTree::skimTree(FND_SIGMC,NT,CATEGORY,"tree_pi0pippimkpkm__B4_PSMC_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  FSModeTree::skimTree(FND_ETAMC,NT,CATEGORY,"tree_pi0pippimkpkm__B4_ETAMC_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  FSModeTree::skimTree(FND_OMEGAMC,NT,CATEGORY,"tree_pi0pippimkpkm__B4_OMEGAMC_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  

  //FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  //FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_BGGEN_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  //FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_PSMC_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_ETAMC_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_OMEGAMC_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  
  //FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  //FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_BGGEN_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_BGGEN_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  //FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_PSMC_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_PSMC_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_ETAMC_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_ETAMC_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_OMEGAMC_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_OMEGAMC_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  
  return;
}
