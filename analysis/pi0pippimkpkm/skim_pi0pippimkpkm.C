
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

void skim_pi0pippimkpkm(bool bggen=false){

  setup();

  TString TREE = "tree_pi0pippimkpkm__B4";
  TString FLAT_BASE = "/sciclone/gluex10/flattened/phi2pi_gx1/tree_pi0pippimkpkm__B4";
  TString FND_DATA = FLAT_BASE+"/data/"+TREE+"_FSROOT_0*.root";
  TString FND_BGGEN = FLAT_BASE+"/*bggen*/"+TREE+"_FSROOT_0*.root";
  TString FND_SIGMC = FLAT_BASE+"/phi2pi_python*/"+TREE+"_FSROOT_0*.root";
  TString FND_ETAMC = FLAT_BASE+"/kpkmeta_mc*/"+TREE+"_FSROOT*.root";
  TString FND_OMEGAMC = FLAT_BASE+"/phiomega_3pi*/"+TREE+"_FSROOT*.root";
  TString CATEGORY = "pi0pippimkpkm";

  FSTree::addFriendTree("Chi2Rank");
  
  FSModeTree::skimTree(FND_DATA,NT,CATEGORY,"tree_pi0pippimkpkm__B4_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  //FSModeTree::skimTree(FND_BGGEN,NT,CATEGORY,"tree_pi0pippimkpkm__B4_BGGEN_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  //FSModeTree::skimTree(FND_SIGMC,NT,CATEGORY,"tree_pi0pippimkpkm__B4_SIGMC_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  //FSModeTree::skimTree(FND_ETAMC,NT,CATEGORY,"tree_pi0pippimkpkm__B4_ETAMC_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  //FSModeTree::skimTree(FND_OMEGAMC,NT,CATEGORY,"tree_pi0pippimkpkm__B4_OMEGAMC_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");
  

  FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  //FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_BGGEN_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  //FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_SIGMC_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  //FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_ETAMC_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  //FSModeTree::createChi2RankingTree("./tree_pi0pippimkpkm__B4_OMEGAMC_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");
  
  FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  //FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_BGGEN_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_BGGEN_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  //FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_SIGMC_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_SIGMC_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  //FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_ETAMC_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_ETAMC_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  //FSModeTree::skimTree("./tree_pi0pippimkpkm__B4_OMEGAMC_GENERAL_SKIM.root",NT,CATEGORY,"tree_pi0pippimkpkm__B4_OMEGAMC_BestChi2_SKIM.root","CUT(rf,chi2rank)");
  
  return;
}
