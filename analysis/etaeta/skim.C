
TString NT("ntFSGlueX_MODECODE");
TString SKIM_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("101_111")->addCategory("pi0pippimeta");

  // FIXED CUTS
  FSCut::defineCut("eBeam","(EnPB>8.0)");
  FSCut::defineCut("chi2","Chi2DOF<20");

  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");
  FSCut::defineCut("eta3pi","MASS([pi0],[pi+],[pi-])>0.35 && MASS([pi0],[pi+],[pi-])<0.85");

  // SET SOME DEFAULT SKIM CUTS
  SKIM_CUTS = "eBeam,chi2,rf,eta3pi";
}

void skim(){

  setup();

  // with mass constraints
  TString FND_DATA = "/volatile/halld/home/jrsteven/flattened/tree_pi0pippimeta__B4_M17/data/tree_pi0pippimeta__B4_M17_FSROOT_*.root";
  TString CATEGORY = "pi0pippimeta";

  FSModeTree::skimTree(FND_DATA,NT,CATEGORY,"tree_pi0pippimeta__B4_M17_GENERAL_SKIM.root","CUT("+SKIM_CUTS+")");

  FSModeTree::createChi2RankingTree("./tree_pi0pippimeta__B4_M17_GENERAL_SKIM.root",NT,CATEGORY,"CUT(rf)");

  return;
}
