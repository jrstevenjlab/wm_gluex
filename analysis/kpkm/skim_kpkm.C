
TString NT("ntFSGlueX_MODECODE");
TString SKIM_CUTS;

void createChi2RankHybridFriendTreesPerFile(const TString& pattern){
  TString cmd = Form("ls -1 %s 2>/dev/null", pattern.Data());
  TString files = gSystem->GetFromPipe(cmd.Data());

  if (files.IsNull()) {
    cout << "No files matched pattern: " << pattern << endl;
    return;
  }

  TObjArray* lines = files.Tokenize("\n");
  for (int i = 0; i < lines->GetEntriesFast(); ++i) {
    TString file = ((TObjString*)lines->At(i))->GetString();
    if (file.IsNull()) continue;
    cout << "Creating Chi2RankHybrid friend tree for: " << file << endl;
    FSModeTree::createRankingTree(file,NT,"","Chi2RankHybrid","Chi2*1000","1==1","Run","Event","EnPB*1000");
  }
  delete lines;
}

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("100_110000")->addCategory("kpkm");

  // FIXED CUTS
  FSCut::defineCut("chi2","Chi2DOF<20");

  // RANKING CUTS
  FSCut::defineCut("chi2rankhybrid","Chi2RankHybrid==1");
}

void skim_period(int period=3){

  TString FND_DATA = "/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/data/tree_kpkm__B4_FSROOT_0506*.root"; //subset of data for testing
  //TString FND_DATA = Form("/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/data/tree_kpkm__B4_FSROOT_%02d*.root", period);
  TString FND_BGGEN = Form("/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/*bggen*/tree_kpkm__B4_FSROOT_%02d*.root", period);
  TString FND_SIGMC = Form("/volatile/halld/home/jrsteven/flattened/tree_kpkm__B4/*phi*/tree_kpkm__B4_FSROOT_%02d*.root", period);
  TString CATEGORY = "kpkm";

  FSTree::addFriendTree("Chi2RankHybrid");

  //createChi2RankHybridFriendTreesPerFile(FND_DATA);
  //createChi2RankHybridFriendTreesPerFile(FND_BGGEN);
  //createChi2RankHybridFriendTreesPerFile(FND_SIGMC);
  //FSModeTree::createChi2RankingTree(Form("./tree_kpkm__B4_PSMC_GENERAL_SKIM_%02d.root",period),NT,"","Chi2RankHybrid","Chi2*1000","1==1","Run","Event","EnPB*1000");
  
  FSModeTree::skimTree(FND_DATA,NT,CATEGORY,Form("tree_kpkm__B4_Hybrid_SKIM_%02d.root",period),"CUT(chi2rankhybrid)");
  //FSModeTree::skimTree(FND_BGGEN,NT,CATEGORY,Form("tree_kpkm__B4_BGGEN_BestChi2_SKIM_%02d.root",period),"CUT(chi2rankhybrid)");
  //FSModeTree::skimTree(FND_SIGMC,NT,CATEGORY,Form("tree_kpkm__B4_PSMC_BestChi2_SKIM_%02d.root",period),"CUT(chi2rankhybrid)");

  return;
}

void skim_kpkm(){
  setup();

  //skim_period(3); // 2017-01
  //skim_period(4); // 2018-01
  skim_period(5); // 2018-08
  //skim_period(7); // 2019-11
}