
TString NT("ntFSGlueX_MODECODE");
TString SELECTION_CUTS;

void setup(){
  if (FSModeCollection::modeVector().size() != 0) return;
  FSHistogram::readHistogramCache();
  FSModeCollection::addModeInfo("100000000_100000")->addCategory("kplamb");

  // STANDARD CUTS
  FSCut::defineCut("unusedE","EnUnusedSh<0.1");
  FSCut::defineCut("unusedTracks","NumUnusedTracks<1");
  FSCut::defineCut("z","ProdVz>=51.2&&ProdVz<=78.8");
  FSCut::defineCut("eBeam","(EnPB>8.2&&EnPB<8.8)");
  FSCut::defineCut("rf","abs(RFDeltaT)<2.004");
  FSCut::defineCut("chi2","Chi2DOF<5");
  FSCut::defineCut("chi2rank","Chi2Rank==1");
  FSCut::defineCut("PARA_0","PolarizationAngle==0");
  FSCut::defineCut("PERP_90","PolarizationAngle==90");
  FSCut::defineCut("PARA_135","PolarizationAngle==135");
  FSCut::defineCut("PERP_45","PolarizationAngle==45");
  FSCut::defineCut("AMO","PolarizationAngle==-1");

  FSCut::defineCut("MCeBeam","(MCEnPB>8.2&&MCEnPB<8.8)");

  // TOPOLOGY SPECIFIC
  FSCut::defineCut("MM2","abs(RMASS2(GLUEXTARGET,B,-1,-2))<0.05");

  // SET SOME DEFAULT CUTS
  SELECTION_CUTS = "eBeam,chi2,rf,unusedE,unusedTracks,z,chi2rank";

  // SET SIGNAL and SIDEBAND
  FSCut::defineCut("Lambda2DSB","abs(MASS(1) - 1.115) < 0.010","abs(MASS(1) - 1.115) > 0.015 && abs(MASS(1) - 1.115) < 0.035",0.5);
}

void skimForAmpTools(bool onlyMC = true){

  setup();

  TString FND_DATA = "tree_kplamb__B4_M18_BestChi2_SKIM*.root";
  TString FND_MC = "tree_kplamb__B4_M18_PSMC_BestChi2_SKIM*.root";
  TString FND_MCGEN = "tree_thrown*PSMCGEN.root";
  TString CATEGORY = "kplamb";

  FSTree::addFriendTree("Chi2Rank");

  string treeName = "tree_kplamb__B4_M18";
  string t = "-1*MASS2(1,-GLUEXTARGET)";
  string MCt = "-1*MCMASS2(1,-GLUEXTARGET)";
  string signalCut = "Lambda2DSB";

  // define kinematic bins
  vector<double> mint = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.3};
  vector<double> maxt = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.3, 2.1};
  
  // loop over t bins
  for(uint it=0; it<mint.size(); it++){
	  
	  //if(it != 1) continue;

	  string cutt = Form("abs(%s)>%f && abs(%s)<%f",t.data(),mint[it],t.data(),maxt[it]);
	  string MCcutt = Form("abs(%s)>%f && abs(%s)<%f",MCt.data(),mint[it],MCt.data(),maxt[it]);
	  FSCut::defineCut("t",cutt);
	  FSCut::defineCut("MCt",MCcutt);
	  
	  string dir = Form("out/t_%0.2f_%0.2f/",mint[it],maxt[it]);
	  system(Form("mkdir -p %s",dir.data()));

	  // phasespace accepted MC for AmpTools
	  FSModeTree::skimTree(FND_MC,NT,CATEGORY,Form("%s/%s_PSMC_SIGNAL_SKIM.root",dir.data(),treeName.data()),Form("CUT(%s,%s,t)",SELECTION_CUTS.Data(),signalCut.data()));
	  
	  // phasespace generated MC for AmpTools 
	  FSModeTree::skimTree(FND_MCGEN,NT,CATEGORY,Form("%s/tree_thrown_PSMCGEN_SKIM.root",dir.data()),"CUT(MCeBeam,MCt)");
	  
	  if(onlyMC) continue;
	  
	  // apply general event selection cuts for signal and background (slowest part...)
	  FSModeTree::skimTree(FND_DATA,NT,CATEGORY,Form("out/%s_BestChi2_SKIM.root",treeName.data()), Form("CUT(%s,t)",SELECTION_CUTS.Data()));
	  
	  // loop over polarization orientations
	  string polLabel[5] = {"PARA_0","PERP_90","PARA_135","PERP_45","AMO"};
	  for(int i=0; i<5; i++) {
		  
		  // sort by polarization
		  FSModeTree::skimTree(Form("out/%s_BestChi2_SKIM.root",treeName.data()), NT,"",Form("out/%s_%s_BestChi2_SKIM.root",treeName.data(),polLabel[i].data()),Form("CUT(%s)",polLabel[i].data()));
		  
		  // signal tree for AmpTools
		  FSModeTree::skimTree(Form("out/%s_%s_BestChi2_SKIM.root",treeName.data(),polLabel[i].data()), NT,"",Form("%s/%s_%s_SIGNAL_SKIM.root",dir.data(),treeName.data(),polLabel[i].data()),Form("CUT(%s)",signalCut.data()));
		  
		  // background tree for AmpTools
		  FSModeTree::skimTree(Form("out/%s_%s_BestChi2_SKIM.root",treeName.data(),polLabel[i].data()), NT,"",Form("%s/%s_%s_SIDEBANDS_SKIM.root",dir.data(),treeName.data(),polLabel[i].data()),Form("CUTSBWT(%s)",signalCut.data()));
		  
		  // background weights for AmpTools
		  vector< pair<TString,TString> > friendTreeContents;
		  friendTreeContents.push_back(pair<TString,TString>("weight",Form("CUTSBWT(%s)",signalCut.data())));
		  FSModeTree::createFriendTree(Form("%s/%s_%s_SIDEBANDS_SKIM.root",dir.data(),treeName.data(),polLabel[i].data()),NT,"","weight",friendTreeContents);
		  
	  } // end polarization loop 
	  
	  // clear files for next bin
	  system("rm out/tree_*");
	  system(Form("rm %s/*.Chi2Rank",dir.data()));
  }

  return;
}
