void process(){

  TTree *tree = (TTree*)_file0->Get("pi0pipmisspim__B1_T1_U1_M7_Effic_Tree");
  tree->Process("DSelector_omega_misspim.C+");
}
