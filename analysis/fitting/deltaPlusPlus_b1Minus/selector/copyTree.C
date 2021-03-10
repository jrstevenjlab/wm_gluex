void copyTree(TString treeName = "anglesOmegaPiAmplitude") {
 gROOT->Reset();

 TString treeFileName = treeName + ".root";
 TFile f(Form("allRecoil/%s", treeFileName.Data()));
 TTree *T = (TTree*)f.Get("kin");

 TFile f2(treeFileName,"recreate");
 TTree *T2 = T->CopyTree("MRecoil<1.35");

 T2->Write();
 f2.Close();
}
