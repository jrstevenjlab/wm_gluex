
void plot_diagnostic(TString dir = "./", TString reac = "") {

	gStyle->SetOptStat(0);

	const int maxFiles = 3;
	const int maxPlots = 4;
	TString fileNames[maxFiles] = {"gen_omegapiAmplitude_b1_mix_1pm_diagnostic.root","gen_omegapiFitResult_b1_real_diagnostic.root","gen_omegapiFitResult_b1_diagnostic.root"};
	TString plotNames[maxPlots] = {"CosTheta_phi","CosThetaH_phiH","CosTheta_Phi_Prod","phi_Phi_Prod"};
    
	// open files 
	TFile *f[maxFiles];
	for(int ifile=0; ifile<maxFiles; ifile++) {
		f[ifile] = TFile::Open(fileNames[ifile]);
		//f[ifile]->ls();
	}

	TCanvas *cc[maxPlots];

	// loop over diagnostic plots to show
	for(int i=0; i<maxPlots; i++) {
		
		cc[i] = new TCanvas(Form("cc%d",i),Form("cc%d",i),800,600);
		cc[i]->Divide(2,2);

		for(int ifile=0; ifile<maxFiles; ifile++) {			
			cc[i]->cd(ifile+1);
			TH2F *h = (TH2F*)f[ifile]->Get(plotNames[i]);
			h->Rebin2D(2,2);
			h->Draw("colz");
			//h->ProjectionX()->Draw();
		}
	}

	//cc->Print(dir+"/fit.pdf");
    
	return;
}
