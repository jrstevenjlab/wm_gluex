
void plot_diagnostic(TString dir = "./", TString reac = "") {

	gStyle->SetOptStat(0);

	const int maxFiles = 4;
	const int maxPlots = 4;
	//TString fileNames[maxFiles] = {"gen_omegapiAmplitude_b1_mix_1pm_diagnostic.root","gen_omegapiFitResult_b1_diagnostic_real.root","gen_omegapiFitResult_b1_diagnostic.root","gen_omegapiFitResult_b1_seed9_diagnostic.root"};
	TString fileNames[maxFiles] = {"gen_omegapi_diagnostic_refl+_real_1pp.root","gen_omegapi_diagnostic_refl-_real_1pp.root","gen_omegapi_diagnostic_refl+_real_1pm.root","gen_omegapi_diagnostic_refl-_real_1pm.root"};
	//TString fileNames[maxFiles] = {"gen_omegapi_diagnostic_refl+_real_1pp.root","gen_omegapi_diagnostic_refl-_real_1pp.root","gen_omegapi_diagnostic_refl+_halfpi_1pp.root","gen_omegapi_diagnostic_refl-_halfpi_1pp.root"};
	TString plotNames[maxPlots] = {"CosTheta_psi","phi_Phi_Prod","PhiH_Psi","PhiH_PsiPrime"};
    
	// open files 
	TFile *f[maxFiles];
	for(int ifile=0; ifile<maxFiles; ifile++) {
		f[ifile] = TFile::Open("testPhase/"+fileNames[ifile]);
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
			h->Rebin2D(3,3);
			if(i==1) {
				h->SetMinimum(10); h->SetMaximum(60);
			}
			h->Draw("colz");

			if(i>1) {
				
				TH1F *h1D = (TH1F*)h->ProjectionX(Form("proj_%d_%d",i,ifile))->Clone();
				h1D->SetMinimum(0);
				h1D->Draw();
			}
				
		}
	}

	//cc->Print(dir+"/fit.pdf");
    
	return;
}
