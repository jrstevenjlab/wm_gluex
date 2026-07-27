
void fit() {

	TCanvas *cc = new TCanvas("cc", "cc", 800, 400);
	cc->Divide(2,1);
	cc->cd(1);

	TFile *f = TFile::Open("/sciclone/gluex10/jrstevens01/gpi0pippim/hist_gpi0pippim_sum.root");
	TH1F *h = (TH1F*)f->Get("Pi0GammaMass");
	//TH2F *h2 = (TH2F*)f->Get("Pi0PiMinusVsPi0GammaDiff");
	//TH1F *h = (TH1F*)h2->ProjectionX();
	h->SetMinimum(0);
	h->Draw();
	double binWidth = h->GetBinWidth(1);

	double low = 0.5;
	double high = 1.1;
	TF1 *fit = new TF1("fit", "gaus(0) + pol2(3)", low, high);
	fit->FixParameter(1, 0.782);
	fit->SetParameter(2, 0.015);
	h->Fit(fit, "", "", low, high);

	TF1 *signal = new TF1("signal", "gaus", low, high);
	for(int i=0; i<3; i++)
		signal->SetParameter(i, fit->GetParameter(i));

	cout<<signal->Integral(0.742, 0.822)/binWidth<<" "<<fit->Integral(0.742, 0.822)/binWidth<<endl;

	cc->cd(2);
	TH1F *hOmegaPiMass = (TH1F*)f->Get("OmegaPiMinusMass");
	hOmegaPiMass->Rebin(75); 
	hOmegaPiMass->SetMinimum(0);
	hOmegaPiMass->Draw();

	return;	

}
