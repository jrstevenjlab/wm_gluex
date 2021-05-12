
void plot_plotter(TString dir = "./", TString reac = "") {

	gStyle->SetOptStat(0);

	const int maxPlots = 9;
	TString plotNames[maxPlots] = {"MOmegaPi","CosTheta","Phi","CosTheta_H","Phi_H","Prod_Ang","MRecoil","MProtonPi","MRecoilPi"};
	const int maxAmps = 6;
	TString ampNames[maxAmps] = {"","_0m","_1p","_1m","_2m","_3m"};
    	TString ampTitles[maxAmps] = {"Fit Result","#[]{0^{#minus}}^{(#pm)} P","#[]{1^{#plus}}^{(#pm)} (S+D)","#[]{1^{#minus}}^{(#pm)} P","#[]{2^{#minus}}^{(#pm)} (P+F)", "#[]{3^{#minus}}^{(#pm)} (F)"};
	TString ampDrawOpt[maxAmps] = {"h","ep","ep","ep","ep","ep"};
    	int ampColors[maxAmps] = {1, 6, 4, 2, 7, 8};
    
    	TFile *f = TFile::Open(dir+"/omegapi_plot.root");

	TCanvas *cc = new TCanvas("cc","cc",800,600);
	cc->Divide(3,3);

    	double textSize = 0.10;
    	TLegend *leg1 = new TLegend(0.1, 0.1, 0.5, 0.9);
    	leg1->SetEntrySeparation(0.01);
	leg1->SetNColumns(2);
	leg1->SetColumnSeparation(1.0);
    	leg1->SetMargin(0.2);
    	leg1->SetFillColor(0);
    	leg1->SetTextSize(textSize);
    	leg1->SetBorderSize(0);
    	leg1->SetLineColor(kWhite);
    
	for(int i=0; i<maxPlots; i++) {
		TH1F *hdat = (TH1F*)f->Get(reac+plotNames[i]+"dat");

		hdat->SetLineColor(kBlack);
		hdat->SetMinimum(0);
		cc->cd(i+1);
	        hdat->SetMarkerStyle(20);
        	hdat->SetMarkerSize(0.5);
	        if(i==0) leg1->AddEntry(hdat, "GlueX Data", "ep");
        	else hdat->Draw();
        
		for(int j=0; j<maxAmps; j++) {
	            TH1F *hacc = (TH1F*)f->Get(reac+plotNames[i]+"acc"+ampNames[j]);
            
       		    if(j==0) {
                	hacc->SetFillColor(31);
	                hacc->SetLineColor(31);
        	        if(i==0) leg1->AddEntry(hacc, ampTitles[j], "f");
                	else hacc->Draw("same" + ampDrawOpt[j]);
	            }
        	    else {
	                hacc->SetLineColor(ampColors[j]);
        	        hacc->SetMarkerColor(ampColors[j]);
                	hacc->SetMarkerSize(0.5);
	                hacc->SetMarkerStyle(20);
        	        if(i==0) leg1->AddEntry(hacc, ampTitles[j], ampDrawOpt[j]);
                	else hacc->Draw("same" + ampDrawOpt[j]);
            	    }
        	}
        
        	if(i!=0) hdat->Draw("same");
	}
    
    	cc->cd(1);
    	leg1->Draw();

	cc->Print(dir+"/fit.pdf");
    
	return;
}
