void plotTopology() {

	gStyle->SetOptStat("");
	gStyle->SetPaintTextFormat("4.3f"); // set precision for histogram text
	gStyle->SetLegendBorderSize(0);

	TFile *f = TFile::Open("/sciclone/gluex10/lelorenti/omega_misspim/hist_bggen_S2018_new_BG_total.root");
	//TFile *f = TFile::Open("/sciclone/home20/lelorenti/wm_gluex/analysis/omega_misspim/selector/hist_omega_misspim_gen_30730.acc.root");

	// plot percentages for different thrown topologies
	TCanvas *bb = new TCanvas("bb","bb",600,400);
	TH1F *hThrownTopologies = (TH1F*)f->Get("hThrownTopologies");
	hThrownTopologies->GetXaxis()->LabelsOption(">"); // order by most common topology
	hThrownTopologies->GetXaxis()->SetRangeUser(0, 10); // only plot first 20 topologies
	hThrownTopologies->Scale(100./hThrownTopologies->GetEntries()); // turn histogram into percentage
	hThrownTopologies->GetYaxis()->SetTitle("Thrown Topology %");
	hThrownTopologies->Draw("htext");

	// draw 3Pi invariant mass histograms for different thrown topologies
	TCanvas *cc = new TCanvas("cc","cc",600,400);
	TCanvas *cc2 = new TCanvas("cc2","cc2",600,400);
	TLegend *leg = new TLegend(0.12,0.65,0.35,0.89);
	TLegend *leg2 = new TLegend(0.12,0.65,0.35,0.89);
	leg->SetNColumns(2);
	leg2->SetNColumns(2);
	THStack *hStack = new THStack("hstack","");

	int locNumTopologies = 11;
	cc->cd();
	for(int i=0; i<locNumTopologies; i++) {
		TH1I *hMassThrownTopology = (TH1I*)f->Get(Form("h3piMass_ThrownTopology_%d", i));
		if(!hMassThrownTopology) break;
		hMassThrownTopology->SetLineColor(kBlack+i);
		hMassThrownTopology->SetFillColor(kBlack+i);
		hMassThrownTopology->Rebin(10);

		TString locLegendTitle = hMassThrownTopology->GetTitle();
		locLegendTitle.ReplaceAll("3pi Mass Topology:","");
		leg->AddEntry(hMassThrownTopology,locLegendTitle,"f");

		if(i==0) {
			hMassThrownTopology->SetTitle("#pi^{#plus}#pi^{#minus}#pi^{0} Invariant Mass; M_{#pi^{#plus}#pi^{#minus}#pi^{0}} (GeV)");
			hMassThrownTopology->Draw();
		}
		else 
			hMassThrownTopology->Draw("same");
		hStack->Add(hMassThrownTopology);
	}

	leg->Draw("same");

	cc2->cd();
	for(int i=0; i<locNumTopologies; i++) {
		TH1I *hMassThrownTopology = (TH1I*)f->Get(Form("hMMOP_ThrownTopology_%d", i));
		if(!hMassThrownTopology) break;
		hMassThrownTopology->SetLineColor(kBlack+i);
		hMassThrownTopology->SetFillColor(kBlack+i);
		hMassThrownTopology->Rebin(10);

		TString locLegendTitle = hMassThrownTopology->GetTitle();
		locLegendTitle.ReplaceAll("MMOP Topology:","");
		leg2->AddEntry(hMassThrownTopology,locLegendTitle,"f");

		if(i==0) {
			hMassThrownTopology->SetTitle("MMOP (GeV)");
			hMassThrownTopology->Draw();
		}
		else 
			hMassThrownTopology->Draw("same");
		//hStack->Add(hMassThrownTopology);
	}

	leg2->Draw("same");

	TCanvas *dd = new TCanvas("dd","dd",600,400);
	hStack->Draw();
	leg->Draw("same");

	// Save pdf images of plots
	bb->Print("topologies/TopologyPercentage2018_01_new.pdf");
	cc->Print("topologies/Topology3PiMass2018_01_new.pdf");
	cc2->Print("topologies/TopologyMMOP2018_01_new.pdf");

	return;
}
